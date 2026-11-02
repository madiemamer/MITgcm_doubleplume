#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_CALC
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_CALC(
     I     futureTime, futureIter,
     I     myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_CALC
C     | o For every plume cell (plumeMask /= 0) solve the double
C     |   plume for that column (DOUBLEPLUME_COLUMN, every
C     |   DP_nUpdate steps) and convert the plume/ocean volume
C     |   exchange into addMass, ramped up over DP_rampTime.
C     | o Fill the plume diagnostics.
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "DOUBLEPLUME.h"

C     !INPUT/OUTPUT PARAMETERS:
      _RL futureTime
      INTEGER futureIter
      INTEGER myThid

#ifdef ALLOW_DOUBLEPLUME
C     !LOCAL VARIABLES:
      INTEGER bi, bj, I, J, K
      LOGICAL doSolve
      _RL rampFac
CEOP

C--   Spin-up ramp of the plume input (0 -> 1 over DP_rampTime)
      IF ( DP_rampTime .GT. 0. _d 0 ) THEN
        rampFac = MIN( 1. _d 0, MAX( 0. _d 0,
     &            DBLE(futureIter - nIter0)*deltaTClock/DP_rampTime ) )
      ELSE
        rampFac = 1. _d 0
      ENDIF
      doSolve = dp_firstCall .OR. ( MOD(futureIter,DP_nUpdate).EQ.0 )

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO J = 1, sNy
         DO I = 1, sNx
          IF ( plumeMask(I,J,bi,bj) .NE. 0. _d 0 ) THEN
           IF ( doSolve )
     &       CALL DOUBLEPLUME_COLUMN( I, J, bi, bj, myThid )
           DO K = 1, Nr
            addMass(I,J,K,bi,bj) = rampFac*rhoConst
     &                           * volFluxDiff3D(I,J,K,bi,bj)
           ENDDO
          ENDIF
         ENDDO
        ENDDO

#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
         CALL DIAGNOSTICS_FILL(QiProf3D,'icefrtQi',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(MiProf3D,'icefrtMi',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(HiProf3D,'icefrtHi',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(BiProf3D,'icefrtBi',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(HgProf3D,'icefrtHg',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(MbProf3D,'icefrtMb',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(CiProf3D,'icefrtCi',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(QoProf3D,'icefrtQo',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(MoProf3D,'icefrtMo',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(HoProf3D,'icefrtHo',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(BoProf3D,'icefrtBo',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(CoProf3D,'icefrtCo',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(mProf3D, 'icefrtMe',0,Nr,1,bi,bj,myThid)
         CALL DIAGNOSTICS_FILL(volFluxDiff3D,'icefrtVX',
     &                         0,Nr,1,bi,bj,myThid)
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */

       ENDDO
      ENDDO

      _BARRIER
      _BEGIN_MASTER(myThid)
      dp_firstCall = .FALSE.
      _END_MASTER(myThid)

#endif /* ALLOW_DOUBLEPLUME */
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_COLUMN
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_COLUMN( I, J, bi, bj, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | Double plume for the column (I,J,bi,bj):
C     | 1. plume type and subglacial discharge
C     | 2. ambient T/S profiles (in-situ T from theta)
C     | 3. source depth (one cell above the deepest cell)
C     | 4. solve the double plume (DOUBLEPLUME_SOLVE)
C     | 5. exchange with the ocean in each model level:
C     |    the net upward flux of the plume system (inner + outer)
C     |    is Qt = Qi + Qo, with Qt = Q_sg below the source and
C     |    Qt = 0 above the surface. In level k (faces h_bot,h_top)
C     |      V_k = -(Qt(h_top) - Qt(h_bot)) + meltwater_k
C     |      H_k = -(Ht(h_top) - Ht(h_bot)) + (ice + gas heat)_k
C     |      S_k = -(Bt(h_top) - Bt(h_bot))
C     |    V_k > 0: plume water with T = H_k/(rho_r cp V_k) - T0,
C     |    S = S_k/V_k is added; V_k < 0: ambient water is removed.
C     |    Sum of V_k = Q_sg + total meltwater (mass conserved);
C     |    heat and salt are conserved where V_k > 0.
C     | 6. layer-mean diagnostics
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "DOUBLEPLUME.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_START.h"
#include "PTRACERS_FIELDS.h"
#endif

C     !INPUT PARAMETERS:
      INTEGER I, J, bi, bj, myThid

#ifdef ALLOW_DOUBLEPLUME
C     !LOCAL VARIABLES:
      INTEGER K, kk, kL, nPass, nF, typ
#ifdef ALLOW_PTRACERS
      INTEGER iTr
#endif
      LOGICAL runPlume
      DOUBLE PRECISION errOut, hMid, depth, dz
      DOUBLE PRECISION Vk(Nr), Hk(Nr), Sk(Nr), mSum(Nr), mLen(Nr)
      DOUBLE PRECISION Ft(3), Fb(3), Y(7), YD(7), xv(20)
      DOUBLE PRECISION mv(DP_maxFine), qh(DP_maxFine), md(DP_maxFine)
      DOUBLE PRECISION Tdet, Sdet, Tlo, Thi, Slo, Shi, Tp, Vmax
      DOUBLE PRECISION sumV, sumM
      _RL tLoc, sLoc, pLoc, zero
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL SW_TEMP, SW_PTMP
      EXTERNAL SW_TEMP, SW_PTMP
CEOP
      zero = 0. _d 0
      dz   = DP_dz

C--   1. Plume type (sign = ice-front orientation):
C     1 ice, no plume; 2 sheet; 3 half-cone; 4 sheet+cone;
C     5 detaching cone
      typ = NINT( ABS(plumeMask(I,J,bi,bj)) )
      useSheetPlume  = ( typ .EQ. 2 .OR. typ .EQ. 4 )
      useConePlume   = ( typ .EQ. 3 .OR. typ .EQ. 4 )
      useDetachPlume = ( typ .EQ. 5 )
      IF ( typ .LT. 1 .OR. typ .GT. 5 ) THEN
        WRITE(msgBuf,'(A)')
     &   'DOUBLEPLUME_COLUMN: plume mask value must be -5..5'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DOUBLEPLUME_COLUMN'
      ENDIF

      delta_x = dxG(I,J,bi,bj)
      delta_y = dyG(I,J,bi,bj)
      w_sg = runoffVel(I,J,bi,bj)
      r_sg = runoffRad(I,J,bi,bj)
      IF ( useSheetPlume ) THEN
        Q_sg = w_sg*r_sg*delta_y
      ELSEIF ( useConePlume .OR. useDetachPlume ) THEN
        Q_sg = w_sg*0.5 _d 0*DP_pi*r_sg**2
      ELSE
        Q_sg = 0. _d 0
      ENDIF

C--   2. Ambient profiles
      DO K = 1, Nr
        prProf(K)  = ABS(rC(K))*rhoConst*gravity*1. _d -4
        saProf(K)  = salt(I,J,K,bi,bj)
        ptProf(K)  = theta(I,J,K,bi,bj)
        taProf(K)  = SW_TEMP( saProf(K), ptProf(K), prProf(K), zero )
        uProf(K)   = ABS(uVel(I,J,K,bi,bj))
        vProf(K)   = ABS(vVel(I,J,K,bi,bj))
        delta_z(K) = drF(K)
      ENDDO
      maxDepth = rF(Nr+1)

C--   3. Plume source: the face above the deepest face of the column
      iceDepth  = 0. _d 0
      iceDepthK = 0
      DO K = 2, Nr+1
        IF ( iceDepthK .EQ. 0 .AND. rF(K) .LE. R_low(I,J,bi,bj) ) THEN
          iceDepth  = rF(K-1)
          iceDepthK = K-1
        ENDIF
      ENDDO
      IF ( iceDepthK .EQ. 0 .OR. iceDepth .EQ. 0. _d 0 ) THEN
        WRITE(msgBuf,'(A,2I6,A)')
     &   'DOUBLEPLUME_COLUMN: plume in cell', I, J,
     &   ' but the water column is too shallow'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DOUBLEPLUME_COLUMN'
      ENDIF
      dp_H        = ABS(iceDepth)
      dp_nAmb     = MAX( 1, iceDepthK-1 )
      plumeLayers = NINT( dp_H/dz ) + 1
      IF ( plumeLayers .GT. DP_maxFine ) THEN
        WRITE(msgBuf,'(A,I8,A)')
     &   'DOUBLEPLUME_COLUMN: need DP_maxFine >=', plumeLayers,
     &   ' (DOUBLEPLUME.h)'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DOUBLEPLUME_COLUMN'
      ENDIF
      nF = plumeLayers

C--   4. Solve the double plume
      runPlume = ( typ .NE. 1 ) .AND.
     &           ( Q_sg .GT. 0. _d 0 .OR. DP_mb0 .GT. 0. _d 0 )
      IF ( runPlume ) THEN
        CALL DOUBLEPLUME_SOLVE( nPass, errOut, myThid )
      ELSE
        DO kk = 1, nF
          QiProf_f(kk) = 0.D0
          MiProf_f(kk) = 0.D0
          HiProf_f(kk) = 0.D0
          HgProf_f(kk) = 0.D0
          BiProf_f(kk) = 0.D0
          MbProf_f(kk) = 0.D0
          CiProf_f(kk) = 0.D0
          QoProf_f(kk) = 0.D0
          MoProf_f(kk) = 0.D0
          HoProf_f(kk) = 0.D0
          BoProf_f(kk) = 0.D0
          CoProf_f(kk) = 0.D0
        ENDDO
        Q_sg = 0. _d 0
      ENDIF

C--   5. Exchange with the ocean
C     meltwater volume, ice + gas heat and melt rate on the fine grid
      DO kk = 1, nF
        mv(kk) = 0.D0
        qh(kk) = 0.D0
        md(kk) = 0.D0
        IF ( QiProf_f(kk) .GT. 0.D0 .AND. MiProf_f(kk) .GT. 0.D0 ) THEN
          Y(1) = QiProf_f(kk)
          Y(2) = MiProf_f(kk)
          Y(3) = HiProf_f(kk)
          Y(4) = HgProf_f(kk)
          Y(5) = BiProf_f(kk)
          Y(6) = MbProf_f(kk)
          Y(7) = CiProf_f(kk)
          CALL DOUBLEPLUME_INNER_RHS( (kk-1)*dz, Y, YD, xv )
          mv(kk) = xv(13)
          qh(kk) = xv(14) + xv(15)
          md(kk) = xv(10)
        ENDIF
      ENDDO

      DO K = 1, Nr
        Vk(K)   = 0.D0
        Hk(K)   = 0.D0
        Sk(K)   = 0.D0
        mSum(K) = 0.D0
        mLen(K) = 0.D0
      ENDDO
      IF ( runPlume ) THEN
C       flux divergence between the faces of each wet level above
C       the source (level K spans heights dp_H-|rF(K+1)| .. dp_H-|rF(K)|)
        DO K = 1, iceDepthK-1
          IF ( K .EQ. 1 ) THEN
C           surface: nothing leaves through the top
            Ft(1) = 0.D0
            Ft(2) = 0.D0
            Ft(3) = 0.D0
          ELSE
            CALL DOUBLEPLUME_TOTFLUX( dp_H - ABS(rF(K)), Ft )
          ENDIF
          IF ( K .EQ. iceDepthK-1 ) THEN
C           source face: subglacial discharge enters from below
            Fb(1) = Q_sg
            Fb(2) = DP_rhoR*DP_cp*(dp_T0src + DP_T0)*Q_sg
            Fb(3) = DP_S0*Q_sg
          ELSE
            CALL DOUBLEPLUME_TOTFLUX( dp_H - ABS(rF(K+1)), Fb )
          ENDIF
          Vk(K) = -( Ft(1) - Fb(1) )
          Hk(K) = -( Ft(2) - Fb(2) )
          Sk(K) = -( Ft(3) - Fb(3) )
        ENDDO
C       meltwater volume and ice/gas heat, integrated per level
        DO kk = 1, nF-1
          hMid  = (kk - 0.5D0)*dz
          depth = dp_H - hMid
          kL = 0
          DO K = 1, iceDepthK-1
            IF ( depth .GE. ABS(rF(K)) .AND. depth .LT. ABS(rF(K+1)) )
     &           kL = K
          ENDDO
          IF ( kL .GT. 0 ) THEN
            Vk(kL) = Vk(kL) + 0.5D0*(mv(kk) + mv(kk+1))*dz
            Hk(kL) = Hk(kL) + 0.5D0*(qh(kk) + qh(kk+1))*dz
            IF ( md(kk) .GT. 0.D0 .AND. md(kk+1) .GT. 0.D0 ) THEN
              mSum(kL) = mSum(kL) + 0.5D0*(md(kk) + md(kk+1))*dz
              mLen(kL) = mLen(kL) + dz
            ENDIF
          ENDIF
        ENDDO
      ENDIF

C     range of plausible temperatures / salinities for the added water
      Tlo =  1.D30
      Thi = -1.D30
      Slo =  1.D30
      Shi = -1.D30
      DO K = 1, dp_nAmb
        Tlo = MIN( Tlo, taProf(K) )
        Thi = MAX( Thi, taProf(K) )
        Slo = MIN( Slo, saProf(K) )
        Shi = MAX( Shi, saProf(K) )
      ENDDO
      IF ( runPlume ) THEN
        Tlo = MIN( Tlo, dp_T0src )
        Slo = MIN( Slo, DP_S0 )
        DO kk = 1, nF
          IF ( QiProf_f(kk) .GT. 0.D0 ) THEN
            Tp  = HiProf_f(kk)/(QiProf_f(kk)*DP_rhoR*DP_cp) - DP_T0
            Tlo = MIN( Tlo, Tp )
            Thi = MAX( Thi, Tp )
            Slo = MIN( Slo, BiProf_f(kk)/QiProf_f(kk) )
            Shi = MAX( Shi, BiProf_f(kk)/QiProf_f(kk) )
          ENDIF
        ENDDO
      ENDIF

      Vmax = 0.D0
      DO K = 1, Nr
        Vmax = MAX( Vmax, ABS(Vk(K)) )
      ENDDO
      sumV = 0.D0
      sumM = 0.D0
      DO K = 1, Nr
        volFluxDiff3D(I,J,K,bi,bj) = Vk(K)
        IF ( Vk(K) .GT. 1.D-9*Vmax .AND. Vk(K) .GT. 0.D0 ) THEN
          Tdet = Hk(K)/(DP_rhoR*DP_cp*Vk(K)) - DP_T0
          Sdet = Sk(K)/Vk(K)
          Tdet = MIN( MAX( Tdet, Tlo ), Thi )
          Sdet = MIN( MAX( Sdet, Slo ), Shi )
          tLoc = Tdet
          sLoc = Sdet
          pLoc = prProf(K)
          temp_addMass3D(I,J,K,bi,bj) =
     &         SW_PTMP( sLoc, tLoc, pLoc, zero )
          salt_addMass3D(I,J,K,bi,bj) = Sdet
        ELSE
C         water removed (or nothing exchanged): no T/S change
          temp_addMass3D(I,J,K,bi,bj) = theta(I,J,K,bi,bj)
          salt_addMass3D(I,J,K,bi,bj) = salt(I,J,K,bi,bj)
        ENDIF
#ifdef ALLOW_PTRACERS
C       the plume does not carry passive tracers: no tracer change
        DO iTr = 1, PTRACERS_num
          ptr_addMass3D(I,J,K,bi,bj,iTr) = pTracer(I,J,K,bi,bj,iTr)
        ENDDO
#endif
        sumV = sumV + Vk(K)
      ENDDO
      DO kk = 1, nF-1
        sumM = sumM + 0.5D0*(mv(kk) + mv(kk+1))*dz
      ENDDO

      IF ( runPlume .AND. debugLevel .GE. debLevB ) THEN
        WRITE(msgBuf,'(A,2I5,A,I3,A,ES10.3,A,I4)')
     &   'DOUBLEPLUME_COLUMN: cell', I, J, ' passes', nPass,
     &   ' err', errOut, ' peels', j_peel
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        WRITE(msgBuf,'(A,ES12.4,A,ES12.4)')
     &   'DOUBLEPLUME_COLUMN: net volume added', sumV,
     &   '  = Q_sg + meltwater', Q_sg + sumM
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

C--   6. Layer-mean diagnostics
      CALL DOUBLEPLUME_LAYERMEAN( QiProf_f, QiProf )
      CALL DOUBLEPLUME_LAYERMEAN( MiProf_f, MiProf )
      CALL DOUBLEPLUME_LAYERMEAN( HiProf_f, HiProf )
      CALL DOUBLEPLUME_LAYERMEAN( HgProf_f, HgProf )
      CALL DOUBLEPLUME_LAYERMEAN( BiProf_f, BiProf )
      CALL DOUBLEPLUME_LAYERMEAN( MbProf_f, MbProf )
      CALL DOUBLEPLUME_LAYERMEAN( CiProf_f, CiProf )
      CALL DOUBLEPLUME_LAYERMEAN( QoProf_f, QoProf )
      CALL DOUBLEPLUME_LAYERMEAN( MoProf_f, MoProf )
      CALL DOUBLEPLUME_LAYERMEAN( HoProf_f, HoProf )
      CALL DOUBLEPLUME_LAYERMEAN( BoProf_f, BoProf )
      CALL DOUBLEPLUME_LAYERMEAN( CoProf_f, CoProf )
      DO K = 1, Nr
        IF ( mLen(K) .GT. 0.D0 ) THEN
C         plume melt rate [m/day], averaged over the plume contact
          mProf(K) = mSum(K)/mLen(K)*86400.D0
        ELSE
          mProf(K) = 0.D0
        ENDIF
        QiProf3D(I,J,K,bi,bj) = QiProf(K)
        MiProf3D(I,J,K,bi,bj) = MiProf(K)
        HiProf3D(I,J,K,bi,bj) = HiProf(K)
        HgProf3D(I,J,K,bi,bj) = HgProf(K)
        BiProf3D(I,J,K,bi,bj) = BiProf(K)
        MbProf3D(I,J,K,bi,bj) = MbProf(K)
        CiProf3D(I,J,K,bi,bj) = CiProf(K)
        QoProf3D(I,J,K,bi,bj) = QoProf(K)
        MoProf3D(I,J,K,bi,bj) = MoProf(K)
        HoProf3D(I,J,K,bi,bj) = HoProf(K)
        BoProf3D(I,J,K,bi,bj) = BoProf(K)
        CoProf3D(I,J,K,bi,bj) = CoProf(K)
        mProf3D (I,J,K,bi,bj) = mProf(K)
      ENDDO

#endif /* ALLOW_DOUBLEPLUME */
      RETURN
      END

C=======================================================================
      SUBROUTINE DOUBLEPLUME_TOTFLUX( h, F )
C     Net upward flux of the plume system (inner + outer) at height h,
C     linearly interpolated on the fine grid:
C     F(1) volume [m^3/s], F(2) heat [J/s], F(3) salt [psu m^3/s]
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION h, F(3)
      INTEGER k1
      DOUBLE PRECISION r, w
      r  = h/DP_dz + 1.D0
      k1 = MAX( 1, MIN( INT(r), plumeLayers-1 ) )
      w  = MAX( 0.D0, MIN( 1.D0, r - k1 ) )
      F(1) = (1.D0-w)*(QiProf_f(k1)   + QoProf_f(k1))
     &     +       w*(QiProf_f(k1+1) + QoProf_f(k1+1))
      F(2) = (1.D0-w)*(HiProf_f(k1)   + HoProf_f(k1))
     &     +       w*(HiProf_f(k1+1) + HoProf_f(k1+1))
      F(3) = (1.D0-w)*(BiProf_f(k1)   + BoProf_f(k1))
     &     +       w*(BiProf_f(k1+1) + BoProf_f(k1+1))
      RETURN
      END

C=======================================================================
      SUBROUTINE DOUBLEPLUME_LAYERMEAN( fine, coarse )
C     Mean of a fine-grid profile over each model level
C     (level k spans depths |rF(k)| <= depth < |rF(k+1)|)
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DOUBLEPLUME.h"
      _RL fine(DP_maxFine), coarse(Nr+1)
      INTEGER k, kk, cnt(Nr)
      DOUBLE PRECISION depth
      DO k = 1, Nr+1
        coarse(k) = 0.D0
      ENDDO
      DO k = 1, Nr
        cnt(k) = 0
      ENDDO
      DO kk = 1, plumeLayers
        depth = dp_H - (kk-1)*DP_dz
        DO k = 1, Nr
          IF ( depth .GE. ABS(rF(k)) .AND. depth .LT. ABS(rF(k+1)) )
     &    THEN
            coarse(k) = coarse(k) + fine(kk)
            cnt(k)    = cnt(k) + 1
          ENDIF
        ENDDO
      ENDDO
      DO k = 1, Nr
        IF ( cnt(k) .GT. 0 ) coarse(k) = coarse(k)/cnt(k)
      ENDDO
      RETURN
      END
