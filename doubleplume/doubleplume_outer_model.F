#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_MODEL
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_MODEL( myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_OUTER_MODEL
C     | o Solves the downward-flowing outer plume of the double
C     |   plume model, given the inner-plume profiles left on the
C     |   fine grid by DOUBLEPLUME_INNER_MODEL.
C     | o Port of doublePlume.py: outer_calc, outer_surf,
C     |   outer_cpic, outer_dpic, outer_fr, outer_derivs,
C     |   outer_stop. Constants: DOUBLEPLUME.h (DOUBLEPLUME_PHYS);
C     |   shared physics: doubleplume_physics.F.
C     *==========================================================*
C
C     Vertical coordinate
C       h = height above the subglacial inlet:
C           h = 0 at the inlet, h = H = |iceDepth| at the surface.
C       Fine-grid index k (1..plumeLayers) sits at h = (k-1)*DP_dz,
C       i.e. depth = H - (k-1)*DP_dz, as in the inner model.
C
C     Outer-plume state (same sign convention as Python):
C       Y(1) = Qo  volume flux   [m^3/s]   (< 0: flows downward)
C       Y(2) = Jo  momentum flux [m^4/s^2] (> 0)
C       Y(3) = Ho  heat flux     [J/s] = rho_r*cp*(To+T0)*Qo
C                  (same convention as the inner plume's Y(3))
C       Y(4) = Bo  salt flux     [psu m^3/s]
C       Y(5) = Co  gas flux      [kg/s]    = co*Qo
C
C     Outputs: QoProf_f, MoProf_f, HoProf_f, BoProf_f, CoProf_f on
C     the fine grid (layer means are formed in DOUBLEPLUME_CALC).
C     Physics routines (DP_RHO, DP_PEEL, ...) are in
C     doubleplume_physics.F.
C     \ev

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DOUBLEPLUME.h"

C     !INPUT PARAMETERS:
      INTEGER myThid

C     !LOCAL VARIABLES:
      INTEGER k, j, jj, jp, nLoop, flag
      DOUBLE PRECISION Y(5), hPeel(DP_maxPeel), htmp
      DOUBLE PRECISION hTop, hElev, hNew, h0, hf, hEnd, depth
      DOUBLE PRECISION Qi, ui, bi, Ti, Si, ci, rho_i, rho_a, Ta, Sa
      DOUBLE PRECISION us, Ep
      DOUBLE PRECISION Qo, To, So, co, uo, rho_o
      LOGICAL hasInner
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      DOUBLE PRECISION DP_RHO, DP_GAS_RHO, DOUBLEPLUME_OUTER_FRSOLVE
      EXTERNAL DP_RHO, DP_GAS_RHO, DOUBLEPLUME_OUTER_FRSOLVE
CEOP

      dpo_nFine = plumeLayers
      jp        = 0

      DO k = 1, dpo_nFine
        QoProf_f(k) = 0.D0
        MoProf_f(k) = 0.D0
        HoProf_f(k) = 0.D0
        BoProf_f(k) = 0.D0
        CoProf_f(k) = 0.D0
      ENDDO

C--   Top of the inner plume = highest fine-grid point where the inner
C     plume exists. The outer plume starts there (Python: zi[-1]),
C     not at a fixed index below the surface.
      dpo_nTop = 0
      DO k = 1, dpo_nFine
        IF ( QiProf_f(k).GT.0.D0 .AND. MiProf_f(k).GT.0.D0 )
     &       dpo_nTop = k
      ENDDO
      IF ( dpo_nTop .LT. 2 ) THEN
        WRITE(msgBuf,'(A)')
     &   'DOUBLEPLUME_OUTER_MODEL: no inner plume, outer plume = 0'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        GOTO 900
      ENDIF
      hTop = (dpo_nTop-1)*DP_dz

C--   Discrete peel heights. The inner model stores peel DEPTHS in
C     z_peel; convert to heights and sort ascending, so hPeel(jp) is
C     the highest peel (Python: zi[i_peel[jp-1]]). A local copy is
C     used so the shared j_peel/z_peel are left untouched.
      jp = MIN( j_peel, DP_maxPeel )
      DO j = 1, jp
        hPeel(j) = dp_H - z_peel(j)
      ENDDO
      DO j = 2, jp
        htmp = hPeel(j)
        jj   = j - 1
   10   CONTINUE
        IF ( jj .GE. 1 ) THEN
          IF ( hPeel(jj) .GT. htmp ) THEN
            hPeel(jj+1) = hPeel(jj)
            jj = jj - 1
            GOTO 10
          ENDIF
        ENDIF
        hPeel(jj+1) = htmp
      ENDDO

C=== 1. Initial condition at the top of the inner plume (outer_surf)
      CALL DOUBLEPLUME_OUTER_INNERVARS( hTop, Qi, ui, bi, Ti, Si,
     &     ci, rho_i, rho_a, Ta, Sa, us, Ep, hasInner )
      depth = dp_H - hTop

      IF ( depth .GT. DP_surfTol + 0.5D0*DP_dz ) THEN
C       Inner plume stopped below the surface: the outer plume is
C       made of inner-plume fluid
        Qo = Qi
        IF ( DP_epsilon .GT. 0.D0 ) Qo = DP_fStall*Qo
        To = Ti
        So = Si
        co = ci
      ELSE
C       Inner plume reached the surface: mix with a fraction fe of
C       ambient fluid
        Qo = (1.D0 + DP_fe)*Qi
        To = (Ti + DP_fe*Ta   )/(1.D0 + DP_fe)
        So = (Si + DP_fe*Sa   )/(1.D0 + DP_fe)
        co = (ci + DP_fe*DP_Ca)/(1.D0 + DP_fe)
      ENDIF
      rho_o = DP_RHO( To, So, depth )
      uo    = DOUBLEPLUME_OUTER_FRSOLVE( Qo, bi, rho_a, rho_o )

      Y(1) = -Qo
      Y(2) = -Y(1)*uo
      Y(3) = DP_rhoR*DP_cp*(To + DP_T0)*Y(1)
      Y(4) = So*Y(1)
      Y(5) = co*Y(1)

      IF ( jp .GT. 0 ) THEN
        hf = hPeel(jp)
      ELSE
        hf = 0.D0
      ENDIF
      CALL DOUBLEPLUME_OUTER_SEGMENT( hTop, hf, Y, hEnd, myThid )
      hElev = hEnd

C=== 2. Walk down to the inlet (while-loop of Python outer_calc)
      nLoop = 0
  100 CONTINUE
      IF ( hElev .GT. DP_bottomTol .AND. nLoop .LT. 100000 ) THEN
        nLoop = nLoop + 1
        flag  = 0
        h0    = hElev
        hf    = 0.D0

        IF ( jp .GT. 0 ) THEN
          IF ( DP_epsilon .EQ. 0.D0 .OR.
     &         ABS(hElev - hPeel(jp)) .LE. DP_peelTol ) THEN
C           At a discrete peel (outer_dpic)
            CALL DOUBLEPLUME_OUTER_DPIC( hPeel(jp), Y )
            h0 = hPeel(jp)
            IF ( jp .GT. 1 ) hf = hPeel(jp-1)
            flag = 1
            jp   = jp - 1
          ELSE
C           Continuous peeling down to the next discrete peel
            CALL DOUBLEPLUME_OUTER_CPIC( hElev, hPeel(jp),
     &                                   h0, Y, flag )
            hf = hPeel(jp)
          ENDIF
        ELSE
          IF ( DP_epsilon .EQ. 0.D0 ) THEN
C           No continuous peeling and no discrete peels left: no more
C           outer plume below (outer_dpic with jp = 0 closes out)
            hElev = 0.D0
            GOTO 100
          ENDIF
          CALL DOUBLEPLUME_OUTER_CPIC( hElev, 0.D0, h0, Y, flag )
        ENDIF

        IF ( flag .EQ. 1 ) THEN
          CALL DOUBLEPLUME_OUTER_SEGMENT( h0, hf, Y, hEnd, myThid )
          hNew = hEnd
        ELSE
C         Not viable: leave zeros and continue below the window
          hNew = h0
        ENDIF

C       Always make downward progress (e.g. zero-width window where
C       the inner plume is absent)
        IF ( hNew .GE. hElev - 1.D-6 ) hNew = hElev - DP_dz
        hElev = hNew
        GOTO 100
      ENDIF
      IF ( nLoop .GE. 100000 ) THEN
        WRITE(msgBuf,'(A)')
     &   'DOUBLEPLUME_OUTER_MODEL: descent loop hit its iteration cap'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

  900 CONTINUE
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_SEGMENT
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_SEGMENT( h0, hf, Y, hEnd, myThid )
C     !DESCRIPTION: \bv
C     Integrates one outer-plume segment downward from h0 to hf with
C     DLSODE, stores the state on the fine grid, and stops early on
C     the Python outer_stop criteria (momentum lost, NaN, bottom).
C     On return hEnd is the height reached and Y the state there.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
C     !INPUT/OUTPUT PARAMETERS:
      DOUBLE PRECISION h0, hf, Y(5), hEnd
      INTEGER myThid
C     !LOCAL VARIABLES:
      EXTERNAL DOUBLEPLUME_OUTER_DERIVS, DP_JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
      INTEGER IWORK(40), i, kg
      DOUBLE PRECISION RWORK(200), T, TOUT, RTOLo, ATOLo, hg
      LOGICAL stopNow, badY
CEOP

      NEQ    = 5
      ITOL   = 1
      RTOLo  = DP_rtolOuter
      ATOLo  = DP_atolOuter
      ITASK  = 1
      ISTATE = 1
      IOPT   = 1
      LRW    = 200
      LIW    = 40
C     MF = 10: non-stiff Adams method, no Jacobian (Python uses the
C     explicit RK45). The old MF = 21 needs a user Jacobian, which the
C     dummy Jacobian never supplied.
      MF     = 10
      DO i = 5, 10
        RWORK(i) = 0.D0
        IWORK(i) = 0
      ENDDO
C     max step (Python max_step), max steps per call, max warnings
      RWORK(6) = DP_hmaxOuter
      IWORK(6) = 20000
      IWORK(7) = 1

      T    = h0
      hEnd = h0
      CALL DOUBLEPLUME_OUTER_STORE( T, Y )
      IF ( hf .GE. h0 ) RETURN
C     Python's outer_stop is also evaluated at the starting point: a
C     plume that starts with Jo below the threshold stops at once.
      IF ( Y(2) .LT. DP_JoMin ) RETURN

C     first fine-grid height strictly below h0
      kg = INT( h0/DP_dz + 1.D-9 )
      hg = kg*DP_dz
      IF ( hg .GE. h0 - 1.D-9 ) hg = hg - DP_dz

      stopNow = .FALSE.
  200 CONTINUE
      IF ( .NOT. stopNow ) THEN
        TOUT = MAX( hg, hf )
        CALL DLSODE( DOUBLEPLUME_OUTER_DERIVS, NEQ, Y, T, TOUT,
     &               ITOL, RTOLo, ATOLo, ITASK, ISTATE, IOPT,
     &               RWORK, LRW, IWORK, LIW, DP_JAC, MF )
        badY = .FALSE.
        DO i = 1, 5
          IF ( Y(i) .NE. Y(i) ) badY = .TRUE.
        ENDDO
        IF ( ISTATE .LT. 0 .OR. badY ) THEN
          stopNow = .TRUE.
        ELSEIF ( Y(2) .LT. DP_JoMin ) THEN
C         outer_stop: momentum lost (round(Jo,3) <= 0). Python's event
C         locates Jo = 0 exactly; here the crossing lies inside the
C         last step, so that point (Jo <= 0) is not stored.
          stopNow = .TRUE.
        ELSE
          CALL DOUBLEPLUME_OUTER_STORE( T, Y )
          IF ( T .LE. 0.D0 )            stopNow = .TRUE.
          IF ( TOUT .LE. hf + 1.D-9 )   stopNow = .TRUE.
          hg = hg - DP_dz
        ENDIF
        hEnd = T
        GOTO 200
      ENDIF

      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_STORE
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_STORE( h, Y )
C     !DESCRIPTION: \bv
C     Stores the outer state at height h if h is a fine-grid point.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION h, Y(5)
      INTEGER k
      DOUBLE PRECISION r
CEOP
      r = h/DP_dz
      k = NINT(r)
      IF ( ABS(r - k) .LT. 1.D-6 ) THEN
        k = k + 1
        IF ( k .GE. 1 .AND. k .LE. dpo_nFine ) THEN
          QoProf_f(k) = Y(1)
          MoProf_f(k) = Y(2)
          HoProf_f(k) = Y(3)
          BoProf_f(k) = Y(4)
          CoProf_f(k) = Y(5)
        ENDIF
      ENDIF
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_DERIVS
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_DERIVS( NEQ, T, Y, YDOT )
C     !DESCRIPTION: \bv
C     Outer plume equations (Python: outer_derivs). T = height h.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
C     !INPUT/OUTPUT PARAMETERS:
      INTEGER NEQ
      DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ)
C     !LOCAL VARIABLES:
      INTEGER i
      DOUBLE PRECISION h, depth, Qi, ui, bi, Ti, Si, ci
      DOUBLE PRECISION rho_i, rho_a, Ta, Sa, us, Ep
      DOUBLE PRECISION uo, bo, To, So, co, rho_o, a1, a2, ent
      LOGICAL hasInner
      DOUBLE PRECISION DP_RHO
      EXTERNAL DP_RHO
CEOP

      h     = MIN( MAX( T, 0.D0 ), dp_H )
      depth = dp_H - h
      CALL DOUBLEPLUME_OUTER_INNERVARS( h, Qi, ui, bi, Ti, Si, ci,
     &     rho_i, rho_a, Ta, Sa, us, Ep, hasInner )

C     Python: if ui == 0 then alpha_1 = alpha_2 = Ep = 0
      IF ( hasInner ) THEN
        a1 = DP_alpha1
        a2 = DP_alpha2
      ELSE
        a1 = 0.D0
        a2 = 0.D0
        Ep = 0.D0
      ENDIF

C     Outer primitive variables (Python getVars, outer part)
      IF ( ABS(Y(1)) .GT. 1.D-12 ) THEN
        uo = Y(2)/Y(1)
        To = Y(3)/(Y(1)*DP_rhoR*DP_cp) - DP_T0
        So = Y(4)/Y(1)
        co = Y(5)/Y(1)
        IF ( uo .NE. 0.D0 ) THEN
          bo = SQRT( ABS( 2.D0*Y(1)/(uo*DP_pi) ) + bi**2 )
        ELSE
          bo = bi
        ENDIF
      ELSE
        uo = 0.D0
        bo = 0.D0
        To = Ta
        So = Sa
        co = DP_Ca
      ENDIF
      rho_o = DP_RHO( To, So, depth )

C     common inner-plume entrainment velocity
      ent = a1*(ui + DP_c1*uo)

C     Volume
      YDOT(1) = DP_pi*bi*( ent + a2*uo )
     &        + DP_pi*bo*DP_alpha3*uo + Ep
C     Momentum
      YDOT(2) = ( -0.5D0*DP_pi*DP_g*(bo**2 - bi**2)/DP_rhoR
     &            *(rho_a - rho_o)
     &          + DP_pi*bi*( ent*uo + a2*uo*ui )
     &          + Ep*ui ) / DP_gammaO
C     Heat (rho_r*cp*(T+T0) convention, as in the inner plume)
      YDOT(3) = DP_rhoR*DP_cp*(
     &            DP_pi*bi*( ent*(To + DP_T0) + a2*uo*(Ti + DP_T0) )
     &          + DP_pi*bo*DP_alpha3*uo*(Ta + DP_T0)
     &          + Ep*(Ti + DP_T0) )
C     Salt
      YDOT(4) = DP_pi*bi*( ent*So + a2*uo*Si )
     &        + DP_pi*bo*DP_alpha3*uo*Sa + Ep*Si
C     Dissolved gas
      YDOT(5) = DP_pi*bi*( ent*co + a2*uo*ci )
     &        + DP_pi*bo*DP_alpha3*uo*DP_Ca + Ep*ci

C     The outer plume is integrated downward: flip the sign (Python)
      DO i = 1, NEQ
        YDOT(i) = -YDOT(i)
      ENDDO

      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_INNERVARS
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_INNERVARS( h, Qi, ui, bi, Ti, Si,
     &     ci, rho_i, rho_a, Ta, Sa, us, Ep, hasInner )
C     !DESCRIPTION: \bv
C     Inner-plume primitive variables, ambient state, bubble slip
C     velocity and peeling flux at height h (the part of Python
C     getVars that the outer plume needs). Inner profiles are
C     linearly interpolated (Python: interp1d) instead of taking
C     the nearest index.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
C     !INPUT/OUTPUT PARAMETERS:
      DOUBLE PRECISION h, Qi, ui, bi, Ti, Si, ci, rho_i, rho_a
      DOUBLE PRECISION Ta, Sa, us, Ep
      LOGICAL hasInner
C     !LOCAL VARIABLES:
      INTEGER k1
      DOUBLE PRECISION r, w, depth, Ji, Hi, Hb, Bsalt, mb, Cg
      DOUBLE PRECISION Tg, P, rho_b, d
      DOUBLE PRECISION DP_RHO, DP_GAS_RHO, DP_SLIP
      DOUBLE PRECISION DP_PEEL
      EXTERNAL DP_RHO, DP_GAS_RHO, DP_SLIP, DP_PEEL
CEOP

      depth = dp_H - h
      CALL DP_AMBIENT( depth, Ta, Sa )
      rho_a = DP_RHO( Ta, Sa, depth )

C     Fine-grid index k sits at h = (k-1)*DP_dz
      r  = h/DP_dz + 1.D0
      k1 = INT(r)
      k1 = MAX( 1, MIN( k1, dpo_nTop-1 ) )
      w  = MAX( 0.D0, MIN( 1.D0, r - k1 ) )
      Qi = (1.D0-w)*QiProf_f(k1) + w*QiProf_f(k1+1)
      Ji = (1.D0-w)*MiProf_f(k1) + w*MiProf_f(k1+1)
      Hi = (1.D0-w)*HiProf_f(k1) + w*HiProf_f(k1+1)
      Hb = (1.D0-w)*HgProf_f(k1) + w*HgProf_f(k1+1)
      Bsalt = (1.D0-w)*BiProf_f(k1) + w*BiProf_f(k1+1)
      mb = (1.D0-w)*MbProf_f(k1) + w*MbProf_f(k1+1)
      Cg = (1.D0-w)*CiProf_f(k1) + w*CiProf_f(k1+1)

      IF ( Qi .LE. 1.D-12 .OR. Ji .LE. 1.D-12 ) THEN
C       No inner plume at this height
        hasInner = .FALSE.
        Qi    = 0.D0
        ui    = 0.D0
        bi    = 0.D0
        Ti    = Ta
        Si    = Sa
        ci    = DP_Ca
        rho_i = rho_a
        us    = 0.D0
        Ep    = 0.D0
        RETURN
      ENDIF

      hasInner = .TRUE.
      ui    = Ji/Qi
      bi    = SQRT( 2.D0*Qi/(DP_pi*ui) )
      Ti    = Hi/(Qi*DP_rhoR*DP_cp) - DP_T0
      Si    = Bsalt/Qi
      ci    = Cg/Qi
      rho_i = DP_RHO( Ti, Si, depth )

C     Bubble slip velocity (Python: gas_density, bubble_diameter,
C     u_slip), or a fixed value if DP_usConst > 0
      IF ( DP_usConst .GT. 0.D0 ) THEN
        us = DP_usConst
      ELSE
        IF ( mb .GT. 0.D0 ) THEN
          Tg = Hb/(mb*DP_cpg) - DP_T0
        ELSE
          Tg = Ti
        ENDIF
        P     = DP_Patm + DP_rhoR*DP_g*depth
        rho_b = DP_GAS_RHO( Tg, P )
        IF ( mb .GT. 0.D0 .AND. dp_nb0 .GT. 0.D0 ) THEN
          d = ( 6.D0*mb/(DP_pi*rho_b*dp_nb0) )**(1.D0/3.D0)
        ELSE
          d = DP_de0/1000.D0
        ENDIF
        us = DP_SLIP( rho_a, rho_b, Ti, d )
      ENDIF

      Ep = DP_PEEL( ui, us, bi, rho_i, rho_a )

      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_DPIC
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_DPIC( hP, Y )
C     !DESCRIPTION: \bv
C     Outer plume initial condition at a discrete peel at height hP
C     (Python: outer_dpic, jp > 0): a fraction fp of the inner
C     volume flux detrains, carrying the inner T, S and gas.
C     The inner model stores the already-reduced (post-peel) state
C     at the peel point, so the pre-peel state is read one fine
C     step below.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION hP, Y(5)
      INTEGER i
      DOUBLE PRECISION hS, depth, Qi, ui, bi, Ti, Si, ci
      DOUBLE PRECISION rho_i, rho_a, Ta, Sa, us, Ep, Qo, uo, rho_o
      LOGICAL hasInner
      DOUBLE PRECISION DP_RHO, DOUBLEPLUME_OUTER_FRSOLVE
      EXTERNAL DP_RHO, DOUBLEPLUME_OUTER_FRSOLVE
CEOP
      hS = MAX( hP - DP_dz, 0.D0 )
      CALL DOUBLEPLUME_OUTER_INNERVARS( hS, Qi, ui, bi, Ti, Si, ci,
     &     rho_i, rho_a, Ta, Sa, us, Ep, hasInner )
      IF ( .NOT. hasInner ) THEN
        DO i = 1, 5
          Y(i) = 0.D0
        ENDDO
        RETURN
      ENDIF

      depth = dp_H - hP
      CALL DP_AMBIENT( depth, Ta, Sa )
      rho_a = DP_RHO( Ta, Sa, depth )
      Qo    = DP_fp*Qi
      rho_o = DP_RHO( Ti, Si, depth )
      uo    = DOUBLEPLUME_OUTER_FRSOLVE( Qo, bi, rho_a, rho_o )

      Y(1) = -Qo
      Y(2) = -Y(1)*uo
      Y(3) = DP_rhoR*DP_cp*(Ti + DP_T0)*Y(1)
      Y(4) = Si*Y(1)
      Y(5) = ci*Y(1)
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_CPIC
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_CPIC( hElev, hFloor, h0, Y, flag )
C     !DESCRIPTION: \bv
C     Outer plume initial condition from continuous peeling
C     (Python: outer_cpic). The peeling flux Ep of the inner plume
C     is integrated (trapezoid rule, like np.trapezoid) over a
C     window of iter inner-plume widths below hElev, for
C     iter = 1..DP_nwidths, never going below hFloor.
C     flag = 1: viable outer plume starting at h0 with state Y.
C     flag = 0: no viable outer plume; continue from h0.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION hElev, hFloor, h0, Y(5)
      INTEGER flag
      INTEGER i, iter, k
      LOGICAL done, hasInner
      DOUBLE PRECISION Qi, ui, bi, Ti, Si, ci, rho_i, rho_a, Ta, Sa
      DOUBLE PRECISION us, Ep, biTop, hUp, hLo, hk, hPrev
      DOUBLE PRECISION fE, fT, fS, fC, sE, sT, sS, sC
      DOUBLE PRECISION ui0, bi0, rhoa0, Ep0
      DOUBLE PRECISION Qo, To, So, co, uo, bo, rho_o, uoS, dQdz
      DOUBLE PRECISION DP_RHO, DOUBLEPLUME_OUTER_FRSOLVE
      EXTERNAL DP_RHO, DOUBLEPLUME_OUTER_FRSOLVE
CEOP
      flag = 0
      h0   = hElev
      DO i = 1, 5
        Y(i) = 0.D0
      ENDDO
      done = .FALSE.
      iter = 1

  300 CONTINUE
      IF ( .NOT. done .AND. iter .LE. DP_nwidths ) THEN

C       Inner plume width at the current elevation
        CALL DOUBLEPLUME_OUTER_INNERVARS( hElev, Qi, ui, biTop, Ti,
     &       Si, ci, rho_i, rho_a, Ta, Sa, us, Ep, hasInner )
        hUp = hElev
        hLo = MAX( hElev - iter*biTop, hFloor )

C       Trapezoid integrals over [hLo, hUp] through hLo, the inner
C       fine-grid points strictly inside, and hUp
        CALL DOUBLEPLUME_OUTER_INNERVARS( hLo, Qi, ui, bi, Ti, Si,
     &       ci, rho_i, rho_a, Ta, Sa, us, Ep, hasInner )
        ui0   = ui
        bi0   = bi
        rhoa0 = rho_a
        Ep0   = Ep
        fE = Ep
        fT = Ep*Ti
        fS = Ep*Si
        fC = Ep*ci
        sE = 0.D0
        sT = 0.D0
        sS = 0.D0
        sC = 0.D0
        hPrev = hLo
        k = INT( hLo/DP_dz + 1.D-9 ) + 2
  310   CONTINUE
        hk = (k-1)*DP_dz
        IF ( hk .GE. hUp - 1.D-9 ) hk = hUp
        IF ( hk .GT. hPrev ) THEN
          CALL DOUBLEPLUME_OUTER_INNERVARS( hk, Qi, ui, bi, Ti, Si,
     &         ci, rho_i, rho_a, Ta, Sa, us, Ep, hasInner )
          sE = sE + 0.5D0*(fE + Ep   )*(hk - hPrev)
          sT = sT + 0.5D0*(fT + Ep*Ti)*(hk - hPrev)
          sS = sS + 0.5D0*(fS + Ep*Si)*(hk - hPrev)
          sC = sC + 0.5D0*(fC + Ep*ci)*(hk - hPrev)
          fE = Ep
          fT = Ep*Ti
          fS = Ep*Si
          fC = Ep*ci
          hPrev = hk
        ENDIF
        IF ( hk .LT. hUp ) THEN
          k = k + 1
          GOTO 310
        ENDIF

        Qo = -sE
        IF ( Qo .LE. 0.D0 ) THEN
C         No net peeling in this window
          h0   = hLo
          flag = 0
          DO i = 1, 5
            Y(i) = 0.D0
          ENDDO
          IF ( hLo .LE. hFloor ) done = .TRUE.
          iter = iter + 1
          GOTO 300
        ENDIF

        To = -sT/Qo
        So = -sS/Qo
        co = -sC/Qo

C       Outer plume starts at the bottom of the window
        rho_o = DP_RHO( To, So, dp_H - hLo )
        uo    = DOUBLEPLUME_OUTER_FRSOLVE( Qo, bi0, rhoa0, rho_o )
        bo    = SQRT( 2.D0*Qo/(DP_pi*uo) + bi0**2 )

        Y(1) = -Qo
        Y(2) = -Y(1)*uo
        Y(3) = DP_rhoR*DP_cp*(To + DP_T0)*Y(1)
        Y(4) = So*Y(1)
        Y(5) = co*Y(1)
        h0   = hLo

C       Viability: the outer plume must gain volume flux as it
C       descends (conditions at the bottom of the window)
        uoS  = -uo
        dQdz = DP_pi*bi0*( DP_alpha1*(ui0 + DP_c1*uoS)
     &                   + DP_alpha2*uoS )
     &       + DP_pi*bo*DP_alpha3*uoS + Ep0
        IF ( dQdz .GT. 0.D0 ) THEN
          flag = 0
          DO i = 1, 5
            Y(i) = 0.D0
          ENDDO
        ELSE
          flag = 1
          done = .TRUE.
          IF ( h0 .LE. 0.D0 ) THEN
            flag = 0
            h0   = 0.D0
            DO i = 1, 5
              Y(i) = 0.D0
            ENDDO
          ENDIF
        ENDIF
        iter = iter + 1
        GOTO 300
      ENDIF

      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_OUTER_FRSOLVE
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DOUBLEPLUME_OUTER_FRSOLVE( Qo, bi,
     &                                              rho_a, rho_o )
C     !DESCRIPTION: \bv
C     Outer plume initial velocity (> 0) from the Froude condition
C     (Python: outer_fr / calculate_outer_fr):
C       uo - Fr*sqrt(|(bo - bi)*g*(rho_a - rho_o)/rho_o|) = 0,
C       bo = sqrt(2*Qo/(pi*uo) + bi^2),  Qo > 0 (magnitude).
C     The residual increases monotonically with uo, so the root is
C     unique and bisection finds the same root as fsolve.
C     \ev
C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION Qo, bi, rho_a, rho_o
      INTEGER it
      DOUBLE PRECISION uLo, uHi, uMid, fMid
      DOUBLE PRECISION DOUBLEPLUME_OUTER_FRRES
      EXTERNAL DOUBLEPLUME_OUTER_FRRES
CEOP
      uLo = 1.D-8
      uHi = 1.D0
      IF ( Qo .LE. 0.D0 ) THEN
        DOUBLEPLUME_OUTER_FRSOLVE = uLo
        RETURN
      ENDIF
      IF ( DOUBLEPLUME_OUTER_FRRES(uLo,Qo,bi,rho_a,rho_o).GE.0.D0 ) THEN
        DOUBLEPLUME_OUTER_FRSOLVE = uLo
        RETURN
      ENDIF
      DO it = 1, 20
        IF ( DOUBLEPLUME_OUTER_FRRES(uHi,Qo,bi,rho_a,rho_o).GT.0.D0 )
     &       GOTO 20
        uHi = 10.D0*uHi
      ENDDO
   20 CONTINUE
      DO it = 1, 200
        uMid = 0.5D0*(uLo + uHi)
        fMid = DOUBLEPLUME_OUTER_FRRES( uMid, Qo, bi, rho_a, rho_o )
        IF ( fMid .GT. 0.D0 ) THEN
          uHi = uMid
        ELSE
          uLo = uMid
        ENDIF
        IF ( uHi - uLo .LT. 1.D-12*uHi ) GOTO 30
      ENDDO
   30 CONTINUE
      DOUBLEPLUME_OUTER_FRSOLVE = 0.5D0*(uLo + uHi)
      RETURN
      END

C=======================================================================
      DOUBLE PRECISION FUNCTION DOUBLEPLUME_OUTER_FRRES( uo, Qo, bi,
     &                                             rho_a, rho_o )
C     Froude-condition residual (Python: outer_fr)
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION uo, Qo, bi, rho_a, rho_o, bo
      bo = SQRT( 2.D0*Qo/(DP_pi*uo) + bi**2 )
      DOUBLEPLUME_OUTER_FRRES = uo - DP_FrOuter
     &     * SQRT( ABS( (bo - bi)*DP_g*(rho_a - rho_o)/rho_o ) )
      RETURN
      END

