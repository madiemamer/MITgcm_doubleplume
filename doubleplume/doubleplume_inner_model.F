#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_INNER_MODEL
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_INNER_MODEL( myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_INNER_MODEL
C     | o Solves the rising inner (bubble) plume, given the outer
C     |   plume left on the fine grid by the previous pass of
C     |   DOUBLEPLUME_OUTER_MODEL (zero on the first pass).
C     | o Port of doublePlume.py: __init__ source conditions,
C     |   main_ic, inner_calc, inner_ic, inner_stop, inner_peel,
C     |   inner_stall, inner_derivs, getVars.
C     *==========================================================*
C
C     Coordinate: h = height above the source (DLSODE's T), from 0
C     to hMax = (plumeLayers-1)*DP_dz; depth = dp_H - h.
C     State (Python y):
C       Y(1) Q  volume flux [m^3/s]    Y(2) J  momentum [m^4/s^2]
C       Y(3) H  heat = rho_r*cp*(T+T0)*Q [J/s]
C       Y(4) Hb bubble heat = mb*cp_g*(Tg+T0) [J/s]
C       Y(5) B  salt [psu m^3/s]       Y(6) mb bubble mass [kg/s]
C       Y(7) C  dissolved gas [kg/s]
C     Events, checked after every DP_dz step and located by linear
C     interpolation (Python locates them exactly with solve_ivp):
C       peel  : plume liquid denser than ambient and u_i < u_s/kappa
C       stall : u_i drops below DP_uStall
C       bubbles no longer buoyant (rho_b > rho_i)
C         -> all three start a new plume segment (Python inner_ic)
C       bubbles dissolved (mb/mb0 < fdis), invalid state
C         -> the inner plume ends
C       surface reached -> the inner plume ends
C     Output: QiProf_f..CiProf_f on the fine grid (0 above the top),
C     j_peel and z_peel (peel depths), dp_nb0, dp_mbSrc, dp_T0src.
C     \ev

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"

C     !INPUT PARAMETERS:
      INTEGER myThid

C     !LOCAL VARIABLES:
      INTEGER NXV
      PARAMETER ( NXV = 20 )
      EXTERNAL DOUBLEPLUME_INNER_DERIVS, DP_JAC
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
      INTEGER IWORK(40), i, k, kLast, nSeg, evType, nStep
      DOUBLE PRECISION RWORK(200), RTOLi, ATOLi
      DOUBLE PRECISION Y(7), Yp(7), Ys(7), Y0(7), YD(7), xv(NXV)
      DOUBLE PRECISION evP(6), evN(6)
      DOUBLE PRECISION T, TOUT, hMax, hPrev, hEv, hc, z0, hg, w, f
      DOUBLE PRECISION depth, Ta, Sa, P, rho_a, rho_b, us, R, u
      DOUBLE PRECISION Qg, Q0, Xi0, Qn, un, bn
      LOGICAL done, dissolved
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      DOUBLE PRECISION DP_RHO, DP_GAS_RHO, DP_SLIP, DP_WUEST
      DOUBLE PRECISION DP_VOIDFRAC
      EXTERNAL DP_RHO, DP_GAS_RHO, DP_SLIP, DP_WUEST, DP_VOIDFRAC
CEOP

      DO k = 1, plumeLayers
        QiProf_f(k) = 0.D0
        MiProf_f(k) = 0.D0
        HiProf_f(k) = 0.D0
        HgProf_f(k) = 0.D0
        BiProf_f(k) = 0.D0
        MbProf_f(k) = 0.D0
        CiProf_f(k) = 0.D0
      ENDDO
      j_peel = 0
      hMax   = (plumeLayers-1)*DP_dz

C=== Source conditions (Python __init__ and main_ic)
      depth = dp_H
      CALL DP_AMBIENT( depth, Ta, Sa )
      rho_a    = DP_RHO( Ta, Sa, depth )
      P        = DP_Patm + DP_rhoR*DP_g*depth
      dp_T0src = DP_lam1*DP_S0 + DP_lam2 + DP_lam3*depth
      dp_mbSrc   = DP_mb0
      rho_b    = DP_GAS_RHO( dp_T0src, P )
      Qg       = dp_mbSrc/rho_b
      dp_nb0   = 6.D0*dp_mbSrc/(rho_b*DP_pi*(DP_de0/1000.D0)**3)
      IF ( DP_usConst .GT. 0.D0 ) THEN
        us = DP_usConst
      ELSE
        us = DP_SLIP( rho_a, rho_b, dp_T0src, DP_de0/1000.D0 )
      ENDIF
      IF ( r_sg .GT. 0.D0 ) THEN
        R = r_sg
      ELSE
        R = DP_R
      ENDIF

      IF ( Q_sg .LE. 0.D0 ) THEN
C       Bubbles only: Wuest et al. (1992) initial velocity
        u = DP_WUEST( us, dp_mbSrc, R, rho_a, rho_b )
      ELSE
C       Discharge (+ bubbles): Lima Neto (2012)
        u = Q_sg/( (1.D0 - Qg/(Qg + Q_sg))*0.5D0*DP_pi*R**2 )
      ENDIF
      IF ( u .LE. 0.D0 ) THEN
        WRITE(msgBuf,'(A)')
     &   'DOUBLEPLUME_INNER_MODEL: no discharge and no bubbles'
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
        RETURN
      ENDIF
      Q0  = 0.5D0*DP_pi*R**2*u
      Xi0 = DP_VOIDFRAC( dp_mbSrc, rho_b, R, us, u )
      IF ( Xi0 .GT. 0.5D0 ) THEN
        WRITE(msgBuf,'(A,F8.3)')
     &   'DOUBLEPLUME_INNER_MODEL: source void fraction > 0.5: ', Xi0
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

      Y(1) = Q0
      Y(2) = Q0*u
      Y(3) = DP_rhoR*DP_cp*(dp_T0src + DP_T0)*Q0
      Y(4) = dp_mbSrc*DP_cpg*(dp_T0src + DP_T0)
      Y(5) = DP_S0*Q0
      Y(6) = dp_mbSrc
      Y(7) = DP_C0*Q0

C=== DLSODE settings. MF = 10: non-stiff Adams (Python: RK45)
      NEQ   = 7
      ITOL  = 1
      RTOLi = DP_rtolInner
      ATOLi = DP_atolInner
      ITASK = 1
      IOPT  = 1
      LRW   = 200
      LIW   = 40
      MF    = 10
      DO i = 5, 10
        RWORK(i) = 0.D0
        IWORK(i) = 0
      ENDDO
C     max step (Python max_step), max steps per call, max warnings
      RWORK(6) = DP_hmaxInner
      IWORK(6) = 20000
      IWORK(7) = 1

C=== Integrate upward (Python inner_calc)
      T      = 0.D0
      kLast  = 1
      CALL DOUBLEPLUME_INNER_STORE( 1, Y )
      CALL DOUBLEPLUME_INNER_EVT( T, Y, evP )
      ISTATE    = 1
      hg        = DP_dz
      done      = .FALSE.
      dissolved = .FALSE.
      nSeg      = 1
      nStep     = 0

  100 CONTINUE
      IF ( .NOT. done ) THEN
        IF ( hg .GT. hMax + 1.D-9 ) THEN
C         reached the surface
          done = .TRUE.
          GOTO 100
        ENDIF
        nStep = nStep + 1
        TOUT  = hg
        hPrev = T
        DO i = 1, 7
          Yp(i) = Y(i)
        ENDDO
        CALL DLSODE( DOUBLEPLUME_INNER_DERIVS, NEQ, Y, T, TOUT,
     &               ITOL, RTOLi, ATOLi, ITASK, ISTATE, IOPT,
     &               RWORK, LRW, IWORK, LIW, DP_JAC, MF )
        IF ( ISTATE .LT. 0 ) THEN
          WRITE(msgBuf,'(A,I4,A,F10.3)')
     &     'DOUBLEPLUME_INNER_MODEL: DLSODE ISTATE=', ISTATE,
     &     ' at h =', hPrev
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          done = .TRUE.
          GOTO 100
        ENDIF
        CALL DOUBLEPLUME_INNER_EVT( TOUT, Y, evN )

C--     Which event (if any) happened first in (hPrev, TOUT] ?
C       evType: 0 none, 1 peel, 2 stall, 3 bubbles not buoyant,
C               4 dissolved, 5 invalid state
        evType = 0
        f      = 2.D0
        IF ( evN(1) .LT. 0.5D0 ) THEN
          evType = 5
          f      = 1.D0
        ELSE
C         peel: dense liquid and u_i < u_s/kappa (Python inner_peel)
          IF ( evP(2) .GE. 0.D0 .AND. evN(2) .LT. 0.D0 ) THEN
            IF ( evP(3) .GE. 0.D0 .AND. evN(3) .LT. 0.D0 ) THEN
C             became dense in this step: locate the density crossing
              hc = evP(3)/(evP(3) - evN(3))
            ELSE
              hc = evP(2)/(evP(2) - evN(2))
            ENDIF
            IF ( hc .LT. f ) THEN
              f      = hc
              evType = 1
            ENDIF
          ENDIF
C         stall (Python inner_stall)
          IF ( evP(4) .GE. 0.D0 .AND. evN(4) .LT. 0.D0 ) THEN
            hc = evP(4)/(evP(4) - evN(4))
            IF ( hc .LT. f ) THEN
              f      = hc
              evType = 2
            ENDIF
          ENDIF
C         bubbles no longer buoyant (Python inner_stop event 2)
          IF ( evP(5) .GT. 0.D0 .AND. evN(5) .LE. 0.D0 ) THEN
            hc = evP(5)/(evP(5) - evN(5))
            IF ( hc .LT. f ) THEN
              f      = hc
              evType = 3
            ENDIF
          ENDIF
C         bubbles dissolved (Python inner_stop event 1)
          IF ( evP(6) .GE. 0.D0 .AND. evN(6) .LT. 0.D0 ) THEN
            hc = evP(6)/(evP(6) - evN(6))
            IF ( hc .LT. f ) THEN
              f      = hc
              evType = 4
            ENDIF
          ENDIF
        ENDIF

        IF ( evType .EQ. 0 ) THEN
C--       Normal step: store and continue
          kLast = NINT(TOUT/DP_dz) + 1
          CALL DOUBLEPLUME_INNER_STORE( kLast, Y )
          DO i = 1, 6
            evP(i) = evN(i)
          ENDDO
          hg = hg + DP_dz
          GOTO 100
        ENDIF

        IF ( evType .GE. 4 ) THEN
C--       Dissolved or invalid: the plume ends at the last stored
C         point (within one DP_dz of the event)
          IF ( evType .EQ. 4 ) dissolved = .TRUE.
          IF ( evType .EQ. 5 ) THEN
            WRITE(msgBuf,'(A,F10.3)')
     &       'DOUBLEPLUME_INNER_MODEL: invalid state near h =', TOUT
            CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                          SQUEEZE_RIGHT, myThid )
          ENDIF
          done = .TRUE.
          GOTO 100
        ENDIF

C--     Peel / stall / non-buoyant bubbles at hEv: re-integrate from
C       the last stored point to hEv to get the state there
        hEv = hPrev + MAX( 0.D0, MIN( 1.D0, f ) )*(TOUT - hPrev)
        DO i = 1, 7
          Ys(i) = Yp(i)
        ENDDO
        IF ( hEv - hPrev .GT. 1.D-9 ) THEN
          T      = hPrev
          ISTATE = 1
          CALL DLSODE( DOUBLEPLUME_INNER_DERIVS, NEQ, Ys, T, hEv,
     &                 ITOL, RTOLi, ATOLi, ITASK, ISTATE, IOPT,
     &                 RWORK, LRW, IWORK, LIW, DP_JAC, MF )
          IF ( ISTATE .LT. 0 ) THEN
            done = .TRUE.
            GOTO 100
          ENDIF
        ENDIF

        IF ( j_peel .GE. DP_maxPeel ) THEN
          WRITE(msgBuf,'(A)')
     &     'DOUBLEPLUME_INNER_MODEL: DP_maxPeel reached, plume ends'
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                        SQUEEZE_RIGHT, myThid )
          done = .TRUE.
          GOTO 100
        ENDIF
        j_peel = j_peel + 1
        z_peel(j_peel) = dp_H - hEv

C--     New plume segment (Python inner_ic)
        CALL DOUBLEPLUME_INNER_RHS( hEv, Ys, YD, xv )
        Qn = (1.D0 - DP_fp)*Ys(1)
C       Wuest velocity with the pre-peel width xv(3)
        un = DP_WUEST( xv(9), Ys(6), xv(3), xv(7), xv(8) )
        IF ( un .LE. 0.D0 ) un = xv(2)
        bn = SQRT( 2.D0*Qn/(DP_pi*un) )
        z0 = hEv + 0.1D0*bn
        Y0(1) = Qn
        Y0(2) = Qn*un
        Y0(3) = DP_rhoR*DP_cp*(xv(4) + DP_T0)*Qn
        IF ( Ys(6) .GT. 0.D0 ) THEN
          Y0(4) = Ys(6)*DP_cpg*(xv(19) + DP_T0)
        ELSE
          Y0(4) = Ys(4)
        ENDIF
        Y0(5) = xv(5)*Qn
        Y0(6) = Ys(6)
        Y0(7) = Ys(7)/Ys(1)*Qn

        IF ( z0 .GE. hMax ) THEN
          done = .TRUE.
          GOTO 100
        ENDIF

C--     Fill the fine grid across the gap (hPrev, z0] by linear
C       interpolation between the pre-peel state (hEv) and the new
C       segment's initial state (z0), as Python's interp1d does
        k = kLast + 1
  200   CONTINUE
        hc = (k-1)*DP_dz
        IF ( hc .LE. z0 + 1.D-9 .AND. k .LE. plumeLayers ) THEN
          IF ( hc .LE. hEv ) THEN
            CALL DOUBLEPLUME_INNER_STORE( k, Ys )
          ELSE
            w = (hc - hEv)/(z0 - hEv)
            DO i = 1, 7
              YD(i) = (1.D0 - w)*Ys(i) + w*Y0(i)
            ENDDO
            CALL DOUBLEPLUME_INNER_STORE( k, YD )
          ENDIF
          kLast = k
          k = k + 1
          GOTO 200
        ENDIF

C--     Restart the integration at z0
        T = z0
        DO i = 1, 7
          Y(i) = Y0(i)
        ENDDO
        ISTATE = 1
        hg     = kLast*DP_dz
        nSeg   = nSeg + 1
        CALL DOUBLEPLUME_INNER_EVT( T, Y, evP )
        IF ( nSeg .GT. 10*DP_maxPeel ) done = .TRUE.
        GOTO 100
      ENDIF

      IF ( debugLevel .GE. debLevB ) THEN
        WRITE(msgBuf,'(A,F9.2,A,I4,A,L2)')
     &   'DOUBLEPLUME_INNER_MODEL: top h =', (kLast-1)*DP_dz,
     &   '  peels =', j_peel, '  dissolved =', dissolved
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                      SQUEEZE_RIGHT, myThid )
      ENDIF

      RETURN
      END

C=======================================================================
      SUBROUTINE DOUBLEPLUME_INNER_STORE( k, Y )
C     Store the inner state Y at fine-grid index k
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      INTEGER k
      DOUBLE PRECISION Y(7)
      IF ( k .GE. 1 .AND. k .LE. plumeLayers ) THEN
        QiProf_f(k) = Y(1)
        MiProf_f(k) = Y(2)
        HiProf_f(k) = Y(3)
        HgProf_f(k) = Y(4)
        BiProf_f(k) = Y(5)
        MbProf_f(k) = Y(6)
        CiProf_f(k) = Y(7)
      ENDIF
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_INNER_EVT
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_INNER_EVT( h, Y, ev )
C     !DESCRIPTION: \bv
C     Event functions of the inner plume at (h, Y):
C       ev(1) = 1 valid state, 0 invalid (Q, J <= 0 or NaN)
C       ev(2) = Python inner_peel: 1 if rho_i <= rho_a,
C               else u_i - u_s/kappa
C       ev(3) = rho_a - rho_i
C       ev(4) = u_i - DP_uStall            (Python inner_stall)
C       ev(5) = rho_i - rho_b              (bubbles still buoyant)
C       ev(6) = mb/mb0 - fdis              (bubbles not dissolved)
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION h, Y(7), ev(6)
      INTEGER i
      DOUBLE PRECISION YD(7), xv(20)
CEOP
      DO i = 2, 6
        ev(i) = 1.D0
      ENDDO
      ev(1) = 1.D0
      DO i = 1, 7
        IF ( Y(i) .NE. Y(i) ) ev(1) = 0.D0
      ENDDO
      IF ( ev(1) .LT. 0.5D0 ) RETURN
      CALL DOUBLEPLUME_INNER_RHS( h, Y, YD, xv )
      ev(1) = xv(1)
      IF ( ev(1) .LT. 0.5D0 ) RETURN
      IF ( xv(6) .LE. xv(7) ) THEN
        ev(2) = 1.D0
      ELSE
        ev(2) = xv(2) - xv(9)/DP_kappaPeel
      ENDIF
      ev(3) = xv(7) - xv(6)
      ev(4) = xv(2) - DP_uStall
      ev(5) = xv(6) - xv(8)
      IF ( dp_mbSrc .GT. 0.D0 ) THEN
        ev(6) = Y(6)/dp_mbSrc - DP_fdis
      ELSE
        ev(6) = 1.D0
      ENDIF
      RETURN
      END

C=======================================================================
C     DLSODE right-hand side wrapper
      SUBROUTINE DOUBLEPLUME_INNER_DERIVS( NEQ, T, Y, YDOT )
      IMPLICIT NONE
      INTEGER NEQ
      DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ)
      DOUBLE PRECISION xv(20)
      CALL DOUBLEPLUME_INNER_RHS( T, Y, YDOT, xv )
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DOUBLEPLUME_INNER_RHS
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_INNER_RHS( h, Y, YDOT, xv )
C     !DESCRIPTION: \bv
C     Inner plume equations (Python: inner_derivs + getVars).
C     Also returns diagnostics in xv:
C       1 valid(1/0)  2 ui   3 bi   4 Ti   5 Si   6 rho_i  7 rho_a
C       8 rho_b  9 us  10 mdot [m/s]  11 Tb  12 Sb
C      13 meltwater volume 2*bi*mdot [m^2/s]
C      14 heat from the ice interface to the plume [W/m]
C      15 heat from the bubbles to the water [W/m]
C      16 void fraction  17 Ta  18 Sa  19 gas temperature
C      20 bubble diameter [m]
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION h, Y(7), YDOT(7), xv(20)
      INTEGER i, k1
      DOUBLE PRECISION depth, Ta, Sa, P, rho_a, rho_i, rho_b
      DOUBLE PRECISION ui, bi, Ti, Si, mb, ci, mdot, Sb, Tb
      DOUBLE PRECISION r, w, Qout, Jout, Hout, Bout, Cout
      DOUBLE PRECISION uo, bwo, To, So, co
      DOUBLE PRECISION Tg, d, nb, us, Ep, betaM, betaT, Cs, Xi, Fpd
      DOUBLE PRECISION ent, bub, heat, sqCd
      DOUBLE PRECISION DP_RHO, DP_GAS_RHO, DP_SLIP, DP_PEEL
      DOUBLE PRECISION DP_SHERWOOD, DP_GAS_DIFF, DP_GAS_SOL
      DOUBLE PRECISION DP_VOIDFRAC
      EXTERNAL DP_RHO, DP_GAS_RHO, DP_SLIP, DP_PEEL
      EXTERNAL DP_SHERWOOD, DP_GAS_DIFF, DP_GAS_SOL, DP_VOIDFRAC
CEOP
      DO i = 1, 20
        xv(i) = 0.D0
      ENDDO
      DO i = 1, 7
        YDOT(i) = 0.D0
      ENDDO

      depth = MAX( dp_H - h, 0.D0 )
      CALL DP_AMBIENT( depth, Ta, Sa )
      P     = DP_Patm + DP_rhoR*DP_g*depth
      rho_a = DP_RHO( Ta, Sa, depth )
      xv(7)  = rho_a
      xv(17) = Ta
      xv(18) = Sa
      IF ( .NOT.( Y(1) .GT. 0.D0 .AND. Y(2) .GT. 0.D0 ) ) RETURN
      xv(1) = 1.D0

C--   Inner primitive variables
      ui = Y(2)/Y(1)
      bi = SQRT( 2.D0*Y(1)/(ui*DP_pi) )
      Ti = Y(3)/(Y(1)*DP_rhoR*DP_cp) - DP_T0
      Si = Y(5)/Y(1)
      mb = Y(6)
      ci = Y(7)/Y(1)

C--   Melting (Holland & Jenkins 1999)
      CALL DP_MELT( Y(1), Y(2), Y(3), Y(5), depth, mdot, Sb, Tb )

C--   Outer plume at h: linear interpolation on the fine grid
      r  = h/DP_dz + 1.D0
      k1 = MAX( 1, MIN( INT(r), plumeLayers-1 ) )
      w  = MAX( 0.D0, MIN( 1.D0, r - k1 ) )
      Qout = (1.D0-w)*QoProf_f(k1) + w*QoProf_f(k1+1)
      Jout = (1.D0-w)*MoProf_f(k1) + w*MoProf_f(k1+1)
      Hout = (1.D0-w)*HoProf_f(k1) + w*HoProf_f(k1+1)
      Bout = (1.D0-w)*BoProf_f(k1) + w*BoProf_f(k1+1)
      Cout = (1.D0-w)*CoProf_f(k1) + w*CoProf_f(k1+1)
      IF ( ABS(Qout) .LT. 1.D-12 ) THEN
        uo  = 0.D0
        bwo = 0.D0
        To  = Ta
        So  = Sa
        co  = DP_Ca
      ELSE
        uo  = Jout/Qout
        To  = Hout/(Qout*DP_rhoR*DP_cp) - DP_T0
        So  = Bout/Qout
        co  = Cout/Qout
        IF ( uo .NE. 0.D0 ) THEN
          bwo = SQRT( ABS( 2.D0*Qout/(uo*DP_pi) ) + bi**2 )
        ELSE
          bwo = bi
        ENDIF
      ENDIF

C--   Densities, bubbles
      IF ( mb .GT. 0.D0 ) THEN
        Tg = Y(4)/(mb*DP_cpg) - DP_T0
        nb = dp_nb0
      ELSE
        Tg = Ti
        nb = 0.D0
      ENDIF
      rho_i = DP_RHO( Ti, Si, depth )
      rho_b = DP_GAS_RHO( Tg, P )
      IF ( mb .GT. 0.D0 .AND. dp_nb0 .GT. 0.D0 ) THEN
        d = ( 6.D0*mb/(DP_pi*rho_b*dp_nb0) )**(1.D0/3.D0)
      ELSE
        d = DP_de0/1000.D0
      ENDIF
      IF ( DP_usConst .GT. 0.D0 ) THEN
        us = DP_usConst
      ELSE
        us = DP_SLIP( rho_a, rho_b, Ti, d )
      ENDIF
      Ep    = DP_PEEL( ui, us, bi, rho_i, rho_a )
      betaM = DP_SHERWOOD( Ti, rho_a, us, d, DP_GAS_DIFF(Ti) )
      betaT = DP_SHERWOOD( Ti, rho_a, us, d, 1.46D-7 )
      Cs    = DP_GAS_SOL( Ti, Si, P )
      Xi    = DP_VOIDFRAC( mb, rho_b, bi, us, ui )
      Fpd   = DP_lambdaR**2*Xi*(rho_a - rho_b)
      sqCd  = SQRT(DP_Cd)

C--   Equations (Python inner_derivs)
      ent = DP_alpha1*(ui + DP_c1*uo)
C     volume
      YDOT(1) = DP_pi*bi*( ent + DP_alpha2*uo ) + Ep
     &        + 2.D0*bi*mdot
C     momentum
      YDOT(2) = ( 0.5D0*DP_pi*DP_g*bi**2/DP_rhoR
     &            *( Fpd + DP_lambdaS**2*(1.D0 - Xi)*(rho_a - rho_i) )
     &          + DP_pi*bi*( ent*uo + DP_alpha2*uo*ui )
     &          + Ep*ui - DP_Cd*2.D0*bi*ui**2 ) / DP_gammaI
C     heat (water)
      heat = DP_rhoR*DP_cp*(
     &          DP_pi*bi*( ent*(To + DP_T0)
     &                   + DP_alpha2*uo*(Ti + DP_T0) )
     &        + Ep*(Ti + DP_T0)
     &        - 2.D0*DP_gammaT*sqCd*bi*ui*(Ti - Tb) )
     &     + DP_rhoR*DP_cp*2.D0*bi*mdot*(Tb + DP_T0)
C     bubble mass (dissolution)
      bub = DP_pi*nb*d**2/(ui + us)
      YDOT(6) = -( bub*betaM*(Cs - ci) )
C     bubble heat, incl. enthalpy carried off by the dissolved gas
      YDOT(4) = -bub*rho_b*DP_cpg*betaT*(Tg - Ti)
     &        + DP_cpg*(Tg + DP_T0)*YDOT(6)
C     the water receives both
      YDOT(3) = heat - YDOT(4)
C     salt
      YDOT(5) = DP_pi*bi*( ent*So + DP_alpha2*uo*Si ) + Ep*Si
     &        + 2.D0*bi*mdot*Sb
     &        - 2.D0*DP_gammaS*sqCd*bi*ui*(Si - Sb)
C     dissolved gas
      YDOT(7) = DP_pi*bi*( ent*co + DP_alpha2*uo*ci ) + Ep*ci
     &        - YDOT(6)

      xv(2)  = ui
      xv(3)  = bi
      xv(4)  = Ti
      xv(5)  = Si
      xv(6)  = rho_i
      xv(8)  = rho_b
      xv(9)  = us
      xv(10) = mdot
      xv(11) = Tb
      xv(12) = Sb
      xv(13) = 2.D0*bi*mdot
      xv(14) = DP_rhoR*DP_cp*( 2.D0*bi*mdot*(Tb + DP_T0)
     &       - 2.D0*DP_gammaT*sqCd*bi*ui*(Ti - Tb) )
      xv(15) = -YDOT(4)
      xv(16) = Xi
      xv(19) = Tg
      xv(20) = d
      RETURN
      END
