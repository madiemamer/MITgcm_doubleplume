#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_PLUME_MODEL
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_OUTER_MODEL (mythid)
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_PLUME_MODEL:
C     | o Uses ODEPACK to solve one of 3 models describing subglacial
C     |   runoff driven plumes
C     | 
C     *===========================================================
c	\ev
C     !USES:
      IMPLICIT NONE
C     == Global variables ==

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"
#include "DOUBLEPLUME_OPTIONS.h"

C     !INPUT/OUTPUT PARAMETERS:

C     ==Local variables==

C     Plume models
      EXTERNAL DOUBLEPLUME_OUTER_DERIVS, JEX_OP

C     For ODEPACK
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(200), LIW, LRW,
     &         MF, NEQ, JT, NG, JROOT
      _RL T, TOUT
      _RL  RWORK(200), Y(5)

C For Brentq root solver.
      INTEGER MAXITER, IFLAG
      DOUBLE PRECISION ATOL_BQ, RTOL_BQ, FTOL_BQ
      DOUBLE PRECISION AX, BX, FAX, FBX, XZERO, FZERO
      PARAMETER (MAXITER=100)

C     Local Variables
      INTEGER  k, mythid, i, j
      INTEGER Nz, z0, zf, temp, qi_idx_upper
      INTEGER end_i, start_i
      DOUBLE PRECISION dz
      DOUBLE PRECISION tambient, sambient, cambient, depth, P
      DOUBLE PRECISION rho_i, rho_a, rho_o, rho_t, rho_b
      DOUBLE PRECISION Qi, Mi, Hi, Bi, Ci
      DOUBLE PRECISION Qi_loc(plumeLayers), Mi_loc(plumeLayers) 
      DOUBLE PRECISION Hi_loc(plumeLayers), Bi_loc(plumeLayers)
      DOUBLE PRECISION Ci_loc(plumeLayers) 
      DOUBLE PRECISION Qo, Mo, Ho, Bo, Co
      DOUBLE PRECISION ui, binner, Ti, Si
      DOUBLE PRECISION uo, bouter, To, So
      DOUBLE PRECISION us, Qg, R, u0
      DOUBLE PRECISION diff_l, smallDiff_l
      DOUBLE PRECISION sum_vol, sum_mom, sum_heat, sum_gasHeat
      DOUBLE PRECISION sum_salt, sum_gasFlux, sum_gasCon
      LOGICAL DISSOLVED, IMAG, done
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C     Temporary holders for continuous peel model    
      INTEGER i_lower, i_upper, iter, qi_idx, KK
      DOUBLE PRECISION z_upper, z_lower
      DOUBLE PRECISION To_t, So_t, P_t
      DOUBLE PRECISION Qi_t, Mi_t, Hi_t, Ci_t, Bi_t
      DOUBLE PRECISION Ti_t, Si_t, binner_t, ui_t, Ep_t, dqdz
      DOUBLE PRECISION Ep_sum, Vol_sum, Mom_sum, Heat_sum
      DOUBLE PRECISION Salt_sum, Gas_sum
      DOUBLE PRECISION Qi_upper, Mi_upper, Hi_upper, Bi_upper, ci_upper
      DOUBLE PRECISION ui_upper, binner_upper, Ti_upper, Si_upper
      DOUBLE PRECISION rho_a_upper, rho_i_upper, Ep_upper
      DOUBLE PRECISION Qi_lower, Mi_lower, Hi_lower, Bi_lower, ci_lower
      DOUBLE PRECISION ui_lower, binner_lower, Ti_lower, Si_lower 
      DOUBLE PRECISION rho_a_lower, rho_i_lower, Ep_lower
      DOUBLE PRECISION epsQ
      PARAMETER (epsQ = 1d-12)

C     External models for gas bubbles
      DOUBLE PRECISION GAS_RHO_OP, SLIP_VELOCITY_OP, WUEST_OP
      DOUBLE PRECISION CP_MODEL_OP, RHO_OP

C     Y is input/output vector for DLSODE
C       Y(1) = outer plume volume flux
C       Y(2) = outer plume momentum flux
C       Y(3) = outer plume heat flux
C       Y(4) = outer plume salt flux
C       Y(5) = outer gas concentration
C     For ODEPACK solver. See ODEPACK documentation for definitions.
      NEQ = 5
      LRW = 400
      LIW = 200
      ITOL = 1
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      MF = 21
      NG = 1
      IWORK(7) = 2 ! To limit number of times repeat error messages are printed
      JT = 2

C     For BrentQ Solver.
      ATOL_BQ = 1.0D-10
      RTOL_BQ = 1.0D-6
      FTOL_BQ = 1.0D-10

C      RTOL = 1.0D-6
C      ATOL = 1.0D-8

      ! set positive dz and derive Nz
      dz = 0.1D0                        ! physical vertical spacing in meters
      zf = 0
      z0 = ABS(iceDepth)               ! e.g. 490.0
      Nz = INT((z0 - zf) / dz) + 1
      iter = 1
      us = 0.20D0

      DO K = 1,Nr+1
        zProf(k) = rf(k)
        zProfAbs(k) = abs(rf(k)) ! Consistent with ta/saProf
      END DO
      
      DO K = 1, Nz
      zProf_f(k)      = -z0 + dz * (K-1) ! Consistent with QiProf
      zProfAbs_f(k)   = abs(-z0 + dz * (K-1)) 
      QoProf_f(K) = 0.D0
      MoProf_f(K) = 0.D0
      HoProf_f(K) = 0.D0
      BoProf_f(K) = 0.D0
      CoProf_f(K) = 0.D0
      END DO

      ! Use surface inner plume at QiProf_f(Nz) since QiProf_f(Nz)=surface
      Qi = QiProf_f(plumeLayers-20)
      Mi = MiProf_f(plumeLayers-20)
      Hi = HiProf_f(plumeLayers-20)
      Bi = BiProf_f(plumeLayers-20)
      Ci = CiProf_f(plumeLayers-20)

      ui = Mi / Qi
      Ti = Hi / (Qi * 4185.D0 * 1031.D0) - 273.15D0
      Si = Bi / Qi
      ci = Ci / Qi
      binner = SQRT(2.D0 * Qi / (PI * ui))

      print *, "Initial inner: ", ui, Ti, Si, binner

      Qo = (1.D0 + 0.1D0) * Qi 
      To = (Ti + taProf(1) * 0.1D0) * Qi / Qo
      So = (Si + saProf(1) * 0.1D0) * Qi / Qo
      Co = (Ci + 0.D0 * 0.1D0) * Qi / Qo

      ! initial densities at top surface depth
      depth = zProfAbs_f(plumeLayers)   ! 0
      rho_a = RHO_OP(taProf(1), saProf(1), depth)
      rho_i = RHO_OP(Ti, Si, depth)
      rho_o = RHO_OP(To, So, depth)

      print *, "in outer plume, depth: ", depth 

C     Use a Fr Condition to set initial velocity
      AX = 0.00001D0
      BX = 1.D0 

      FAX = WUEST_OP(AX, Qo, binner, rho_a, rho_o)
      FBX = WUEST_OP(BX, Qo, binner, rho_a, rho_o)
      print *, "in surf, FAX, FBX: ", FAX, FBX

      CALL BRENTQ_OP(WUEST_OP, AX, BX, FAX, FBX, XZERO, 
     &              FZERO, IFLAG, MAXITER, ATOL_BQ, RTOL_BQ, FTOL_BQ,
     &              Qo, binner, rho_a, rho_o)
        
      IF (IFLAG .EQ. 0) THEN
C        PRINT *, 'Root found: U = ', XZERO
C       PRINT *, 'Residual = ', FZERO
        uo = XZERO
      ELSE
C        PRINT *, 'Root not found, IFLAG = ', IFLAG
        uo = 0.D0
      END IF

      Qo = -Qo
      uo = -uo

      Y(1) = Qo 
      Y(2) = -Qo * 0.02563317!uo
      Y(3) = 0.15701733 * Qo
      Y(4) = 32.816 * Qo 
      Y(5) = Co * Qo
      print *, "Mo: ", Mo
      print *, "Qo, uo, To, So, Co, rho_o: ", Qo, uo, To, So, Co, rho_o 
      ! --------------------------------------------------------------------
      ! Main solve loop: integrate outer plume downward in depth (from surface to bottom)
      ! We step T from z0 down to 0 in increments dz: start at surface depth = z0
      ! In this code we keep T as current depth (meters)
      ! --------------------------------------------------------------------
      T = z0
      TOUT = T - dz    ! next depth to integrate to (one step down)
      DO IOUT = 1, Nz
        T = z0 - (IOUT-1) * dz
        TOUT = T - dz

        IF (ISTATE .GT. 0) THEN
        
          CALL DLSODE(DOUBLEPLUME_OUTER_DERIVS, NEQ, Y, T, 
     &            TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT,
     &            RWORK, LRW, IWORK, LIW, JEX_OP, MF)
        END IF ! ISTATE .GT. 0

      ! If solver succeeded, update T and prepare next step
        IF (ISTATE .GT. 0) THEN
        IF (Y(1) .LT. 0) THEN 
      ! Save profile
          QoProf_f(Nz - IOUT + 1) = Y(1)   ! store reversed: index Nz->surface mapping
          MoProf_f(Nz - IOUT + 1) = Y(2)
          HoProf_f(Nz - IOUT + 1) = Y(3)
          BoProf_f(Nz - IOUT + 1) = Y(4)
          CoProf_f(Nz - IOUT + 1) = Y(5)
        ELSE 
          QoProf_f(Nz - IOUT + 1) = 0.d0   ! store reversed: index Nz->surface mapping
          MoProf_f(Nz - IOUT + 1) = 0.d0
          HoProf_f(Nz - IOUT + 1) = 0.d0
          BoProf_f(Nz - IOUT + 1) = 0.d0
          CoProf_f(Nz - IOUT + 1) = 0.d0
        END IF
        END IF ! ISTATE .GT. 0

        IF (Y(2) .LT. 0.005) THEN
          ISTATE = -1 
        END IF

        depth = ABS(iceDepth) - T 
        IF (j_peel .GT. 0) THEN
        IF (NINT(depth) .EQ. NINT(z_peel(j_peel))) THEN 
          print *, "Reached depth of peeling event."
          print *, "j_peel: ", j_peel 
          print *, "z_peel, depth: ", z_peel(j_peel), depth 
          j_peel = j_peel - 1
          ISTATE = -1
        END IF 
        END IF

      ! If solver failed / outer plume stopped, try continuous-peeling (CP) search
      ! (ISTATE <= 0 indicates the ODE solver stopped / signaled)
        IF (ISTATE .LE. 0) THEN
          done = .FALSE.
          iter = 1
          DO WHILE (.NOT. done .AND. iter .LE. 1)
            IF (T .LT. 1) THEN
              T = 1.d0 
            END IF 
            depth = ABS(iceDepth) - T

            smallDiff_l = 1.0D30
            i_upper = 1
            DO K = 1, Nz
              diff_l = ABS(zProfAbs_f(K) - depth)
              IF (diff_l .LT. smallDiff_l) THEN
                smallDiff_l = diff_l
                i_upper = K
              END IF
            END DO

            Qi = QiProf_f(i_upper)
            Mi = QiProf_f(i_upper)
            Hi = QiProf_f(i_upper)
            Bi = BiProf_f(i_upper)

            ui = Mi/Qi 
            binner = SQRT(2 * Qi / (PI * ui))
            Ti = Hi/(Qi * 4000 * 1031) - 273.15 
            Si = Bi/Qi

      ! Determine Search Window
            z_upper = depth 
            z_lower = z_upper + iter * binner   !1.D0 = iter

            IF (z_lower .LT. 1.D0) z_lower = 1.D0 ! Hit diffuser

            smallDiff_l = 1.0D30
            i_upper = 1
            DO K = 1, Nz
              diff_l = ABS(zProfAbs_f(K) - z_upper)
              IF (diff_l .LT. smallDiff_l) THEN
                smallDiff_l = diff_l
                i_upper = K
              END IF
            END DO

            smallDiff_l = 1.0D30
            i_lower = 1
            DO K = 1, Nz
              diff_l = ABS(zProfAbs_f(K) - z_lower)
              IF (diff_l .LT. smallDiff_l) THEN
                smallDiff_l = diff_l
                i_lower = K
              END IF
            END DO

      ! Ensure indices are in correct order (i_lower <= i_upper)
            IF (i_lower .GT. i_upper) THEN
              temp = i_lower
              i_lower = i_upper
              i_upper = temp
            END IF

            IF (i_lower .LT. 1) THEN 
              i_lower = 1
            END IF 

            CALL linint_op(Nr, zProfAbs, taProf, z_upper, Tambient)
            CALL linint_op(Nr, zProfAbs, saProf, z_upper, Sambient)

            Qi_upper = QiProf_f(i_upper)
            Mi_upper = MiProf_f(i_upper)
            Hi_upper = HiProf_f(i_upper)
            Bi_upper = BiProf_f(i_upper)

            ui_upper = Mi_upper/Qi_upper 
            binner_upper = SQRT(2 * Qi_upper / (PI * ui_upper))
            Ti_upper = Hi_upper/(Qi_upper * 1031 * 4000) - 273.15
            Si_upper = Bi_upper/Qi_upper

            rho_a_upper = RHO_OP(Tambient, Sambient, z_upper)
            rho_i_upper = RHO_OP(Ti_upper, Si_upper, z_upper)
            Ep_upper = CP_MODEL_OP(ui_upper, us, binner_upper, 
     &                        rho_i_upper, rho_a_upper)

      ! Lower Data

            CALL linint_op(Nr, zProfAbs, taProf, z_lower, Tambient)
            CALL linint_op(Nr, zProfAbs, saProf, z_lower, Sambient)

            Qi_lower = QiProf_f(i_lower)
            Mi_lower = MiProf_f(i_lower)
            Hi_lower = HiProf_f(i_lower)
            Bi_lower = BiProf_f(i_lower)

            ui_lower = Mi_lower/Qi_lower
            binner_lower = SQRT(2 * Qi_lower / (PI * ui_lower))
            Ti_lower = Hi_lower/(Qi_lower * 1031 * 4000) - 273.15
            Si_lower = Bi_lower/Qi_lower

            rho_a_lower = RHO_OP(Tambient, Sambient, z_lower)
            rho_i_lower = RHO_OP(Ti_lower, Si_lower, z_lower)
            Ep_lower = CP_MODEL_OP(ui_lower, us, binner_lower, 
     &                        rho_i_lower, rho_a_lower)
     
            Ep_sum = Ep_lower * dz
            Vol_sum = Ep_lower * dz
            Mom_sum = ui_lower * Ep_lower * dz
            Heat_sum = Ti_lower * Ep_lower * dz
            Salt_sum = Si_lower * Ep_lower * dz
            Gas_sum = Ci_lower * Ep_lower * dz

            DO K = i_lower+1, i_upper-1
              Qi_t = QiProf_f(K)
              Mi_t = MiProf_f(K)
              Hi_t = HiProf_f(K)
              Bi_t = BiProf_f(K)
              Ci_t = CiProf_f(K)
        ! compute primitives (protect against zero Q)
              ui_t = Mi_t / Qi_t
              binner_t = SQRT(2.D0 * Qi_t /(PI * ui_t) )
              Ti_t = Hi_t / (Qi_t * 1031.D0 * 4000.D0) - 273.15D0
              Si_t = Bi_t / Qi_t
              ci_t = Ci_t / Qi_t

        ! ambient at this K
              depth = zProfAbs_f(K)
              IF (depth .GT. ABS(iceDepth)) THEN
                Tambient = taProf(Nr)   ! bottom ambient (if taProf indexed by depth)
                Sambient = saProf(Nr)
              ELSE
                CALL linint_op(Nr, zProfAbs, taProf, depth, Tambient)
                CALL linint_op(Nr, zProfAbs, saProf, depth, Sambient)
              END IF
        
              rho_a = RHO_OP(Tambient, Sambient, depth)
              rho_i = RHO_OP(Ti_t, Si_t, depth)

              Ep_t = CP_MODEL_OP(ui_t, us, binner_t, rho_i, rho_a)
        
              Ep_sum = Ep_sum + Ep_t * dz
              Vol_sum = Vol_sum + Ep_t * dz
              Mom_sum = Mom_sum + ui_t * Ep_t* dz
              Heat_sum = Heat_sum + Ti_t * Ep_t * dz
              Salt_sum = Salt_sum + Si_t * Ep_t * dz
              Gas_sum = Gas_sum + Ci_t * Ep_t * dz
            END DO 
            Ep_sum = Ep_sum + Ep_lower * dz
            Vol_sum = Vol_sum + Ep_lower * dz
            Mom_sum = Mom_sum + ui_lower * Ep_lower * dz
            Heat_sum = Heat_sum + Ti_lower * Ep_lower * dz
            Salt_sum = Salt_sum + Si_lower * Ep_lower * dz
            Gas_sum = Gas_sum + Ci_lower * Ep_lower * dz

            Qo = -Ep_sum 
            uo = -Mom_sum
            To = -Heat_sum
            So = -Salt_sum
            co = -Gas_sum 
            IF (ABS(Qo). LT. 0.000001) THEN 
              To = Tambient
              So = Sambient 
              Co = 0.d0
            ELSE
              To = To / (Qo)
              So = So / Qo
              co = Co / Qo
            END IF 
            rho_o = RHO_OP(To, So, z_lower)


      ! find outer-plume velocity magnitude with BrentQ using the magnitude of Qo
      ! pass positive Qo magnitude to root solver, then make uo negative (descending)
            FAX = WUEST_OP(AX, Qo, binner_upper, 
     &                  rho_a_upper, rho_o)
            FBX = WUEST_OP(BX, Qo, binner_upper, 
     &                  rho_a_upper, rho_o)

            CALL BRENTQ_OP(WUEST_OP, AX, BX, FAX, FBX, XZERO, 
     &            FZERO, IFLAG, MAXITER, ATOL_BQ, RTOL_BQ, FTOL_BQ, 
     &            Qo, binner_upper, rho_a_upper, rho_o)

            IF (IFLAG .EQ. 0) THEN
              uo = - ABS(XZERO)
            ELSE
              uo = 0.D0
            END IF

            bouter = SQRT(2 * ABS(Qo) / (ABS(uo) * PI) 
     &                  + binner_upper**2)

            IF (Qo .GT. 0.D0) THEN
              Qo = -Qo   ! ensure negative
            END IF
            Mo = ABS(Qo * uo) 
            Ho = To * Qo
            Bo = So * Qo 
            Co = co * Qo 

            Ep_t = Ep_lower
            ui = ui_lower 

            dqdz = PI * binner_upper * (0.055 * (ui + 0.0 * uo) 
     &                   + 0.11 * uo)+ PI * bouter * 0.11 * uo + Ep_t

            IF (dqdz .LT. 0.D0) THEN 
              done = .TRUE. 
              Y(1) = Qo
              Y(2) = Mo 
              Y(3) = Ho 
              Y(4) = Bo 
              Y(5) = Co 
              ISTATE = 1
              IWORK(11) = 0
            ELSE 
              iter = iter + 1
              IF (iter .GT. 3) THEN
                done = .TRUE.
              END IF 
              Y(1) = 0.D0 
              Y(2) = 0.D0 
              Y(3) = 0.D0 
              Y(4) = 0.D0 
              Y(5) = 0.D0 
              ISTATE = -1
            END IF 
            T = TOUT 
            TOUT = T - dz
          END DO ! WHILE DO

        END IF ! ISTATE .LT. 0

        QoProf_f(Nz - IOUT + 1) = Y(1)
        MoProf_f(Nz - IOUT + 1) = Y(2)
        HoProf_f(Nz - IOUT + 1) = Y(3)
        BoProf_f(Nz - IOUT + 1) = Y(4)
        CoProf_f(Nz - IOUT + 1) = Y(5)

      END DO  ! IOUT loop

C      J = 1
C      sum_vol = 0.D0
C      sum_mom = 0.D0
C      sum_heat = 0.D0
C      sum_gasHeat = 0.D0
C      sum_salt = 0.D0
C      sum_gasFlux = 0.D0
C      sum_gasCon = 0.D0
      DO K = 1, Nr+1 
        QoProf(K) = 0.d0
        MoProf(K) = 0.d0 
        HoProf(K) = 0.d0
        BoProf(K) = 0.d0 
        CoProf(K) = 0.d0 
      END DO

      end_i = 100
      start_i = 0
      DO K = 1, Nr
      DO KK = start_i, end_i
      QoProf(K) = QoProf_f(KK) + QoProf(K)
      MoProf(K) = MoProf_f(KK) + MoProf(K)
      HoProf(K) = HoProf_f(KK) + HoProf(K)
      BoProf(K) = BoProf_f(KK) + BoProf(K)
      END DO
      start_i = end_i 
      end_i = start_i + 100
      END DO

C      DO K = 1, Nz+1
C        IF (MOD(K, 100) .EQ. 0) THEN
C          QoProf(J) = sum_vol
C          MoProf(J) = sum_mom
C          HoProf(J) = sum_heat
C          BoProf(J) = sum_salt
C          CoProf(J) = sum_gasCon
C          J = J + 1
C          sum_vol = 0.D0
C          sum_mom = 0.D0
C          sum_heat = 0.D0
C          sum_gasHeat = 0.D0
C          sum_salt = 0.D0
C          sum_gasFlux = 0.D0
C          sum_gasCon = 0.D0
          
C        ELSE
C          sum_vol = sum_vol + QoProf_f(K)
C          sum_mom = sum_mom + MoProf_f(K)
C          sum_heat = sum_heat + HoProf_f(K)
C          sum_salt = sum_salt + BoProf_f(K)
C          sum_gasCon = sum_gasCon + CoProf_f(K)  
C        ENDIF
C      ENDDO

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE DOUBLEPLUME_OUTER_DERIVS (NEQ, T, Y, YDOT)
      
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

C      EXTERNAL gas_model ! Not sure if this is correct


      INTEGER NEQ
      DOUBLE PRECISION T, Y(5), YDOT(5)

      DOUBLE PRECISION Tambient, Sambient, Cambient, Tin
      DOUBLE PRECISION m_dot, Sb, Tb
      DOUBLE PRECISION To, So, Uo, bouter
      DOUBLE PRECISION Ti, Si, ui, binner, ci 
      DOUBLE PRECISION Qi, Mi, Hi, Bi
      DOUBLE PRECISION rho_i, rho_a, rho_b 
      DOUBLE PRECISION us, Ep, beta_M, beta_T, nb, Cg, Xi, Fpd
      DOUBLE PRECISION P, co, delDiss, depth


      DOUBLE PRECISION RHO_OP, CP_MODEL_OP
      EXTERNAL RHO_OP, CP_MODEL_OP
      IF (T .LT. 0) THEN 
        depth = ABS(iceDepth) - 1
      ELSE 
        depth = ABS(iceDepth) - T
      END IF 

C Interpolate from imposed ambient profiles 
      IF ( depth .LE. 3 ) THEN    ! If we are at the surface
        Tambient = taProf(1)
        Sambient = saProf(1)
        Qi = QiProf_f(plumeLayers-20)
        Mi = MiProf_f(plumeLayers-20)
        Hi = HiProf_f(plumeLayers-20)
        Bi = BiProf_f(plumeLayers-20)
      ELSEIF ( depth .GE. ABS(iceDepth) - 2 ) THEN ! If we are at diffuser depth
        Tambient = taProf(iceDepthK)
        Sambient = saProf(iceDepthK)
        Qi = QiProf_f(1)
        Mi = MiProf_f(1)
        Hi = HiProf_f(1)
        Bi = BiProf_f(1)
      ELSE
        CALL linint_op(Nr, zProfAbs, taProf, depth, Tambient)
        CALL linint_op(Nr, zProfAbs, saProf, depth, Sambient)
        Qi = QiProf_f(NINT(T/0.1))
        Mi = MiProf_f(NINT(T/0.1))
        Hi = HiProf_f(NINT(T/0.1))
        Bi = BiProf_f(NINT(T/0.1))
        Ci = CiProf_f(NINT(T/0.1))
C        CALL linint_op(plumeLayer, zProfAbs_f, 
C     &            QiProf_f(plumeLayers:1:-1), depth, Qi)
C        CALL linint_op(plumeLayers, zProfAbs_f,
C     &            MiProf_f(plumeLayers:1:-1), depth, Mi)
C        CALL linint_op(plumeLayers, zProfAbs_f, 
C     &            HiProf_f(plumeLayers:1:-1), depth, Hi)
C        CALL linint_op(plumeLayers, zProfAbs_f, 
C     &            BiProf_f(plumeLayers:1:-1), depth, Bi)

      ENDIF
      
      
C Get primative variables from integral quantities
      ui = Mi / Qi
      IF (ui .LE. 0.01) THEN
        ui = 0.01
      ENDIF
      binner = SQRT(2 * Qi / (PI * ui))
      Ti = Hi / (4000.d0 * 1031.d0 * Qi) - 273.15
      Si = Bi / Qi
      ci = Ci / Qi
      
      IF (ABS(Y(1)) .GE. 0.0001D0) THEN
C        print *, "Outer plume exists"
        To = Y(3) / (Y(1))
        So = Y(4) / Y(1)
        uo = Y(2) / Y(1)
        co = Y(5) / Y(1)
        bouter = SQRT(2 * ABS(Y(1)) / (ABS(uo) * PI) + binner**2)
      ELSE 
        To = Tambient
        So = Sambient
        uo = 0
        co = 0
        bouter = 0
      ENDIF

      rho_i = RHO_OP(Ti, Si, depth)
      rho_a = RHO_OP(Tambient, Sambient, depth)
      rho_o = RHO_OP(To, So, depth)

      Cambient = 0.D0

      us = 0.2D0 !SLIP_VELOCITY(rho_a, rho_b, Ti, d_gas) 
      IF (ui .EQ. 0.d0) THEN 
      EP = 0.d0
      ELSE
      Ep = CP_MODEL_OP(ui, us, binner, rho_i, rho_a)
      END IF

C      IF (ABS(Ep) .GT. 5) Ep = -0.00001d0
          
C Mass Conservation (water)
      YDOT(1) = PI * binner * (0.055 * (ui + 0.0*uo) + 0.11 * uo)
     &             + PI * bouter * 0.11 * uo + Ep

C Momentum Conservation
      YDOT(2) = 1/1.10 * (-0.5 * PI * 9.81 * (bouter**2 - binner**2) 
     &            / 1031 * (rho_a - rho_o) + PI * binner * 
     &            (0.055 * (ui + 0.0*uo)*uo + 0.11 * uo * ui) 
     &            + Ep*ui)

C Heat Conservation (water)
      YDOT(3) = (PI * binner * (0.055 * (ui + 0.0*uo) * 
     &            (To) + 0.11 * uo * (Ti)) 
     &            + PI * bouter * 0.11 * uo * (Tambient) 
     &            + Ep * (Ti)); 

C Salinity Conservation
      YDOT(4) = PI * binner * (0.055 * (ui + 0.0*uo) * So 
     &            + 0.11 * uo * Si) 
     &            + PI * bouter * 0.11 * uo * Sambient 
     &            + Ep*Si;

C Gas Concentration
      YDOT(5) = PI * binner * (0.055 * (ui + 0.0*uo) * co 
     &            + 0.11 * uo * ci) 
     &            + PI * bouter * 0.11 * uo * cambient 
     &            + Ep * ci 

      ydot = -ydot

C      print *, "YDOT:"
C      print *, "YDOT(1): ", YDOT(1)
C      print *, "YDOT(2): ", YDOT(2)
C      print *, "YDOT(3): ", YDOT(3)
C      print *, "YDOT(4): ", YDOT(4)
C      print *, "YDOT(5): ", YDOT(5)
C      print *, "YDOT(6): ", YDOT(6)
C      print *, "YDOT(7): ", YDOT(7)

      RETURN 
      END 

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION WUEST_OP(uo, Qo, binner, rho_a, 
     &                                  rho_o)
 
C     #include "DOUBLEPLUME_OPTIONS.h"
C     #include "SIZE.h"
C     #include "EEPARAMS.h"
C     #include "PARAMS.h"
C     #include "GRID.h"
C     #include "DYNVARS.h"
C     #include "DOUBLEPLUME.h"

      DOUBLE PRECISION uo, Qo, binner, rho_a, rho_o
      DOUBLE PRECISION Fr, bouter
      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

C      g = 9.81
C      Fr0 = 1.6
C      rho_r = 1025.0D0

      Fr = 0.1
      bouter = SQRT(2 * Qo/(PI * uo) + binner**2)

C Residual: target Froude number
      WUEST_OP = uo - Fr * SQRT(abs((bouter - binner) * 9.81 
     &            *(rho_a - rho_o)/rho_o))

      
      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION CP_MODEL_OP(ui, us, binner, 
     &                                  rho_i, rho_a)

#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      DOUBLE PRECISION ui, us, binner, rho_i, rho_a
      
C      c2 = 1
C      epsilon = 0.0175
C      rho_r = 1025
C      print *, "in Cp model, us, ui, binner: ", us, ui, binner
C      print *, "rho_a, rho_i: ", rho_a, rho_i 
C      print *, "PI: ", PI 

      CP_MODEL_OP = 1 * 0.075 * (us/ui)**2 * 9.81 * 
     &          (rho_a - rho_i)/1031 * PI * binner**2 / (2 * ui)
      
C      print *, "CP_MODEL_OP: ", CP_MODEL_OP
      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RHO_OP(T,S,z)

C Equation of state (UNESCO 1983)

C     T = temperature (deg C)
C     S = salinity (PSU)
C     z = depth (m)

      DOUBLE PRECISION T,S,z
      DOUBLE PRECISION rho_0, g, P
      DOUBLE PRECISION kw, Aw, Bw, k0
      DOUBLE PRECISION bulk_modulus
      DOUBLE PRECISION A, B, C, rho_w,rho_zero

      PARAMETER(rho_0=1027)
      PARAMETER(g=9.81)

      P= rho_0*g*abs(z)*1.0E-5
   
C ... RHO_1 (in situ) (saved in Y(7))
      kw= 19652.21+ 148.4206*T- 2.327105*T**2+ 
     &    1.360477e-2*(T**3)-5.155288e-5*(T**4)
      Aw= 3.239908+ 1.43713e-3*T+ 1.16092e-4*T**2-
     &    5.77905e-7*T**3
      Bw= 8.50935e-5- 6.12293e-6*T + 5.2787e-8*(T**2)
      k0= kw + (54.6746- 0.603459*T+ 1.09987e-2*(T**2)
     &    -6.1670e-5*(T**3))*S +(7.944e-2 + 1.6483e-2*
     &    T- 5.3009e-4*(T**2))*(S**1.5)
      A=  Aw+ (2.2838e-3- 1.0981e-5*T- 1.6078e-6*(T**2))
     &    *S+ 1.91075e-4*(S**1.5)
      B= Bw+ (-9.9348e-7+ 2.0816e-8*T+ 9.1697e-10*T**2)*S
      bulk_modulus= k0+ A*P+ B*P**2

      A= 8.24493e-1- 4.0899e-3*T+ 7.6438e-5*T**2-
     &   8.2467e-7*T**3+5.3875e-9*T**4
      B= -5.72466e-3 + 1.0227e-4*T- 1.6546e-6*T**2
      C= 4.8314e-4
      rho_w= 999.842594 + 6.793952e-2*T- 9.095290e-3*T**2+
     &       1.001685e-4*T**3-1.120083e-6*T*
     &       *4+ 6.536336e-9*T**5
      rho_zero= rho_w+ A*S + B*(S**1.5)+ C*(S**2)

      RHO_OP = rho_zero/(1- (P/bulk_modulus))

      RETURN
      END

C-------------------------------------------------------------

C Dummy routine for ODEPACK. Necessary for Jacobian matrix if stiff ODEs.

      SUBROUTINE jex_op()
      RETURN
      END
C----------------------------------------------------------------------

       SUBROUTINE BRENTQ_OP(FUNC, AX, BX, FAX, FBX, XZERO, FZERO, IFLAG,
     &                  MAXITER, ATOL_BQ, RTOL, FTOL,
     &                  Qo, binner, rho_a, rho_o)
C=======================================================================
C  BRENTQ — Find a root of FUNC between AX and BX.
C  Classic Brent’s method, F77-compatible (no modules or types).
C=======================================================================
C  Arguments:
C     FUNC     - external function providing f(x)
C     AX, BX   - initial bracketing interval
C     FAX, FBX - function values at AX, BX
C     XZERO    - returned root approximation
C     FZERO    - function value at XZERO
C     IFLAG    - status flag: 0=root found, -2=max iterations reached
C     MAXITER  - maximum number of iterations
C     ATOL     - absolute tolerance
C     RTOL     - relative tolerance
C     FTOL     - function tolerance
C=======================================================================

      IMPLICIT NONE
      INTEGER MAXITER, IFLAG, I
      DOUBLE PRECISION AX, BX, FAX, FBX, XZERO, FZERO
      DOUBLE PRECISION ATOL_BQ, RTOL, FTOL
      DOUBLE PRECISION XPRE, XCUR, XBLK, FPRE, FCUR, FBLK
      DOUBLE PRECISION SPRE, SCUR, SBIS, DELTA, STRY, DPRE, DBLK, XDELTA
      DOUBLE PRECISION Qo, binner, rho_a, rho_o
      LOGICAL cond1, cond2
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC

C-----------------------------------------------------------------------
C  Initialize
C-----------------------------------------------------------------------
      IFLAG = 0
      XPRE = AX
      XCUR = BX
      FPRE = FAX
      FCUR = FBX

C-----------------------------------------------------------------------
C  Main iteration loop
C-----------------------------------------------------------------------
      DO 100 I = 1, MAXITER
C  Ensure a bracketing interval
         IF (FPRE*FCUR .LT. 0.0D0) THEN
            XBLK = XPRE
            FBLK = FPRE
            SCUR = XCUR - XPRE
            SPRE = SCUR
         END IF

C  Make sure |f(xblk)| < |f(xcur)|
         IF (DABS(FBLK) .LT. DABS(FCUR)) THEN
            XPRE = XCUR
            XCUR = XBLK
            XBLK = XPRE
            FPRE = FCUR
            FCUR = FBLK
            FBLK = FPRE
         END IF

         DELTA = 0.5D0*(ATOL_BQ + RTOL*DABS(XCUR))
         SBIS = 0.5D0*(XBLK - XCUR)

C  Convergence test
         IF (DABS(FCUR) .LE. FTOL .OR. DABS(SBIS) .LT. DELTA) GOTO 200

C-----------------------------------------------------------------------
C  Attempt interpolation or bisection
C-----------------------------------------------------------------------
          cond1 = DABS(SPRE) .GT. DELTA
          cond2 = DABS(FCUR) .LT. DABS(FPRE)

          IF (cond1 .AND. cond2) THEN
            IF (XPRE .EQ. XBLK) THEN
C              Linear interpolation
               STRY = -FCUR*(XCUR - XPRE)/(FCUR - FPRE)
            ELSE
C              Inverse quadratic interpolation
               DPRE = (FPRE - FCUR)/(XPRE - XCUR)
               DBLK = (FBLK - FCUR)/(XBLK - XCUR)
               STRY = -FCUR*(FBLK*DBLK - FPRE*DPRE) /
     &                (DBLK*DPRE*(FBLK - FPRE))
            END IF

C           Check step size
            IF (2.0D0*DABS(STRY) .LT.
     &          DMIN1(DABS(SPRE), 3.0D0*DABS(SBIS) - DELTA)) THEN
               SPRE = SCUR
               SCUR = STRY
            ELSE
               SPRE = SBIS
               SCUR = SBIS
            END IF
         ELSE
C           Bisect
            SPRE = SBIS
            SCUR = SBIS
         END IF

C-----------------------------------------------------------------------
C  Update
C-----------------------------------------------------------------------
         XPRE = XCUR
         FPRE = FCUR

         IF (DABS(SCUR) .GT. DELTA) THEN
            XCUR = XCUR + SCUR
         ELSE
            IF (SBIS .GT. 0.0D0) THEN
               XDELTA = DELTA
            ELSE
               XDELTA = -DELTA
            END IF
            XCUR = XCUR + XDELTA
         END IF
         
         FCUR = FUNC(XCUR, Qo, binner, rho_a, rho_o)

         IF (DABS(FCUR) .LE. FTOL) GOTO 200
         IF (I .EQ. MAXITER) IFLAG = -2

  100 CONTINUE

C-----------------------------------------------------------------------
C  Converged or exited
C-----------------------------------------------------------------------
  200 CONTINUE
      XZERO = XCUR
      FZERO = FCUR
      RETURN
      END

C---------------------------------------------------------

      DOUBLE PRECISION FUNCTION GAS_RHO_OP(T_local, P)
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"      
      DOUBLE PRECISION T_local, P
C Now declare local variables
      DOUBLE PRECISION Ti_k, Z, b, ac, mu, alpha, a, AD, BD
      DOUBLE PRECISION coeffs(4), zr(3), zi(3), Vm
      EXTERNAL rroots_chebyshev_cubic

C        DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

C        Ru = 8.314472           ! Universal gas constant in J/(mol K)
C        Tcr_O = 154.58          ! Critical Temperature Oxygen in K
C        Pcr_O = 5042825.270     ! Critical Pressure Oxygen in Pa
C        omega_g = 0.0216          ! Bubble eccentricity factor
C        M_O = 2*15.9994/1000    ! Molar weight kg/mol
        
        Ti_k = T_local + 273.15      ! Convert to Kelvin

C Peng-Robinson EOS parameters
        b = 0.07780 * Ru * Tcr_O / Pcr_O
        ac = 0.45724 * Ru**2 * Tcr_O**2 / Pcr_O
        mu = 0.37464 + 1.54226*omega_g - 0.26992*omega_g**2
        alpha = (1 + mu * (1 - SQRT(Ti_k/Tcr_O)))**2
        a = ac * alpha

C Dimensionless parameters
        AD = a * P / (Ru**2 * Ti_k**2)
        BD = b * P / (Ru * Ti_k)

C Cubic coefficients for Z-factor: Z^3 + c2 Z^2 + c1 Z + c0 = 0
        coeffs(1) = 1
        coeffs(2) = -(1 - BD)
        coeffs(3) = AD - 3*BD**2 - 2*BD
        coeffs(4) = -(AD*BD - BD**2 - BD**3)

C Solve cubic for Z (compressibility factor)

      CALL rroots_chebyshev_cubic(coeffs, zr, zi)
        
C        if SIZE(zr) == 0:
C            raise RuntimeError("No real Z roots found for PR EOS.")
        Z = MAXVAL(zr)          ! choose vapor root

C Compute molar volume
        Vm = Z * Ru * Ti_k / P  ! m^3/mol

C Gas density
        GAS_RHO_OP = M_O / Vm      ! kg/m^3

        RETURN
        END
C---------------------------------------------------------

      DOUBLE PRECISION FUNCTION SLIP_VELOCITY_OP(rho_a, rho_b, 
     &                  T_local, d_gas)
      
      
C Might not need to include any of these but PARAMS.h? 10/23/25
C d_gas might be in params file and doesn't need to be passed 10/23/25
C #include "DOUBLEPLUME_OPTIONS.h"
C #include "SIZE.h"
C #include "EEPARAMS.h"
C #include "PARAMS.h"
C #include "GRID.h"
C #include "DYNVARS.h"
C #include "DOUBLEPLUME.h"
      DOUBLE PRECISION rho_a, rho_b, T_local, g_local
      DOUBLE PRECISION mu, mu_w, sigma, M
      DOUBLE PRECISION Eo, H, J, Re, Nd
C     DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0
      g_local = 9.81 ! gravitational accelaration


C Viscosity of continuous phase
      mu = 1721.20 * EXP(-0.0251 * T_local) / 0.000001 ! Pa s

C Viscosity in Briada's expts
      mu_w = 0.0009

C Calculate surface tension for water
      sigma = (75.7060 - 0.1511 * T_local) / 1000

C Calculate non-dimensionals
      M = g_local * mu**4 * (rho_a - rho_b) / (rho_a**2 * sigma**3)
      Eo = g_local * (rho_a - rho_b) * d_gas**2 / sigma


C Get reynolds number at terminal velocity
      H = 4.0D0/3.0D0 * Eo * M**(-0.149) * (mu/mu_w)**(-0.14)
C      print *, "H: ", H

      IF (H .GT. 2) THEN
C Elliptical wobbling regimes
        IF (H .GT. 59.3) THEN
            J = 3.42 * H**(0.441)
        ELSEIF (H .GT. 2) THEN
            J = 0.94 * H**(0.757)
        ENDIF

        Re = M**(-0.149) * (J-0.857)

      ELSE
C Spherical regime (Clift p.113)
        Nd = 4 * rho_a * (rho_a - rho_b) * g * d_gas**3 / (3 * mu**2)
C        print *, "Nd: ", H
        IF (Nd .LE. 73) THEN
            Re = Nd/24 - 1.7569e-4 * Nd**2 + 6.9252e-7 * Nd**3 
     &          - 2.3027e-10*Nd**4
        ELSEIF (Nd .LE. 580) THEN
            Re = 10**(-1.7095 + 1.33438 * LOG(Nd) 
     &          - 0.11591 * LOG(Nd)**2)
        ELSEIF (Nd .LE. 1.55e7) THEN
            Re = 10**(-1.81391 + 1.34671*LOG(Nd) - 0.12427*LOG(Nd)**2 
     &          + 0.006344*LOG(Nd)**3)
        ELSE
            Re = 10**(5.33283 - 1.21728*LOG(Nd) 
     &          + 0.19007*LOG(Nd)**2 - 0.007005*LOG(Nd)**3)
        ENDIF
      ENDIF

      SLIP_VELOCITY_OP = mu / (rho_a * d_gas) * Re
    
      IF (SLIP_VELOCITY_OP .LT. 0) THEN
        SLIP_VELOCITY_OP = 0.001
      ENDIF
    
      RETURN 
      END

C-------------------------------------------------------------
      subroutine linint_op(nx,xtab,ytab,x,y)
c
c     Given a value of x return a value of y based on interpolation
c     within a table of y values (ytab) corresponding to the x values
c     contained in the array xtab.  The subroutine assumes that the
c     values in xtab increase monotonically
c
c    John Mahaffy 2/12/95
c    Modified slightly TRC 2014

      integer nx
      double precision xtab(nx), ytab(nx), x, y

c     local variables
      integer i, i1
      double precision  wx

      if (x.lt.(xtab(1)).or.x.GT.(xtab(nx))) then
         write(6,*) 'x = ', x, '  is out of table range'
         stop
      endif
      do 100 i=2,nx
         if (x.le.xtab(i)) go to 200
  100    continue
  200 i1=i-1
      wx=(x-xtab(i1))/(xtab(i1+1)-xtab(i1))
      y=(1-wx)*ytab(i1)+wx*ytab(i1+1)
      return
      end