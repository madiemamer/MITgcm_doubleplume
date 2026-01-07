#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_PLUME_MODEL
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_INNER_MODEL (mythid)
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_PLUME_MODEL:
C     | o Uses ODEPACK to solve inner plume model
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
      EXTERNAL DOUBLEPLUME_INNER_DERIVS, JEX, JAC, EVENTFUNC

C     For ODEPACK
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(200), LIW, LRW,
     &         MF, NEQ, JT, NG, JROOT
      _RL T, TOUT

      _RL  RWORK(200), Y(7)

C     Y is input/output vector for DLSODE
C       Y(1) = inner plume volume flux
C       Y(2) = inner plume momentum flux
C       Y(3) = inner plume heat flux
C       Y(4) = inner gas heat flux
C       Y(5) = inner plume salt flux
C       Y(6) = inner gas flux
C       Y(7) = inner gas concentration

C     Other local variables
      INTEGER  k, mythid, i, j, KK
      INTEGER end_i, start_i
      DOUBLE PRECISION Tambient, Sambient, depth, P
      DOUBLE PRECISION rho_plume, rho_a, rho_b
      DOUBLE PRECISION mb, us, Qsg, Qg
      DOUBLE PRECISION R, u0
      DOUBLE PRECISION Ti_t, Si_t, P_t 
      DOUBLE PRECISION z0, zf, dz, Nz
      DOUBLE PRECISION Qi, Ji, ui, binner, Ti, Si, Ci
      DOUBLE PRECISION C0, S0, T0, q_m3s
      DOUBLE PRECISION sum_vol, sum_mom, sum_heat, sum_gasHeat
      DOUBLE PRECISION sum_salt, sum_gasFlux, sum_gasCon

      DOUBLE PRECISION GAS_RHO, SLIP_VELOCITY, WUEST, RHO, LIMANETO

      CHARACTER*(MAX_LEN_MBUF) msgBuf
      LOGICAL DISSOLVED, IMAG

      INTEGER MAXITER, IFLAG
      DOUBLE PRECISION ATOL_BQ, RTOL_BQ, FTOL_BQ
      DOUBLE PRECISION AX, BX, FAX, FBX, XZERO, FZERO
      PARAMETER (MAXITER=100)

C     For Brentq root solver.
      ATOL_BQ = 1.0D-10
      RTOL_BQ = 1.0D-6
      FTOL_BQ = 1.0D-10

C     For ODEPACK solver. See ODEPACK documentation for definitions.
      NEQ = 7
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

      DISSOLVED = .FALSE.

      Tambient = taProf(iceDepthK)
      Sambient = saProf(iceDepthK)

C     Set up local numerical grid, coordinate system is reversed
      dz = 0.1 ! m
      z0 = 0 ! 0 = depth of inlet
      zf = ABS(iceDepth) ! Surface
      Nz = INT((zf - z0) / dz)

      mb    = 1.0D-5        ! initial gas flux
      C0    = 0.005D0       ! initial concentration
      S0    = 0.01D0        ! initial salinity
      T0    = 0.293D0       ! initial temperature
      R     = 0.05D0        ! plume radius
      Qsg   = 0.01D0        ! p.Qsg
      us    = 0.2D0         ! slip velocity
      P     = ABS(iceDepth) * 9.81 * 1031 + 101000 ! pressure
      rho_b = (2 * 16/1000.D0 * P) / (8.314472 * (T0 + 273.15)) ! gas density
      q_m3s = mb / rho_b 
      Qg    = 1e-6          ! Initial gas flux

      IF (Qsg .LE. 0.D0) THEN
        AX = 0.00001D0
        BX = 0.01D0 
        FAX = WUEST(AX, us, mb, R, rho_a, rho_b)
        FBX = WUEST(BX, us, mb, R, rho_a, rho_b)
        CALL BRENTQ(WUEST, AX, BX, FAX, FBX, XZERO, FZERO, IFLAG, 
     &              MAXITER, ATOL_BQ, RTOL_BQ, FTOL_BQ,
     &              us, mb, R, 1031.d0, rho_b)
        IF (IFLAG .EQ. 0) THEN
C            PRINT *, 'Root found: U = ', XZERO
C            PRINT *, 'Residual = ', FZERO
            u0 = XZERO
        ELSE
C            PRINT *, 'Root not found, IFLAG = ', IFLAG
            u0 = 0.D0
        END IF
      ELSE
        u0 = LIMANETO(Qsg, q_m3s, R)
      ENDIF

      Qi    = PI * R**2 * u0 /2 ! Initial total flux
      Ji    = Qi * u0           ! Initial total momentum flux

C Initial conditions
      Y(1)  = Qi                ! Initial volume flux
      Y(2)  = Ji                ! Initial momentum flux
      Y(3)  = 4000 * 1031 * (T0 + 273.15) * Qi  ! Initial heat flux
      Y(4)  = mb * (0.5 * 4000) * (T0 + 273.15) ! Initial gas heat flux
      Y(5)  = S0 * Qi           ! Initial brine flux
      Y(6)  = Qi * rho_b        ! Initial gas volume flux
      Y(7)  = C0 * Qi           ! Initial gas concentration

      DO K = 1,Nr+1
        zProf(k) = rf(k)
        zProfAbs(k) = abs(rf(k)) 
      END DO

      DO K = 1,Nz+1
        zProf_f(k)      = -zf + dz * (K-1)
        zProfAbs_f(k)   = abs(-zf + dz * (K-1)) 
        QiProf_f(k)     = 0.D0
        MiProf_f(k)     = 0.D0
        HiProf_f(k)     = 0.D0
        HgProf_f(k)     = 0.D0
        BiProf_f(k)     = 0.D0
        MbProf_f(k)     = 0.D0
        CiProf_f(k)     = 0.D0
      END DO

C Start at depth of subglacial inlet
      T = z0 

C Next point at which to retrieve values
      TOUT = z0 + dz

C Set intial conditions
      QiProf_f(1) = Y(1)
      MiProf_f(1) = Y(2)
      HiProf_f(1) = Y(3)
      HgProf_f(1) = Y(4)
      BiProf_f(1) = Y(5)
      MbProf_f(1) = Y(6)
      CiProf_f(1) = Y(7)

      DO IOUT = 1, Nz+1     
       
        IF (ISTATE .GT. -1) THEN
          CALL DLSODE(DOUBLEPLUME_INNER_DERIVS, NEQ,Y,T,
     &               TOUT,ITOL,RTOL, ATOL, ITASK,ISTATE,
     &               IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)

C Checking Stop Criteria
C Stop Event 1: No longer buoyant
          Ti_t = Y(3) / (4000 * 1031 * Y(1)) - 273.15
          Si_t = Y(5) / Y(1)
          depth = ABS(iceDepth) - T 
          rho_plume = RHO(Ti_t, Si_t, depth)

          IF (depth .LE. 1) THEN                      ! At surface 
            tambient = taProf(1)
            sambient = saProf(1)
          ELSE IF (depth .GE. ABS(iceDepth) - 1) THEN ! At inlet
            tambient = taProf(iceDepthK)
            sambient = saProf(iceDepthK)
          ELSE                                        ! Elsewhere
            CALL linint(Nr,zProfAbs,taProf, 
     &                  depth, Tambient)
            CALL linint(Nr,zProfAbs,saProf, 
     &                  depth, Sambient)
          ENDIF 

          rho_a = RHO(tambient, sambient, depth)

C Stop Event 2: Bubbles Dissolved
          IF (Y(6)/mb .LT. 0.001) THEN
            print *, "Bubbles have dissolved."
            print *, "Y(6), mb: ", Y(6), mb
            ISTATE = -1
            DISSOLVED = .TRUE.
          ENDIF

C Stop Event 3: Are the bubbles still buoyant?
          P_t = zProfAbs_f(IOUT) * 9.81 * 1031 + 101000 
          rho_b = GAS_RHO(Ti_t, P_t)

          IF (rho_b .GT. rho_a) THEN 
            print *, "Bubbles have lost buoyancy."
            ISTATE = -1
          ENDIF

C Stop Event 4: Momentum < = 0
          IF (Y(2)/Y(1) .LE. 0.05) THEN
            print *, "Momentum too low. Vel, Mom: ", Y(2)/Y(1), Y(2)
            ISTATE = -1
          ENDIF

C Stop Event 5: Plume has reaced surface
          IF (depth .GE. 1) THEN
            T = TOUT           ! Define present Depth
            TOUT = T + dz      ! Define next Depth
          ELSE 
            ISTATE = -1        ! Stop integrating if reached surface
          END IF

C This section is entered once the plume has hit a stop event
        ELSE ! (ISTATE .GT. -1)
          IF (.NOT. DISSOLVED) THEN
            IF (IOUT .LE. Nz+1) THEN
              depth = ABS(iceDepth) - T
              IF (j_peel .GT. 1) THEN 
                IF (depth .LT. z_peel(j_peel) - 5 .AND. 
     &             depth .GT. 5) THEN
                j_peel = j_peel + 1 ! update peel counter
                z_peel(j_peel) = depth
                print *, "Peeling event has occured, T: ", T
                print *, "Depth: ", depth
                print *, "j_peel = ", j_peel
                END IF
              ELSE 
                j_peel = j_peel + 1 ! update peel counter
                z_peel(j_peel) = depth
                print *, "Peeling event has occured, T: ", T
                print *, "Depth: ", depth
                print *, "j_peel = ", j_peel
              END IF 

              ISTATE = 1
              Qi = (1 - 0.95) * Y(1)
              binner = SQRT(2 * Y(1)**2 / (PI * Y(2))) ! initial guess of b

              FAX = WUEST(0.001d0, us, Y(6), binner, rho_a, rho_b)
              FBX = WUEST(1.d0, us, Y(6), binner, rho_a, rho_b)

              CALL BRENTQ(WUEST, 0.001d0, 0.3d0, FAX, FBX, XZERO, FZERO,  
     &              IFLAG, MAXITER, ATOL_BQ, RTOL_BQ, FTOL_BQ,
     &              us, Y(6), binner, rho_a, rho_b)
              IF (IFLAG .EQ. 0) THEN
C                PRINT *, 'Root found: U = ', XZERO
C                PRINT *, 'Residual = ', FZERO
                ui = XZERO
              ELSE
C                PRINT *, 'Root not found, IFLAG = ', IFLAG
                ui = 0.D0
              END IF

              IF (ui .LE. 0.01) THEN
                ui = 0.01
              ENDIF

              print *, "ui: ", ui

              binner = SQRT(2 * ABS(Qi) / (PI * ui))
              Ti = Y(3) / (4000 * 1031 * Y(1)) - 273.15
              Si = Y(5) / Y(1)
              Ci = Y(7) / Y(1)

              T = T + 0.1 * binner
              TOUT = TOUT + 0.1 * binner
              
              Y(1) = Qi
              Y(2) = Qi * ui
              Y(3) = 1031 * 4000 * (Ti + 273.15) * Qi
              Y(4) = Y(6) * (4000 * 0.5) * (Ti + 273.15)               ! Bubble Heat Flux J/s, Tg = Ti
              Y(5) = Si * Qi                                  ! Salt Flux m^3 psu /s
              Y(6) = Y(6)                                     ! Bubble mass Flux kg/s
              Y(7) = Ci * Qi  

            ELSE ! If we have reached the surface
              Y(1) = 0
              Y(2) = 0
              Y(3) = 0
              Y(4) = 0
              Y(5) = 0
              Y(6) = 0
              Y(7) = 0
            ENDIF
          ELSE  ! If all bubbles have dissolved
            Y(1) = 0
            Y(2) = 0
            Y(3) = 0
            Y(4) = 0
            Y(5) = 0
            Y(6) = 0
            Y(7) = 0
            
          ENDIF
        ENDIF ! (ISTATE .GT. -1)

        QiProf_f(IOUT) = Y(1)
        MiProf_f(IOUT) = Y(2)
        HiProf_f(IOUT) = Y(3)
        HgProf_f(IOUT) = Y(4)
        BiProf_f(IOUT) = Y(5)
        MbProf_f(IOUT) = Y(6)
        CiProf_f(IOUT) = Y(7)
      
      ENDDO ! IOUT = iceDepthK-1, 1, -1

C       TRANSLATE LOCAL GRID TO MAIN GRID HERE 10/29/31
C      J = 1
C      sum_vol = 0.D0
C      sum_mom = 0.D0
C      sum_heat = 0.D0
C      sum_gasHeat = 0.D0
C      sum_salt = 0.D0
C      sum_gasFlux = 0.D0
C      sum_gasCon = 0.D0
      DO K = 1, Nr+1 
        QiProf(K) = 0.d0
        MiProf(K) = 0.d0 
        HiProf(K) = 0.d0
        HgProf(K) = 0.d0
        BiProf(K) = 0.d0 
        MbProf(K) = 0.d0
        CiProf(K) = 0.d0 
      END DO

C      end_i = 100
C      start_i = 0
C      DO K = 1, Nr
C      DO KK = start_i, end_i
C      QiProf(K) = QiProf_f(KK) + QiProf(K)
C      MiProf(K) = MiProf_f(KK) + MiProf(K)
C      HiProf(K) = HiProf_f(KK) + HiProf(K)
C      BiProf(K) = BiProf_f(KK) + BiProf(K)
C      HgProf(K) = HgProf_f(KK) + HgProf(K)
C      MbProf(K) = MbProf_f(KK) + MbProf(K)
C      CiProf(K) = CiProf_f(KK) + CiProf(K)
C      END DO
C      start_i = end_i 
C      end_i = start_i + 100
C      END DO


C Need to make sure this stops at ice depth and fill the rest with zeros
C      DO K = 1, Nz
C        IF (MOD(K, 100) .EQ. 0) THEN
C          QiProf(J) = sum_vol
C          MiProf(J) = sum_mom
C          HiProf(J) = sum_heat
C          HgProf(J) = sum_gasHeat
C          BiProf(J) = sum_salt
C          MbProf(J) = sum_gasFlux
C          CiProf(J) = sum_gasCon
C          J = J + 1
C          sum_vol = 0.D0
C          sum_mom = 0.D0
C          sum_heat = 0.D0
C          sum_gasHeat = 0.D0
C          sum_salt = 0.D0
C          sum_gasFlux = 0.D0
C          sum_gasCon = 0.D0
C        ELSE
C          sum_vol = sum_vol + QiProf_f(K)
C          sum_mom = sum_mom + MiProf_f(K)
C          sum_heat = sum_heat + HiProf_f(K)
C          sum_gasHeat = sum_gasHeat + HgProf_f(K)
C          sum_salt = sum_salt + BiProf_f(K)
C          sum_gasFlux = sum_gasFlux + MbProf_f(K)
C          sum_gasCon = sum_gasCon + CiProf_f(K)  
C        ENDIF
C      ENDDO
C      
      RETURN
      END 

C---------------------------------------------------------

      SUBROUTINE DOUBLEPLUME_INNER_DERIVS (NEQ, T, Y, YDOT)
      
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

C      EXTERNAL gas_model ! Not sure if this is correct


      INTEGER NEQ
      DOUBLE PRECISION T, Y(7), YDOT(7)

      DOUBLE PRECISION Tambient, Sambient
      DOUBLE PRECISION m_dot, Sb, Tb
      DOUBLE PRECISION To, So, Uo, bouter
      DOUBLE PRECISION Qo, Mo, Ho, Bo, Co
      DOUBLE PRECISION Ti, Si, ui, binner, mb, ci 
      DOUBLE PRECISION rho_i, rho_a, rho_b 
      DOUBLE PRECISION us, Ep, beta_M, beta_T, nb, Cg, Xi, Fpd
      DOUBLE PRECISION P, delDiss, depth_local

      DOUBLE PRECISION RHO
      DOUBLE PRECISION MDOT, BOUNDARY_SALINITY, BOUNDARY_TEMPERATURE
      DOUBLE PRECISION GAS_RHO, SLIP_VELOCITY, CP_MODEL
      DOUBLE PRECISION VOID_FRACTION, GAS_MT, GAS_HT, GAS_SOL

      EXTERNAL RHO
      EXTERNAL MDOT, BOUNDARY_SALINITY, BOUNDARY_TEMPERATURE
      EXTERNAL GAS_RHO, SLIP_VELOCITY, CP_MODEL
      EXTERNAL VOID_FRACTION, GAS_MT, GAS_HT, GAS_SOL

      depth_local = ABS(iceDepth) - T

C Interpolate from imposed ambient profiles
      IF (depth_local .LE. 1) THEN
       Tambient = taProf(1)
       Sambient = saProf(1)
       Qo = QoProf_f(1)
       Mo = MoProf_f(1)
       Ho = HoProf_f(1)
       Bo = BoProf_f(1)
       Co = CoProf_f(1)
      ELSEIF (depth_local .GE. ABS(iceDepth)-1) THEN
       Tambient = taProf(iceDepthK)
       Sambient = saProf(iceDepthK)
       Qo = QoProf_f(plumeLayers)
       Mo = MoProf_f(plumeLayers)
       Ho = HoProf_f(plumeLayers)
       Bo = BoProf_f(plumeLayers)
       Co = CoProf_f(plumeLayers)
      ELSE
       CALL linint(Nr,zProfAbs, taProf,
     &              depth_local, Tambient)
       CALL linint(Nr,zProfAbs, saProf,
     &              depth_local, Sambient)
        Qo = QoProf_f(NINT(T/0.1))
        Mo = MoProf_f(NINT(T/0.1))
        Ho = HoProf_f(NINT(T/0.1))
        Bo = BoProf_f(NINT(T/0.1))
        Co = CoProf_f(NINT(T/0.1))
      ENDIF

C Get primative variables from integral quantities
      ui = Y(2) / Y(1)
      binner = SQRT(2 * Y(1) / (PI * Ui))
      Ti = Y(3) / (4000 * 1031 * Y(1)) - 273.15
      mb = Y(6)
      Si = Y(5) / Y(1)
      ci = Y(7) / Y(1)

      IF (Qo .EQ. 0.D0) THEN
      To = Tambient
      So = Sambient
      uo = 0.D0
      co = 0.D0
      bouter = 0.D0
      ELSE 
      To = Ho/Qo 
      So = Bo/Qo 
      uo = Mo/Qo 
      co = Co/Qo 
      bouter = SQRT(2 * ABS(Qo) / (PI * ABS(uo)) + binner**2)
      ENDIF

      rho_i = RHO(Ti, Si, depth_local)
      rho_a = RHO(To, So, depth_local)
      P = depth_local  * 9.81 * rho_a  + 101000
      rho_b = GAS_RHO(Ti, P)

      m_dot = MDOT(Y, depth_local) 
      Sb = BOUNDARY_SALINITY(Y, m_dot) 
      Tb = BOUNDARY_TEMPERATURE(Sb, depth_local)

      us = 0.2D0
      Ep = CP_MODEL(ui, us, binner, rho_i, rho_a)
      beta_M = GAS_MT(Ti, rho_a, us, 0.01)
      beta_T = GAS_HT(Ti, rho_a, us, 0.01)
      nb = 6 * (mb/rho_b) / (PI * 0.01**3)

      Cg = GAS_SOL(rho_a, Ti, Si, P, mb) 

      Xi = VOID_FRACTION(mb, rho_b, binner, us, ui)
      Fpd = 0.64D0 * Xi * (rho_a - rho_b)

C Mass Conservation (water)
      YDOT(1) = PI * binner * (0.055 * (ui + 0.0*uo) + 0.11*uo) 
     &          + Ep + 2 * binner * m_dot

C Momentum Conservation
      YDOT(2) = (1/1.1) * (0.5 * PI * 9.81 * binner**2 / 1031 * 
     &          (Fpd + 1.10**2 * (1 - Xi) * (rho_a - rho_i)) +
     &          PI * binner * (0.055 * (ui + 0.0 * uo) * uo + 
     &          0.11 * uo * ui) + Ep*ui - 0.001 * 2 * binner * ui**2)

C Heat Conservation (water)
      YDOT(3) = 1031 * 4000 * (PI * binner * (0.055 * (ui + 0.0 * uo) 
     &          * (To + 273.15) + 0.110 * uo * (Ti + 273.15)) + 
     &          Ep*(Ti + 273.15) - 
     &          4 * 0.022 * SQRT(0.001) * binner * ui * (Ti - Tb))

C Heat Conservation (gas)
      YDOT(4) = -PI * nb * (0.01)**2/(ui + us) * rho_b *
     &          (0.5 * 4000) * beta_T * (Tg - Ti)
      YDOT(3) = YDOT(3) - YDOT(4) ! Subtract gas content

C Salinity Conservation
      YDOT(5) = PI * binner * (0.055 * (ui + 0.0*uo) * So 
     &          + 0.110 * uo * Si) 
     &          + Ep * Si 
     &          - 4 * 0.00062 * SQRT(0.001) * binner * ui * (Si - Sb)

C Mass Conservation (gas)
      YDOT(6) = -(PI * nb * (0.01)**2/(ui + us) * beta_M * (Cg - ci))
      delDiss = YDOT(6)

C Gas Concentration
      YDOT(7) = PI * binner * (0.055 * (ui + 0.0*uo) * co 
     &          + 0.110 * uo * ci) + Ep * ci - delDiss 

      RETURN 
      END 

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION WUEST(u, us, mb, binner, rho_a, 
     &                                  rho_b)
 
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      DOUBLE PRECISION u, us, mb, binner, rho_a, rho_b
      DOUBLE PRECISION xi_gas, rho_p
C      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

C      g = 9.81
C      Fr0 = 1.6
C      rho_r = 1025.0D0
      xi_gas = mb / ( 0.5 * PI * binner**2 * rho_b * 0.8D0**2 *
     &      (us + 2*u/(1+0.8D0**2)))

      rho_p = ((xi_gas * rho_b) +
     &       (1 - (xi_gas)) * rho_a)

C Residual: target Froude number
      WUEST = 1.6 - u / SQRT(2 * 0.8D0 * binner * 9.81 * 
     &        (rho_a - rho_p) / 1031)

      
      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION LIMANETO(Qg, Qsg, R)
 
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      DOUBLE PRECISION Qg, Qsg, R
C      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

      A = PI * R**2 / 2

      epsilon_o = Qg / (Qg + Qsg)
      LIMANETO = Qsg / ((1 - epsilon_o) * A)

      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION MDOT(Y, z)
 
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

C      EXTERNAL gas_model ! Not sure if this is correct

      DOUBLE PRECISION Y(7)
      DOUBLE PRECISION z
      DOUBLE PRECISION Ti, J, Q, B 
      DOUBLE PRECISION beta_loc, a, bb, c_loc 

C      cp = 4185.5
C      rho_r = 1025
C      gammaS = 6.2e-4
C      gammaT = 2.2e-2
C      Cd = 0.001
C      L_i = 334000
C      lam1 = -5.73e-2
C      lam2 = 0.08
C      lam3 = -7.61e-4
C z = zProfAbs(Nr) - T 
C      print *, "In MDOT function"
C      print *, "Depth: ", z

      Ti = Y(3) / (4000 * 1031 * Y(1)) - 273.15
      J = Y(2)
      Q = Y(1)
      B = Y(5)

C      print *, "Ti: ", Ti

      beta_loc = 4000 * 0.022 * SQRT(0.001) / 334000
      a = 1
      bb = 6.2e-4 * SQRT(0.001) * J/Q - beta_loc * 
     &          J/Q * (Ti - 0.08 - -7.61e-4 * z)
      c_loc = beta_loc * 6.2e-4 * SQRT(0.001) * J**2/Q**2 * (B/Q  
     &     * -5.73e-2 - Ti + 0.08 + -7.61e-4 * z)

      MDOT = (-bb + SQRT(bb**2 - 4 * a * c_loc)) / (2 * a)

      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BOUNDARY_SALINITY(Y, m_dot)
      
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      DOUBLE PRECISION Y(7), m_dot
      DOUBLE PRECISION Q, J, B

C      gammaS = 6.2e-4
C      Cd = 0.001

      Q = Y(1)
      J = Y(2)
      B = Y(5)

      BOUNDARY_SALINITY = 6.2e-4 * SQRT(Cd) * J/Q * B/Q / 
     &                  (m_dot + 6.2e-4 * SQRT(0.001) * J/Q)

      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BOUNDARY_TEMPERATURE(Sb, z)

#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      DOUBLE PRECISION Sb
      DOUBLE PRECISION z

C      lam1 = -5.73e-2
C      lam2 = 0.08
C      lam3 = -7.61e-4
    
      BOUNDARY_TEMPERATURE = -5.73e-2 * Sb + 08 + -7.61e-4 * z

      RETURN
      END
C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION CP_MODEL(ui, us, binner, rho_i, rho_a)

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

      CP_MODEL = 1 * 0.0175 * (us/ui)**2 * 9.81 * (rho_a - rho_i)/1031
     &          * PI * binner**2 / (2 * ui)
      
      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RHO(T,S,z)

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

      RHO= rho_zero/(1- (P/bulk_modulus))

      RETURN
      END

C-------------------------------------------------------------

C Dummy routine for ODEPACK. Necessary for Jacobian matrix if stiff ODEs.

      SUBROUTINE jex()
      RETURN
      END
C----------------------------------------------------------------------

       SUBROUTINE BRENTQ(FUNC, AX, BX, FAX, FBX, XZERO, FZERO, IFLAG,
     &                  MAXITER, ATOL_BQ, RTOL, FTOL,
     &                  us, mb, R, rho_a, rho_b)
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
      DOUBLE PRECISION us, mb, R, rho_a, rho_b 
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
      ! Initialize all working variables
      XBLK = XPRE
      FBLK = FPRE
      SPRE = BX - AX
      SCUR = BX - AX

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

         FCUR = FUNC(XCUR, us, mb, R, rho_a, rho_b)
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

      DOUBLE PRECISION FUNCTION GAS_RHO(T_local, P)
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
        GAS_RHO = M_O / Vm      ! kg/m^3

        RETURN
        END
C---------------------------------------------------------

      DOUBLE PRECISION FUNCTION SLIP_VELOCITY(rho_a, rho_b, 
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
      print *, "H: ", H

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

      SLIP_VELOCITY = mu / (rho_a * d_gas) * Re
    
      IF (SLIP_VELOCITY .LT. 0) THEN
        SLIP_VELOCITY = 0.001
      ENDIF
    
      RETURN 
      END


C---------------------------------------------------------

      DOUBLE PRECISION FUNCTION VOID_FRACTION(mb, rho_b, binner, us, ui)
      
      
C Might not need to include any of these but PARAMS.h? 10/23/25
C d_gas might be in params file and doesn't need to be passed 10/23/25
C #include "DOUBLEPLUME_OPTIONS.h"
C #include "SIZE.h"
C #include "EEPARAMS.h"
C #include "PARAMS.h"
C #include "GRID.h"
C #include "DYNVARS.h"
C #include "DOUBLEPLUME.h"

      DOUBLE PRECISION mb, rho_b, binner, us, ui
      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

C      gas_lam1 = 0.8
      
C      print *, "In VOID_FRACTION"
C      print *, "mb: ", mb 
C      print *, "rho_b: ", rho_b 
C      print *, "binner: ", binner 
C      print *, "us: ", us 
C      print *, "ui: ", ui 


      VOID_FRACTION = mb /(rho_b * (0.5 * PI * binner**2 
     &                * 0.8D0 * ( us + 2 * ui/(1 + 0.8D0**2))))

      RETURN
      END

C---------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION GAS_MT(T_loc, rho_a, us)
      
      
C Might not need to include any of these but PARAMS.h? 10/23/25
C d_gas might be in params file and doesn't need to be passed 10/23/25
C #include "DOUBLEPLUME_OPTIONS.h"
C #include "SIZE.h"
C #include "EEPARAMS.h"
C #include "PARAMS.h"
C #include "GRID.h"
C #include "DYNVARS.h"
C #include "DOUBLEPLUME.h"
      
      DOUBLE PRECISION mu_w, rho_a, T_loc, us
      DOUBLE PRECISION B, dE, D, Sc, Pe, Re, Sh
C      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

C      Ru = 8.314472 ! Universal gas constant in J/(mol K)
C      K = 1.0                                             ! Solubility of bubbles in [--]
      mu_w = 1721.2 * EXP(-0.0251 * T_loc)/1e6                              ! Pa s
        
C Calcuate the diffusion coefficient for each gas in m^2/s from Wise
C and Houghton (1966). These are for dilute, binary systems. 
      B = 4.2 * 1e-2        ! [cm^2/s]
      dE = 4390/0.23886     ! [J/mol]
      D = B * EXP(-dE/(8.314472 * (T_loc + 273.15)))/100.D0**2                      ! [m^2/s]
        
C Calculate the Schmidt number 
        Sc = mu_w/(rho_a * D)
    
C Calculate the Peclet number
        Pe = us * (0.01)/D
        
C Calculate the Reynolds number 
        Re = rho_a * us * (0.01) / mu_w
        
C Calculate the Sherwood number for a rigid particle
        Sh = 0 
        IF (Re .LT. 100) THEN
            Sh = 1 + (1 + 1/Pe)**(1.0D0/3.0D0)*Re**0.41
     &          *Sc**(1.0D0/3.0D0)
        ELSEIF (Re .LT. 2000) THEN
            Sh = 1 + 0.724*Re**0.48*Sc**(1.0D0/3.0D0)
        ELSE
            Sh = 1 + 0.425*Re**0.55*Sc**(1.0D0/3.0D0)
        ENDIF

        GAS_MT = 1.D0 * Sh * D/(0.01)                                 ! [m/s]

        RETURN 
        END

C---------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION GAS_HT(T_loc, rho_a, us)
      
C Might not need to include any of these but PARAMS.h? 10/23/25
C d_gas might be in params file and doesn't need to be passed 10/23/25
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      DOUBLE PRECISION mu_w, k_loc, Sc, Pe, Re, Sh
      DOUBLE PRECISION T_loc, rho_a, us
C      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0


C Calculate the heat transfer coeff. of each droplet in a liquid mixture
C Enter the physical constant for seawater 
      mu_w = 1721.2 * EXP(-0.0251*T_loc)/1e6 ! Pa s

C Calculate the thermal conductivity of seawater
      k_loc = 1.46e-7 ! m^2/s

C Calculate the Schmidt number 
      Sc = mu_w/(rho_a * k_loc); 

C Calculate the Peclet number
      Pe = us * 0.01 / k_loc 

C Calculate the Reynolds number
      Re = rho_a * us * (0.01)/mu_w 

      IF (Re .LT. 100) THEN
        Sh = 1 + (1 + 1/Pe)**(1.0D0/3.0D0)*Re**0.41*Sc**(1.0D0/3.0D0) 
      ELSEIF (Re .LT. 2000) THEN
        Sh = 1 + 0.724*Re**0.48*Sc**(1.0D0/3.0D0) 
      ELSE
        Sh = 1 + 0.425*Re**0.55*Sc**(1.0D0/3.0D0) 
      ENDIF

C Use the Sherwood number to obtain the mass transfer coeff. and store the
C value in beta
C beta = gas.K_t.*Sh(1).*k./gas.de(1); %m/s

      GAS_HT = 1.0*Sh*k_loc/(0.01) ! m/s Is it K or K_t? (9/11/2025)


      RETURN
      END

C---------------------------------------------------------
      
      DOUBLE PRECISION FUNCTION GAS_SOL(rho_a, T_loc, S_loc, P, mb)
      
C Might not need to include any of these but PARAMS.h? 10/23/25
C d_gas might be in params file and doesn't need to be passed 10/23/25
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      EXTERNAL rroots_chebyshev_cubic 

      DOUBLE PRECISION rho_a, T_loc, Ti_k, S_loc, P, mb
      DOUBLE PRECISION H0, nu_bar, dH_solR, T_0, moles, H 
      DOUBLE PRECISION b, ac, mu, alpha, a, AD, BD, coeffs(4)
      DOUBLE PRECISION zr(3), zi(3), Z, lnphi, phi, fug
      DOUBLE PRECISION GAS_HT
C      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846D0

C      M_O = 2*15.9994/1000! molar weight kg/mol
C      dhSol_gas = 428733
C      Ru = 8.314472 ! Universal gas constant in J/(mol K)
C      Tcr_O = 154.58
C      Pcr_O = 5042825.270
C      omega_g = 0.0216

C Convert pressure and Temp to absolute
      P = P
      Ti_k = T_loc + 273.15

C Henry's law coef
      H0 = 1.3 * 1e-3*1000/101325 * M_O ![kg/Pa m^3]

C Partial molar volume at infinite dilution
      nu_bar = 32/1e6              ! m^3/mol              

C Heat of Solution [K]
      dH_solR = dhSol_gas / Ru * M_O

C Compute the temperature effect
      T_0 = 298.15 ! May need to change?
      H0 = H0 * EXP(dH_solR * (1/Ti_k - 1/T_0))

C Compute the salting-out effect. These numbers are for CO2 in water from Weiss 1974
C H0 = H0 * np.exp(S * (0.027766 - 0.025888 * T/100 + 0.0050578 * (T/100)**2)); 
      moles = mb/M_O                                          ! Molar flux [mol/s]
C Pk = P * moles                                          

C Compute the pressure effect
      H = H0 * EXP(P * nu_bar/(Ru * Ti_k))                             ! [kg/Pa m^3]

C Compute Peng-Robinson mixture coeffs.
      b = 0.07780 * Ru * Tcr_O/Pcr_O 
      ac = 0.45724 * Ru**2 * Tcr_O**2/Pcr_O 
      mu = 0.37464 + 1.54226 * omega_g - 0.26992 * omega_g**2 
      alpha = (1 + mu * (1 - SQRT(Ti_k/Tcr_O)))**2 
      a = alpha * ac

      AD = a * P / (Ru**2 * Ti_k**2)
      BD = b * P / (Ru * Ti_k)

C Compressibility
C coefficients for cubic in Z: Z^3 + c2 Z^2 + c1 Z + c0 = 0
      coeffs(1) = 1
      coeffs(2) = -(1.0 - BD)
      coeffs(3) = (AD - 3.0*BD**2 - 2.0*BD)
      coeffs(4) = -(AD*BD - BD**2 - BD**3)

      CALL rroots_chebyshev_cubic(coeffs, zr, zi)
        
C        if SIZE(zr) == 0:
C            raise RuntimeError("No real Z roots found for PR EOS.")
        
      Z = MAXVAL(zr)          ! choose vapor root

      lnphi = (Z - 1.0) - LOG(Z - BD) - (AD / (2.0*SQRT(2.0D0)*BD)) * 
     &          LOG((Z + (1.0 + SQRT(2.0D0))*BD) 
     &           / (Z + (1.0 - SQRT(2.0D0))*BD))

      phi = EXP(lnphi)                                                 ! Fugacity coeff. []
      fug = phi * P                                                    ! [Pa]

C Solubility [kg/m^3]
      GAS_SOL = fug * H

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE rroots_chebyshev_cubic(coeffs, zr, zi)

      DOUBLE PRECISION, DIMENSION(4), INTENT(IN)  :: coeffs
      DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: zr, zi

      INTEGER :: l
      DOUBLE PRECISION :: a, b, c, d, p, t1, t2, t3, t4, t, x1, x2, x3
      DOUBLE PRECISION, PARAMETER :: sqrt3 = SQRT(3.0D0)
      DOUBLE PRECISION, PARAMETER :: s = 1.0D0 / 3.0D0
      DOUBLE PRECISION, PARAMETER :: small = 1.0D-16   ! or 10.0D0**INT(LOG(EPSILON(1.0D0)))

      ! coefficients
      a = coeffs(1)
      b = coeffs(2)
      c = coeffs(3)
      d = coeffs(4)

      main: BLOCK
      t = sqrt3
      t2 = b*b
      t3 = 3.0D0*a
      t4 = t3*c
      p = t2 - t4
      x3 = ABS(p)
      x3 = SQRT(x3)
      x1 = b*(t4 - p - p) - 3.0D0*t3*t3*d
      x2 = ABS(x1)
      x2 = x2**s
      t2 = 1.0D0/t3
      t3 = b*t2

      IF (x3 > small*x2) THEN
            t1 = 0.5D0*x1/(p*x3)
            x2 = ABS(t1)
            t2 = x3*t2
            t = t*t2
            t4 = x2*x2
            IF (p < 0.0D0) THEN
            p = x2 + SQRT(t4 + 1.0D0)
            p = p**s
            t4 = -1.0D0/p
            IF (t1 >= 0.0D0) t2 = -t2
            x1 = (p + t4)*t2
            x2 = -0.5D0*x1
            x3 = 0.5D0*t*(p - t4)
            l = 2
            EXIT main
            ELSE
            x3 = ABS(1.0D0 - t4)
            x3 = SQRT(x3)
            IF (t4 > 1.0D0) THEN
            p = (x2 + x3)**s
            t4 = 1.0D0/p
            IF (t1 < 0.0D0) t2 = -t2
            x1 = (p + t4)*t2
            x2 = -0.5D0*x1
            x3 = 0.5D0*t*(p - t4)
            l = 2
            EXIT main
            ELSE
            t4 = ATAN2(x3, t1)*s
            x3 = COS(t4)
            t4 = SQRT(1.0D0 - x3*x3)*t
            x3 = x3*t2
            x1 = x3 + x3
            x2 = t4 - x3
            x3 = -(t4 + x3)
            IF (x2 <= x3) THEN
                  t2 = x2
                  x2 = x3
                  x3 = t2
            ENDIF
            ENDIF
            ENDIF
      ELSE
            IF (x1 < 0.0D0) x2 = -x2
            x1 = x2*t2
            x2 = -0.5D0*x1
            x3 = -t*x2
            IF (ABS(x3) > small) THEN
            l = 2
            EXIT main
            ENDIF
            x3 = x2
      ENDIF
      l = 0
      IF (x1 <= x2) THEN
            t2 = x1
            x1 = x2
            x2 = t2
            IF (t2 <= x3) THEN
            x2 = x3
            x3 = t2
            ENDIF
      ENDIF
      x3 = x3 - t3
      END BLOCK main

      x1 = x1 - t3
      x2 = x2 - t3

      SELECT CASE (l)
      CASE (0)
      zr = [x1, x2, x3]
      zi = 0.0D0
      CASE (2)
      zr = [x1, x2, x2]
      zi = [0.0D0, x3, -x3]
      END SELECT

      RETURN
      END SUBROUTINE rroots_chebyshev_cubic

C-------------------------------------------------------------
      subroutine linint(nx,xtab,ytab,x,y)
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