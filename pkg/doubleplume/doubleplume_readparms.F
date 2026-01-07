
#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C !ROUTINE: DOUBLEPLUME_READPARMS

      SUBROUTINE DOUBLEPLUME_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize DOUBLEPLUME parameters, read in data.DOUBLEPLUME

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
#include "PARAMS.h"
#include "GRID.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
C  iUnit                :: unit number for I/O
C  msgBuf               :: message buffer
C  bi,bj,i,j,k		:: loop indices

      INTEGER iUnit
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER bi, bj, j, i, k

CEOP
#ifdef ALLOW_DOUBLEPLUME

C    PARMO1 - DOUBLEPLUME inputs
      NAMELIST /DOUBLEPLUME_PARM01/
     &     runoffVelfile, runoffRadfile, plumeMaskFile,
     &     T_sg, S_sg, E_0, w_sg, r_sg, theta_sg,
     &     RTOL, ATOL,
     &     IceTemp,
     &     usePlumeDiagnostics,
     &     conserveMass,
     &     rho_r,
     &     g,
     &     cp,
     &     c_i,
     &     L,
     &     lam1,
     &     lam2,
     &     lam3,
     &     gammaT,
     &     gammaS,
     &     Cd,
     &     useSheetPlume,
     &     useConePlume,
     &     maxDepth,
     &     backgroundVel,
     &     c1, c2, epsilon,
     &     alpha_1, alpha_2,
     &     gamma_i,
     &     Fr0, 
     &     gas_lam1, gas_lam2, d_gas,
     &     fp_g,
     &     cp_g,
     &     T_g,
     &     Ru,
     &     Tcr_O, Pcr_O, M_O, dhSol_gas,
     &     omega_g, K_g, plumeLayers

C    PARMO2 - icefront inputs required for thermodynamics
      NAMELIST /DOUBLEPLUME_PARM02/
     &     ICEFRONTlatentHeat, ICEFRONTHeatCapacity_Cp,
     &     applyIcefrontTendT,applyIcefrontTendS

#ifdef ALLOW_PTRACERS
C    PARM03 - inputs relating to passive tracers
      NAMELIST /DOUBLEPLUME_PARM03/
     &     ptracerMaskFile,
     &     useInputPtracers
#endif /* ALLOW_PTRACERS */ 

      _BEGIN_MASTER(myThid)

C     Default values for PARM01
      E_0      = 0.1D0
      T_sg     = 1.0D-3
      S_sg     = 1.0D-3
      w_sg     = 1.0D-2
      r_sg     = 1.0D-2
      theta_sg = PI/2.
      RTOL     = 1.0D-5
      ATOL     = 1.0D-5
      iceTemp  = 0.
      runoffVelfile = ' '
      runoffRadfile = ' '
      plumeMaskFile = ' '
      usePlumeDiagnostics = .FALSE.
      conserveMass    = .FALSE.
      useSheetPlume   = .FALSE.
      useConePlume    = .FALSE.
      rho_r  = 1025.D0
      g      = 9.81D0
      cp     = 3994.D0
      c_i    = 2009.D0
      L      = 334000.D0
      lam1   = -0.0573D0
      lam2   = 0.0832D0
      lam3   = 0.000761D0
      gammaT = 0.022
      gammaS = 0.00062
      Cd     = 0.001
      backgroundVel = 0.001
      plumeLayers = 4901

C Default values for PARMO2
      ICEFRONTlatentHeat       =  334.0 _d 3
      ICEFRONTHeatCapacity_Cp  = 2000.0 _d 0
      applyIcefrontTendT       = .TRUE.
      applyIcefrontTendS       = .TRUE.

C Default values for PARM03
      c1 = 0.0
      c2 = 1.0
      epsilon   = 0.0175
      alpha_1   = 0.055
      alpha_2   = 0.110
      gamma_i   = 1.10
      Fr0       = 1.6
      gas_lam1  = 0.8
      gas_lam2  = 1.10
      d_gas     = 0.01
      fp_g        = 0.95
      cp_g      = 0.53 * cp
      T_g       = T_sg
      Ru        = 8.314472
      Tcr_O     = 154.58
      Pcr_O     = 5042825.270
      omega_g   = 0.0216
      M_O       = 2*15.9994/1000
      K_g         = 1.0
      dhSol_gas = 428733.0

#ifdef ALLOW_PTRACERS
C Default values for PARM03
       ptracerMaskFile= ' '
       useInputPtracers = .FALSE.
#endif /* ALLOW_PTRACERS */

C Open and read the data.icefront file
      WRITE(msgBuf,'(A)') 'DOUBLEPLUME_READPARMS: 
     &.         opening data.DOUBLEPLUME'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT, myThid )
      CALL OPEN_COPY_DATA_FILE(
     I                   'data.DOUBLEPLUME', 'DOUBLEPLUME_READPARMS',
     O                   iUnit,
     I                   myThid )
      READ(UNIT=iUnit,NML=DOUBLEPLUME_PARM01)
      READ(UNIT=iUnit,NML=DOUBLEPLUME_PARM02)
#ifdef ALLOW_PTRACERS
      READ(UNIT=iUnit,NML=DOUBLEPLUME_PARM03)
#endif /* ALLOW_PTRACERS */
      WRITE(msgBuf,'(A)')
     &  ' DOUBLEPLUME_READPARMS: finished reading data.DOUBLEPLUME'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT, myThid )

C Close the open data file
      CLOSE(iUnit)

C Zero values for T_sg and S_sg can cause the routine to crash
      IF ( T_sg.EQ.0.D0 ) T_sg = 1.0D-3
      IF ( S_sg.EQ.0.D0 ) S_sg = 1.0D-3

      _END_MASTER(myThid)
C--   Everyone else must wait for the parameters to be loaded
C     and tiled OB indices to be set.
      _BARRIER

#endif /* ALLOW_DOUBLEPLUME */
      RETURN
      END
