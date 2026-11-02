#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C !ROUTINE: DOUBLEPLUME_READPARMS

      SUBROUTINE DOUBLEPLUME_READPARMS( myThid )

C !DESCRIPTION:
C     Initialize DOUBLEPLUME parameters, read in data.DOUBLEPLUME.
C     Namelists (in this order in data.DOUBLEPLUME):
C       DOUBLEPLUME_PARM01 : files, switches, background melt
C       DOUBLEPLUME_PARM02 : icefront thermodynamics
C       DOUBLEPLUME_PARM03 : ptracers (only with ALLOW_PTRACERS)
C       DOUBLEPLUME_PHYS   : plume physics (optional; defaults are
C                            the values of const.py/gas.py/params.py)

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !LOCAL VARIABLES: ====================================================
      INTEGER iUnit, ios
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

#ifdef ALLOW_DOUBLEPLUME

      NAMELIST /DOUBLEPLUME_PARM01/
     &     runoffVelfile, runoffRadfile, plumeMaskFile,
     &     iceTemp, c_i, backgroundVel,
     &     usePlumeDiagnostics,
     &     conserveMass,
     &     useSheetPlume,
     &     useConePlume

      NAMELIST /DOUBLEPLUME_PARM02/
     &     ICEFRONTlatentHeat, ICEFRONTHeatCapacity_Cp,
     &     applyIcefrontTendT,applyIcefrontTendS

#ifdef ALLOW_PTRACERS
      NAMELIST /DOUBLEPLUME_PARM03/
     &     ptracerMaskFile,
     &     useInputPtracers
#endif /* ALLOW_PTRACERS */

      NAMELIST /DOUBLEPLUME_PHYS/
     &  DP_g, DP_Patm, DP_T0, DP_cp, DP_cpgFrac, DP_Li,
     &  DP_lam1, DP_lam2, DP_lam3,
     &  DP_Ru, DP_Mgas, DP_PcrGas, DP_TcrGas, DP_omegaGas, DP_dhSol,
     &  DP_de0, DP_Kgas, DP_lambdaR, DP_fdis,
     &  DP_rhoR, DP_C0, DP_S0, DP_R, DP_mb0, DP_u0,
     &  DP_Cd, DP_gammaT, DP_gammaS,
     &  DP_alpha1, DP_alpha2, DP_alpha3, DP_lambdaS,
     &  DP_epsilon, DP_c1, DP_c2, DP_gammaI, DP_gammaO,
     &  DP_kappaPeel, DP_fe, DP_fp, DP_nwidths,
     &  DP_dz, DP_FrOuter, DP_FrWuest, DP_uStall, DP_JoMin,
     &  DP_Ca, DP_usConst,
     &  DP_rtolInner, DP_atolInner, DP_hmaxInner,
     &  DP_rtolOuter, DP_atolOuter, DP_hmaxOuter,
     &  DP_surfTol, DP_bottomTol, DP_peelTol, DP_fStall,
     &  DP_maxIter, DP_convTol, DP_rampTime, DP_nUpdate

      _BEGIN_MASTER(myThid)

C--   Defaults: PARM01
      runoffVelfile = ' '
      runoffRadfile = ' '
      plumeMaskFile = ' '
      iceTemp       = 0. _d 0
      c_i           = 2009. _d 0
      backgroundVel = 0.001 _d 0
      usePlumeDiagnostics = .FALSE.
      conserveMass  = .FALSE.
      useSheetPlume = .FALSE.
      useConePlume  = .FALSE.

C--   Defaults: PARM02
      ICEFRONTlatentHeat       =  334.0 _d 3
      ICEFRONTHeatCapacity_Cp  = 2000.0 _d 0
      applyIcefrontTendT       = .TRUE.
      applyIcefrontTendS       = .TRUE.

#ifdef ALLOW_PTRACERS
C--   Defaults: PARM03
      ptracerMaskFile  = ' '
      useInputPtracers = .FALSE.
#endif /* ALLOW_PTRACERS */

C--   Defaults: DOUBLEPLUME_PHYS -- model/const.py
      DP_g        = 9.81 _d 0
      DP_Patm     = 101325. _d 0
      DP_T0       = 273.15 _d 0
      DP_cp       = 4000. _d 0
      DP_cpgFrac  = 0.53 _d 0
      DP_Li       = 334000. _d 0
      DP_lam1     = -5.73 _d -2
      DP_lam2     = 0.08 _d 0
      DP_lam3     = -7.61 _d -4
      DP_Ru       = 8.314472 _d 0
      DP_Mgas     = 2. _d 0*15.9994 _d 0/1000. _d 0
      DP_PcrGas   = 5042825.270 _d 0
      DP_TcrGas   = 154.58 _d 0
      DP_omegaGas = 0.0216 _d 0
      DP_dhSol    = 428733. _d 0
C--   model/gas.py
      DP_de0      = 10.0 _d 0
      DP_Kgas     = 1.0 _d 0
      DP_lambdaR  = 0.80 _d 0
      DP_fdis     = 0.001 _d 0
C--   params.py
      DP_rhoR     = 1031. _d 0
      DP_C0       = 0.005 _d 0
      DP_S0       = 0. _d 0
      DP_R        = 0.05 _d 0
      DP_mb0      = 0.002 _d 0
      DP_u0       = 0.11 _d 0
      DP_Cd       = 0.001 _d 0
      DP_gammaT   = 2.20 _d -2
      DP_gammaS   = 6.2 _d -4
      DP_alpha1   = 0.055 _d 0
      DP_alpha2   = 0.110 _d 0
      DP_alpha3   = 0.110 _d 0
      DP_lambdaS  = 1.10 _d 0
      DP_epsilon  = 0.0175 _d 0
      DP_c1       = 0. _d 0
      DP_c2       = 1. _d 0
      DP_gammaI   = 1.10 _d 0
      DP_gammaO   = 1.10 _d 0
      DP_kappaPeel= 2. _d 0
      DP_fe       = 0.1 _d 0
      DP_fp       = 0.95 _d 0
      DP_nwidths  = 1
C--   values hard-coded in doublePlume.py, and numerics
      DP_dz        = 0.1 _d 0
      DP_FrOuter   = 0.1 _d 0
      DP_FrWuest   = 0.3 _d 0
      DP_uStall    = 1. _d -3
      DP_JoMin     = 5. _d -4
      DP_Ca        = 0. _d 0
      DP_usConst   = 0. _d 0
      DP_rtolInner = 1. _d -8
      DP_atolInner = 1. _d -10
      DP_hmaxInner = 0.1 _d 0
      DP_rtolOuter = 1. _d -6
      DP_atolOuter = 1. _d -9
      DP_hmaxOuter = 0.1 _d 0
      DP_surfTol   = 0. _d 0
      DP_bottomTol = 1.0 _d 0
      DP_peelTol   = 0.01 _d 0
      DP_fStall    = 0.1 _d 0
      DP_maxIter   = 10
      DP_convTol   = 0.01 _d 0
      DP_rampTime  = 30. _d 0*86400. _d 0
      DP_nUpdate   = 1

      WRITE(msgBuf,'(A)')
     &  'DOUBLEPLUME_READPARMS: opening data.DOUBLEPLUME'
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
      READ(UNIT=iUnit,NML=DOUBLEPLUME_PHYS,IOSTAT=ios)
      IF ( ios .GT. 0 ) THEN
        WRITE(msgBuf,'(A,I6)')
     &   'DOUBLEPLUME_READPARMS: error reading DOUBLEPLUME_PHYS,'
     &   //' iostat =', ios
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DOUBLEPLUME_READPARMS'
      ELSEIF ( ios .LT. 0 ) THEN
        WRITE(msgBuf,'(A)')
     &   'DOUBLEPLUME_READPARMS: no DOUBLEPLUME_PHYS group,'
     &   //' using defaults'
        CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                     SQUEEZE_RIGHT, myThid )
      ENDIF
      WRITE(msgBuf,'(A)')
     &  ' DOUBLEPLUME_READPARMS: finished reading data.DOUBLEPLUME'
      CALL PRINT_MESSAGE(msgBuf, standardMessageUnit,
     &                   SQUEEZE_RIGHT, myThid )
#ifdef SINGLE_DISK_IO
      CLOSE(iUnit)
#else
      CLOSE(iUnit,STATUS='DELETE')
#endif /* SINGLE_DISK_IO */

C--   Derived quantities and checks
      DP_pi  = 4. _d 0*ATAN(1. _d 0)
      DP_cpg = DP_cpgFrac*DP_cp
      IF ( DP_lam3 .GT. 0. _d 0 ) THEN
        WRITE(msgBuf,'(A)')
     &   'DOUBLEPLUME_READPARMS: DP_lam3 must be < 0 (depth > 0)'
        CALL PRINT_ERROR( msgBuf, myThid )
        STOP 'ABNORMAL END: S/R DOUBLEPLUME_READPARMS'
      ENDIF
      IF ( DP_nwidths .LT. 1 ) DP_nwidths = 1
      IF ( DP_nUpdate .LT. 1 ) DP_nUpdate = 1
      IF ( DP_maxIter .LT. 0 ) DP_maxIter = 0

      _END_MASTER(myThid)
C--   Everyone else must wait for the parameters to be loaded
      _BARRIER

#endif /* ALLOW_DOUBLEPLUME */
      RETURN
      END
