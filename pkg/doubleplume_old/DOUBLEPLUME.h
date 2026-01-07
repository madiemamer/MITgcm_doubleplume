
#include "DOUBLEPLUME_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#endif /* ALLOW_PTRACERS */

#ifdef ALLOW_DOUBLEPLUME

C Header file pkg/DOUBLEPLUME

      COMMON /DOUBLEPLUME_INTEGERS/
     &       iceDepthK,
     &       runoffLocation, plumeInputLayer

      INTEGER iceDepthK
      INTEGER runoffLocation  (100,2)
      INTEGER plumeInputLayer (100)

      COMMON /DOUBLEPLUME_LOGICALS/
     &      usePlumeDiagnostics,
     &      dispersePlumeOutput,
     &      useSheetPlume,
     &      useConePlume,
     &      useDetachPlume,
     &      conserveMass,
     &      useInputPtracers
      LOGICAL usePlumeDiagnostics
      LOGICAL dispersePlumeOutput
      LOGICAL useSheetPlume
      LOGICAL useConePlume
      LOGICAL useDetachPlume
      LOGICAL conserveMass
      LOGICAL useInputPtracers

      COMMON /DOUBLEPLUME_FIELDS/
     &     runoffVel, runoffVel0, runoffVel1,
     &     runoffRad, runoffRad0, runoffRad1,
     &     plumeMask, 
     &     HeatFlux3D, FwFlux3D,
     &     QiProf3D, MiProf3D,
     &     HiProf3D, BiProf3D,
     &     mProf3D, mProfAv3D,
     &     QiProf3DLocal, MiProf3DLocal,
     &     HiProf3DLocal, BiProf3DLocal,
     &     mProf3DLocal, mProfAv3DLocal,
     &     temp_addMass3D, salt_addMass3D,
     &     volFluxDiff3D,
     &     Qin
      _RL runoffVel  (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL runoffVel0 (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL runoffVel1 (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL runoffRad  (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL runoffRad0 (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL runoffRad1 (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RL plumeMask  (1-Olx:sNx+Olx,1-Oly:sNy+Oly,nSx,nSy)
      _RS HeatFlux3D  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS FwFlux3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS QiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS MiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS HiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS BiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS mProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS mProfAv3D       (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL temp_addMass3D  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL salt_addMass3D  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL volFluxDiff3D   (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS QiProf3DLocal    (sNx,sNy,Nr)
      _RS MiProf3DLocal    (sNx,sNy,Nr)
      _RS HiProf3DLocal    (sNx,sNy,Nr)
      _RS BiProf3DLocal    (sNx,sNy,Nr)
      _RS mProf3DLocal    (sNx,sNy,Nr)
      _RS mProfAv3DLocal       (sNx,sNy,Nr)
      _RS Qin (100)

      COMMON /DOUBLEPLUME_FIELDS_PROFILES/
     &     saProf, taProf, ptProf, prProf, uProf, vProf,
     &     QiProf, MiProf, HiProf, BiProf,
     &     uProfPlume, mProfPlume, mProfAv, mProf, zProf, zProfAbs,
     &     volFLux, volFluxDiff,
     &     FwFlux, HeatFlux,
     &     delta_z, aProf, mIntProf
      _RL saProf  (Nr)
      _RL taProf  (Nr)
      _RL ptProf (Nr)
      _RL prProf (Nr)
      _RL uProf  (Nr)
      _RL vProf  (Nr)
      _RL QiProf (Nr+1)
      _RL MiProf (Nr+1)
      _RL HiProf (Nr+1)
      _RL BiProf (Nr+1)
      _RL aProf (Nr+1)   ! integrated area of plume from bottom
      _RL mIntProf (Nr+1)   ! integrated melt of plume from bottom
      _RL uProfPlume (Nr)     
      _RL mProfPlume (Nr)
      _RL mProfAv    (Nr)
      _RL mProf      (Nr)
      _RL zProf      (Nr+1)
      _RL zProfAbs   (Nr+1)
      _RL volFlux     (Nr+1)
      _RL volFluxDiff (Nr)
      _RL HeatFlux(Nr)
      _RL FwFlux(Nr)
      _RL delta_z(Nr)
      
#ifdef DOUBLEPLUME_ALLOW_DETACHED_PLUME      
      COMMON /DOUBLEPLUME_FIELDS_DETACHED/
     & thetaProf, distanceProf,
     & thetaProf3D, distanceProf3D,
     & thetaProf3DLocal, distanceProf3DLocal
      _RL thetaProf (Nr+1)
      _RL distanceProf (Nr+1)
      _RL thetaProf3D
     &  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL distanceProf3D
     &  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)      
      _RL thetaProf3DLocal (sNx,sNy,Nr)
      _RL distanceProf3DLocal (sNx,sNy,Nr)
#endif      

      COMMON /DOUBLEPLUME_PARM_R/
     &     E_0,
     &     Q_sg, T_sg, S_sg, r_sg, w_sg,
     &     theta_sg,
     &     delta_x, delta_y,
     &     RTOL, ATOL,
     &     iceTemp,
     &     outputThickness,
     &     kap, rho_ref, g, c_w, c_i, L,
     &     gam, lambda1, lambda2, lambda3,
     &     GamT, GamS, Cd,
     &     sOutPlume, tOutPlume,
     &     iceDepth,
     &     backgroundVel,
     &     ptracerfluxSum,
     &     maxDepth
      _RS E_0
      _RS Q_sg
      _RL T_sg
      _RL S_sg
      _RL w_sg
      _RL r_sg
      _RL theta_sg
      _RL delta_x
      _RL delta_y
      _RL RTOL
      _RL ATOL
      _RL iceTemp
      _RL outputThickness
      _RL kap
      _RL rho_ref
      _RL g
      _RL c_w
      _RL c_i
      _RL L
      _RL gam
      _RL lambda1
      _RL lambda2
      _RL lambda3
      _RL sOutPlume
      _RL tOutPlume
      _RL gamS
      _RL gamT
      _RL cd
      _RL iceDepth
      _RL backgroundVel
      _RL ptracerFluxSum
      _RL maxDepth

      COMMON /DOUBLEPLUME_FILES/
     &	    runoffVelfile,
     &	    runoffRadfile,
     &      plumeMaskFile
C Again, doesn't like MAX_LEN_FNAM
      CHARACTER*(512)
     &	    runoffVelfile,
     &	    runoffRadfile,
     &      plumeMaskFile

C-----------------------------------------
C Parameters relating to the icefront package

      COMMON /DOUBLEPLUME_ICEFRONT_PARMS_R/
     &     ICEFRONTlatentHeat, 
     &     ICEFRONTheatCapacity_Cp
      _RL ICEFRONTlatentHeat
      _RL ICEFRONTheatCapacity_Cp

      COMMON /DOUBLEPLUME_ICEFRONT_FIELDS_RL/
     &     icefront_TendT,
     &     icefront_TendS
      _RL icefront_TendT (1:sNx,1:sNy,Nr,nSx,nSy)
      _RL icefront_TendS (1:sNx,1:sNy,Nr,nSx,nSy)

      COMMON /DOUBLEPLUME_ICEFRONT_FIELDS_RS/
     &     R_icefront,
     &     icefrontlength
      _RS R_icefront     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS icefrontlength (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C      LOGICAL ICEFRONTisOn
      LOGICAL applyIcefrontTendT
      LOGICAL applyIcefrontTendS
      COMMON /DOUBLEPLUME_ICEFRONT_PARMS_L/
C     &     ICEFRONTisOn,
     &     applyIcefrontTendT,
     &     applyIcefrontTendS

C---------------------------------------------

#ifdef ALLOW_PTRACERS
C Parameters relating to PTRACERS

      COMMON /DOUBLEPLUME_PTRACERS_RL/
     &     ptr_addMass3D,
     &     ptracerMask
      _RL ptr_addMass3d
     &     (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy,PTRACERS_num)
      _RL ptracerMask
     &     (1-Olx:sNx+Olx,1-Oly:sNy+Oly,PTRACERS_num,nSx,nSy)
   
      COMMON /DOUBLEPLUME_PTRACERS_FILES/
     &     ptracerMaskFile
      CHARACTER*(512) 
     &     ptracerMaskFile

#endif /* ALLOW_PTRACERS */

# endif /* ALLOW_DOUBLEPLUME */ 
