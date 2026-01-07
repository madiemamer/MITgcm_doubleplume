
#include "DOUBLEPLUME_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#endif /* ALLOW_PTRACERS */

#ifdef ALLOW_DOUBLEPLUME

C Header file pkg/DOUBLEPLUME

      COMMON /DOUBLEPLUME_INTEGERS/
     &       iceDepthK,
     &       runoffLocation, plumeInputLayer, j_peel

      INTEGER iceDepthK
      INTEGER runoffLocation  (100,2)
      INTEGER plumeInputLayer (100)
      INTEGER j_peel

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
     &     volFluxDiff3D, QoProf3DLocal,
     &     MoProf3DLocal, HoProf3DLocal,
     &     BoProf3DLocal, QoProf3D, MoProf3D,
     &     HoProf3D, BoProf3D,
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
      _RS HgProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS CiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS MbProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS HiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS BiProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS mProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS QoProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS MoProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS HoProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS BoProf3D    (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS mProfAv3D       (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL temp_addMass3D  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL salt_addMass3D  (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RL volFluxDiff3D   (1-OLx:sNx+OLx,1-Oly:sNy+Oly,Nr,nSx,nSy)
      _RS QiProf3DLocal    (sNx,sNy,Nr)
      _RS MiProf3DLocal    (sNx,sNy,Nr)
      _RS HiProf3DLocal    (sNx,sNy,Nr)
      _RS CiProf3DLocal    (sNx,sNy,Nr)
      _RS HgProf3DLocal    (sNx,sNy,Nr)
      _RS MbProf3DLocal    (sNx,sNy,Nr)
      _RS BiProf3DLocal    (sNx,sNy,Nr)
      _RS QoProf3DLocal    (sNx,sNy,Nr)
      _RS MoProf3DLocal    (sNx,sNy,Nr)
      _RS HoProf3DLocal    (sNx,sNy,Nr)
      _RS BoProf3DLocal    (sNx,sNy,Nr)
      _RS mProf3DLocal    (sNx,sNy,Nr)
      _RS mProfAv3DLocal       (sNx,sNy,Nr)
      _RS Qin (100)

      COMMON /DOUBLEPLUME_FIELDS_PROFILES/
     &     saProf, taProf, ptProf, prProf, uProf, vProf,
     &     QiProf, MiProf, HiProf, BiProf,
     &     uProfPlume, mProfPlume, mProfAv, mProf, zProf, zProfAbs,
     &     volFLux, volFluxDiff, zProfAbs_f,
     &     FwFlux, HeatFlux,
     &     delta_z, aProf, mIntProf, QoProf, MoProf,
     &     HoProf, BoProf, QiProf_f, MiProf_f, HiProf_f, HgProf_f,
     &     BiProf_f, MbProf_f, CiProf_f, QoProf_f, MoProf_f, HoProf_f,
     &     BoProf_f, CoProf_f, z_peel, iter_l
      _RL saProf  (Nr)
      _RL taProf  (Nr)
      _RL ptProf (Nr)
      _RL prProf (Nr)
      _RL uProf  (Nr)
      _RL vProf  (Nr)
      _RL QiProf_f (INT(ABS(490) / 0.1D0)+ 1)
      _RL MiProf_f (INT(ABS(490) / 0.1D0)+ 1)
      _RL HiProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL HgProf_f(INT(ABS(490) / 0.1D0) + 1)
      _RL BiProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL MbProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL CiProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL QiProf (Nr+1)
      _RL MiProf (Nr+1)
      _RL HiProf (Nr+1)
      _RL HgProf (Nr+1)
      _RL BiProf (Nr+1)
      _RL MbProf (Nr+1)
      _RL CiProf (Nr+1)
      _RL aProf  (Nr+1)   ! integrated area of plume from bottom
      _RL mIntProf (Nr+1)   ! integrated melt of plume from bottom
      _RL uProfPlume (Nr)     
      _RL mProfPlume (Nr)
      _RL mProfAv    (Nr)
      _RL mProf      (Nr)
      _RL zProf      (Nr+1)
      _RL zProfAbs   (Nr+1)
      _RL zProf_f      (INT(ABS(490) / 0.1D0) + 1)
      _RL zProfAbs_f   (INT(ABS(490) / 0.1D0) + 1)
      _RL volFlux     (Nr+1)
      _RL volFluxDiff (Nr+1)
      _RL HeatFlux(Nr)
      _RL FwFlux(Nr)
      _RL delta_z(Nr)
      _RL QoProf_f (INT(ABS(490) / 0.1D0)+ 1)
      _RL MoProf_f (INT(ABS(490) / 0.1D0)+ 1)
      _RL HoProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL BoProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL CoProf_f (INT(ABS(490) / 0.1D0) + 1)
      _RL QoProf (Nr+1)
      _RL MoProf (Nr+1)
      _RL HoProf (Nr+1)
      _RL BoProf (Nr+1)
      _RL CoProf (Nr+1)
      DOUBLE PRECISION iter_l
      DOUBLE PRECISION z_peel (Nr)
      
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
     &     kap, rho_r, g, cp, c_i, L,
     &     gam, lam1, lam2, lam3,
     &     gammaT, gammaS, Cd,
     &     sOutPlume, tOutPlume,
     &     iceDepth,
     &     backgroundVel,
     &     ptracerfluxSum,
     &     maxDepth, plumeLayers
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
      _RL rho_r
      _RL g
      _RL cp
      _RL c_i
      _RL L
      _RL gam
      _RL lam1
      _RL lam2
      _RL lam3
      _RL sOutPlume
      _RL tOutPlume
      _RL gammaS
      _RL gammaT
      _RL Cd
      _RL iceDepth
      _RL backgroundVel
      _RL ptracerFluxSum
      _RL maxDepth
      INTEGER plumeLayers

      COMMON /DOUBLEPLUME_FILES/
     &	    runoffVelfile,
     &	    runoffRadfile,
     &      plumeMaskFile
C Again, doesnt like MAX_LEN_FNAM
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
C-----------------------------------------
C Parameters relating to the Bubble DoublePlume package

      COMMON /DOUBLEPLUME_FILES/
     &     c1, c2, epsilon,
     &     alpha_1, alpha_2, gamma_i,
     &     Fr0, gas_lam1, gas_lam2, d_gas,
     &     fp_g, cp_g, T_g, Ru, Tcr_O, Pcr_O, omega_g,
     &     M_O, K_g, dhSol_gas
      _RL c1
      _RL c2
      _RL epsilon
      _RL alpha_1
      _RL alpha_2
      _RL gamma_i
      _RL Fr0
      _RL gas_lam1
      _RL gas_lam2
      _RL d_gas
      _RL fp_g
      _RL cp_g
      _RL T_g
      _RL Ru
      _RL Tcr_O
      _RL Pcr_O
      _RL omega_g
      _RL M_O
      _RL K_g 
      _RL dhSol_gas

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
