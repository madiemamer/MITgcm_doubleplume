
#include "DOUBLEPLUME_OPTIONS.h"

#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#endif /* ALLOW_PTRACERS */

#ifdef ALLOW_DOUBLEPLUME

C     *==========================================================*
C     | DOUBLEPLUME.h : header file of pkg/doubleplume
C     | Double (inner + outer) bubble plume model, a port of
C     | doublePlume.py. Physical constants and model parameters
C     | (DP_*) mirror model/const.py, model/gas.py and params.py
C     | and are read from namelist DOUBLEPLUME_PHYS in
C     | data.DOUBLEPLUME (see DOUBLEPLUME_READPARMS).
C     *==========================================================*

C--   Compile-time array limits
C     DP_maxFine :: max points on the fine plume grid; must be
C                   >= |grounding-line depth|/DP_dz + 1
C                   (10001 allows 1000 m at DP_dz = 0.1 m)
C     DP_maxPeel :: max number of discrete peels per plume
      INTEGER DP_maxFine, DP_maxPeel
      PARAMETER ( DP_maxFine = 10001, DP_maxPeel = 200 )

C--   Integers
C     iceDepthK   :: index of the model face at the plume source
C     plumeLayers :: number of fine-grid points in the current column
C     j_peel      :: number of discrete peels of the inner plume
      COMMON /DOUBLEPLUME_INTEGERS/
     &       iceDepthK, plumeLayers,
     &       runoffLocation, plumeInputLayer, j_peel
      INTEGER iceDepthK
      INTEGER plumeLayers
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
     &      useInputPtracers,
     &      dp_firstCall
      LOGICAL usePlumeDiagnostics
      LOGICAL dispersePlumeOutput
      LOGICAL useSheetPlume
      LOGICAL useConePlume
      LOGICAL useDetachPlume
      LOGICAL conserveMass
      LOGICAL useInputPtracers
      LOGICAL dp_firstCall

C--   2-D / 3-D fields
C     volFluxDiff3D :: volume exchanged with the ambient ocean by the
C                      plume in each cell, before the spin-up ramp
C                      [m^3/s], > 0 = plume water added
      COMMON /DOUBLEPLUME_FIELDS/
     &     runoffVel, runoffVel0, runoffVel1,
     &     runoffRad, runoffRad0, runoffRad1,
     &     plumeMask,
     &     HeatFlux3D, FwFlux3D,
     &     QiProf3D, MiProf3D, HiProf3D, BiProf3D,
     &     HgProf3D, MbProf3D, CiProf3D,
     &     QoProf3D, MoProf3D, HoProf3D, BoProf3D, CoProf3D,
     &     mProf3D, mProfAv3D,
     &     temp_addMass3D, salt_addMass3D,
     &     volFluxDiff3D,
     &     Qin
      _RL runoffVel  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL runoffVel0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL runoffVel1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL runoffRad  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL runoffRad0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL runoffRad1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL plumeMask  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL HeatFlux3D (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL FwFlux3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL QiProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL MiProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HiProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL BiProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HgProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL MbProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL CiProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL QoProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL MoProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL HoProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL BoProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL CoProf3D   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL mProf3D    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL mProfAv3D  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL temp_addMass3D (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL salt_addMass3D (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL volFluxDiff3D  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL Qin (100)

C--   Column profiles (model levels) and fine-grid plume profiles.
C     Fine-grid index k sits at height h = (k-1)*DP_dz above the
C     plume source (depth = |iceDepth| - h).
C     Inner plume: QiProf_f volume [m^3/s], MiProf_f momentum
C       [m^4/s^2], HiProf_f heat = rho_r*cp*(T+T0)*Q [J/s], HgProf_f
C       bubble heat [J/s], BiProf_f salt [psu m^3/s], MbProf_f
C       bubble mass [kg/s], CiProf_f dissolved gas [kg/s].
C     Outer plume: QoProf_f .. CoProf_f, same conventions (Qo < 0).
      COMMON /DOUBLEPLUME_FIELDS_PROFILES/
     &     saProf, taProf, ptProf, prProf, uProf, vProf,
     &     QiProf, MiProf, HiProf, HgProf, BiProf, MbProf, CiProf,
     &     QoProf, MoProf, HoProf, BoProf, CoProf,
     &     mProf, mProfAv, FwFlux, HeatFlux, delta_z,
     &     QiProf_f, MiProf_f, HiProf_f, HgProf_f,
     &     BiProf_f, MbProf_f, CiProf_f,
     &     QoProf_f, MoProf_f, HoProf_f, BoProf_f, CoProf_f,
     &     z_peel
      _RL saProf  (Nr)
      _RL taProf  (Nr)
      _RL ptProf  (Nr)
      _RL prProf  (Nr)
      _RL uProf   (Nr)
      _RL vProf   (Nr)
      _RL QiProf  (Nr+1)
      _RL MiProf  (Nr+1)
      _RL HiProf  (Nr+1)
      _RL HgProf  (Nr+1)
      _RL BiProf  (Nr+1)
      _RL MbProf  (Nr+1)
      _RL CiProf  (Nr+1)
      _RL QoProf  (Nr+1)
      _RL MoProf  (Nr+1)
      _RL HoProf  (Nr+1)
      _RL BoProf  (Nr+1)
      _RL CoProf  (Nr+1)
      _RL mProf   (Nr)
      _RL mProfAv (Nr)
      _RL FwFlux  (Nr)
      _RL HeatFlux(Nr)
      _RL delta_z (Nr)
      _RL QiProf_f (DP_maxFine)
      _RL MiProf_f (DP_maxFine)
      _RL HiProf_f (DP_maxFine)
      _RL HgProf_f (DP_maxFine)
      _RL BiProf_f (DP_maxFine)
      _RL MbProf_f (DP_maxFine)
      _RL CiProf_f (DP_maxFine)
      _RL QoProf_f (DP_maxFine)
      _RL MoProf_f (DP_maxFine)
      _RL HoProf_f (DP_maxFine)
      _RL BoProf_f (DP_maxFine)
      _RL CoProf_f (DP_maxFine)
C     z_peel :: DEPTHS [m, > 0] of the discrete peels
      _RL z_peel  (DP_maxPeel)

#ifdef DOUBLEPLUME_ALLOW_DETACHED_PLUME
      COMMON /DOUBLEPLUME_FIELDS_DETACHED/
     & thetaProf, distanceProf,
     & thetaProf3D, distanceProf3D
      _RL thetaProf (Nr+1)
      _RL distanceProf (Nr+1)
      _RL thetaProf3D
     &  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL distanceProf3D
     &  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#endif

C--   Column / source quantities (set in DOUBLEPLUME_CALC)
C     Q_sg, w_sg, r_sg :: subglacial discharge [m^3/s], velocity
C                         [m/s], radius [m]
C     iceTemp, c_i     :: ice temperature and heat capacity, used by
C                         the background melt (DOUBLEPLUME_MELTRATE)
C     backgroundVel    :: minimum velocity for background melt [m/s]
C     iceDepth         :: depth of the plume source (< 0)
      COMMON /DOUBLEPLUME_PARM_R/
     &     Q_sg, w_sg, r_sg,
     &     delta_x, delta_y,
     &     iceTemp, c_i, backgroundVel,
     &     iceDepth, maxDepth
      _RL Q_sg
      _RL w_sg
      _RL r_sg
      _RL delta_x
      _RL delta_y
      _RL iceTemp
      _RL c_i
      _RL backgroundVel
      _RL iceDepth
      _RL maxDepth

      COMMON /DOUBLEPLUME_FILES/
     &      runoffVelfile,
     &      runoffRadfile,
     &      plumeMaskFile
      CHARACTER*(512) runoffVelfile
      CHARACTER*(512) runoffRadfile
      CHARACTER*(512) plumeMaskFile

C-----------------------------------------
C     Parameters relating to the icefront package
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

      LOGICAL applyIcefrontTendT
      LOGICAL applyIcefrontTendS
      COMMON /DOUBLEPLUME_ICEFRONT_PARMS_L/
     &     applyIcefrontTendT,
     &     applyIcefrontTendS

C-----------------------------------------
C     Double plume constants / parameters (namelist DOUBLEPLUME_PHYS)
C
C  model/const.py
C     DP_pi       :: pi (derived)
C     DP_g        :: gravitational acceleration        [m/s^2]
C     DP_Patm     :: atmospheric pressure               [Pa]
C     DP_T0       :: 0 degC in Kelvin                   [K]
C     DP_cp       :: heat capacity of seawater          [J/kg/K]
C     DP_cpgFrac  :: gas heat capacity / DP_cp          [-]
C     DP_cpg      :: heat capacity of gas (derived)     [J/kg/K]
C     DP_Li       :: latent heat of ice                 [J/kg]
C     DP_lam1/2/3 :: liquidus Tb = lam1*S + lam2 + lam3*depth,
C                    depth > 0, so DP_lam3 < 0
C     DP_Ru       :: universal gas constant             [J/mol/K]
C     DP_Mgas     :: molar mass of the gas (O2)         [kg/mol]
C     DP_PcrGas, DP_TcrGas, DP_omegaGas :: critical pressure [Pa],
C                    critical temperature [K], acentric factor
C     DP_dhSol    :: heat of dissolution                [J/kg]
C  model/gas.py
C     DP_de0      :: initial bubble diameter            [mm]
C     DP_Kgas     :: bubble mass-transfer factor K      [-]
C     DP_lambdaR  :: bubble / velocity spreading ratio  [-]
C     DP_fdis     :: bubbles dissolved when mb/mb0 < fdis
C  params.py
C     DP_rhoR     :: reference density                  [kg/m^3]
C     DP_C0       :: initial O2 concentration           [kg/m^3]
C     DP_S0       :: salinity of the discharge          [psu]
C     DP_R        :: source radius if runoffRad = 0     [m]
C     DP_mb0      :: bubble mass flux at the source     [kg/s]
C     DP_u0       :: initial velocity guess             [m/s]
C     DP_Cd, DP_gammaT, DP_gammaS :: melt parameterisation
C     DP_alpha1/2/3 :: entrainment: inner, outer->inner,
C                    ambient->outer                    [-]
C     DP_lambdaS  :: conc./velocity spreading ratio     [-]
C     DP_epsilon  :: continuous peeling parameter       [-]
C     DP_c1, DP_c2:: recirculation / peeling constants  [-]
C     DP_gammaI/O :: momentum amplification inner/outer [-]
C     DP_kappaPeel:: discrete peel when u_i < u_s/kappa [-]
C     DP_fe       :: ambient fraction engulfed at surface[-]
C     DP_fp       :: inner flux fraction detrained at a peel
C     DP_nwidths  :: plume widths searched per CP window
C  Values hard-coded inside doublePlume.py, and numerics
C     DP_dz       :: fine plume grid spacing            [m]
C     DP_FrOuter  :: outer plume initial Froude no.     [-]
C     DP_FrWuest  :: inner (Wuest) initial Froude no.   [-]
C     DP_uStall   :: inner plume stall velocity         [m/s]
C     DP_JoMin    :: outer momentum stop threshold  [m^4/s^2]
C     DP_Ca       :: ambient dissolved gas          [kg/m^3]
C     DP_usConst  :: > 0: fixed bubble slip velocity    [m/s]
C                    <= 0: computed from the bubble state
C     DP_rtolInner, DP_atolInner, DP_hmaxInner :: inner DLSODE
C     DP_rtolOuter, DP_atolOuter, DP_hmaxOuter :: outer DLSODE
C     DP_surfTol  :: inner plume "reached the surface" tolerance [m]
C     DP_bottomTol:: outer plume stops this far above source [m]
C     DP_peelTol  :: "at a discrete peel" tolerance     [m]
C     DP_fStall   :: outer/inner flux ratio at the top when the
C                    inner plume stops below the surface
C     DP_maxIter  :: max inner/outer iterations (Python: 10)
C     DP_convTol  :: convergence tolerance (Python: 0.01)
C     DP_rampTime :: addMass spin-up ramp time (<= 0: none) [s]
C     DP_nUpdate  :: re-solve the plume every DP_nUpdate steps
      _RL DP_pi, DP_g, DP_Patm, DP_T0, DP_cp, DP_cpgFrac, DP_cpg
      _RL DP_Li, DP_lam1, DP_lam2, DP_lam3
      _RL DP_Ru, DP_Mgas, DP_PcrGas, DP_TcrGas, DP_omegaGas, DP_dhSol
      _RL DP_de0, DP_Kgas, DP_lambdaR, DP_fdis
      _RL DP_rhoR, DP_C0, DP_S0, DP_R, DP_mb0, DP_u0
      _RL DP_Cd, DP_gammaT, DP_gammaS
      _RL DP_alpha1, DP_alpha2, DP_alpha3, DP_lambdaS
      _RL DP_epsilon, DP_c1, DP_c2, DP_gammaI, DP_gammaO
      _RL DP_kappaPeel, DP_fe, DP_fp
      _RL DP_dz, DP_FrOuter, DP_FrWuest, DP_uStall, DP_JoMin
      _RL DP_Ca, DP_usConst
      _RL DP_rtolInner, DP_atolInner, DP_hmaxInner
      _RL DP_rtolOuter, DP_atolOuter, DP_hmaxOuter
      _RL DP_surfTol, DP_bottomTol, DP_peelTol, DP_fStall
      _RL DP_convTol, DP_rampTime
      INTEGER DP_nwidths, DP_maxIter, DP_nUpdate
      COMMON /DOUBLEPLUME_PHYS_R/
     &  DP_pi, DP_g, DP_Patm, DP_T0, DP_cp, DP_cpgFrac, DP_cpg,
     &  DP_Li, DP_lam1, DP_lam2, DP_lam3,
     &  DP_Ru, DP_Mgas, DP_PcrGas, DP_TcrGas, DP_omegaGas, DP_dhSol,
     &  DP_de0, DP_Kgas, DP_lambdaR, DP_fdis,
     &  DP_rhoR, DP_C0, DP_S0, DP_R, DP_mb0, DP_u0,
     &  DP_Cd, DP_gammaT, DP_gammaS,
     &  DP_alpha1, DP_alpha2, DP_alpha3, DP_lambdaS,
     &  DP_epsilon, DP_c1, DP_c2, DP_gammaI, DP_gammaO,
     &  DP_kappaPeel, DP_fe, DP_fp,
     &  DP_dz, DP_FrOuter, DP_FrWuest, DP_uStall, DP_JoMin,
     &  DP_Ca, DP_usConst,
     &  DP_rtolInner, DP_atolInner, DP_hmaxInner,
     &  DP_rtolOuter, DP_atolOuter, DP_hmaxOuter,
     &  DP_surfTol, DP_bottomTol, DP_peelTol, DP_fStall,
     &  DP_convTol, DP_rampTime
      COMMON /DOUBLEPLUME_PHYS_I/ DP_nwidths, DP_maxIter, DP_nUpdate

C--   Work variables of the plume solver (current column)
C     dp_H     :: height of the column above the source = |iceDepth|
C     dp_nb0   :: source bubble number flux (Python self.nb0) [1/s]
C     dp_mbSrc   :: source bubble mass flux (dissolution test) [kg/s]
C     dp_T0src :: source (liquidus) temperature  [degC]
C     dp_nAmb  :: deepest wet level above the source (ambient T/S)
C     dpo_nFine:: number of fine-grid points (= plumeLayers)
C     dpo_nTop :: highest fine-grid index with an inner plume
      _RL dp_H, dp_nb0, dp_mbSrc, dp_T0src
      INTEGER dp_nAmb, dpo_nFine, dpo_nTop
      COMMON /DOUBLEPLUME_WORK_R/ dp_H, dp_nb0, dp_mbSrc, dp_T0src
      COMMON /DOUBLEPLUME_WORK_I/ dp_nAmb, dpo_nFine, dpo_nTop

C---------------------------------------------

#ifdef ALLOW_PTRACERS
C Parameters relating to PTRACERS

      COMMON /DOUBLEPLUME_PTRACERS_RL/
     &     ptr_addMass3D,
     &     ptracerMask
      _RL ptr_addMass3d
     &     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,PTRACERS_num)
      _RL ptracerMask
     &     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,PTRACERS_num,nSx,nSy)

      COMMON /DOUBLEPLUME_PTRACERS_FILES/
     &     ptracerMaskFile
      CHARACTER*(512)
     &     ptracerMaskFile

#endif /* ALLOW_PTRACERS */

#endif /* ALLOW_DOUBLEPLUME */
