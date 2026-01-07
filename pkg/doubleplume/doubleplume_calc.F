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
C     | o Send ambient conditions to plume model
C     | o Calculate source/sink terms to parameterise movement of
C     |   water and tracers between vertical layers by resulting plume
C     | o Calcualte melt rates and tendencies due to melting of
C     |   ice front in none-plume locations
C     |   
C     *==========================================================*

C     \ev
C     !USES:
      IMPLICIT NONE
C     === Global variables ===
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

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
      _RL futureTime
      INTEGER futureIter
      INTEGER myThid
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifdef ALLOW_DOUBLEPLUME

C     !LOCAL VARIABLES:
C     == Local variables ==
C     I,J,K,bi,bj  :: loop indices
C     msgBuf       :: Informational/error message buffer
C     BiProf, HiProf, uProf, vProf :: salt, pot. temperature and
C                            uVel and vVel profiles
C     ptrProf  :: ambient ptracer profile
C     ptrPlumeCum :: cumulative quantity of ptracer in plume
C     ptrPlume :: quantity of ptracer in plume outflow
C     eps5     :: for thermodynamics (see pkg icefront)
C     maxDepth :: vertical extent of domain (m)
C     plumeAreaInCell :: surface area of plume in contact with ice in that cell (m^2)
C     negSum, posSum :: sum of negative and positive contributions to the plume volume
C     posNegRatio    :: ratio of the above
C     wVelCell, vVelCell, uVelCell :: velocities at cell centres
C     meanVel :: ice tangental velocity
C     rho_0 :: average density of seawater

      INTEGER bi, bj
      INTEGER J, K, I, iter
      LOGICAL DONE 
      DOUBLE PRECISION error, QiProf_temp(plumeLayers)
C      DOUBLE PRECISION QiProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION MiProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION HiProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION HgProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION BiProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION MbProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION CiProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION QoProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION MoProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION HoProf_f_3D(sNx, sNy, Nr)
C      DOUBLE PRECISION BoProf_f_3D(sNx, sNy, Nr)
      DOUBLE PRECISION heatFluxDiff(Nr+1), saltFluxDiff(Nr+1)
      DOUBLE PRECISION Qt(plumeLayers), Bt(plumeLayers)
      DOUBLE PRECISION Ht(plumeLayers), volDiff(plumeLayers)
      DOUBLE PRECISION Qt_c(Nr), Ht_c(Nr), Bt_c(Nr)
      DOUBLE PRECISION T_c(Nr), S_c(Nr), volDiff_c(Nr)
      DOUBLE PRECISION ramp
      INTEGER end_i, start_i, KK

      _RL eps5
      _RL plumeAreaInCell
      _RL negSum, posSum, posNegRatio
      _RL wVelCell, vVelCell, uVelCell, meanVel
      _RL sw_temp, sw_ptmp
      _RL rho_0
      external SW_TEMP
      external SW_PTMP

      PARAMETER(rho_0=1027.0D0)
      print *, "iter_l: n", iter_l
      iter_l = iter_l + 1.d0
      print *, "iter_l: n", iter_l
      print *, "In doubleplume_calc"
      print *, "I,x range: ", 1-OLx,sNx+Olx
      print *, "J,y range: ", 1-OLy,sNy+OLy
      print *, "plumeMask shape: ", SHAPE(plumeMask)
      print *, "current iteration: ", iter_l
      ramp = 30.d0 * 24.d0 * 60.d0 * 60.d0 ! ramp input over 30 days
C -----------------------------------------
C Enter into loops
      DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
          DO I = 1-OLx,sNx+Olx
            DO J = 1-OLy,sNy+OLy              
              IF (J .EQ. 5 .AND. I .EQ. 3) THEN
                plumeMask(I,J,bi,bj) = 1 
              ELSE 
                plumeMask(I,J,bi,bj) = 0 
              END IF
C Check to see if there is ice in that cell. If not, skip to end.
              IF ( plumeMask(I,J,bi,bj) .NE. 0 ) THEN
              print *, "plume mask at: I,J", I, J

              DO K = 1, Nr
              QiProf(K) = 0.d0
              MiProf(K) = 0.d0
              HiProf(K) = 0.d0
              BiProf(K) = 0.d0
              QoProf(K) = 0.d0
              MoProf(K) = 0.d0
              HoProf(K) = 0.d0
              BoProf(K) = 0.d0
              END DO
C Read from the plume mask which type of plume should be used in this cell.
C 1 = ice but no plume (melting only)
C 2 = sheet plume (Jenkins)
C 3 = half-conical plume (Morton/Slater)
C 4 = both sheet plume and half-conical plume (NOT YET IMPLEMENTED)
C 5 = detaching conical plume (Goldberg)

C POSITIVE values indicate ice front is orientated north-south
C NEGATIVE values indicate ice front is orientated east-west

C If there is subglacial discharge but no plume type defined, there will be no
C plume.
              IF ( plumeMask(I,J,bi,bj) .EQ. -1
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 1) THEN
                useSheetPlume = .FALSE.
                useConePlume = .FALSE.
                useDetachPlume = .FALSE.

              ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -2
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 2) THEN
                useSheetPlume = .TRUE.
                useConePlume = .FALSE.
                useDetachPlume = .FALSE.

              ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -3
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 3) THEN
                useSheetPlume = .FALSE.
                useConePlume = .TRUE.
                useDetachPlume = .FALSE.

              ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -4
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 4) THEN
                useSheetPlume = .TRUE.
                useConePlume = .TRUE.
                useDetachPlume = .FALSE.

              ELSEIF ( plumeMask(I,J,bi,bj) .EQ. -5
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 5) THEN

              ELSE
                WRITE(msgBuf,'(2A)')
     &         'Plume mask value must be between -5 and 5'
                CALL PRINT_ERROR( msgBuf, myThid )
                STOP 'ABNORMAL END: S/R DOUBLEPLUME_CALC'

              ENDIF

C Read in the subglacial discharge for this place and time
C (only the lowermost cell in column is read)
C             w_sg = runoffVel(I,J,bi,bj)
C             r_sg = runoffRad(I,J,bi,bj)

C             IF ( useSheetPlume ) THEN
C               Q_sg = w_sg*r_sg*delta_y
C             ELSEIF ( useConePlume ) THEN
C               Q_sg = w_sg*((pi*r_sg**2.)/2.)
C             ELSEIF ( useDetachPlume ) THEN
C               Q_sg = w_sg*((pi*r_sg**2.)/2.)
C             ELSE
C               Q_sg = 0.D0
C             ENDIF
              print *, "setting ambient profiles"
              DO K = 1,Nr
C           Tracers
                prProf(k) = ABS(rC(k))*rho_0*9.81*1.0E-6 ! Pressure (dbar)
                saProf(K)  = salt(I,J,K,bi,bj)         ! Salinity
                ptProf(K) = theta(I,J,K,bi,bj)        ! Potential Temperature
                taProf(k)  = 
     &           SW_TEMP(saProf(k),ptProf(k),prProf(k),0. _d 0) ! Temperature
                vProf(k) = ABS(vVel(I,J,K,bi,bj))          ! v velocity
                uProf(K) = ABS(uVel(I,J,K,bi,bj))          ! u Velocity
                delta_z(k) = drF(K)

              ENDDO
              Q_sg = 1.D0

C         Cell resolution
              delta_x = dxG(I,J,bi,bj)
              delta_y = dyG(I,J,bi,bj)

C Create variables with temperature, salinity
C and velocity profiles for that column

C Vertical extent of domain
              maxDepth = rF(Nr+1)

C Grounding line is in cell above deepest cell
C to avoid interpolation issues in plume model   
              print *, "finding iceDepth"

              DO K = 1,Nr+1
                IF ( rF(K) .LE. R_low(I,J,bi,bj) ) THEN
                  iceDepth = rF(k-1)
                  EXIT
                ENDIF
              ENDDO

              IF ( iceDepth .EQ. 0 ) THEN
                WRITE(msgBuf,'(2A)')
     &          'Plume specified in cell I = ', I, ', J = ', J,
     &          ', but depth of this cell = 0'
                CALL PRINT_ERROR( msgBuf, myThid )
                STOP 'ABNORMAL END: S/R DOUBLEPLUME_CALC III'
              ENDIF
              
C     Update Plume layers 
              plumeLayers = NINT(ABS(iceDepth) / 0.1) + 1

C Find grid layer at depth of ice face (move this to initialisation)
              icedepthK = 0
              DO K=1,Nr+1
                IF ( rF(K) .EQ. iceDepth ) iceDepthK = K
              ENDDO

C If a matching grid layer is not found, this may be because the bottom layer is a partial cell
C In this case, start in cell above partial cell
              IF ( iceDepthK .EQ. 0 ) THEN
                DO K=1,Nr+1
                  IF ( rF(K) .GT. iceDepth ) THEN
                    IF ( rF(K+1) .LT. iceDepth ) THEN
                      iceDepthK = K
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF

C If we still cannot find the bottom cell

              IF ( iceDepthK .EQ. 0 ) THEN
                WRITE(msgBuf,'(2A)')
     &          'Unable to identify index of cell',
     &          'at grounding line.',
     &          'This may be because this is a partial cell.'
                CALL PRINT_ERROR( msgBuf, myThid )
                STOP 'ABNORMAL END: S/R DOUBLEPLUME_CALC IV'
              ENDIF

C --- If there is subglacial outflow in that column, then parameterise plume ---
              
C              Q_sg = 20.D0
              j_peel = 0
              iter = 1
              error = 1000.D0 
              DONE = .FALSE.

              DO K = 1, plumeLayers
              QiProf_temp(K) = 0.D0 
              QoProf_f(K) = 0.D0
              QiProf_f(K) = 0.D0
              END DO
              IF ( Q_sg.GT.0 ) THEN
C                DO WHILE (.NOT. DONE)
C This routine calculates T, S and W and r profiles for plume
                CALL DOUBLEPLUME_INNER_MODEL (mythid)
                print *, "Finished running inner plume."
                print *, "j_peel = ", j_peel
                print *, "z_peel : ", z_peel(1:j_peel)

                print *, "Now running outer plume."
                CALL DOUBLEPLUME_OUTER_MODEL (mythid)
                IF (iter .EQ. 1) THEN 
                error = error
                ELSE
                error = MAXVAL((QiProf_temp(10:plumeLayers-10) 
     &                          - QiProf_f(10:plumeLayers - 10)) 
     &                / QiProf_temp(10:plumeLayers-10))
                END IF 

                print *, "Percent error, iter: ", error, iter 
                iter = iter + 1
                IF (iter .GT. 1) DONE = .TRUE.
                IF (error .LT. 0.1) DONE = .TRUE.
C                END DO

C Calculate vertical plume volume flux...
C Calculate coarse grid profiles
      
                Qt = QoProf_f + QiProf_f
                Bt = BoProf_f + BiProf_f
                DO K = 1, plumeLayers
                IF (QoProf_f(K) .EQ. 0) THEN
                Ht(K) = HiProf_f(K)/(1031.d0 * 4000.d0 * QiProf_f(K)) 
     &                  - 273.15d0
                ELSE IF (QiProf_f(K) .EQ. 0) THEN
                Ht(K) = 0.d0
                ELSE
                Ht(K) = HoProf_f(K) / QoProf_f(K) 
                END IF 
                END DO
                Ht = Ht * Qt

                DO K = 1, plumeLayers-1
                IF (ABS(Qt(K)) . LT. 0.1) THEN
                Qt(K) = 0.d0
                Bt(K) = 0.d0
                Ht(K) = 0.d0
                END IF
                volDiff(K) = Qt(K+1) - Qt(K)
                END DO
                
                DO K = 1, Nr
                volDiff_c(K) = 0.d0 
                Qt_c(K) = 0.d0 
                Ht_c(K) = 0.d0
                Bt_c(K) = 0.d0
                END DO

                end_i = 100
                start_i = 0
                DO K = 1, Nr
                DO KK = start_i, end_i
                Qt_c(K) = Qt(KK) + Qt_c(K)
                Ht_c(K) = Ht(KK) + Ht_c(K)
                Bt_c(K) = Bt(KK) + Bt_c(K)
                volDiff_c(K) = volDiff(KK) + volDiff_c(K)
                END DO
                start_i = end_i 
                end_i = start_i + 100
                END DO

                end_i = 100
                start_i = 0
                DO K = 1, Nr
                DO KK = start_i, end_i
                QiProf(K) = QiProf_f(KK) + QiProf(K)
                MiProf(K) = MiProf_f(KK) + MiProf(K)
                HiProf(K) = HiProf_f(KK) + HiProf(K)
                BiProf(K) = BiProf_f(KK) + BiProf(K)
                HgProf(K) = HgProf_f(KK) + HgProf(K)
                MbProf(K) = MbProf_f(KK) + MbProf(K)
                CiProf(K) = CiProf_f(KK) + CiProf(K)
                END DO
                start_i = end_i 
                end_i = start_i + 100
                END DO

                T_c = Ht_c / Qt_c
                S_c = Bt_c / Qt_c
                                
C A couple of corrections:
C - even if plume is still buoyant, it cannot flow through the fjord surface
C                volFlux(1) = 0.D0
C - the initial volume flux is equal to runoff
C                volflux(iceDepthK) = Q_sg
              T_c = T_c(SIZE(T_c):1:-1)
              S_c = S_c(SIZE(S_c):1:-1)
              volDiff_c = -volDiff_c(SIZE(volDiff_c):1:-1)
              print *, "entering add mass"
              print *, "T_c:", T_c 
              print *, "S_c: ", S_c 
              print *, "volDiff_c: ", volDiff_c
C Find temperature and salinity of plume outflow
                DO k = 1,Nr
                  IF ( volFluxDiff(k) .GT. 0 ) THEN
                    temp_AddMass3D(I,J,K,bi,bj) = 
     &             SW_PTMP(S_c(k),T_c(k),
     &                    ptProf(k),0. _d 0)
                    salt_AddMass3D(I,J,K,bi,bj) = S_c(k)

                  ELSE 
                    temp_AddMass3D(I,J,K,bi,bj) = 0.D0 
                    salt_AddMass3D(I,J,K,bi,bj) = 0.D0 

                  ENDIF /* ( volFluxDiff(k) .GT. 0 ) */
C Convert m3/s into kg/s
                  addMass(I,j,k,bi,bj) = MIN(1.d0, iter_l * deltaT
     &                                  /ramp) * volDiff_c(k) *1000 
                ENDDO /* k = 1,Nr line 494 */

              print *, "exit add mass."



              ELSE /* ( Q_sg .EQ. 0 ) */

C If no subglacial output, then there is no plume
                DO K = 1, plumeLayers
                  QiProf_f(K) = 0.D0
                  MiProf_f(K) = 0.D0
                  HiProf_f(K) = 0.D0
                  BiProf_f(K) = 0.D0
                  MbProf_f(K) = 0.D0
                  CiProf_f(K) = 0.D0
                  HgProf_f(K) = 0.D0

                  QoProf_f(K) = 0.D0
                  MoProf_f(K) = 0.D0
                  HoProf_f(K) = 0.D0
                  BoProf_f(K) = 0.D0

                  temp_AddMass3D(I,J,K,bi,bj) = 0.D0 
                  salt_AddMass3D(I,J,K,bi,bj) = 0.D0
                  addMass(I,j,k,bi,bj) = 0.D0

                ENDDO

              ENDIF /* ( Q_sg .EQ. 0 ) */

C Send outputs to 3D grid for diagnostics
C              IF ( usePlumeDiagnostics ) THEN
C                print *, "Doing Diagnostics."
C                DO K = 1, Nr                  
C                  QiProf3D(I,J,K,bi,bj) = QiProf(K)
C                  MiProf3D(I,J,K,bi,bj) = MiProf(K)
C                  HiProf3D(I,J,K,bi,bj) = HiProf(K)
C                  BiProf3D(I,J,K,bi,bj) = BiProf(K)
C                  CiProf3D(I,J,K,bi,bj) = CiProf(K)
C                  MbProf3D(I,J,K,bi,bj) = MbProf(K)
C                  HgProf3D(I,J,K,bi,bj) = HgProf(K)
C                ENDDO
                
C              ENDIF /* ( usePlumeDiagnostics ) */

C-------- Calculate melt rates ----------------------------
C              print *, "In calculate melt rates"
C              DO K = 1,Nr

C Check to see if we are above the sea bed
C                IF ( K .GE. iceDepthK ) THEN

C If not then there is no melting             
C                  mProfAv(k) = 0.D0
C                  mProf(k) = 0.D0
C                  mProf(k)      = 0.D0

C                ELSE /* ( K .GE. iceDepthK ) */

C If there is a plume in that cell, then need to calculate plume melt rate
C distinct to background melt rate. Plume melt rate is already encorporated in 
C the plrume model, and taken into account in the temperature and salinity of the
C plume outflow. It is useful though to have it available as a diagnostic.

C                  plumeAreaInCell = 0.0
            
C                  IF ( ( Q_sg .NE. 0 ) .AND. 
C     &            (useConePlume .OR. useSheetPlume.OR.useDetachPlume) )
C     &                                                             THEN

C                    plumeAreaInCell = aProf(k) - aProf(k+1)
              
C                    IF (plumeAreaInCell .gt. 0.0) THEN

C                      mProf(k) =(mIntProf(k) 
C     &                           -mIntProf(k+1))/ plumeAreaInCell
C     &                           * 86400.

C                    ELSE /* (plumeAreaInCell .gt. 0.0) */
C                      mProf (k) = 0.0

C                    ENDIF /* (plumeAreaInCell .gt. 0.0) */              

C                  ELSE /* ( Q_sg .NE. 0 ) */

C If there is no plume in that cell, set plume melt rate to zero
C                    mProf(k) = 0.D0

C                  ENDIF /* ( Q_sg .NE. 0 ) */

C Calculate the background melt rate (i.e. not generated by plumes). This will then be 
C used to update the temperature and salinity in the adjacent cells.
C Velocities are calculated at cell faces - find averages for cell centres.
C Does not include velocity perpendicular to ice - this differs depending on 
C orientation of ice front

C                  vVelCell = (ABS(vVel(I,J,K,bi,bj))
C     &                 +ABS(vVel(I,J+1,K,bi,bj))) / 2.

C                  uVelCell = (ABS(uVel(I,J,K,bi,bj))
C     &                 +ABS(uVel(I+1,J,K,bi,bj))) / 2.

C                  IF ( K .LT. Nr ) THEN
C                    wVelCell = (ABS(wVel(I,J,K,bi,bj))
C     &                 +ABS(wVel(I,J,K+1,bi,bj))) / 2.
C                  ELSE
C                    wVelCell = ABS(wVel(I,J,K,bi,bj)) / 2.
C                  ENDIF /* ( K .LT. Nr ) */

C                  IF ( plumeMask(I,J,bi,bj) .LT. 0 ) THEN
C Negative mask values indicate north-south ice front orientation
C                    meanVel = ((wVelCell**2.)+(uVelCell**2.))**0.5
C                  ELSEIF ( plumeMask(I,J,bi,bj) .GT. 0 ) THEN
C Positive mask values indicate east-west ice front orientation
C                    meanVel = ((wVelCell**2.)+(vVelCell**2.))**0.5
C                  ENDIF /* ( plumeMask(I,J,bi,bj) .LT. 0 ) */

C                  CALL DOUBLEPLUME_MELTRATE(
C     I              (HiProf(k)/(4185 * 1025 * QiProf(k)) - 273.15), 
C     O              BiProf(k)/QiProf(k), meanVel,rC(k), mProf(k) )

C Get average melt rate. This is useful for visualling melt patterns and 
C assessing overall melt rate of glacier.

C the following should apply to both conical and sheet plume models
C                  IF ( ( Q_sg .NE. 0 ) ) THEN

C                    plumeAreaInCell = aProf(k) - aProf(k+1)

C                    IF ( plumeAreaInCell .LE. delta_y*delta_z(k) ) THEN
C                      IF ( plumeAreaInCell .LE. 0 ) THEN
C If there is no plume in cell, then the melt rate is equal to the background melt rate.
C                        mprofAv(k) = mProf(K)
C                      ELSE
C If there is a plume in cell, calculate average melt rate
C                        mProfAv(k) = (mProf(k)*plumeAreaInCell
C     &                 +mProf(k)*(delta_y*delta_z(k)-plumeAreaInCell)) /
C     &                 (delta_y * delta_z(k))

C Scale down background melt rate to account for area occupied by plume
C (necessary so that tendency terms aren't over estimated)

C                       mProf(k) = mProf(k)*(1-plumeAreaInCell/
C     &                  (delta_y*delta_z(k)))

C                      ENDIF /* ( plumeAreaInCell .LE. 0 ) */

C                    ELSE ! ( plumeAreaInCell .LE. delta_y*delta_z(k) ) 

C if the plume contact area is larger than the cell area, we assume there is
C no background melting

C                      mProfAv(k) = mProf(k)*plumeAreaInCell / 
C     &                     (delta_y*delta_z(k))
C                      mProf(k) = 0.

C                    ENDIF ! ( plumeAreaInCell .LE. delta_y*delta_z(k) ) 


C                  ELSE ! ( ( Q_sg .NE. 0 ) ) 

C If it is not a plume cell, then no plume melting.
C                    mProf(k) = 0.D0
C                    mProfAv(k) = mProf(k)

C                  ENDIF ! ( ( Q_sg .NE. 0 ) ) 
C                ENDIF /* ( K .GE. iceDepthK ) */

C Send outputs to 3D grid for diagnostics
C                IF ( usePlumeDiagnostics ) THEN
C                  mProf3D(I,J,K,bi,bj) = mProf(k)
C                  mProfAv3D(I,J,K,bi,bj) = mProfAv(k)
C                ENDIF /* ( usePlumeDiagnostics ) */


C ------------Tendencies------------------------------
C These are applied in cells where there is no plume. The idea is that in these areas there are likely to be local subgrid convection cells. As such, it is most realistic to apply changes in T and S to ice edge cell. Where there is a plume, products of melt are transported upwards so no local changes applied.
C The thermodynamics in this section are taken from pkg/icefront
C Tendencies are applied in S/R EXTERNAL_FORCING
                print *, "In tendencies"
C To convert from melt rate (m d^-1) to freshwater flux (kg m^-2 s^-1): mProf (m/d) * rho_shelfice_ref (mg/m^3) / 86400 (s/d) = mProf/94.22
C                FwFlux(k) = -mProf(k)/94.22

C Heat required to melt that much ice (W m^-2): FwFlux / latent heat of fusion (333.55E+3)
C                heatflux(k) = -FwFlux(k)*333.55E+3

C Create local (no overlap) arrays of heat and freshwater flux from background
C melting

C                IF ( ( J .GT. 0 ) .AND. ( J .LT. sNy+1 ) ) THEN
C                  Fwflux3D(I,J,k,bi,bj) = FwFlux(k)
C                  HeatFlux3D(I,J,k,bi,bj) = HeatFlux(k)
C                ENDIF
C                print *, "set fw/heatflux3d"
C     Compute tendencies (as for pkg/icefront)
C      eps5 = mass2rUnit/HeatCapacity_Cp

C                icefront_TendT(I,J,K,bi,bj) =
C     &               - HeatFlux3D(I,J,K,bi,bj)*eps5
C                icefront_TendS(I,J,K,bi,bj) = 
C    &                FWFlux3D(I,J,K,bi,bj)*
C     &                mass2rUnit * (BiProf(k)/QiProf(k)) ! UNSURE

C     Scale by icefrontlength, which is the ratio of the horizontal length
C     of the ice front in each model grid cell divided by the grid cell area.
C     (icefrontlength = dy / dxdy = 1 / dx)
C                print *, "set icefronttendt/icefronttends"
C                icefront_TendT(I,j,K,bi,bj) = 
C     &            icefront_TendT(I,j,K,bi,bj)
C     &            * 1./delta_x
C                icefront_TendS(I,j,K,bi,bj) = 
C     &            icefront_TendS(I,j,K,bi,bj)
C     &            * 1./delta_x
C                print *, "modified icefronttendt/icefronttends"
C              ENDDO /* DO K = 1,Nr line 321 */


C The plume transport is undertaken using the addMass terms.
C addMass terms for volume (kg/s), pot. temperature, salinity and ptracers are input 
C into the correct locations in 3D arrays

              

            ENDIF ! plumeMask .NE. 0

	    ENDDO ! DO J loop

      ENDDO ! DO I loop

C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

C Save plume values
#ifdef ALLOW_DIAGNOSTICS
C        IF ( usePlumeDiagnostics ) THEN
          print *, "QiProf: ", QiProf 
          print *, "MiProf: ", MiProf
         DO k = 1, Nr
          DO j = 1, sNy
            DO i = 1, sNx
            IF ( plumeMask(i,j,bi,bj) .NE. 0 ) THEN
              QiProf3D(i,j,k,bi,bj) = QiProf(k)
              MiProf3D(i,j,k,bi,bj) = MiProf(k)
              HiProf3D(i,j,k,bi,bj) = HiProf(k)
              HgProf3D(i,j,k,bi,bj) = HgProf(k)
              BiProf3D(i,j,k,bi,bj) = BiProf(k)
              MbProf3D(i,j,k,bi,bj) = MbProf(k)
              CiProf3D(i,j,k,bi,bj) = CiProf(k)
              QoProf3D(i,j,k,bi,bj) = QoProf(k)
              MoProf3D(i,j,k,bi,bj) = MoProf(k)
              HoProf3D(i,j,k,bi,bj) = HoProf(k)
              BoProf3D(i,j,k,bi,bj) = BoProf(k)

            ELSE
              QiProf3D(i,j,k,bi,bj) = 0.D0
              MiProf3D(i,j,k,bi,bj) = 0.D0
              HiProf3D(i,j,k,bi,bj) = 0.D0
              HgProf3D(i,j,k,bi,bj) = 0.D0
              BiProf3D(i,j,k,bi,bj) = 0.D0
              MbProf3D(i,j,k,bi,bj) = 0.D0
              CiProf3D(i,j,k,bi,bj) = 0.D0
              QoProf3D(i,j,k,bi,bj) = 0.D0
              MoProf3D(i,j,k,bi,bj) = 0.D0
              HoProf3D(i,j,k,bi,bj) = 0.D0
              BoProf3D(i,j,k,bi,bj) = 0.D0
            END IF
            END DO
          END DO
        END DO

C Transfer to local (no bi,bj indices) and interior only arrays
C (can't seem to get the diagnostics to output properly with the
C other options!)
C          DO I = 1,sNx
C            DO J = 1,sNy
C              DO K = 1,Nr
C                MiProf3dLocal(I,J,K) =
C     &            MiProf3D(I,J,K,bi,bj)
C                HiProf3dLocal(I,J,K) =
C     &            HiProf3D(I,J,K,bi,bj)
C                BiProf3dLocal(I,J,K) =
C     &            BiProf3D(I,J,K,bi,bj)
C                QiProf3dLocal(I,J,K) =
C     &            QiProf3D(I,J,K,bi,bj)
C                HgProf3dLocal(I,J,K) =
C     &            HgProf3D(I,J,K,bi,bj)
C                MbProf3dLocal(I,J,K) =
C     &            MbProf3D(I,J,K,bi,bj)
C                CiProf3dLocal(I,J,K) =
C     &            CiProf3D(I,J,K,bi,bj)
C                MoProf3dLocal(I,J,K) =
C     &            MoProf3D(I,J,K,bi,bj)
C                HoProf3dLocal(I,J,K) =
C     &            HoProf3D(I,J,K,bi,bj)
C                BoProf3dLocal(I,J,K) =
C     &            BoProf3D(I,J,K,bi,bj)
C                QoProf3dLocal(I,J,K) =
C     &            QoProf3D(I,J,K,bi,bj)
C              ENDDO
C            ENDDO
C          ENDDO

        CALL DIAGNOSTICS_FILL(MiProf3d, 'icefrtMi',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(QiProf3d, 'icefrtQi',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(HiProf3d, 'icefrtHi',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(BiProf3d, 'icefrtBi',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(CiProf3d, 'icefrtCi',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(QoProf3d, 'icefrtQo',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(MoProf3d, 'icefrtMo',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(HoProf3d, 'icefrtHo',
     &          0, Nr, 1, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(BoProf3d, 'icefrtBo',
     &          0, Nr, 1, bi, bj, myThid)

C Clear local arrays otherwise results replicate on other tiles
C          DO I = 1,sNx
C            DO J = 1,sNy
C              DO K = 1,Nr
C                MiProf3dLocal(I,J,K) = 0.d0
C                HiProf3dLocal(I,J,K) = 0.d0
C                BiProf3dLocal(I,J,K) = 0.d0
C                QiProf3dLocal(I,J,K) = 0.d0
C                CiProf3dLocal(I,J,K) = 0.d0
C                MbProf3dLocal(I,J,K) = 0.d0
C                HgProf3dLocal(I,J,K) = 0.d0
C                MoProf3dLocal(I,J,K) = 0.d0
C                HoProf3dLocal(I,J,K) = 0.d0
C                BoProf3dLocal(I,J,K) = 0.d0
C                QoProf3dLocal(I,J,K) = 0.d0
C              ENDDO
C            ENDDO
C          ENDDO

C        ENDIF /* ( useDiagnostics .AND. usePlumeDiagnostics ) */

#endif /* ALLOW_DIAGNOSTICS */ 

       ENDDO ! End bi loop
      ENDDO ! End bj loop

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END



