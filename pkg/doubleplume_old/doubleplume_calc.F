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
      INTEGER J, K, I

      _RL eps5
      _RL plumeAreaInCell
      _RL negSum, posSum, posNegRatio
      _RL wVelCell, vVelCell, uVelCell, meanVel
      _RL sw_temp, sw_ptmp
      _RL rho_0
      external SW_TEMP
      external SW_PTMP

      PARAMETER(rho_0=1027.0D0)


C -----------------------------------------
C Enter into loops
      DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
          DO I = 1-OLx,sNx+Olx
            DO J = 1-OLy,sNy+OLy


C Check to see if there is ice in that cell. If not, skip to end.
              IF ( plumeMask(I,J,bi,bj) .NE. 0 ) THEN
C                print *, "plume exists, mask = ", plumeMask(I,J,bi,bj)

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
              Q_sg = 0.D0
              IF ( plumeMask(I,J,bi,bj) .EQ. -1
     &           .OR. plumeMask(I,J,bi,bj) .EQ. 1) THEN
                useSheetPlume = .FALSE.
                useConePlume = .FALSE.
                useDetachPlume = .FALSE.
                Q_sg = 0.D0

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
                Q_sg = 20.D0

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

C         Cell resolution
              delta_x = dxG(I,J,bi,bj)
              delta_y = dyG(I,J,bi,bj)

C Create variables with temperature, salinity
C and velocity profiles for that column

	        DO K = 1,Nr
C           Tracers
                prProf(k) = ABS(rC(k))*rho_0*9.81*1.0E-6 ! Pressure (dbar)
                BiProf(K)  = salt(I,J,K,bi,bj)         ! Salinity
                ptProf(K) = theta(I,J,K,bi,bj)        ! Potential Temperature
                HiProf(k)  = 
     &          SW_TEMP(BiProf(k),ptProf(k),prProf(k),0. _d 0) ! Temperature

C           Velocities
                vProf(k) = ABS(vVel(I,J,K,bi,bj))          ! v velocity
                uProf(K) = ABS(uVel(I,J,K,bi,bj))          ! u Velocity
                delta_z(k) = drF(K)

	        ENDDO

C Vertical extent of domain
              maxDepth = rF(Nr+1)

C Grounding line is in cell above deepest cell
C to avoid interpolation issues in plume model

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
C              print *, "Q_sg = ", Q_sg
              IF ( Q_sg.GT.0 ) THEN

C This routine calculates T, S and W and r profiles for plume
C                print *, "jumping into plume model"
                CALL DOUBLEPLUME_PLUME_MODEL (mythid)

C Calculate vertical plume volume flux...

                DO k=1,Nr
C ... after checking to see if we are above the base of the ice face...
                  IF ( K .LT. iceDepthK ) THEN
                    volFlux(k) = QiProf(k)
                  ELSE
                    volFlux(k) = 0.D0
                  ENDIF
                ENDDO
C                print *, "VolFlux = ", volFlux
                                
C A couple of corrections:
C - even if plume is still buoyant, it cannot flow through the fjord surface
                volFlux(1) = 0.D0
C - the initial volume flux is equal to runoff
                volflux(iceDepthK) = Q_sg


C Calculate volume flux differential to give entrainment / extrainment
C First clear volfluxdiff

                DO K = 1,Nr
                  volfluxdiff(K) = 0.D0
                ENDDO

                DO k=1,iceDepthK-1
                  volFluxDiff(k) = volFlux(k+1) - volFlux(k)
                ENDDO

C                print *, "volFluxDiff = ", volFluxDiff

                IF ( conserveMass ) THEN
C Scale output to compensate for entrainment lost in expanding of output layer
C i.e. so that there is no net flow over boundary

                  negSum = 0.D0
                  posSum = 0.D0

                  DO K = 1,Nr
                    IF ( volFluxDiff(K) .LT. 0 ) THEN
                      negSum = volFluxDiff(K) + negSum
                    ELSE
                      posSum = volFluxDiff(K) + posSum
                    ENDIF
                  ENDDO

                  IF ( posSum .NE. 0 ) THEN
                    posNegRatio = -negSum / posSum

                    DO K = 1,Nr
                      IF ( volFluxDiff(K) .GT. 0 )
     &                  volFluxDiff(K) = volFluxDiff(K) * posNegRatio
                    ENDDO
                  ENDIF /* ( posSum .NE. 0 ) */

                ENDIF /* ( conserveMass ) */


              ELSE /* ( Q_sg .EQ. 0 ) */

C If no subglacial output, then there is no plume
                DO k = 1,Nr
                  QiProf(K) = 0.D0
                  MiProf(K) = 0.D0
                  HiProf(K) = 0.D0
                  BiProf(K) = 0.D0
                  aProf(K) = 0.D0
                  mIntProf(K) = 0.D0
                ENDDO

              ENDIF /* ( Q_sg .EQ. 0 ) */

C Send outputs to 3D grid for diagnostics
              IF ( usePlumeDiagnostics ) THEN
C               print *, "Doing Diagnostics."
                IF ( Q_sg .GT. 0 ) THEN
C                  print *, "Diagnostics QiProf = ", QiProf
                ENDIF
                DO K = 1,Nr
                  QiProf3D(I,J,K,bi,bj) = QiProf(k)
                  MiProf3D(I,J,K,bi,bj) = MiProf(k)
                  HiProf3D(I,J,K,bi,bj) = HiProf(k)
                  BiProf3D(I,J,K,bi,bj) = BiProf(k)
                  IF ( Q_sg .GT. 0 ) THEN
C                    print *, "QiProf3D = ", K, QiProf3D(I,J,K,bi,bj)
                  ENDIF
                ENDDO
                
              ENDIF /* ( usePlumeDiagnostics ) */

C-------- Calculate melt rates ----------------------------

              DO K = 1,Nr

C Check to see if we are above the sea bed
                IF ( K .GE. iceDepthK ) THEN

C If not then there is no melting             
                  mProfAv(k) = 0.D0
                  mProf(k) = 0.D0
                  mProf(k)      = 0.D0

                ELSE /* ( K .GE. iceDepthK ) */

C If there is a plume in that cell, then need to calculate plume melt rate
C distinct to background melt rate. Plume melt rate is already encorporated in 
C the plrume model, and taken into account in the temperature and salinity of the
C plume outflow. It is useful though to have it available as a diagnostic.

                  plumeAreaInCell = 0.0
            
                  IF ( ( Q_sg .NE. 0 ) .AND. 
     &            (useConePlume .OR. useSheetPlume.OR.useDetachPlume) )
     &                                                             THEN

                    plumeAreaInCell = aProf(k) - aProf(k+1)
              
                    IF (plumeAreaInCell .gt. 0.0) THEN

                      mProf(k) =(mIntProf(k) 
     &                           -mIntProf(k+1))/ plumeAreaInCell
     &                           * 86400.

                    ELSE /* (plumeAreaInCell .gt. 0.0) */
                      mProf (k) = 0.0

                    ENDIF /* (plumeAreaInCell .gt. 0.0) */              

                  ELSE /* ( Q_sg .NE. 0 ) */

C If there is no plume in that cell, set plume melt rate to zero
                    mProf(k) = 0.D0

                  ENDIF /* ( Q_sg .NE. 0 ) */

C Calculate the background melt rate (i.e. not generated by plumes). This will then be 
C used to update the temperature and salinity in the adjacent cells.
C Velocities are calculated at cell faces - find averages for cell centres.
C Does not include velocity perpendicular to ice - this differs depending on 
C orientation of ice front

                  vVelCell = (ABS(vVel(I,J,K,bi,bj))
     &                 +ABS(vVel(I,J+1,K,bi,bj))) / 2.

                  uVelCell = (ABS(uVel(I,J,K,bi,bj))
     &                 +ABS(uVel(I+1,J,K,bi,bj))) / 2.

                  IF ( K .LT. Nr ) THEN
                    wVelCell = (ABS(wVel(I,J,K,bi,bj))
     &                 +ABS(wVel(I,J,K+1,bi,bj))) / 2.
                  ELSE
                    wVelCell = ABS(wVel(I,J,K,bi,bj)) / 2.
                  ENDIF /* ( K .LT. Nr ) */

                  IF ( plumeMask(I,J,bi,bj) .LT. 0 ) THEN
C Negative mask values indicate north-south ice front orientation
                    meanVel = ((wVelCell**2.)+(uVelCell**2.))**0.5
                  ELSEIF ( plumeMask(I,J,bi,bj) .GT. 0 ) THEN
C Positive mask values indicate east-west ice front orientation
                    meanVel = ((wVelCell**2.)+(vVelCell**2.))**0.5
                  ENDIF /* ( plumeMask(I,J,bi,bj) .LT. 0 ) */

                  CALL DOUBLEPLUME_MELTRATE(
     I              HiProf(k)/QiProf(k),BiProf(k)/QiProf(k),
     O              meanVel,rC(k), mProf(k) )

C Get average melt rate. This is useful for visualling melt patterns and 
C assessing overall melt rate of glacier.

C the following should apply to both conical and sheet plume models
                  IF ( ( Q_sg .NE. 0 ) ) THEN

                    plumeAreaInCell = aProf(k) - aProf(k+1)

                    IF ( plumeAreaInCell .LE. delta_y*delta_z(k) ) THEN
                      IF ( plumeAreaInCell .LE. 0 ) THEN
C If there is no plume in cell, then the melt rate is equal to the background melt rate.
                        mprofAv(k) = mProf(K)
                      ELSE
C If there is a plume in cell, calculate average melt rate
                        mProfAv(k) = (mProf(k)*plumeAreaInCell
     &                 +mProf(k)*(delta_y*delta_z(k)-plumeAreaInCell)) /
     &                 (delta_y * delta_z(k))

C Scale down background melt rate to account for area occupied by plume
C (necessary so that tendency terms aren't over estimated)

                       mProf(k) = mProf(k)*(1-plumeAreaInCell/
     &                  (delta_y*delta_z(k)))

                      ENDIF /* ( plumeAreaInCell .LE. 0 ) */

                    ELSE ! ( plumeAreaInCell .LE. delta_y*delta_z(k) ) 

C if the plume contact area is larger than the cell area, we assume there is
C no background melting

                      mProfAv(k) = mProf(k)*plumeAreaInCell / 
     &                     (delta_y*delta_z(k))
                      mProf(k) = 0.

                    ENDIF ! ( plumeAreaInCell .LE. delta_y*delta_z(k) ) 


                  ELSE ! ( ( Q_sg .NE. 0 ) ) 

C If it is not a plume cell, then no plume melting.
                    mProf(k) = 0.D0
                    mProfAv(k) = mProf(k)

                  ENDIF ! ( ( Q_sg .NE. 0 ) ) 
                ENDIF /* ( K .GE. iceDepthK ) */

C Send outputs to 3D grid for diagnostics
                IF ( usePlumeDiagnostics ) THEN
                  mProf3D(I,J,K,bi,bj) = mProf(k)
                  mProfAv3D(I,J,K,bi,bj) = mProfAv(k)
                ENDIF /* ( usePlumeDiagnostics ) */


C ------------Tendencies------------------------------
C These are applied in cells where there is no plume. The idea is that in these areas there are likely to be local subgrid convection cells. As such, it is most realistic to apply changes in T and S to ice edge cell. Where there is a plume, products of melt are transported upwards so no local changes applied.
C The thermodynamics in this section are taken from pkg/icefront
C Tendencies are applied in S/R EXTERNAL_FORCING

C To convert from melt rate (m d^-1) to freshwater flux (kg m^-2 s^-1): mProf (m/d) * rho_shelfice_ref (mg/m^3) / 86400 (s/d) = mProf/94.22
                FwFlux(k) = -mProf(k)/94.22

C Heat required to melt that much ice (W m^-2): FwFlux / latent heat of fusion (333.55E+3)
                heatflux(k) = -FwFlux(k)*333.55E+3

C Create local (no overlap) arrays of heat and freshwater flux from background
C melting

                IF ( ( J .GT. 0 ) .AND. ( J .LT. sNy+1 ) ) THEN
                  Fwflux3D(I,J,k,bi,bj) = FwFlux(k)
                  HeatFlux3D(I,J,k,bi,bj) = HeatFlux(k)
                ENDIF
    
C     Compute tendencies (as for pkg/icefront)
      eps5 = mass2rUnit/HeatCapacity_Cp

                icefront_TendT(I,J,K,bi,bj) =
     &               - HeatFlux3D(I,J,K,bi,bj)*eps5
                icefront_TendS(I,J,K,bi,bj) = 
     &                FWFlux3D(I,J,K,bi,bj)*
     &                mass2rUnit * BiProf(k)

C     Scale by icefrontlength, which is the ratio of the horizontal length
C     of the ice front in each model grid cell divided by the grid cell area.
C     (icefrontlength = dy / dxdy = 1 / dx)
 
                icefront_TendT(I,j,K,bi,bj) = 
     &            icefront_TendT(I,j,K,bi,bj)
     &            * 1./delta_x
                icefront_TendS(I,j,K,bi,bj) = 
     &            icefront_TendS(I,j,K,bi,bj)
     &            * 1./delta_x
              ENDDO /* DO K = 1,Nr line 321 */


C The plume transport is undertaken using the addMass terms.
C addMass terms for volume (kg/s), pot. temperature, salinity and ptracers are input 
C into the correct locations in 3D arrays
              IF ( Q_sg.NE.0 ) THEN
C Find temperature and salinity of plume outflow
                DO k = 1,Nr
                  IF ( volFluxDiff(k) .GT. 0 ) THEN
                    temp_AddMass3D(I,J,K,bi,bj) = 0.D0 ! convert to potential temp
C     &             SW_PTMP(BiProf(k),HiProf(k),ptProf(k),0. _d 0)
                    salt_AddMass3D(I,J,K,bi,bj) = 0.D0 !BiProf(k)

                  ELSE 
                    temp_AddMass3D(I,J,K,bi,bj) = 0.D0 !HiProf(k)
                    salt_AddMass3D(I,J,K,bi,bj) = 0.D0 !BiProf(k)

                  ENDIF /* ( volFluxDiff(k) .GT. 0 ) */
C Convert m3/s into kg/s
                  addMass(I,j,k,bi,bj) = 0.D0 !volFLuxDiff(k)*1000 
                ENDDO /* k = 1,Nr line 494 */

              ELSE ! No Plume
                DO K = 1,Nr
                  temp_AddMass3D(I,J,K,bi,bj) = 2.D0 !taProf(k) !HiProf(k)
                  salt_AddMass3D(I,J,K,bi,bj) = 34.D0 !saProf(k) !BiProf(k)
                  addMass(I,j,k,bi,bj) = 0.D0
                ENDDO

              ENDIF /* ( Q_sg.NE.0 )  */

            ENDIF ! plumeMask .NE. 0

	    ENDDO ! DO J loop

        ENDDO ! DO I loop

C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

C Save plume values
#ifdef ALLOW_DIAGNOSTICS

        IF ( usePlumeDiagnostics ) THEN

C Transfer to local (no bi,bj indices) and interior only arrays
C (can't seem to get the diagnostics to output properly with the
C other options!)
          DO I = 1,sNx
            DO J = 1,sNy
              DO K = 1,Nr
                MiProf3dLocal(I,J,K) =
     &            MiProf3D(I,J,K,bi,bj)
                HiProf3dLocal(I,J,K) =
     &            HiProf3D(I,J,K,bi,bj)
                BiProf3dLocal(I,J,K) =
     &            BiProf3D(I,J,K,bi,bj)
                QiProf3dLocal(I,J,K) =
     &            QiProf3D(I,J,K,bi,bj)
                mProf3dLocal(I,J,K) =
     &            mProf3D(I,J,K,bi,bj)
                mProfAv3dLocal(I,J,K) =
     &            mProfAv3D(I,J,K,bi,bj)
                IF ( Q_sg .GT. 0 ) THEN
C                    print *, "QiProfLocal = ", K, QiProf3dLocal(I,J,K)
                ENDIF

              ENDDO
            ENDDO
          ENDDO


C Output diagnostics
          DO K=1,Nr
C Here, the first 'k' is the layer in the output field in which to save the data
C and the second k is layer from which this data is taken in the original field
            CALL DIAGNOSTICS_FILL_RS(MiProf3DLocal,'icefrntW',
     &         k,k,3,bi,bj,myThid)
            CALL DIAGNOSTICS_FILL_RS(HiProf3DLocal,'icefrntT',
     &         k,k,3,bi,bj,myThid)
            CALL DIAGNOSTICS_FILL_RS(BiProf3DLocal,'icefrntS',
     &         k,k,3,bi,bj,myThid)
            CALL DIAGNOSTICS_FILL_RS(QiProf3DLocal,'icefrntR',
     &         k,k,3,bi,bj,myThid)
            CALL DIAGNOSTICS_FILL_RS(mProf3DLocal,'icefrntM',
     &         k,k,3,bi,bj,myThid)
            CALL DIAGNOSTICS_FILL_RS(mProfAv3DLocal,'icefrntA',
     &         k,k,3,bi,bj,myThid)

          ENDDO

C Clear local arrays otherwise results replicate on other tiles
          DO I = 1,sNx
            DO J = 1,sNy
              DO K = 1,Nr
                MiProf3dLocal(I,J,K) = 0.d0
                HiProf3dLocal(I,J,K) = 0.d0
                BiProf3dLocal(I,J,K) = 0.d0
                QiProf3dLocal(I,J,K) = 0.d0
                mProf3dLocal(I,J,K) = 0.d0
                mProfAv3dLocal(I,J,K) = 0.d0

              ENDDO
            ENDDO
          ENDDO

        ENDIF /* ( useDiagnostics .AND. usePlumeDiagnostics ) */

#endif /* ALLOW_DIAGNOSTICS */ 

       ENDDO ! End bi loop
      ENDDO ! End bj loop

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END



