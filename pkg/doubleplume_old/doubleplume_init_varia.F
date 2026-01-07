
#include "DOUBLEPLUME_OPTIONS.h"

CBOP
      SUBROUTINE DOUBLEPLUME_INIT_VARIA( myThid )
C     *============================================================*
C     | SUBROUTINE DOUBLEPLUME_INIT_VARIA                             |
C     | o Routine to initialize DOUBLEPLUME variables.                |
C     *============================================================*
C     | Initialize DOUBLEPLUME parameters and variables.              |
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "FFIELDS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "DOUBLEPLUME.h"
#ifdef ALLOW_PTRACERS
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#endif /* ALLOW_PTRACERS */

C     === Routine arguments ===
C     myThid -  Number of this instance of ICEFRONT_INIT_VARIA
      INTEGER myThid
      INTEGER myIter
      INTEGER fp, iRec0
      INTEGER iTracer
      _RL aWght, bWght
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      _RL     myTime
CEndOfInterface

#ifdef ALLOW_DOUBLEPLUME
C     === Local variables ===
C     I,J,K,bi,bj - Loop counters
      INTEGER I, J, K, bi, bj
CEOP
      fp = readBinaryPrec
      iRec0 = 1

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO K = 1, Nr
         DO J = 1, sNy
          DO I = 1, sNx
           icefront_TendT (i,j,k,bi,bj) = 0. _d 0
           icefront_TendS (i,j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      DO k = 1,Nr
      BiProf  (K) = 0. _d 0
      HiProf  (K) = 0. _d 0
      ptProf (K) = 0. _d 0
      prProf (K) = 0. _d 0
      uProf  (K) = 0. _d 0
      vProf  (K) = 0. _d 0
      mProf  (K) = 0. _d 0
      mProfAv (K) =0. _d 0
      QiProf (K) = 0. _d 0
      MiProf (K) = 0. _d 0
      HiProf (K) = 0. _d 0
      BiProf (K) = 0. _d 0
      uProfPlume (K) = 0. _d 0
      volFlux     (K) = 0. _d 0
      volFluxDiff (K) = 0. _d 0
      HeatFlux(K) = 0. _d 0
      FwFlux  (K) = 0. _d 0
      delta_z     (K) = 0. _d 0
      ENDDO

      DO bi = 1,nSx
      DO bj = 1,nSy
       DO I = 1-OLx,sNx+OLx
        DO J = 1-Oly,sNy+Oly
         DO K = 1,Nr
          HeatFlux3D  (I,J,K,bi,bj) = 0.0
          FwFlux3D    (I,J,K,bi,bj) = 0.0
          QiProf3D (I,J,K,bi,bj) = 0.0
          MiProf3D (I,J,K,bi,bj) = 0.0
          HiProf3D (I,J,K,bi,bj) = 0.0
          BiProf3D (I,J,K,bi,bj) = 0.0
          mProf3D (I,J,K,bi,bj) = 0.0
          mProfAv3D    (I,J,K,bi,bj) = 0.0
          temp_addMass3D (I,J,K,bi,bj) = 0.D0
          salt_addMass3D (I,J,K,bi,bj) = 0.D0
          volFluxDiff3D (I,J,K,bi,bj) = 0.D0
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      ENDDO

       DO I = 1,sNx
        DO J = 1,sNy
         DO K = 1,Nr
          QiProf3DLocal (I,J,K) = 0.0
          MiProf3DLocal (I,J,K) = 0.0
          HiProf3DLocal (I,J,K) = 0.0
          BiProf3DLocal (I,J,K) = 0.0
          mProf3DLocal (I,J,K) = 0.0
          mProfAv3DLocal    (I,J,K) = 0.0
         ENDDO
        ENDDO
       ENDDO

       DO I = 1,100
        runoffLocation(I,1) = 0.0
        runoffLocation(I,2) = 0.0
        plumeInputLayer(I)  = 0.0
        Qin(I)              = 0.0
       ENDDO

C Clear addMass from previous timestep

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J = 1-OLy,sNy+OLy
         DO I = 1-OLx,sNx+OLx
          DO K = 1,Nr
            temp_AddMass3D(I,J,K,bi,bj) = theta(I,J,K,bi,bj)
            salt_AddMass3D(I,J,K,bi,bj) = salt(I,J,K,bi,bj)
#ifdef ALLOW_PTRACERS
            DO iTracer = 1,PTRACERS_num
             ptr_Addmass3D(I,J,K,bi,bj,iTracer)
     &                 = pTracer(I,J,K,bi,bj,iTracer)
            ENDDO
#endif /* ALLOW_PTRACERS */

            addMass(I,j,k,bi,bj) = 0.D0

          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO


#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END
