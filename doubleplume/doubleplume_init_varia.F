
#include "DOUBLEPLUME_OPTIONS.h"

CBOP
      SUBROUTINE DOUBLEPLUME_INIT_VARIA( myThid )
C     *============================================================*
C     | SUBROUTINE DOUBLEPLUME_INIT_VARIA
C     | o Initialize DOUBLEPLUME variables.
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
      INTEGER myThid
CEOP

#ifdef ALLOW_DOUBLEPLUME
C     === Local variables ===
      INTEGER I, J, K, bi, bj
#ifdef ALLOW_PTRACERS
      INTEGER iTracer
#endif

      _BEGIN_MASTER(myThid)
      dp_firstCall = .TRUE.
      j_peel = 0
      DO K = 1, DP_maxPeel
        z_peel(K) = 0. _d 0
      ENDDO
      DO K = 1, Nr
        saProf (K) = 0. _d 0
        taProf (K) = 0. _d 0
        ptProf (K) = 0. _d 0
        prProf (K) = 0. _d 0
        uProf  (K) = 0. _d 0
        vProf  (K) = 0. _d 0
        mProf  (K) = 0. _d 0
        mProfAv(K) = 0. _d 0
        HeatFlux(K)= 0. _d 0
        FwFlux (K) = 0. _d 0
        delta_z(K) = 0. _d 0
      ENDDO
      DO K = 1, Nr+1
        QiProf (K) = 0. _d 0
        MiProf (K) = 0. _d 0
        HiProf (K) = 0. _d 0
        HgProf (K) = 0. _d 0
        BiProf (K) = 0. _d 0
        MbProf (K) = 0. _d 0
        CiProf (K) = 0. _d 0
        QoProf (K) = 0. _d 0
        MoProf (K) = 0. _d 0
        HoProf (K) = 0. _d 0
        BoProf (K) = 0. _d 0
        CoProf (K) = 0. _d 0
      ENDDO
      DO K = 1, DP_maxFine
        QiProf_f(K) = 0. _d 0
        MiProf_f(K) = 0. _d 0
        HiProf_f(K) = 0. _d 0
        HgProf_f(K) = 0. _d 0
        BiProf_f(K) = 0. _d 0
        MbProf_f(K) = 0. _d 0
        CiProf_f(K) = 0. _d 0
        QoProf_f(K) = 0. _d 0
        MoProf_f(K) = 0. _d 0
        HoProf_f(K) = 0. _d 0
        BoProf_f(K) = 0. _d 0
        CoProf_f(K) = 0. _d 0
      ENDDO
      DO I = 1,100
        runoffLocation(I,1) = 0
        runoffLocation(I,2) = 0
        plumeInputLayer(I)  = 0
        Qin(I)              = 0. _d 0
      ENDDO
      _END_MASTER(myThid)

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO K = 1, Nr
         DO J = 1, sNy
          DO I = 1, sNx
           icefront_TendT(I,J,K,bi,bj) = 0. _d 0
           icefront_TendS(I,J,K,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
         DO J = 1-OLy, sNy+OLy
          DO I = 1-OLx, sNx+OLx
           HeatFlux3D(I,J,K,bi,bj) = 0. _d 0
           FwFlux3D  (I,J,K,bi,bj) = 0. _d 0
           QiProf3D  (I,J,K,bi,bj) = 0. _d 0
           MiProf3D  (I,J,K,bi,bj) = 0. _d 0
           HiProf3D  (I,J,K,bi,bj) = 0. _d 0
           BiProf3D  (I,J,K,bi,bj) = 0. _d 0
           HgProf3D  (I,J,K,bi,bj) = 0. _d 0
           MbProf3D  (I,J,K,bi,bj) = 0. _d 0
           CiProf3D  (I,J,K,bi,bj) = 0. _d 0
           QoProf3D  (I,J,K,bi,bj) = 0. _d 0
           MoProf3D  (I,J,K,bi,bj) = 0. _d 0
           HoProf3D  (I,J,K,bi,bj) = 0. _d 0
           BoProf3D  (I,J,K,bi,bj) = 0. _d 0
           CoProf3D  (I,J,K,bi,bj) = 0. _d 0
           mProf3D   (I,J,K,bi,bj) = 0. _d 0
           mProfAv3D (I,J,K,bi,bj) = 0. _d 0
           volFluxDiff3D (I,J,K,bi,bj) = 0. _d 0
C          no T/S change from addMass until the plume is computed
           temp_addMass3D(I,J,K,bi,bj) = theta(I,J,K,bi,bj)
           salt_addMass3D(I,J,K,bi,bj) = salt(I,J,K,bi,bj)
#ifdef ALLOW_PTRACERS
           DO iTracer = 1,PTRACERS_num
            ptr_addMass3D(I,J,K,bi,bj,iTracer)
     &                 = pTracer(I,J,K,bi,bj,iTracer)
           ENDDO
#endif /* ALLOW_PTRACERS */
           addMass(I,J,K,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END
