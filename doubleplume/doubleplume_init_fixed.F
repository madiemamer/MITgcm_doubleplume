
#include "DOUBLEPLUME_OPTIONS.h"

      SUBROUTINE DOUBLEPLUME_INIT_FIXED( myThid )
C     *============================================================*
C     | SUBROUTINE DOUBLEPLUME_INIT_FIXED
C     | o Read the plume mask, the runoff and the ptracer mask;
C     |   register the diagnostics.
C     *============================================================*
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DOUBLEPLUME.h"

C     === Routine arguments ===
      INTEGER myThid

#ifdef ALLOW_DOUBLEPLUME
C     === Local variables ===
      INTEGER I, J, bi, bj, iRec0
#ifdef ALLOW_PTRACERS
      INTEGER iTracer
#endif
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      iRec0 = 1

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO J = 1-OLy,sNy+OLy
         DO I = 1-OLx,sNx+OLx
          plumeMask(I,J,bi,bj) = 0. _d 0
          runoffVel(I,J,bi,bj) = 0. _d 0
          runoffRad(I,J,bi,bj) = 0. _d 0
#ifdef ALLOW_PTRACERS
          DO iTracer = 1,PTRACERS_num
           ptracerMask(I,J,iTracer,bi,bj) = 0. _d 0
          ENDDO
#endif /* ALLOW_PTRACERS */
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C--   Plume mask: standard MITgcm read into the tiled array
C     (the old direct-access read of "plumeMask.bin" into the whole
C     tiled array, overlaps included, misplaced the plume cells)
      IF ( plumeMaskFile .NE. ' ' ) THEN
        CALL READ_REC_XY_RL( plumeMaskFile, plumeMask,
     &                       iRec0, nIter0, myThid )
      ENDIF

C--   Runoff (if periodic forcing is used it is read elsewhere)
      IF ( .NOT. periodicExternalForcing ) THEN
        WRITE(msgBuf,'(1X,A)')
     &       'DOUBLEPLUME_INIT_FIXED: reading initial runoff'
        CALL PRINT_MESSAGE(msgBuf,standardMessageUnit,
     &       SQUEEZE_RIGHT,myThid)
        IF ( runoffVelFile .NE. ' ') THEN
          CALL READ_REC_XY_RL(
     &         runoffVelFile,runoffVel,iRec0,nIter0,myThid)
        ENDIF
        IF ( runoffRadFile .NE. ' ') THEN
          CALL READ_REC_XY_RL(
     &         runoffRadFile,runoffRad,iRec0,nIter0,myThid)
        ENDIF
      ENDIF

#ifdef ALLOW_PTRACERS
      IF ( ptracerMaskFile .NE. ' ' ) THEN
        CALL READ_REC_3D_RL
     &     (ptracerMaskFile,readBinaryPrec,PTRACERS_num,ptracerMask,
     &      iRec0,nIter0,myThid)
      ENDIF
#endif /* ALLOW_PTRACERS */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       CALL DOUBLEPLUME_DIAGNOSTICS_INIT(myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END
