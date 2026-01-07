
#include "DOUBLEPLUME_OPTIONS.h"

      SUBROUTINE DOUBLEPLUME_CHECK( myThid )
C     /==========================================================\
C     | SUBROUTINE DOUBLEPLUME_CHECK                                |
C     | o Validate basic package setup and inter-package         |
C     | dependencies.                                            |
C     \==========================================================/
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DOUBLEPLUME.h"

C     === Routine arguments ===
C     myThid -  Number of this instance of DOUBLEPLUME_CHECK
      INTEGER myThid

#ifdef ALLOW_DOUBLEPLUME

C     === Local variables ===
C     msgBuf      - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf

      _BEGIN_MASTER(myThid)

      WRITE(msgBuf,'(A)') 'DOUBLEPLUME_CHECK: #define ALLOW_DOUBLEPLUME'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &     SQUEEZE_RIGHT, myThid )

C    INSERT OTHER CHECKS HERE AS NECESSARY

       IF ( selectAddFluid .NE. 1 ) THEN
       WRITE(msgBuf,'(2A)')
     &  'Cannot use DOUBLEPLUME if selectAddFluid not equal to 1'
                CALL PRINT_ERROR( msgBuf, 1 )
                STOP 'ABNORMAL END: S/R DOUBLEPLUME_CHECK'
      ENDIF

      _END_MASTER(myThid)

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END
