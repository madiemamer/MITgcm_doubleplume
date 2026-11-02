#include "DOUBLEPLUME_OPTIONS.h"

C !INTERFACE:
      SUBROUTINE DOUBLEPLUME_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "DOUBLEPLUME.h"



C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP
#ifdef ALLOW_DOUBLEPLUME
#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      diagName  = 'icefrntW'
      diagTitle = 'plume vertical velocity'
      diagUnits = 'm/s        '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrntT'
      diagTitle = 'plume temperature'
      diagUnits = 'Deg. C           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrntS'
      diagTitle = 'plume salinity'
      diagUnits = 'PSU        '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrntR'
      diagTitle = 'plume thickness'
      diagUnits = 'm           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtQi'
      diagTitle = 'inner plume vol flux'
      diagUnits = 'm^3 s^-1        '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtMi'
      diagTitle = 'inner plume momentum flux'
      diagUnits = 'm^4 s^-2         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtHi'
      diagTitle = 'inner plume heat flux'
      diagUnits = 'J s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtBi'
      diagTitle = 'inner plume salt flux'
      diagUnits = 'm^3 ppt s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtCi'
      diagTitle = 'inner plume dissolved gas flux'
      diagUnits = 'kg s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtQo'
      diagTitle = 'outer plume vol flux'
      diagUnits = 'm^3 s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtMo'
      diagTitle = 'outer plume momentum flux'
      diagUnits = 'm^4 s^-2         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtHo'
      diagTitle = 'outer plume heat flux'
      diagUnits = 'J s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtBo'
      diagTitle = 'outer plume salt flux'
      diagUnits = 'm^3 ppt s^-1         '
      diagCode  = 'SM      MR      '            ! Change R to X for custom levels
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C      CALL DIAGNOSTICS_SETKLEV( diagName, plumeLayers , myThid )

      diagName  = 'icefrtHg'
      diagTitle = 'inner plume bubble heat flux'
      diagUnits = 'J s^-1          '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'icefrtMb'
      diagTitle = 'inner plume bubble mass flux'
      diagUnits = 'kg s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'icefrtCo'
      diagTitle = 'outer plume dissolved gas flux'
      diagUnits = 'kg s^-1         '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'icefrtMe'
      diagTitle = 'plume melt rate of the ice face'
      diagUnits = 'm d^-1          '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'icefrtVX'
      diagTitle = 'plume-ocean volume exchange (+ = added, unramped)'
      diagUnits = 'm^3 s^-1        '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END
