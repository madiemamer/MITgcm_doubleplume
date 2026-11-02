
      SUBROUTINE DOUBLEPLUME_MELTRATE(
     I           temperature, salinity, velocity, depth,
     O           meltrate )

      IMPLICIT NONE

C Global variables
#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DOUBLEPLUME.h"

C Local variables

      _RL temperature
      _RL salinity
      _RL velocity
      _RL meltrate
      _RL depth
      _RL absVelocity
      _RL a,b,c, tb, sb, mdot

#ifdef ALLOW_DOUBLEPLUME

C Routine can't cope with zero velocity. Unlikely to occur anyway with currents,
C waves, convection etc. This isn't very physical, but will do for now.
      IF ( velocity .LT. backgroundVel ) velocity = backgroundVel

      absVelocity = abs(velocity)

C Calculate melt rate from 3 equation formualtion (as for plume models)
C Constants from DOUBLEPLUME_PHYS; depth is POSITIVE (m) since
C DP_lam3 < 0 (pass ABS(rC(k)), not rC(k)).

C     Equations for Sb, Tb and mdot

      a = DP_lam1*(DP_gammaT*DP_cp-DP_gammaS*c_i)

      b = DP_gammaS*c_i*(DP_lam1*salinity-DP_lam2-DP_lam3*depth+
     &         iceTemp-(DP_Li/c_i))
     &        -DP_gammaT*DP_cp*(temperature-DP_lam2-DP_lam3*depth)

      c = DP_gammaS*salinity*(c_i*(DP_lam2+DP_lam3*depth-iceTemp)
     &    +DP_Li)

      Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
      Tb   = DP_lam1*Sb+DP_lam2+DP_lam3*depth !Tb
      mdot = DP_gammaS*(DP_Cd**0.5)*absVelocity*(salinity-Sb)/Sb ! mdot

C Meltrate is in m d^-1
       meltrate = 86400.*mdot

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END
