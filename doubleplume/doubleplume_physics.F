#include "DOUBLEPLUME_OPTIONS.h"

C     *==========================================================*
C     | doubleplume_physics.F
C     | Physics shared by the inner and outer plume models. Each
C     | routine is a direct port of the doublePlume.py method named
C     | in its header; all constants come from DOUBLEPLUME.h.
C     *==========================================================*

C=======================================================================
CBOP
C     !ROUTINE: DP_AMBIENT
C     !INTERFACE:
      SUBROUTINE DP_AMBIENT( depth, Ta, Sa )
C     !DESCRIPTION: \bv
C     Ambient in-situ temperature and salinity at depth [m, > 0],
C     linearly interpolated between cell centres |rC(k)| (taProf and
C     saProf are cell-centre values), clamped to the top and to the
C     deepest wet level above the source (dp_nAmb).
C     Python: ambient_interp in getVars.
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION depth, Ta, Sa
      INTEGER k, n
      DOUBLE PRECISION w, z1, z2
CEOP
      n = MAX( 1, dp_nAmb )
      IF ( n .EQ. 1 .OR. depth .LE. ABS(rC(1)) ) THEN
        Ta = taProf(1)
        Sa = saProf(1)
        RETURN
      ENDIF
      IF ( depth .GE. ABS(rC(n)) ) THEN
        Ta = taProf(n)
        Sa = saProf(n)
        RETURN
      ENDIF
      DO k = 1, n-1
        z1 = ABS(rC(k))
        z2 = ABS(rC(k+1))
        IF ( depth .GE. z1 .AND. depth .LE. z2 ) THEN
          w  = (depth - z1)/(z2 - z1)
          Ta = (1.D0-w)*taProf(k) + w*taProf(k+1)
          Sa = (1.D0-w)*saProf(k) + w*saProf(k+1)
          RETURN
        ENDIF
      ENDDO
      Ta = taProf(n)
      Sa = saProf(n)
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_RHO
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_RHO( T, S, z )
C     !DESCRIPTION: \bv
C     In-situ density, UNESCO 1983 (T in-situ degC, S psu, z depth m).
C     (Python uses gsw.rho_t_exact; see notes in the package README.)
C     \ev
      IMPLICIT NONE
      DOUBLE PRECISION T, S, z
      DOUBLE PRECISION rho_0, gg, P
      DOUBLE PRECISION kw, Aw, Bw, k0, bulk_modulus
      DOUBLE PRECISION A, B, C, rho_w, rho_zero
      PARAMETER ( rho_0 = 1027.D0 )
      PARAMETER ( gg    = 9.81D0 )
CEOP
C     pressure in bar
      P  = rho_0*gg*ABS(z)*1.D-5
      kw = 19652.21D0 + 148.4206D0*T - 2.327105D0*T**2
     &   + 1.360477D-2*T**3 - 5.155288D-5*T**4
      Aw = 3.239908D0 + 1.43713D-3*T + 1.16092D-4*T**2
     &   - 5.77905D-7*T**3
      Bw = 8.50935D-5 - 6.12293D-6*T + 5.2787D-8*T**2
      k0 = kw + (54.6746D0 - 0.603459D0*T + 1.09987D-2*T**2
     &   - 6.1670D-5*T**3)*S
     &   + (7.944D-2 + 1.6483D-2*T - 5.3009D-4*T**2)*S**1.5D0
      A  = Aw + (2.2838D-3 - 1.0981D-5*T - 1.6078D-6*T**2)*S
     &   + 1.91075D-4*S**1.5D0
      B  = Bw + (-9.9348D-7 + 2.0816D-8*T + 9.1697D-10*T**2)*S
      bulk_modulus = k0 + A*P + B*P**2

      A  = 8.24493D-1 - 4.0899D-3*T + 7.6438D-5*T**2
     &   - 8.2467D-7*T**3 + 5.3875D-9*T**4
      B  = -5.72466D-3 + 1.0227D-4*T - 1.6546D-6*T**2
      C  = 4.8314D-4
      rho_w = 999.842594D0 + 6.793952D-2*T - 9.095290D-3*T**2
     &      + 1.001685D-4*T**3 - 1.120083D-6*T**4
     &      + 6.536336D-9*T**5
      rho_zero = rho_w + A*S + B*S**1.5D0 + C*S**2

      DP_RHO = rho_zero/(1.D0 - P/bulk_modulus)
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_PR_Z
C     !INTERFACE:
      SUBROUTINE DP_PR_Z( Tk, P, Z, AD, BD )
C     !DESCRIPTION: \bv
C     Peng-Robinson compressibility factor Z (vapour = largest root)
C     and the dimensionless A, B, for the gas at Tk [K], P [Pa].
C     Newton's method from above all roots (Cauchy bound); O2 is
C     supercritical at ocean temperatures, so the cubic has a single
C     real root. Python: gas_density / gas_SOL (np.roots).
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION Tk, P, Z, AD, BD
      INTEGER it
      DOUBLE PRECISION bp, ac, mu, alpha, ap
      DOUBLE PRECISION cc0, cc1, cc2, f, fp, dZ
CEOP
      bp    = 0.07780D0*DP_Ru*DP_TcrGas/DP_PcrGas
      ac    = 0.45724D0*DP_Ru**2*DP_TcrGas**2/DP_PcrGas
      mu    = 0.37464D0 + 1.54226D0*DP_omegaGas
     &      - 0.26992D0*DP_omegaGas**2
      alpha = ( 1.D0 + mu*(1.D0 - SQRT(Tk/DP_TcrGas)) )**2
      ap    = ac*alpha
      AD    = ap*P/(DP_Ru**2*Tk**2)
      BD    = bp*P/(DP_Ru*Tk)
C     Z^3 + cc2 Z^2 + cc1 Z + cc0 = 0
      cc2   = -(1.D0 - BD)
      cc1   = AD - 3.D0*BD**2 - 2.D0*BD
      cc0   = -(AD*BD - BD**2 - BD**3)
      Z     = 1.D0 + MAX( ABS(cc2), ABS(cc1), ABS(cc0) )
      DO it = 1, 100
        f  = ((Z + cc2)*Z + cc1)*Z + cc0
        fp = (3.D0*Z + 2.D0*cc2)*Z + cc1
        IF ( fp .EQ. 0.D0 ) GOTO 20
        dZ = f/fp
        Z  = Z - dZ
        IF ( ABS(dZ) .LT. 1.D-13*ABS(Z) ) GOTO 20
      ENDDO
   20 CONTINUE
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_GAS_RHO
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_GAS_RHO( T_loc, P )
C     !DESCRIPTION: \bv
C     Gas density [kg/m^3], Peng-Robinson EOS (Python: gas_density).
C     T_loc [degC], P [Pa].
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION T_loc, P
      DOUBLE PRECISION Tk, Z, AD, BD
CEOP
      Tk = T_loc + DP_T0
      CALL DP_PR_Z( Tk, P, Z, AD, BD )
      DP_GAS_RHO = DP_Mgas*P/(Z*DP_Ru*Tk)
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_GAS_SOL
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_GAS_SOL( T_loc, S_loc, P )
C     !DESCRIPTION: \bv
C     Gas solubility [kg/m^3] (Python: gas_SOL): Henry's law with
C     temperature, salting-out (O2, Weiss 1970) and pressure effects,
C     times the Peng-Robinson fugacity. T_loc [degC], P [Pa].
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION T_loc, S_loc, P
      DOUBLE PRECISION Tk, H0, nu_bar, dH_solR, Hc, Z, AD, BD
      DOUBLE PRECISION sq2, lnphi
CEOP
      Tk      = T_loc + DP_T0
      H0      = 1.3D-3*1000.D0/101325.D0*DP_Mgas
      nu_bar  = 32.D-6
      dH_solR = DP_dhSol/DP_Ru*DP_Mgas
      H0 = H0*EXP( dH_solR*(1.D0/Tk - 1.D0/298.15D0) )
      H0 = H0*EXP( S_loc*( -0.033096D0 + 0.014259D0*Tk/100.D0
     &                    - 0.0017000D0*(Tk/100.D0)**2 ) )
      Hc = H0*EXP( -(P - 101325.D0)*nu_bar/(DP_Ru*Tk) )
      CALL DP_PR_Z( Tk, P, Z, AD, BD )
      sq2   = SQRT(2.D0)
      lnphi = (Z - 1.D0) - LOG(Z - BD) - AD/(2.D0*sq2*BD)
     &      * LOG( (Z + (1.D0 + sq2)*BD)/(Z + (1.D0 - sq2)*BD) )
      DP_GAS_SOL = EXP(lnphi)*P*Hc
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_SLIP
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_SLIP( rho_a, rho_b, T_loc, de )
C     !DESCRIPTION: \bv
C     Bubble slip velocity [m/s] (Python: u_slip, Clift et al.).
C     de = bubble diameter [m], T_loc [degC].
C     NOTE: in the spherical regime Clift et al. use log10; the Python
C     model uses the natural log. Kept identical to Python.
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION rho_a, rho_b, T_loc, de
      DOUBLE PRECISION mu, mu_w, sigma, Mo, Eo, Hn, Jn, Re, Nd, lnNd
      DOUBLE PRECISION us
CEOP
      mu    = 1721.2D0*EXP(-0.0251D0*T_loc)/1.D6
      mu_w  = 0.0009D0
      sigma = (75.7060D0 - 0.1511D0*T_loc)/1000.D0
      Mo    = DP_g*mu**4*(rho_a - rho_b)/(rho_a**2*sigma**3)
      Eo    = DP_g*(rho_a - rho_b)*de**2/sigma
      Hn    = 4.D0/3.D0*Eo*Mo**(-0.149D0)*(mu/mu_w)**(-0.14D0)
      IF ( Hn .GT. 2.D0 ) THEN
        IF ( Hn .GT. 59.3D0 ) THEN
          Jn = 3.42D0*Hn**0.441D0
        ELSE
          Jn = 0.94D0*Hn**0.757D0
        ENDIF
        Re = Mo**(-0.149D0)*(Jn - 0.857D0)
      ELSE
        Nd   = 4.D0*rho_a*(rho_a - rho_b)*DP_g*de**3/(3.D0*mu**2)
        lnNd = LOG(Nd)
        IF ( Nd .LE. 73.D0 ) THEN
          Re = Nd/24.D0 - 1.7569D-4*Nd**2 + 6.9252D-7*Nd**3
     &       - 2.3027D-10*Nd**4
        ELSEIF ( Nd .LE. 580.D0 ) THEN
          Re = 10.D0**( -1.7095D0 + 1.33438D0*lnNd
     &                - 0.11591D0*lnNd**2 )
        ELSEIF ( Nd .LE. 1.55D7 ) THEN
          Re = 10.D0**( -1.81391D0 + 1.34671D0*lnNd
     &                - 0.12427D0*lnNd**2 + 0.006344D0*lnNd**3 )
        ELSE
          Re = 10.D0**( 5.33283D0 - 1.21728D0*lnNd
     &                + 0.19007D0*lnNd**2 - 0.007005D0*lnNd**3 )
        ENDIF
      ENDIF
      us = mu/(rho_a*de)*Re
      IF ( us .LT. 0.D0 ) us = 0.001D0
      DP_SLIP = us
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_PEEL
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_PEEL( ui, us, bi, rho_i, rho_a )
C     !DESCRIPTION: \bv
C     Continuous peeling flux Ep [m^2/s] (Python: cp_model,
C     Socolofsky et al. 2008, half cone). Only Ep < 0 is kept.
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION ui, us, bi, rho_i, rho_a, Ep
CEOP
      IF ( ui .LE. 0.D0 ) THEN
        DP_PEEL = 0.D0
        RETURN
      ENDIF
      Ep = DP_c2*DP_epsilon*(us/ui)**2*DP_g*(rho_a - rho_i)/DP_rhoR
     &   * DP_pi*bi**2/(2.D0*ui)
      DP_PEEL = MIN( Ep, 0.D0 )
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_SHERWOOD
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_SHERWOOD( T_loc, rho_a, us, de, Df )
C     !DESCRIPTION: \bv
C     Transfer coefficient K*Sh*Df/de [m/s] for a rigid sphere (Clift
C     et al. 1979); Df = diffusivity [m^2/s]. Used for mass transfer
C     (Python: gas_MT, Df = gas diffusivity) and heat transfer
C     (Python: gas_HT, Df = 1.46e-7).
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION T_loc, rho_a, us, de, Df
      DOUBLE PRECISION mu_w, Sc, Pe, Re, Sh
CEOP
      mu_w = 1721.2D0*EXP(-0.0251D0*T_loc)/1.D6
      Sc   = mu_w/(rho_a*Df)
      Pe   = us*de/Df
      Re   = rho_a*us*de/mu_w
      IF ( Re .LT. 100.D0 ) THEN
        Sh = 1.D0 + (1.D0 + 1.D0/Pe)**(1.D0/3.D0)*Re**0.41D0
     &     * Sc**(1.D0/3.D0)
      ELSEIF ( Re .LT. 2000.D0 ) THEN
        Sh = 1.D0 + 0.724D0*Re**0.48D0*Sc**(1.D0/3.D0)
      ELSE
        Sh = 1.D0 + 0.425D0*Re**0.55D0*Sc**(1.D0/3.D0)
      ENDIF
      DP_SHERWOOD = DP_Kgas*Sh*Df/de
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_GAS_DIFF
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_GAS_DIFF( T_loc )
C     !DESCRIPTION: \bv
C     O2 diffusivity in water [m^2/s], Wise & Houghton (1966)
C     (Python: gas_MT).
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION T_loc
CEOP
      DP_GAS_DIFF = 4.2D-2*EXP( -(4390.D0/0.23886D0)
     &            /(DP_Ru*(DP_T0 + T_loc)) )/100.D0**2
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_VOIDFRAC
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_VOIDFRAC( mb, rho_b, bi, us, ui )
C     !DESCRIPTION: \bv
C     Gas void fraction (Python: void_fraction).
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION mb, rho_b, bi, us, ui
CEOP
      DP_VOIDFRAC = mb/( rho_b*0.5D0*DP_pi*bi**2*DP_lambdaR**2
     &            *( us + 2.D0*ui/(1.D0 + DP_lambdaR**2) ) )
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_MELT
C     !INTERFACE:
      SUBROUTINE DP_MELT( Q, J, Hh, B, depth, mdot, Sb, Tb )
C     !DESCRIPTION: \bv
C     Ice-wall melt rate mdot [m/s], boundary salinity Sb and
C     temperature Tb (Python: melting, boundary_salinity,
C     boundary_temperature; Holland & Jenkins 1999).
C     Q, J, Hh, B = inner volume, momentum, heat, salt fluxes.
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION Q, J, Hh, B, depth, mdot, Sb, Tb
      DOUBLE PRECISION T, u, beta, bb, cc, disc, gs
CEOP
      T    = Hh/(DP_rhoR*DP_cp*Q) - DP_T0
      u    = J/Q
      gs   = DP_gammaS*SQRT(DP_Cd)*u
      beta = DP_cp*DP_gammaT*SQRT(DP_Cd)/DP_Li
      bb   = gs - beta*u*(T - DP_lam2 - DP_lam3*depth)
      cc   = beta*gs*u*( B/Q*DP_lam1 - T + DP_lam2 + DP_lam3*depth )
      disc = MAX( bb**2 - 4.D0*cc, 0.D0 )
      mdot = ( -bb + SQRT(disc) )/2.D0
      Sb   = gs*(B/Q)/(mdot + gs)
      Tb   = DP_lam1*Sb + DP_lam2 + DP_lam3*depth
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_WUEST_RES
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_WUEST_RES( u, us, mb, bi,
     &                                        rho_a, rho_b )
C     !DESCRIPTION: \bv
C     Froude-number residual for the Wuest et al. (1992) inner plume
C     initial velocity (Python: wuest_residual).
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION u, us, mb, bi, rho_a, rho_b
      DOUBLE PRECISION xi, rho_p, lr
CEOP
      lr    = DP_lambdaR
      xi    = mb/( 0.5D0*DP_pi*bi**2*rho_b*lr**2
     &      *( us + 2.D0*u/(1.D0 + lr**2) ) )
      rho_p = xi*rho_b + (1.D0 - xi)*rho_a
      DP_WUEST_RES = DP_FrWuest - u
     &      / SQRT( 2.D0*lr*bi*DP_g*(rho_a - rho_p)/DP_rhoR )
      RETURN
      END

C=======================================================================
CBOP
C     !ROUTINE: DP_WUEST
C     !INTERFACE:
      DOUBLE PRECISION FUNCTION DP_WUEST( us, mb, bi, rho_a, rho_b )
C     !DESCRIPTION: \bv
C     Wuest et al. (1992) initial velocity [m/s] (Python: wuest).
C     The residual decreases monotonically with u (positive as
C     u -> 0, -> -inf as u -> inf), so the root is unique and found
C     by bisection. Returns -1 if there is no bubble buoyancy.
C     \ev
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "DOUBLEPLUME.h"
      DOUBLE PRECISION us, mb, bi, rho_a, rho_b
      INTEGER it
      DOUBLE PRECISION uLo, uHi, uMid
      DOUBLE PRECISION DP_WUEST_RES
      EXTERNAL DP_WUEST_RES
CEOP
      IF ( mb .LE. 0.D0 .OR. bi .LE. 0.D0
     &     .OR. rho_a .LE. rho_b ) THEN
        DP_WUEST = -1.D0
        RETURN
      ENDIF
      uLo = 1.D-10
      uHi = 1.D0
      DO it = 1, 30
        IF ( DP_WUEST_RES(uHi,us,mb,bi,rho_a,rho_b) .LT. 0.D0 )
     &       GOTO 10
        uHi = 10.D0*uHi
      ENDDO
   10 CONTINUE
      DO it = 1, 200
        uMid = 0.5D0*(uLo + uHi)
        IF ( DP_WUEST_RES(uMid,us,mb,bi,rho_a,rho_b) .GT. 0.D0 ) THEN
          uLo = uMid
        ELSE
          uHi = uMid
        ENDIF
        IF ( uHi - uLo .LT. 1.D-13*uHi ) GOTO 20
      ENDDO
   20 CONTINUE
      DP_WUEST = 0.5D0*(uLo + uHi)
      RETURN
      END

C=======================================================================
C     Dummy Jacobian for DLSODE (not called with MF = 10, but DLSODE
C     needs an EXTERNAL in its argument list).
      SUBROUTINE DP_JAC( NEQ, T, Y, ML, MU, PD, NRPD )
      IMPLICIT NONE
      INTEGER NEQ, ML, MU, NRPD
      DOUBLE PRECISION T, Y(*), PD(NRPD,*)
      RETURN
      END
