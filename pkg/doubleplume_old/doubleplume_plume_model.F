#include "DOUBLEPLUME_OPTIONS.h"

CBOP
C     !ROUTINE: DOUBLEPLUME_PLUME_MODEL
C     !INTERFACE:
      SUBROUTINE DOUBLEPLUME_PLUME_MODEL (mythid)
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE DOUBLEPLUME_PLUME_MODEL:
C     | o Uses ODEPACK to solve one of 3 models describing subglacial
C     |   runoff driven plumes
C     | 
C     *===========================================================
c	\ev
C     !USES:
      IMPLICIT NONE
C     == Global variables ==

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"
#include "DOUBLEPLUME_OPTIONS.h"


C     !INPUT/OUTPUT PARAMETERS:

C     ==Local variables==

C     Plume models
      EXTERNAL DOUBLEPLUME_INNER, JENKINS, HALFCONE, DETACHED_CONE, JEX

C     For ODEPACK
      INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(20), LIW, LRW,
     &         MF, NEQ
      _RL T, TOUT

      _RL  RWORK(116), Y(4)
#endif /* DOUBLEPLUME_ALLOW_DETACHED_PLUME */

C     Y is input/output vector for DLSODE
C       Y(1) = plume volume flux
C       Y(2) = plume momentum flux
C       Y(3) = plume heat flux
C       Y(4) = plume salt flux
 

C     Other local variables
      INTEGER  k, mythid
      DOUBLE PRECISION RHO
      DOUBLE PRECISION temperature, salinity, depth
      DOUBLE PRECISION tambient, sambient
      DOUBLE PRECISION rho_plume, rho_ambient
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifdef ALLOW_DOUBLEPLUME

C For ODEPACK solver. See ODEPACK documentation for definitions.
      NEQ = 4
      LRW = 116
      LIW = 116
      ITOL = 1
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      MF = 10
      IWORK(7) = 2 ! To limit number of times repeat error messages are printed

C Initial conditions
      Y(1) = 20.D0
      Y(2) = 10.D0          ! initial vertical velocity
      Y(3) = 0.D0          ! initial temperature
      Y(4) = 0.D0          ! initial salinity

C Prepare profiles
      DO K = 1,Nr+1
        zProf(k) = rF(k)
        zProfAbs(k) = abs(rF(k))
        QiProf(k) = 0.D0
        MiProf(k) = 0.D0
        HiProf(k) = 0.D0
        BiProf(k) = 0.D0
C        aProf(k) = 0.D0
C        mIntProf(k) = 0.D0
      ENDDO

C Start at bottom of ice face
      T = iceDepth

C Next point at which to retrieve values
      TOUT = zProf(iceDepthK-1)

C Clear cells below depth of ice face
      DO k = Nr+1,iceDepthK+1
        QiProf(k) = 0.D0
        MiProf(k) = 0.D0
        HiProf(k) = 0.D0
        BiProf(k) = 0.D0
C        aProf(k)  = 0.D0
C        mIntProf(k) = 0.D0
      ENDDO
   
C Set intial conditions
      QiProf(iceDepthK) = Y(1)
      MiProf(iceDepthK) = Y(2)
      HiProf(iceDepthK) = Y(3)
      BiProf(iceDepthK) = Y(4)

C      print *, "In plume model."
C      print *, "Initial conditions = ", Y(1), Y(2), Y(3), Y(4)

C Move up through water column from lowest layer

      DO IOUT = iceDepthK-1,1,-1
      
C Check to make sure plume hasn't reached neutral buoyancy in a lower layer
        IF ( ISTATE .GT. -1 ) THEN 


            CALL DLSODE (DOUBLEPLUME_INNER,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,
     &               ITASK,ISTATE, IOPT, RWORK, LRW, IWORK, LIW, 
     &               JEX, MF)


C Test to see if neutral buoyancy has now been reached.
C If solver returns ISTATE = -1, then it has been unable to meet required tolerances
C at this level. This generally occurs because plume has reached neutral buoyancy and
C run out of momentum, and so is no longer rising. At this point, we therefore end
C the call to the plume model.
C Our aim is to catch the plume at the point of neutral buoyancy. We therefore perform
C a manual comparrison of ambient and plume density. If plume density >= ambient density
C we assign ISTATE = -1, again ending the call to the plume model.

          temperature = Y(3) / Y(1) ! Hi / Qi
          salinity = Y(4) / Y(1) ! Bi / Qi
          depth = T

C Calculate plume density
          rho_plume = RHO(temperature,salinity,depth)
C          print *, "rho_plume =", rho_plume

          IF (IOUT.EQ.1) THEN
            tambient = taprof(1)
            sambient = saprof(1)
          ELSE
            tambient = .5*(taprof(IOUT-1)+taprof(IOUT))
            sambient = .5*(saprof(IOUT-1)+saprof(IOUT))
          ENDIF ! (IOUT.EQ.1)

C Calculate ambient density
          rho_ambient = RHO(2.D0,34.D0,depth)
C          print *, "rho_ambient = ", rho_ambient

          IF ( rho_plume .GT. rho_ambient ) ISTATE = -1

C If ISTATE is now < 0, then plume has reached neutral buoyancy 
          IF ( ISTATE .LT. 0 ) THEN

C If we have reached neutral buoyancy then there is no volume flux out of this cell,
C so plume area and velocity equal zero.
C Other values are kept for use in determining plume outflow properties.

            Y(1) = 0.D0
            Y(2) = 0.D0

          ELSE ! ( ISTATE .LT. 0 )

C If the plume has not reached neutral buoyancy, then we assign a depth at which to calculate the next value and loop round to call the plume model again.
C Make sure we're not at the surface

            IF ( iout .NE. 1 ) THEN

              T=TOUT                ! Define present depth
              TOUT = zProf(IOUT-1)  ! Define next depth

            ENDIF ! ( iout .NE. 1 )
          ENDIF ! ( ISTATE .LT. 0 )
          

C This section is entered once the plume has reached neutral buoyancy
        ELSE ! ( ISTATE .GT. -1 )

C Once plume has reached neutral buoyancy, no plume values
          Y(1) = 0.0
          Y(2) = 0.0
          Y(3) = 0.0
          Y(4) = 0.0

        ENDIF ! ( ISTATE .GT. -1 )

C Save results
C       print *, "Y(2)/Y(1) = ", Y(2) / Y(1)
C       print *, "Y(4)/Y(1) = ", Y(4) / Y(1)
        QiProf(iout) = Y(1)
        MiProf(iout) = Y(2)
        HiProf(iout) = Y(3)
        BiProf(iout) = Y(4)
C        aProf(iout) = Y(5)
C        mIntProf(iout) = Y(6)

      ENDDO ! IOUT = iceDepthK-1,1,-1
C      print *, "Succesful inner plume run."
C      print *, "Qi = ", QiProf

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE DOUBLEPLUME_INNER (NEQ, T, Y, YDOT)
      
      implicit none

#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"


      INTEGER NEQ
      DOUBLE PRECISION T, Y(4), YDOT(4)

      _RL Tambient, Sambient, Tin, rho_a
      DOUBLE PRECISION m_dot, Sb, Tb
      DOUBLE PRECISION To, So, Uo, bo
      DOUBLE PRECISION Ti, Si, Ui, bi, rho_i, rho_m
      DOUBLE PRECISION a, b, c, beta_loc
      DOUBLE PRECISION E_i, E_o, E_p
      DOUBLE PRECISION RHO
C      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846

C Y(1) = Qi
C Y(2) = Mi
C Y(3) = Hi
C Y(4) = Bi
C RPAR(1) = Qo
C RPAR(2) = Mo
C RPAR(3) = Ho
C RPAR(4) = Bo

C Interpolate from imposed ambient profiles
      IF ( abs(T) .LE. zProfAbs(1) ) THEN
       Tambient = taProf(1)
       Sambient = saProf(1)
      ELSEIF ( abs(T) .GE. zProfAbs(Nr) ) THEN
       Tambient = taProf(Nr)
       Sambient = saProf(Nr)
      ELSE
       CALL linint(Nr,zProfAbs,taProf,abs(T), Tambient)
       CALL linint(Nr,zProfAbs,saProf,abs(T), Sambient)
      ENDIF
      
	Tin = T

      Ti = Y(3) / Y(1)
      Si = Y(4) / Y(1)
      Ui = Y(2) / Y(1)
      bi = Y(1) / SQRT(PI * Y(2))

C      To = RPAR(3) / RPAR(1)
C      So = RPAR(4) / RPAR(1)
C      Uo = RPAR(2) / RPAR(1)
C      bo = RPAR(1) / SQRT(PI * RPAR(2))


C      IF (Uo .EQ. 0.D0) THEN
      To = 2.D0 !Tambient
      So = 34.D0 !Sambient
      Uo = 0.D0
      bo = bi
C      ENDIF

      rho_i = RHO(Ti, Si, abs(Tin))
      rho_a = RHO(To, So, abs(Tin))
      rho_m = 0.002 * 1.204 + (1-0.002) * rho_i

      beta_loc = 4000 * 0.02 * SQRT(0.001) / 334000
      a = 1
      
      b = 0.00062 * SQRT(0.001) * Ui - beta_loc * Ui * 
     &         (Ti - 0.08 - -0.000761 * abs(Tin))
    
      c = beta_loc * 0.00062 * SQRT(0.001) * Ui**2 * 
     &         (Si * -0.0573 - Ti + 0.08 + -0.000761 * abs(Tin))
    
      m_dot = (-b + SQRT(b**2 - 4 * a * c)) / (2 * a)

      Sb = 0.00062 * SQRT(0.001) * Ui * Si                  
     &        / (m_dot + 0.00062 * SQRT(0.001) * Ui)
      
      Tb = -0.0573 * Sb + 0.08 + -0.000761 * abs(Tin)

      E_i = 2 * PI * bi * 0.055 * (Ui)
      E_o = 0.D0 !- 2 * PI * bi * 0.110 * Uo
      E_p = - 0.015 * (0.1 / Ui)**2 * 
     &       (9.81 * ((rho_a - rho_i) / 1025.0) * Ui * PI * bi**2)
     &         / Ui**2
      
      IF (E_p .LT. 0) THEN
       E_p = 0.D0
      ELSEIF (Ui .LT. 0.25) THEN
       E_p = 0.D0
      ELSE
       E_p = - 0.015 * (0.1 / Ui)**2 * 
     &       (9.81 * ((rho_a - rho_i) / 1025.0) * Ui * PI * bi**2)
     &         / Ui**2
      ENDIF


      YDOT(1) = E_i - E_o - E_p + b * m_dot

      YDOT(2) = PI * 9.81 * bi**2 * (1**2 * (rho_a - rho_m)) /
     &         (1.1 * 1025) + E_i * Uo - E_o * Ui - 
     &          E_p * Ui - 0.001 * bi * Ui**2
      
      YDOT(3) = E_i * To - E_o * Ti - E_p * Ti + bi * m_dot * Tb 
     &        - 4 * 0.02 * SQRT(0.001) * bi * Ui * (Ti - Tb)

      YDOT(4) = E_i * So - E_o * Si - E_p * Si + bi * m_dot * Sb 
     &        - 4 * 0.00062 * SQRT(0.001) * bi * Ui * (Si - Sb)

      RETURN
      END

C---------------------------------------------------

      SUBROUTINE  JENKINS (NEQ, T, Y, YDOT)

#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      INTEGER  NEQ
      DOUBLE PRECISION T, Y(7), YDOT(6)
      _RL Tambient, Sambient, rho_0, rho_1, Tin
      DOUBLE PRECISION mdot, Sb, Tb
      DOUBLE PRECISION a,b,c
      DOUBLE PRECISION RHO

#ifdef ALLOW_DOUBLEPLUME

C     Y(1) = r
C     Y(2) = w
C     Y(3) = T
C     Y(4) = S

      
C Interpolate from imposed ambient profiles

      IF ( abs(T) .LE. zProfAbs(1) ) THEN
       Tambient = taProf(1)
       Sambient = saprof(1)
      ELSEIF ( abs(T) .GE. zProfAbs(Nr) ) THEN
       Tambient = taProf(Nr)
       Sambient = saprof(Nr)
      ELSE
       CALL linint(Nr,zProfAbs,taProf,abs(T),Tambient)
       CALL linint(Nr,zProfAbs,saprof,abs(T),Sambient)
      ENDIF
           
	  Tin = T
	  
C Plume density
      rho_1   = RHO(Y(3),Y(4),Tin)
C Ambient density
      rho_0 = RHO(Tambient,Sambient,Tin)

C Equations for Sb, Tb and mdot

      a = lambda1*(GamT*c_w-GamS*c_i)

      b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*Tin+
     &         iceTemp-(L/c_i))
     &        -GamT*c_w*(Y(3)-lambda2-lambda3*Tin)

      c = GamS*Y(4)*(c_i*(lambda2+lambda3*Tin-iceTemp)+L)

      Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
      Tb   = lambda1*Sb+lambda2+lambda3*Tin !Tb
      mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb ! mdot

C Differential equations
C Plume thickness
      YDOT(1)=2*E_0+Cd-(g*Y(1)/(Y(2)**2))*(rho_0-rho_1)
     &  /rho_ref+2*mdot/Y(2)

C Plume vertical velocity
      YDOT(2)=-(Y(2)/Y(1))*(E_0+Cd+mdot/Y(2))
     &  +(g/Y(2))*(rho_0-rho_1)/rho_ref

C Plume temperature
       YDOT(3)=E_0*Tambient/Y(1)-(Y(3)/Y(1))
     &   *(E_0+mdot/Y(2))+(mdot/(Y(1)*Y(2)))
     &   *(Tb-(L/c_w)-(c_i/c_w)*(Tb-iceTemp))

C Plume salinity
       YDOT(4)=E_0*Sambient/Y(1)-(Y(4)/Y(1))
     &  *(E_0+mdot/Y(2));

C along-plume integrated melt rate and contact area  - 
      YDOT(5) = Y(7)
      YDOT(6) = Y(7) * mdot

      T = Tin

#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END

C - ------------------------------------------------------------------

      SUBROUTINE  HALFCONE (NEQ, T, Y, YDOT)

#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      INTEGER  NEQ
      DOUBLE PRECISION T, Y(7), YDOT(6), Tin
      _RL Tambient, Sambient, rho_0, rho_1
      DOUBLE PRECISION mdot, Sb, Tb
      DOUBLE PRECISION a,b,c
      DOUBLE PRECISION RHO

#ifdef ALLOW_DOUBLEPLUME

C     Y(1) = r
C     Y(2) = w
C     Y(3) = T
C     Y(4) = S

     
C     Interpolate from imposed ambient profiles

      IF ( abs(T) .LE. zProfAbs(1) ) THEN
       Tambient = taProf(1)
       Sambient = saProf(1)
      ELSEIF ( abs(T) .GE. zProfAbs(Nr) ) THEN
       Tambient = taProf(Nr)
       Sambient = saProf(Nr)
      ELSE
       CALL linint(Nr,zProfAbs,taProf,abs(T),Tambient)
       CALL linint(Nr,zProfAbs,saProf,abs(T),Sambient)
      ENDIF

      Tambient = 2.D0
      Sambient = 34.D0

      Tin = T
      rho_1   = RHO(Y(3),Y(4),Tin)
      rho_0 = RHO(Tambient,Sambient,Tin)
      
C      print *, "In Halfcone. rho = ", rho_1, rho_0

C     Equations for Sb, Tb and mdot

      a = lambda1*(GamT*c_w-GamS*c_i)

      b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*Tin+
     &         iceTemp-(L/c_i))
     &        -GamT*c_w*(Y(3)-lambda2-lambda3*Tin)

      c = GamS*Y(4)*(c_i*(lambda2+lambda3*Tin-iceTemp)+L)

      Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
      Tb   = lambda1*Sb+lambda2+lambda3*Tin !Tb
      mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb ! mdot
      
      print *, "GamS = ", GamS
      print *, "Tin = ", Tin 
      print *, "iceTemp = ", iceTemp
      print *, "Y(4) = ", Y(4)
      print *, "L = ", L 
      print *, "lambda3 = ", lambda3
      YDOT(1) = 2.*E_0+4.*mdot/(pi*Y(2))-Y(1)*g*
     &          (rho_0-rho_1)/(2.*Y(2)*Y(2)*rho_ref)+2.*Cd/pi

      YDOT(2) = -2.*E_0*Y(2)/Y(1)-4.*mdot/(pi*Y(1))+g*
     &          (rho_0-rho_1)/(Y(2)*rho_ref)-4.*Cd*Y(2)/(pi*Y(1))

      YDOT(3) = 2.*E_0*(TAMBIENT-Y(3))/Y(1)+4.*mdot*
     &           (Tb-Y(3))/(pi*Y(1)*Y(2))-4.*
     &           GamT*(Cd**0.5)*(Y(3)-Tb)/(pi*Y(1))

      YDOT(4) = 2.*E_0*(Sambient-Y(4))/Y(1)+4.*mdot*
     &           (Sb-Y(4))/(pi*Y(1)*Y(2))-4.*
     &           GamS*(Cd**0.5)*(Y(4)-Sb)/(pi*Y(1))

      YDOT(5) = 2.*Y(1)
      YDOT(6) = 2.*Y(1)*mdot

      print *, "YDOT(1) = ", YDOT(1)
      print *, "YDOT(3) = ", YDOT(3)
      T = Tin
#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END

C - ------------------------------------------------------------------

#ifdef DOUBLEPLUME_ALLOW_DETACHED_PLUME

      SUBROUTINE  DETACHED_CONE (NEQ, T, Y, YDOT)

#include "DOUBLEPLUME_OPTIONS.h"
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "DOUBLEPLUME.h"

      INTEGER  NEQ
      DOUBLE PRECISION T, Y(9), YDOT(8), Tin
      _RL Tambient, Sambient, rho_0, rho_1
      _RL Wplume, Cplume
      DOUBLE PRECISION mdot, Sb, Tb
      DOUBLE PRECISION a,b,c
      DOUBLE PRECISION RHO

      EXTERNAL PLUME_GEOMETRY

#ifdef ALLOW_DOUBLEPLUME

C     Y(1) = r
C     Y(2) = w
C     Y(3) = T
C     Y(4) = S

     
C     Interpolate from imposed ambient profiles

      IF ( abs(T) .LE. zProfAbs(1) ) THEN
       Tambient = taProf(1)
       Sambient = saProf(1)
      ELSEIF ( abs(T) .GE. zProfAbs(Nr) ) THEN
       Tambient = taProf(Nr)
       Sambient = saProf(Nr)
      ELSE
       CALL linint(Nr,zProfAbs,taProf,abs(T),Tambient)
       CALL linint(Nr,zProfAbs,saProf,abs(T),Sambient)
      ENDIF

      Tin = T
      rho_1   = RHO(Y(3),Y(4),Tin)
      rho_0 = RHO(Tambient,Sambient,Tin)

C     Equations for Sb, Tb and mdot

      a = lambda1*(GamT*c_w-GamS*c_i)

      b = GamS*c_i*(lambda1*Y(4)-lambda2-lambda3*Tin+
     &         iceTemp-(L/c_i))
     &        -GamT*c_w*(Y(3)-lambda2-lambda3*Tin)

      c = GamS*Y(4)*(c_i*(lambda2+lambda3*Tin-iceTemp)+L)

      Sb   = (1./(2.*a))*(-b-((b**2.-4.*a*c)**0.5)) !Sb
      Tb   = lambda1*Sb+lambda2+lambda3*Tin !Tb
      mdot = GamS*(Cd**0.5)*Y(2)*(Y(4)-Sb)/Sb ! mdot


      CALL PLUME_GEOMETRY(Y(1),Y(8),Cplume,Wplume)


      YDOT(1) = 2*E_0*Cplume/sin(Y(7))+
     &          2*Wplume*mdot/(sin(Y(7))*Y(2))-
     &          Y(1)*gravity*(rho_0-rho_1)/(Y(2)*Y(2)*rho_ref)+
     &          Cd*Wplume

      YDOT(2) = -E_0*Cplume*Y(2)/Y(1)/sin(Y(7))-
     &          Wplume*mdot/(Y(1)*sin(Y(7)))+
     &          gravity*(rho_0-rho_1)/(Y(2)*rho_ref)-
     &          Cd*Wplume*Y(2)/Y(1)

      YDOT(3) = (E_0*Cplume*(Tambient-Y(3))/Y(1)+
     &           Wplume*mdot*(Tb-Y(3))/(Y(1)*Y(2))-
     &           GamT*sqrt(Cd)*Wplume*(Y(3)-Tb)/(Y(1)))
     &          /sin(Y(7))


      YDOT(4) = (E_0*Cplume*(Sambient-Y(4))/Y(1)+
     &           Wplume*mdot*(Sb-Y(4))/(Y(1)*Y(2))-
     &           GamS*sqrt(Cd)*Wplume*(Y(4)-Sb)/(Y(1)))
     &          /sin(Y(7))

      YDOT(5) = Wplume
      YDOT(6) = Wplume*mdot

      YDOT(7) = (Y(2)*Y(2)*YDOT(1)+2*Y(2)*Y(1)*YDOT(2)) /
     &          (Y(2)*Y(2)*Y(1))/TAN(Y(7))
      YDOT(8) = 1.0/(TAN(Y(7)))

      T = Tin


#endif /* ALLOW_DOUBLEPLUME */

      RETURN
      END

#endif /* DOUBLEPLUME_ALLOW_DETACHED_PLUME */

C----------------------------------------------------------------------
C Dummy routine for ODEPACK. Necessary for Jacobian matrix if stiff ODEs.

      SUBROUTINE jex()
      RETURN
      END

C ---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RHO(t,S,z)

C Equation of state (UNESCO 1983)

C     T = temperature (deg C)
C     S = salinity (PSU)
C     z = depth (m)

      DOUBLE PRECISION T,S,z
      DOUBLE PRECISION rho_0, g, P
      DOUBLE PRECISION kw, Aw, Bw, k0
      DOUBLE PRECISION bulk_modulus
      DOUBLE PRECISION A, B, C, rho_w,rho_zero

      PARAMETER(rho_0=1027)
      PARAMETER(g=9.81)

      P= rho_0*g*abs(z)*1.0E-5
   
C ... RHO_1 (in situ) (saved in Y(7))
      kw= 19652.21+ 148.4206*T- 2.327105*T**2+ 
     &    1.360477e-2*(T**3)-5.155288e-5*(T**4)
      Aw= 3.239908+ 1.43713e-3*T+ 1.16092e-4*T**2-
     &    5.77905e-7*T**3
      Bw= 8.50935e-5- 6.12293e-6*T + 5.2787e-8*(T**2)
      k0= kw + (54.6746- 0.603459*T+ 1.09987e-2*(T**2)
     &    -6.1670e-5*(T**3))*S +(7.944e-2 + 1.6483e-2*
     &    T- 5.3009e-4*(T**2))*(S**1.5)
      A=  Aw+ (2.2838e-3- 1.0981e-5*T- 1.6078e-6*(T**2))
     &    *S+ 1.91075e-4*(S**1.5)
      B= Bw+ (-9.9348e-7+ 2.0816e-8*T+ 9.1697e-10*T**2)*S
      bulk_modulus= k0+ A*P+ B*P**2

      A= 8.24493e-1- 4.0899e-3*T+ 7.6438e-5*T**2-
     &   8.2467e-7*T**3+5.3875e-9*T**4
      B= -5.72466e-3 + 1.0227e-4*T- 1.6546e-6*T**2
      C= 4.8314e-4
      rho_w= 999.842594 + 6.793952e-2*T- 9.095290e-3*T**2+
     &       1.001685e-4*T**3-1.120083e-6*T*
     &       *4+ 6.536336e-9*T**5
      rho_zero= rho_w+ A*S + B*(S**1.5)+ C*(S**2)

      RHO= rho_zero/(1- (P/bulk_modulus))

      RETURN
      END

C-------------------------------------------------------------


      subroutine linint(nx,xtab,ytab,x,y)
c
c     Given a value of x return a value of y based on interpolation
c     within a table of y values (ytab) corresponding to the x values
c     contained in the array xtab.  The subroutine assumes that the
c     values in xtab increase monotonically
c
c    John Mahaffy 2/12/95
c    Modified slightly TRC 2014

      integer nx
      double precision xtab(nx), ytab(nx), x, y

c     local variables
      integer i, i1
      double precision  wx

      if (x.lt.(xtab(1)).or.x.GT.(xtab(nx))) then
         write(6,*) 'x = ', x, '  is out of table range'
         stop
      endif
      do 100 i=2,nx
         if (x.le.xtab(i)) go to 200
  100    continue
  200 i1=i-1
      wx=(x-xtab(i1))/(xtab(i1+1)-xtab(i1))
      y=(1-wx)*ytab(i1)+wx*ytab(i1+1)
      return
      end

C-------------------------------------------------------------

#ifdef DOUBLEPLUME_ALLOW_DETACHED_PLUME
      SUBROUTINE PLUME_GEOMETRY(A,X,C,W)

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

      _RL A   ! input: plume cross-sectional area
      _RL X   ! input: distance of center from wall
      _RL C   ! output: perimeter exposed to ambient
      _RL W   ! output: perimeter exposed to wall

      _RL B, f, fp, tol, adj
      INTEGER ITERMAX, I
      CHARACTER*(MAX_LEN_MBUF) msgBuf

C ---------------------------------------------------

      B = sqrt(A/pi)
    
      IF (B.le.X) THEN

!      PLUME IS FULLY DETACHED    

       C = 2*PI*B
       W = 0

      ELSE
     
!      PLUME IS PARTLY OR NON-DETACHED (half-conical)
!      FIND SHAPE CONSISTENT WITH X AND A

       TOL = 1.e-12
       ITERMAX = 300

       DO I=1,ITERMAX
        F = B*B * (PI-acos(X/B)) + X*(B*B-X*X)**0.5 - A
        FP = 2*B*(PI - acos(X/B)) - X/(1 - X*X/(B*B))**0.5 + 
     &   (X*B)/(B*B - X*X)**0.5
        ADJ = f/fp
        B = B - ADJ
        if (B .LT. X) then
            B=X+.001
        endif
       if ( ABS(ADJ) .LT. TOL) GOTO 20
       enddo

   20  CONTINUE

       if (I.eq.ITERMAX) then
                WRITE(msgBuf,'(2A)')
     &          'cannot use detaching plume without ',
     &          '#define DOUBLEPLUME_ALLOW_DETACHED_PLUME'
                CALL PRINT_ERROR( msgBuf, 1 )
                STOP 'ABNORMAL END: S/R DOUBLEPLUME_PLUME_MODEL'
       endif
    
       W = 2 * sqrt(B*B-X*X)
       C = 2 * (pi-acos(X/b)) * B
     
      ENDIF

      return
      end
#endif /* DOUBLEPLUME_ALLOW_DETACHED_PLUME */

