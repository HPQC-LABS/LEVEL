c***********************************************************************
c Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
c***********************************************************************
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	SUBROUTINE SCATLENGHT_MAPTAN 
     &           (Mes, U, Scale, Rmin, C4, Rbar, alpha, SL, SL_err)
C----------------------------------------------------------------------
C                            Input variables 
C----------------------------------------------------------------------
	INTEGER              Mes 
C     Mes is the munber of grid points; it must be an even integer number
	DOUBLE PRECISION     U, Scale, Rmin, C4
C     U is the potential defined as external subroutine-function U(r) where
C     r is the radial coordinate defined on the interval [Rmin,infinity] 
C     it is assumed that U(r->infinity)->0
C---------------------------------------------------------------------
C     Rmin is the left boundary point
C     Scale is the scaled factor of the energy: hbar^2/2mu
C     C4 is the long range coefficient for the C4/r^4 potential
C     C4 should zero for the Cn/r^n potential where n>4
C----------------------------------------------------------------------
C            Mapping parameters of the function
C            y(r) = 2*arctan[alpha*(r/Rbar-1)]/Pi 
C----------------------------------------------------------------------
	DOUBLE PRECISION     Rbar, alpha
C     Rbar>Rmin should be close to eqvilibrium distance
C     alpha should be close to the magic number n/2-1
C----------------------------------------------------------------------
C							Output variables  
C----------------------------------------------------------------------
	DOUBLE PRECISION     SL, SL_err
C     SL is the scattering lenght
C     SL_err is the estimate of absolute error by Richardson method
C----------------------------------------------------------------------
	DOUBLE PRECISION     Pi, Pi2, Pi4
	PARAMETER           ( Pi=3.1415926535897932384626433832795D0 )
	PARAMETER           ( Pi2 = Pi/2.d0, Pi4 = Pi2**2 )
C
	INTEGER              Mes_M, Mes_M2, j, k
	PARAMETER            (Mes_M = 10000000, Mes_M2 = Mes_M/2)
	DOUBLE PRECISION     Q(0:Mes_M), Q2(0:Mes_M2)
C
	DOUBLE PRECISION     H, Ymin, drdy, r, A, B 
	DOUBLE PRECISION     w, H3, H6, z
	DOUBLE PRECISION     Yl, Yl2, SL2
C  ************************** Initialization *************************
	      A = 0.d0 ! n>4 case
		  if(C4.NE.0.d0) A = Scale*C4/(Rbar/Pi2/alpha)**2 ! n=4 case
	      B = (Pi4 + A)/3.D0
C----------------------------------------------------------------------
C	           y(r) = 2*arctan[alpha*(r/rbar-1)]/Pi 
C                r(y) = rbar*[1+tan(Pi*y/2)/alpha]
C----------------------------------------------------------------------
           Ymin = DATAN(alpha*(Rmin/Rbar-1.d0))/Pi2 
	     H = (1.d0 - Ymin)/DFLOAT(Mes) 
C----------------------------------------------------------------------
		 DO j=0, Mes-1
              y    = Ymin + H*DFLOAT(j)
 	        r    = rbar*(1.d0 + DTAN(Pi2*y)/alpha)
		    drdy = rbar*Pi2/DCOS(Pi2*y)**2/alpha
	        Q(j) = Pi4 - Scale*U(r)*(drdy**2)
	      END DO
C  ************************** Initialization **************************
C  *********  Log derivative Johnson's method with step h   ***********
C  ********************************************************************
	        H3=H*H/3.D0
		    H6=H3/2.D0
C  *********** Estimation of initial z1-value by WKB method ***********
	        z=H*(DSQRT(-Q(0))-H*Q(0)/3.D0)
C  *********************** Outward propagation ************************
              DO j=1,Mes-1
	        k = j/2
	        w = 2.d0*Q(j)
	        if(j.NE.2*k) then
			w = 2.d0*w/(1.d0+H6*Q(j)) 
			else
			Q2(k)=Q(j)
			endif    
              z=z/(1.d0+z)-H3*w
              END DO
C  ******************* final step at Ymax point ***********************
              Yl = z/(1.d0+z)/H-H*B
C  ************************** Initialization **************************
C  *****  Log derivative Johnson's method with double step 2*H   ******
C  ********************************************************************
	        Mes=Mes/2
	        H=H*2.d0
	        H3=H3*4.D0
		    H6=H6*4.D0
C  *********** Estimation of initial z-value by WKB method ***********
	        z=H*(DSQRT(-Q(0))-H*Q(0)/3.D0)
C  *********************** Outward propagation ************************
              DO j=1,Mes-1
	        k = j/2
	        w = 2.d0*Q2(j)
	        if(j.NE.2*k) w = 2.d0*w/(1.d0+H6*Q2(j))    
              z=z/(1.d0+z)-H3*w
              END DO
C  ******************* final step at Ymax point ***********************
              Yl2 = z/(1.d0+z)/H-H*B
C======================================================================
c             Richardson's extrapolation to the zero step h->0
C======================================================================
	        SL     = Rbar*(Yl/Pi2/alpha+1.d0)
		    SL2    = Rbar*(Yl2/Pi2/alpha+1.d0)
	        SL_err = ( SL - SL2 )/15.d0
	        SL     = SL + SL_err
C======================================================================
	return
	END
C----------------------------------------------------------------------
	SUBROUTINE SCATLENGHT_MAPR1
     &           (Mes, U, Scale, Rmin, C4, Rbar, SL, SL_err)
C----------------------------------------------------------------------
C                            Input variables 
C----------------------------------------------------------------------
	INTEGER   Mes 
C     Mes is the munber of grid points; it must be an even integer number
	DOUBLE PRECISION     U, Scale, Rmin, C4
C     U is the potential defined as external subroutine-function U(r) where
C     r is the radial coordinate defined on the interval [Rmin,infinity] 
C     it is assumed that U(r->infinity)->0
C---------------------------------------------------------------------
C     Rmin is the left boundary point
C     Scale is the scaled factor of the energy: hbar^2/2mu
C     C4 is the long range coefficient for the C4/r^4 potential
C     C4 should zero for the Cn/r^n potential where n>4
C----------------------------------------------------------------------
C            Mapping parameter of the function
C	              y(r) = 1/(1+Rbar/r) 
C----------------------------------------------------------------------
	DOUBLE PRECISION     Rbar
C     Rbar>Rmin should be close to eqvilibrium distance
C----------------------------------------------------------------------
C							Output variables  
C----------------------------------------------------------------------
	DOUBLE PRECISION     SL, SL_err
C     SL is the scattering lenght
C     SL_err is the estimate of absolute error by Richardson method
C----------------------------------------------------------------------
C
	INTEGER              Mes_M, Mes_M2, j, k
	PARAMETER            (Mes_M = 10000000, Mes_M2 = Mes_M/2)
	DOUBLE PRECISION     Q(0:Mes_M), Q2(0:Mes_M2)
C
	DOUBLE PRECISION     H, Ymin, drdy, r, A, B 
	DOUBLE PRECISION     w, H3, H6, z
	DOUBLE PRECISION     Yl, Yl2, SL2
C  ************************** Initialization **************************
	     A = 0.d0                           ! n>4 case
		 if(C4.NE.0d0) A = Scale*C4/Rbar**2 ! n=4 case
C======================================================================           
C	     y(r) = 1/(1+Rbar/r);  y is defined on [Ymin,1] 
C          r(y) = Rbar*y/(1-y);  r is defined on [Rmin,infinity] 
C----------------------------------------------------------------------
           Ymin = 1.d0/(1.d0+Rbar/Rmin)
	     H  = (1.d0 - Ymin)/DFLOAT(Mes)
C
		 DO j = 0, Mes-1
              y    = Ymin + H*DFLOAT(j)
		    r    = Rbar*y /(1.d0-y)
		    drdy = Rbar / (1.d0-y)**2
	        Q(j) = -Scale * U(r) * (drdy**2)
	     ENDDO
C  ************************** Initialization **************************
C  *********  Log derivative Johnson's method with step H   ***********
C  ********************************************************************
		    H3 = H*H / 3.D0
		    H6 = H3 / 2.D0
C  *********** Estimation of initial z-value by WKB method ************
	        z = H*(DSQRT(-Q(0))-H*Q(0)/3.D0)
C  *********************** Outward propagation ************************
              DO j = 1, Mes-1
	           k = j/2
	           w = 2.d0*Q(j)
	           if(j.NE.2*k) then
		       w = 2.d0*w/(1.d0 + H6*Q(j)) 
			   else
			   Q2(k)=Q(j)
			   endif    
                 z = z/(1.d0+z)-H3*w
              ENDDO
C  ******************* final step at Ymax point ***********************
              Yl = z/(1.d0+z)/H-H*A/3.D0
C  ************************** Initialization **************************
C  *****  Log derivative Johnson's method with double step 2*H   ******
C  ********************************************************************
	        Mes = Mes/2
	        H = H*2.d0
	        H3 = H3*4.D0
		    H6 = H6*4.D0
C  *********** Estimation of initial z-value by WKB method ************
	        z = H*(DSQRT(-Q(0))-H*Q(0)/3.D0)
C  *********************** Outward propagation ************************
              DO j = 1, Mes-1
	           k = j/2
	           w = 2.d0*Q2(j)
	           if(j.NE.2*k) w = 2.d0*w/(1.d0+h6*Q2(j))    
                 z = z/(1.d0+z)-h3*w
              END DO
C  ******************* final step at Ymax point ***********************
              Yl2 = z/(1.d0+z)/H-H*A/3.D0
C======================================================================
c             Richardson's extrapolation to the zero step h->0
C======================================================================
		    SL      = Rbar * (Yl-1.d0)
			SL2     = Rbar * (Yl2-1.d0)
	        SL_err  = (SL-SL2)/15.d0
			SL      =  SL + SL_err
C======================================================================
	return
	END
