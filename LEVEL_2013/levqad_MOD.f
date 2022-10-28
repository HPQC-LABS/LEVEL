c!!
c!! Partially fixed version of levqad for unequally spaced points.
c!! some slippage merging this with  WIDTH
c!!
c!!           22 February 2009
c**********************************************************************
      SUBROUTINE LEVQAD(Y1,Y2,Y3,X1,X2,X3,H,RT,ANS1,ANS2)
c** Subroutine "LEVQAD" fits quadratic  Y = A + B*X + C*X**2  through
c  function values  Y1, Y2, Y3  at points  X1, X2 & X3, where  Y1 < 0
c  and (Y2,Y3 .ge.0), locates the function zero (at RT, relative to 
c  between points X1 & X2, and evaluates the integral from RT to X3
c  of   1/sqrt(Y)  , returned as  ANS1, and the integral (same range)
c  of  sqrt(Y) , which is ANS2
c** Alternately, if Y1 & Y3 both  < 0  and only the middle point
c  Y2.ge.0 ,   fit the points to:  Y = A - B*(X-X0)**2 , locate the
c  turning points between which  Y(X) > 0  and evaluate these integrals
c  on this interval.  *************************************************
c----------------------------------------------------------------------
      REAL*8  A,ANS1,ANS2,B,C,CQ,H,HPI,R1,R2,RCQ,RR,RT,SL3,SLT,
     1        X0,Y1,Y2,Y3,X1,X2,X3,ZT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DATA HPI/1.570796326794896D0/
      IF((Y1.GE.0.d0).OR.(Y2.LT.0.d0)) GO TO 99
c     IF(Y3.LT.0.d0) GO TO 50
c** Here treat case where both 'Y2' & 'Y3' are positive
      IF(DABS((Y2-Y1)/(Y3-Y2) -1.D0).LT.1.d-10) THEN
c ... special case of true (to 1/10^10) linearity ...
          RT= -H*Y2/(Y2-Y1)
          ANS1= 2.d0*(H-RT)/DSQRT(Y3)
          ANS2= ANS1*Y3/3.D0
          RETURN
          ENDIF
      B= (Y2-Y3)/(X2-X3)
      C= ((Y1-Y2)/(X1-X2) - B)/(X1-X3)
      B= B - C*(X2+X3)
      A= Y2 - X2*(B + C*X2)
      CQ= B**2- 4.d0*A*C
      RCQ= DSQRT(CQ)
      R1= (-B-RCQ)/(2.d0*C)
      R2= R1+ RCQ/C
      IF(((R2 - X1)*(R2-X2)).LE.0.d0)  RT= R2

      SLT= 2.d0*C*RT + B
      IF(C.GE.0.d0) THEN
          ANS1= DLOG((2.d0*DSQRT(C*Y3)+SL3)/SLT)/DSQRT(C)
        ELSE
          ANS1= -(DASIN(SL3/RCQ)- DSIGN(HPI,SLT))/DSQRT(-C)
        ENDIF
      ANS2= (SL3*DSQRT(Y3)- CQ*ANS1/2.d0)/(4.d0*C)
c     IF(RT.GE.H) WRITE(6,601) H,R1,R2
c 601 FORMAT(' *** CAUTION *** in LEVQAD, turning point not between poin
c    1ts 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)
      RETURN
c!!
c** Here treat case when only 'Y2' is non-negative
c!!  Not yet updated this to the general case of NOT equally spaced points

c  50 RR= (Y2-Y1)/(Y2-Y3)
c     X0= H*(RR-1.d0)/((RR+1.d0)*2.d0)
c     B= (Y2-Y1)/(H*(2.d0*X0+H))
c     A= Y2+ B*X0**2
c     ZT= DSQRT(A/B)
c     RT= X0- ZT
c     ANS1= 2.d0*HPI/DSQRT(B)
c     ANS2= ANS1*A*0.5d0
c     RETURN
c!!
   99 WRITE(6,602) Y1,Y2
  602 FORMAT(' *** ERROR in LEVQAD *** No turning point between 1-st two
     1 points as   Y1=',D10.3,'   Y2=',D10.3)
      ANS1= 0.d0
      ANS2= 0.d0
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

