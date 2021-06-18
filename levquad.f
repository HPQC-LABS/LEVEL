c***********************************************************************
      SUBROUTINE LEVQAD(Y1,Y2,Y3,H,RT,ANS1,ANS2)
c** Subroutine "LEVQAD" fits quadratic  Y = A + B*X + C*X**2  through
c  function values  Y1, Y2, Y3  at equally spaced points separated by
c  distance H, where  Y1 < 0  and (Y2,Y3 .ge.0), locates the function
c  zero (at RT, relative to  X1 < X2 = 0) between points X1 & X2, and
c  evaluates the integral from RT to R3 of   1/sqrt(Y)  , called
c  ANS1, and the integral (same range) of  sqrt(Y) , which is ANS2
c** Alternately, if Y1 & Y3 both  < 0  and only the middle point
c  Y2.ge.0 ,   fit the points to:  Y = A - B*(X-X0)**2 , locate the
c  turning points between which  Y(X) > 0  and evaluate these integrals
c  on this interval.  **************************************************
c-----------------------------------------------------------------------
      REAL*8  A,ANS1,ANS2,B,C,CQ,H,HPI,R1,R2,RCQ,RR,RT,SL3,SLT,
     1        X0,Y1,Y2,Y3,ZT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DATA HPI/1.570796326794896D0/
      IF((Y1.GE.0).OR.(Y2.LT.0)) GO TO 99
      IF(Y3.LT.0.d0) GO TO 50
c** Here treat case where both 'Y2' & 'Y3' are positive
      IF(DABS((Y2-Y1)/(Y3-Y2) -1.D0).LT.1.d-10) THEN
c ... special case of true (to 1/10^10) linearity ...
          RT= -H*Y2/(Y2-Y1)
          ANS1= 2.d0*(H-RT)/DSQRT(Y3)
          ANS2= ANS1*Y3/3.D0
          RETURN
          ENDIF
      C= (Y3-2.d0*Y2+Y1)/(2.d0*H*H)
      B= (Y3-Y2)/H-C*H
      A= Y2
      CQ= B**2- 4.d0*A*C
      RCQ= DSQRT(CQ)
      R1= (-B-RCQ)/(2.d0*C)
      R2= R1+ RCQ/C
      IF((R2.LE.0.d0).AND.(R2.GE.-H)) RT=R2
      IF((R1.LE.0.d0).AND.(R1.GE.-H)) RT=R1
      SL3= 2.d0*C*H+B
      SLT= 2.d0*C*RT+B
      IF(C.LT.0.d0) GO TO 10
      ANS1= DLOG((2.d0*DSQRT(C*Y3)+SL3)/SLT)/DSQRT(C)
      GO TO 20
   10 ANS1= -(DASIN(SL3/RCQ)- DSIGN(HPI,SLT))/DSQRT(-C)
   20 ANS2= (SL3*DSQRT(Y3)- CQ*ANS1/2.d0)/(4.d0*C)
      IF(RT.GE.H) WRITE(6,601) H,R1,R2
  601 FORMAT(' *** CAUTION *** in LEVQAD, turning point not between poin
     1ts 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)
      RETURN
c** Here treat case when only 'Y2' is non-negative
   50 RR= (Y2-Y1)/(Y2-Y3)
      X0= H*(RR-1.d0)/((RR+1.d0)*2.d0)
      B= (Y2-Y1)/(H*(2.d0*X0+H))
      A= Y2+ B*X0**2
      ZT= DSQRT(A/B)
      RT= X0- ZT
      ANS1= 2.d0*HPI/DSQRT(B)
      ANS2= ANS1*A*0.5d0
      RETURN
   99 WRITE(6,602) Y1,Y2
  602 FORMAT(' *** ERROR in LEVQAD *** No turning point between 1-st two
     1 points as   Y1=',D10.3,'   Y2=',D10.3)
      ANS1= 0.d0
      ANS2= 0.d0
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

