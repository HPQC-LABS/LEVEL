c***********************************************************************
c** Program to Perform Interpolation Scheme Tests on the numerical    **
c  potential in a standard LEVEL input deck.  Note that the testing   **
c  subroutine used below can also consider input potentials for which **
c  the first derivatives are known, but this driver program does not  **
c  exploit that capability.   The program works by in turn dropping   **
c  each known point from the input deck and interpolating for it.     **
c  The output consists of the resulting discrepancies yielded by the  **
c  various schemes, in units  (cm-1/10**NPOW) .                       **
c***********************************************************************
      INTEGER I,ISPL,NTP,NPOW,IR2MX,NDGFMX,KDERMX
      REAL*8 EFACT,SCALE,XK(2000),YK(2000),YKD(2000)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP is the number of input turning points;  NPOW is the power of 10
c  by which the differences (in cm-1) are multiplied to facilitate
c  display in a F8.0 format;  test interpolation over the products
c  y*x**(2*IR)  for IR=0 to IR2MX ;  set NDGFMX > 0  to interpolate with
c  fits having up to NDGFMX degrees of freedom (otherwise NDGFMX=0);
c** If input data include potential first derivative at each point, 
c  set KDERMX=1 (otherwise set KDERMX=0).
c** EFACT is a scaling factor which converts read-in energies to  cm-1 .
c-----------------------------------------------------------------------
   10 READ(5,*,END=999) NTP,NPOW,IR2MX,NDGFMX,KDERMX,EFACT
c-----------------------------------------------------------------------
      WRITE(7,601) NPOW
      WRITE(6,601) NPOW
  601 FORMAT(/' Differences displayed in units of  10**(-',i2,')  cm-1')
      SCALE= EFACT*10**NPOW
c** Now read in actual known potential points (XK,YK)
c-----------------------------------------------------------------------
      IF(KDERMX.LE.0) READ(5,*) (XK(I),YK(I),I=1,NTP)
      IF(KDERMX.GT.0) READ(5,*) (XK(I),YK(I),YKD,I=1,NTP)
c-----------------------------------------------------------------------
C** Define various control parameters for INTEST.
      ISPL= 1
      IF(KDERMX.LT.0) KDERMX= 0
      IF(KDERMX.GT.1) KDERMX= 1
c** Now call the actual interpolation testing routine.
      CALL INTEST(NTP,ISPL,KDERMX,IR2MX,NDGFMX,SCALE,XK,YK,YKD)
      GO TO 10
  999 STOP
      END
c***********************************************************************
      SUBROUTINE INTEST(NPT,ISPL,KDERMX,IR2MX,NDGFMX,SCALE,XK,YK,YKD)
c** Subroutine to test various interpolation schemes for a given set of
c  NPT known points {XK(i),YK(i)} which may, for KDERMX>0, include known
c  first derivatives {YKD(i)}, by in turn dropping each known point and
c  interpolating for its value.
c** If ISPL>0 , (only allowed for KDERMX=0) first case considered is a
c  cubic spline.
c** Consider interpolation over YK*XK**(2*IRR) for IRR=1 to IR2MX.
c** NTP is the number of XK's
c** KDERMX > 0  if known derivatives are also read-in with points and 
c             used in interpolation; otherwise  KDERMX.le.0.
c** SCALE is factor multiplied into DIFF's to give good printout in F8.0
c***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,IRR,IR2,IR20,IR2X,IER,IFAIL,IR2M,IR2MX,ISPL,K,KDER,
     1  KDERMX,L,NPT,NPTM,NPTM4,NDGF,NDGFMX,NUSED,NUSE,NUSMAX
      REAL*8 SCALE,XSAV,YSAV,YDSAV
      REAL*8  XK(NPT),YK(NPT),YKD(NPT),XT(2000),YT(2000),YTD(2000),
     1    DIFF(2000),CC(24),ESUM(33),CSPL(8000),XTST(1),YTST(1)
c-----------------------------------------------------------------------
      IF(IR2MX.LT.1) IR2MX=1
      IF(NDGFMX.LT.0) NDGFMX=0
      IF(KDERMX.LT.0) KDERMX=0
c** First loop over inclusion/exclusion of derivatives, as appropriate
      DO 60 KDER=0,KDERMX
c** Loop over various possible numbers of degrees of freedom
          DO 50 NDGF=0,NDGFMX
c** Now loop over possible IR2 powers
              DO 48 IRR=1,IR2MX
                  WRITE(6,603) KDER,NDGF
                  WRITE(7,603) KDER,NDGF
                  IR20= 0
                  IR2X= 2*IRR
                  IF(ISPL.LE.0) THEN
      WRITE(6,601) SCALE,IR20,IR2X,((NUSED,NUSED=4,12,2),IR2=0,1)
      WRITE(7,601) SCALE,IR20,IR2X,((NUSED,NUSED=4,12,2),IR2=0,1)
                      ENDIF
                  IF(ISPL.GT.0) THEN
      WRITE(6,601) SCALE,IR20,IR2X,(-3,(NUSED,NUSED=4,10,2),IR2=0,1)
      WRITE(7,601) SCALE,IR20,IR2X,(-3,(NUSED,NUSED=4,10,2),IR2=0,1)
                      ENDIF
                  DO  I=1,22
                      ESUM(I)= 0.D0
                      ENDDO
c** Now omit each point (XK) in turn & interpolate for it divers ways
                  NPTM=NPT-1
                  NPTM4= 4*NPTM
                  DO 40 K=1,NPT
                      XSAV=XK(K)
                      YSAV=YK(K)
                      IF(KDER.GT.0) YDSAV=YKD(K)
c** Start by preparing reduced array for interpolation
                      DO  I=1,NPTM
                          L=I
                          IF(I.GE.K) L=I+1
                          XT(I)= XK(L)
                          YT(I)= YK(L)
                          IF(KDER.GT.0) YTD(I)= YKD(L)
                          ENDDO
                      L= 0
                      IR2= 0
c** Now loop over possible interpolation schemes
                      NUSMAX= 12
                      IF((KDER.GT.0).OR.(ISPL.GT.0)) NUSMAX= 10
   20                 IF(ISPL.LE.0) GO TO 22
c** Here determine the spline expansion coefficients, and then evaluate
c++ Use the Hutson SPLINE routine.
                      L= L+1
                      XTST(1)= XSAV
                      CALL SPLINT(1,NPTM,XT,YT,1,1,XTST,YTST)
                      DIFF(L)= (YTST(1)-YSAV)*SCALE
                      IF((IR2.GT.0).AND.(XSAV**2.GT.0.))
     1                   DIFF(L)= (YTST(1)/XSAV**IR2 - YSAV)*SCALE
                      IF(IER.NE.0) DIFF(L)= -999999.D0
      IF((K.GT.2).AND.(K.LT.(NPT-1))) ESUM(L)= ESUM(L)+DIFF(L)**2
   22                 DO 24 NUSE=4,NUSMAX,2
                          L=L+1
                          DIFF(L)= 0.D0
                          IF(NUSE.GT.NPTM) GO TO 122
                          IFAIL= 0
                          IF((NDGF.GT.0).OR.(KDER.GT.0))
     1 CALL INTFIT(NPTM,KDER,XT,YT,YTD,XSAV,NUSE,NDGF,CC,KDER,IFAIL)
                          IF((NDGF.LE.0).AND.(KDER.LE.0))
     1 CALL INTP(XT,YT,NPTM,XSAV,CC,NUSE,KDER)
                          DIFF(L)= (CC(1)-YSAV)*SCALE
                          IF((IR2.GT.0).AND.(XSAV**2.GT.0.))
     1                       DIFF(L)= (CC(1)/XSAV**IR2 - YSAV)*SCALE
                          IF(IFAIL.EQ.0) go to 124
  122                     DIFF(L)= -999999.D0
  124 IF((K.GT.2).AND.(K.LT.(NPT-1))) ESUM(L)= ESUM(L)+DIFF(L)**2
   24                     CONTINUE
                      IF((NUSMAX.EQ.12).OR.(ISPL.GT.0)) GO TO 26
                      L= L+1
                      DIFF(L)= -999999.D0
                      ESUM(L)= ESUM(L)+ DIFF(L)**2
   26                 IF(IR2.GT.IR20) GO TO 30
                      IR2=IR2X
                      IR2M=IR2-1
                      DO  I=1,NPTM
      IF(KDER.GT.0) YTD(I)= (IR2*YT(I)+XT(I)*YTD(I))*XT(I)**IR2M
                          YT(I)= YT(I)*XT(I)**IR2
                          ENDDO
                      GO TO 20
   30                 YPR= YSAV*SCALE
                      WRITE(6,602) XSAV,YPR,(DIFF(I),I=1,L)
   40                 CONTINUE
                  DEN= NPT-4
                  DO  I=1,L
                      ESUM(I)= DSQRT(ESUM(I)/DEN)
                      ENDDO
                  WRITE(6,604) (ESUM(I),I=1,10)
                  WRITE(7,604) (ESUM(I),I=1,10)
   48             CONTINUE
              IF(ISPL.GT.0) ISPL=0
   50         CONTINUE
   60     CONTINUE
   99 RETURN
  601 FORMAT(60x,'Differences Multiplied by factor   SCALE =',G10.3//
     1  20X,2(20X,'IR2 =',I2,15X)/4X,'X',6X,'Y',7X,'NUSED =',
     2  10(I3,6X)/2X,55('--'))
  602 FORMAT(1X,F7.3,F13.0,10(1X,F8.0))
  603 FORMAT(/' For   KDER=',i2,'  allow for   NDGF =',i2,
     1  '  degrees of freedom, with:')
  604 FORMAT(2X,55('--')/21X,10(1X,F8.0)/2X,55('--'))
      END
c***********************************************************************
      SUBROUTINE INTP(XI,YI,NPT,RR,C,N,IDER)
c** From the NPT known mesh points (XI,YI) , given in order of
c  increasing XI(I), select the N points (XJ,YJ) surrounding
c  the given point RR and by fitting an (N-1)-th degree polynomial
c  through them, interpolate to find the function value CC(1) and its
c  first IDER derivatives  (CC(I+1),I=1,IDER)  evaluated at RR.
c** Adapted from algorithm #416, Comm.A.C.M.
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,I1,I2,IDER,IFC,IM,J,J1,K,N,NH,NPT 
      DIMENSION  XI(NPT),YI(NPT),C(N),XJ(20),YJ(20)
      IF((N.GT.20).OR.(N.GT.NPT)) GO TO 101
      NH=N/2
c** First locate the known mesh points (XJ,YJ) bracketing RR
      I1=1
      I2=N
      IF(N.EQ.NPT) GO TO 16
      DO  I=1,NPT
          IM=I
          IF(XI(I).GT.RR) GO TO 12
          ENDDO
   12 I1=IM-NH
      IF(I1.LE.0) I1=1
      I2=I1+N-1
      IF(I2.GT.NPT) THEN
          I2=NPT
          I1=I2-N+1
          ENDIF
   16 J=0
      DO  I=I1,I2
          J=J+1
          XJ(J)=XI(I)-RR
          YJ(J)=YI(I)
          ENDDO
c** Now determine polynomial coefficients C(I) .
      DO  I=2,N
          I1=I-1
          K=I1+1
          DO  J=1,I1
              K=K-1
              YJ(K)=(YJ(K+1)-YJ(K))/(XJ(I)-XJ(K))
              ENDDO
          ENDDO
      C(1)=YJ(1)
      DO 40 I=2,N
          XX=XJ(I)
          C(I)=C(I-1)
          IF(I.EQ.2) GO TO 40
          I1=I-1
          K=I1+1
          DO  J=2,I1
              K=K-1
              C(K)=-XX*C(K)+C(K-1)
              ENDDO
   40     C(1)=YJ(I)-XX*C(1)
      IF(IDER.LE.1) GO TO 99
c** Finally, convert polynomial coefficients to derivatives at RR.
      IFC=1
      IF(IDER.GE.N) IDER=N-1
      DO  I=2,IDER
          J=I+1
          IFC=IFC*I
          C(J)=C(J)*IFC
          ENDDO
      IF(J.GE.N) GO TO 99
      J1=J+1
      DO  I=J1,N
          C(I)=0.D+0
          ENDDO
   99 RETURN
  101 WRITE(6,601) N,N,NPT
      STOP
  601 FORMAT(/' ***** Dimensioning error in  INTP :  Either   (N=',
     1  I2,' .GT. 20)   or   (N=',I2,' .GT. NPT=',I3,')')
      END
c***********************************************************************
      SUBROUTINE INTFIT(NPT,KDER,XI,YI,YD,RR,NUSE,NDGF,CC,NDER,IFAIL)
c** Given a set of known points YI (and for KDER>0 also their first
c  derivatives YD) at the NPT values of the independent variable XI (the
c  latter arranged in increasing order), select the NUSE points (XJ,YJ,
c  YDJ) bracketing the desired point RR and fit a polynomial through
c  them defined so as to have NDGF degrees of freedom [e.g., for KDER=0,
c  the polynomial order is (N-1-NDGF); for KDER=1, the order is
c  (2*N-1-NDGF)].
c** Returns values of the function CC(1) and its forst NDER derivatives
c  (CC(I+1),I=1,NDER) evaluated at RR.
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER  I,I1,I2,IM,ID,IFC,IORD,J,J1,NPT,NH,NUSE,NPAR,NPAR1,NDGF,
     1  MXDATA,MXPARM,KDER,NDER,IFAIL
      DIMENSION  XI(NPT),YI(NPT),YD(NPT),XJ(20),YJ(20),YDJ(20),
     1  XX(40,25),YO(40),UY(40),DY(40),CC(25),UCC(25),SCC(25),CM(25,25)
      DATA Z0/0.D0/,Z1/1.D0/,MXDATA/40/,MXPARM/25/
      IF((NUSE.GT.20).OR.(NUSE.GT.NPT)) GO TO 101
      NH=NUSE/2
c** First locate the (known) mesh points {XJ,YJ,YDJ} bracketing RR
      I1=1
      I2=NUSE
      IF(NUSE.EQ.NPT) GO TO 16
      DO  I=1,NPT
          IM=I
          IF(XI(I).GT.RR) GO TO 12
          ENDDO
   12 I1=IM-NH
      IF(I1.LE.0) I1=1
      I2=I1+NUSE-1
      IF(I2.LE.NPT) GO TO 16
      I2=NPT
      I1=I2-NUSE+1
   16 J=0
      DO  I=I1,I2
          J=J+1
          XJ(J)= XI(I)-RR
          YJ(J)= YI(I)
          YDJ(J)= YD(I)
          ENDDO
c** Now prepare XX array for call to LLSQF
      IORD= NUSE-1-NDGF
      IF(KDER.GT.0) IORD= IORD+NUSE
      NPAR= IORD+1
      NPAR1= NPAR+1
      ID= NUSE
      DO  I=1,NUSE
          XX(I,1)= Z1
          UY(I)= Z1
          DO  J=1,IORD
              J1=J+1
              XX(I,J1)= XX(I,J)*XJ(I)
              ENDDO
          YO(I)= YJ(I)
          IF(KDER.GT.0) THEN
              ID= NUSE+I
              XX(ID,1)= Z0
              UY(ID)= Z1
              DO  J=1,IORD
                  J1= J+1
                  XX(ID,J1)= J*XX(I,J)
                  ENDDO
              YO(ID)= YDJ(I)
              ENDIF
          ENDDO
      CALL LLSQF(ID,NPAR,MXDATA,MXPARM,YO,UY,XX,DY,CC,UCC,SCC,CM,SERR)
   38 IF(NDER.LE.0) GO TO 99
c** If required, convert polynomial coefficients CC(J) to derivatives.
      IFC=1
      IF(NDER.GE.NUSE) NDER=NUSE-1
      DO  I=2,NDER
          J=I+1
          IFC=IFC*I
          CC(J)=CC(J)*IFC
          ENDDO
      IF(J.GE.NUSE) GO TO 99
      J1=J+1
      DO 60 I=J1,NUSE
   60 CC(I)=0.D+0
   99 RETURN
  101 WRITE(6,601) NUSE,NUSE,NPT
      STOP
  601 FORMAT(/' *** Dimensioning error in INTFIT :  either   (NUSE=',
     1  I2,' .GT. 20)   or   (NUSE=',I2,' .GT. NPT=',I3,')')
      END
c***********************************************************************
      SUBROUTINE SPLINT(LNPT,NTP,R1,V1,MBEG,MEND,XX,YY)
c** Subroutine to generate (if LNPT.ge.0) 4*NTP coefficients CSP(J)
c  of a cubic spline passing through the NTP points (R1(J),V1(J))
c  and to then calculate values of the resulting function YY(I) at the
c  entering abscissae values XX(I) for  I=MBEG to MEND.
c** If LNPT < 0 , generate function values at the given XX(I) using
c  the coefficients CSP(J) obtained and SAVEd on a preceding call.
c** Assumes both R1(J) & XX(I) are monotonic increasing.
c+++++ Calls only subroutine SPLINE +++++++++++++++++++++++++++++++++++
c======================================================================
      INTEGER MAXSP
      PARAMETER (MAXSP=2400)
      INTEGER  I,IER,I1ST,IDER,JK,K,KK,LNPT,N2,N3,NIPT,NTP,MBEG,MEND
      REAL*8 EPS,R2,RI,RRR,TTMP,R1(NTP),V1(NTP),CSP(MAXSP),
     1  YY(MEND),XX(MEND)
      SAVE CSP
c
      IF(4*NTP.GT.MAXSP) THEN
          WRITE(6,602) MAXSP,NTP
          STOP
          ENDIF
      EPS= 1.D-6*(R1(2)-R1(1))
      N2= 2*NTP
      N3= 3*NTP
      IF(LNPT.GT.0) THEN
c** On first pass for a given data set, generate spline function
c  coefficients in subroutine SPLINE
c** Start by using a cubic polynomial at each end of the range to get
c  the first derivative at each end for use in defining the spline.
          IDER= 1
          NIPT= 4
          I1ST= NTP-3
          CALL INTP(R1(I1ST),V1(I1ST),NIPT,R1(NTP),CSP,NIPT,IDER)
          TTMP= CSP(2)
          CALL INTP(R1,V1,NIPT,R1(1),CSP,NIPT,IDER)
          CSP(1)= CSP(2)
          CSP(2)= TTMP
c** Now call routine to actually generate spline coefficients
          CALL SPLINE(R1,V1,NTP,3,CSP,MAXSP,IER)
          IF(IER .NE. 0) THEN
              WRITE(6,604)
              STOP
              ENDIF
          ENDIF
      IF(MEND.LT.MBEG) GO TO 99
c** Now, use spline to generate function at desired points XX(I)
      DO  I= MBEG,MEND
          RI= XX(I)
          RRR= RI-EPS
          KK= 1
c** For a monotonic increasing distance array XX(I),  this statement 
c  speeds up the search for which set of cubic coefficients to use.
          IF(I.GT.MBEG) THEN
              IF(XX(I).GT.XX(I-1)) KK= JK
              ENDIF
          DO  K= KK,NTP
              JK= K
              IF(R1(K).GE.RRR) GO TO 64
              ENDDO
   64     CONTINUE
          JK= JK-1
          IF(JK.LT.1) JK= 1
          R2= RI-R1(JK)
          YY(I)= CSP(JK)+R2*(CSP(NTP+JK)+R2*(CSP(N2+JK)+R2*CSP(N3+JK)))
          ENDDO
   99 RETURN
  602 FORMAT(' *** ERROR in SPLINT ***  Array dimension  MAXSP=',I4,
     1  ' cannot contain spline coefficients for  NTP=',I4)
  604 FORMAT(' *** ERROR in generating spline coefficients in SPLINE')
      END
c**********************************************************************
      SUBROUTINE SPLINE(X,Y,N,IOPT,C,N4,IER)
c** Subroutine for generating cubic spline coefficients
c  C(J), (J=1,N4=4*N) through the N points X(I), Y(I).
c** C(I+M*N), M=0-3  are the coefficients of order  0-3  of cubic
c  polynomial expanded about X(I) so as to describe the interval:
c             -  X(I) to X(I+1)  , if  X(I)  in increasing order
c             -  X(I-1) to X(I)  , if  X(I)  in decreasing order.
c** IOPT indicates boundary conditions used in creating the  spline .
c*  If (IOPT=0)  second derivatives = zero at both ends of range.
c*  If (IOPT=1)  1st derivative at first point X(1) fixed at C(1),
c                and 2nd derivative at X(N) = zero.
c*  If (IOPT=2)  1st derivative at last point X(N) fixed at C(2),
c                and 2nd derivative at X(1) = zero.
c*  If (IOPT=3)  constrain first derivatives at end points to have
c                (read in) values  C(1)  at  X(1)  &  C(2)  at  X(N)
c** IER is the error flag.  IER=0  on return if routine successful.
c-----------------------------------------------------------------------
      INTEGER I,II,IER,IOH,IOL,IOPT,J,J1,J2,J3,NER,N,N4,JMP
      REAL*8  A,H,R,DY2,DYA,DYB,XB,XC,YA,YB, X(N),Y(N),C(N4)
c
      JMP= 1
      NER= 1000
      IF(N.LE.1) GO TO 250
c** Initialization
      XC= X(1)
      YB= Y(1)
      H= 0.D0
      A= 0.D0
      R= 0.D0
      DYB= 0.D0
      NER= 2000
c
c  IOL=0 - given derivative at firstpoint
c  IOH=0 - given derivative at last point
c
      IOL= IOPT-1
      IOH= IOPT-2
      IF(IOH.EQ.1) THEN
          IOL= 0
          IOH= 0
          ENDIF
      DY2= C(2)
c
c  Form the system of linear equations
c  and eliminate subsequentially
c
      J= 1
      DO  I= 1,N
          J2= N+I
          J3= J2+N
          A= H*(2.D0-A)
          DYA= DYB+H*R
          IF(I.GE.N) THEN
c
c  set derivative dy2 at last point
c
              DYB= DY2
              H= 0.D0
              IF(IOH.EQ.0) GOTO 200
              DYB= DYA
              GOTO 220
              ENDIF
          J= J+JMP
          XB= XC
          XC= X(J)
          H= XC-XB
c
c  II= 0 - increasing abscissae
c  II= 1 - decreasing abscissae
c
          II= 0
          IF(H.LT.0) II= 1
          IF(H.EQ.0) GO TO 250
          YA= YB
          YB= Y(J)
          DYB= (YB-YA)/H
          IF(I.LE.1) THEN
              J1= II
              IF(IOL.NE.0) GO TO 220
              DYA= C(1)
              ENDIF
200       IF(J1.NE.II) GO TO 250
          A= 1.D0/(H+H+A)
220       R= A*(DYB-DYA)
          C(J3)= R
          A= H*A
          C(J2)= A
          C(I)= DYB
          ENDDO
c
c  back substitution of the system of linear equations
c     and computation of the other coefficients
c
      A= 1.D0
      J1= J3+N+II-II*N
      I= N
      DO  IOL= 1,N
          XB= X(J)
          H= XC-XB
          XC= XB
          A= A+H
          YB= R
          R= C(J3)-R*C(J2)
          YA= R+R
          C(J3)= YA+R
          C(J2)= C(I)-H*(YA+YB)
          C(J1)= (YB-R)/A
          C(I)= Y(J)
          A= 0.D0
          J= J-JMP
          I= I-1
          J2= J2-1
          J3= J3-1
          J1= J3+N+II
          ENDDO
      IER= 0
      RETURN
250   IER= NER
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE LLSQF(NDATA,NPARM,MXDATA,MXPARM,YO,YU,DYDP,YD,PV,PU,PS,
     1                 CM,DSE)
c**  Program for performing linear least squares fits using orthogonal 
c  decomposition of the Design (partial derivative) matrix.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1997  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** This version of the program is designed for the data sets of modest
c  size where it is convenient to generate and store the complete 
c  partial derivative matrix prior to calling LLSQF.  If this is not the
c  case, subroutine version LLSQFVL, which generates this partial 
c  derivative array one row at a time through calls to a user-supplied
c  subroutine, should be used.
c
c** On entry: NDATA  is the number of data to be fitted (.le.MXDATA)
c             NPARM  the number of parameters to be varied (.le.MXPARM)
c                 If NPARM.le.0 , assume  YD(i)=YO(i)  and calculate the
c                 (RMS dimensionless deviation)=DSE  from them & YU(i)
c             MXDATA & MXPARM are array dimension parameters (see below)
c                 Internal array sizes currently assume  MXPARM .le. 60
c             YO(i)  are the NDATA 'observed' data;  for iterative 
c                  non-linear fits these are:  [Y(obs,i) - Y(trial,i)]
c             YU(i)  are the uncertainties in these YO(i) values
c             DYDP(i,j)  is the partial derivative array  dYO(i)/dPV(j)
c
c** On Exit: PV(j)  are the fitted parameter values;  for iterative
c                  non-linear fits these are the parameter changes
c            PU(j) are 95% confidence limit uncertainties in the PV(j)'s
c            PS(j) are 'parameter sensitivities' for the PV(j)'s, defined
c               such that the RMS displacement of predicted data  due to
c               rounding off parameter-j by PS(j) is .le. DSE/10*NPARM
c            DSE  is predicted (dimensionless) standard error of the fit
c            YD(i) is the array of differences  [YO(i) - Ycalc(i)]
c            CM(j,k)  is the correlation matrix obtained by normalizing
c   variance/covariance matrix:  CM(j,k) = CM(j,k)/SQRT[CM(j,j)*CM(k,k)]
c** The squared 95% confidence limit uncertainty in a property F({PV(j)})
c  defined in terms of the fitted parameters {PV(j)} is (where the
c  L.H.S. involves  [row]*[matrix]*[column]  multiplication):
c  [D(F)]^2 = [PU(1)*dF/dPV(1), PU(2)*dF/dPV(2), ...]*[CM(j,k)]*
c                              [PU(2)*dF/dPV(1), PU(2)*dF/dPV(2), ...]
c** Externally dimension:  YO, YU and YD  .ge. NDATA (say as MXDATA),
c             PV, PU  and  PS  .ge.  NPARM (say as MXPARM), 
c             DYDP  with column length MXDATA and row length .ge. NPARM
c             CM   as a square matrix with column & row length  MXPARM
c  Authors: Robert J. Le Roy  &  Michael Dulick, Department of Chemistry
c    U. of Waterloo, Waterloo, Ontario  N2L 3G1.    Version of: 29/10/97
c***********************************************************************
      INTEGER I,J,K,L,M,IDF,NDATA,MXDATA,NPARM,MXPARM
      REAL*8  YO(NDATA), YU(NDATA), YD(NDATA), PV(NPARM), PU(NPARM), 
     1   PS(NPARM), DYDP(MXDATA,NPARM), CM(MXPARM,MXPARM), DSE,
     2   PX(60), F95(10), TFACT, S, U
      DATA F95/12.7062D0,4.3027D0,3.1824D0,2.7764D0,2.5706D0,2.4469D0,
     1  2.3646D0,2.3060D0,2.2622D0,2.2281D0/
c
      IF((NDATA.GT.MXDATA).OR.(NPARM.GT.MXPARM).OR.(NPARM.GT.60)
     1                    .OR.(NPARM.GT.NDATA)) THEN
c** If array dimensioning inadequate, print warning & then STOP
          WRITE(6,601) NDATA,MXDATA,NPARM,MXPARM
          STOP
          ENDIF
      IF(NPARM.LE.0) THEN
c** If no parameters varied - simply calculate RMS deviation = DSE
          DSE= 0.D0
          DO 2 I= 1,NDATA
              YD(I)= YO(I)
    2         DSE= DSE+ (YD(I)/YU(I))**2
          DSE= DSQRT(DSE/DFLOAT(NDATA))
          RETURN
          ENDIF
c** TFACT  is 95% student t-value for (NDATA-NPARM) degrees of freedom.
c [Approximate expression for (NDATA-NPARM).GT.10 accurate to ca. 0.002]
      TFACT= 0.D0
      IF(NDATA.GT.NPARM) THEN
          IDF= NDATA-NPARM
          IF(IDF.GT.10) TFACT= 1.960D0*DEXP(1.265D0/DFLOAT(IDF))
          IF(IDF.LE.10) TFACT= F95(IDF) 
        ELSE
          TFACT= 0.D0
        ENDIF
      DO 10 I = 1,NPARM
          PS(I) = 0.D0
          PU(I) = 0.D0
          PX(I) = 0.D0
          DO 8 J = 1,NPARM
    8         CM(I,J) = 0.D0
   10     CONTINUE
c
c** Begin by forming the Jacobian Matrix from the input partial 
c  derivative matrix DYDP.  For VERY large data sets, these partial 
c  derivatives may be generated inside this loop (see version LLSQFVL).
      DO 14 I = 1,NDATA
          S = 1.D0 / YU(I)
          U = YO(I) * S
          DO 12 J = 1,NPARM
   12         PV(J) = DYDP(I,J) * S
          CALL QROD(NPARM,MXPARM,MXPARM,CM,PV,PX,U,PS,PU)
   14     CONTINUE
c
c** Compute the inverse of  CM 
      CM(1,1) = 1.D0 / CM(1,1)
      DO 20 I = 2,NPARM
          L = I - 1
          DO 18 J = 1,L
              S = 0.D0
              DO 16 K = J,L
   16             S = S + CM(K,I) * CM(J,K)
   18         CM(J,I) = -S / CM(I,I)
   20     CM(I,I) = 1.D0 / CM(I,I)
c
c** Solve for parameter values  PV(j)
      DO 26 I = 1,NPARM
          J = NPARM - I + 1
          PV(J) = 0.D0
          DO 24 K = J,NPARM
   24         PV(J) = PV(J) + CM(J,K) * PX(K)
   26     CONTINUE
c
c** Get (upper triangular) "dispersion Matrix" [variance-covarience
c  matrix  without the sigma^2 factor].
      DO 30 I = 1,NPARM
          DO 30 J = I,NPARM
              U = 0.D0
              DO 28 K = J,NPARM
   28             U = U + CM(I,K) * CM(J,K)
   30         CM(I,J) = U
c** Generate core of Parameter Uncertainties  PU(j) and (symmetric)
c   correlation matrix  CM
      DO 36 J = 1,NPARM
          PU(J) = DSQRT(CM(J,J))
          DO 32 K= J,NPARM
   32         CM(J,K)= CM(J,K)/PU(J)
          DO 34 K= 1,J
              CM(K,J)= CM(K,J)/PU(J)
   34         CM(J,K)= CM(K,J)
   36     PX(J)= 0.d0
c
c** Generate differences:   YD(i) = [YO(i) - Ycalc(i)] , standard error
c  DSE = sigma^2,  and prepare to calculate Parameter Sensitivities PS
      DSE= 0.D0
      DO 40 I = 1,NDATA
          S = 1.D0 / YU(I)
          U = 0.D0
          DO 38 J = 1,NPARM
              PX(J)= PX(J)+ (DYDP(I,J)*S)**2
   38         U = U + DYDP(I,J) * PV(J)
          YD(I) = YO(I) - U
   40     DSE= DSE+ (S*YD(I))**2
      IF(NDATA.GT.NPARM) THEN
          DSE= DSQRT(DSE/(NDATA-NPARM))
        ELSE
          DSE= 0.d0
        ENDIF
c** Use DSE to get final (95% confid. limit) parameter uncertainties PU
c** Calculate 'parameter sensitivities', changes in PV(j) which would
c  change predictions of input data by an RMS average of  DSE*0.1/NPARM
      U= DSE*0.1d0/DFLOAT(NPARM)
      S= DSE*TFACT
      DO 44 J = 1,NPARM
          PU(J)= S* PU(J)
   44     PS(J)= U*DSQRT(NDATA/PX(J))
c
      RETURN
  601 FORMAT(/' *** Dimensioning problems in LLSQF *** (NDATA, MXDATA, N
     1PARM, MXPARM)  =  (',I5,4(' ,',I5),' )')
      END
c***********************************************************************
      SUBROUTINE QROD(N,NR,NC,A,R,F,B,GC,GS)
C** Performs ORTHOGONAL DECOMPOSITION OF THE LINEAR LEAST-SQUARES    
C            EQUATION J * X = F TO A * X = B(TRANSPOSE) * F WHERE   
C            J IS THE JACOBIAN IN WHICH THE FIRST N ROWS AND COLUMNS
C            ARE TRANSFORMED TO THE UPPER TRIANGULAR MATRIX A      
C            (J = B * A), X IS THE INDEPENDENT VARIABLE VECTOR, AND
C            F IS THE DEPENDENT VARIABLE VECTOR. THE TRANSFORMATION
C            IS APPLIED TO ONE ROW OF THE JACOBIAN MATRIX AT A TIME.
C  PARAMETERS :                                                   
C      N   -  (INTEGER) DIMENSION OF A TO BE TRANSFORMED.        
C      NR  -  (INTEGER) ROW DIMENSION OF A DECLARED IN CALLING PROGRAM.
C      NC  -  (INTEGER) Column DIMENSION OF F DECLARED IN CALLING PROGRAM.
C      A   -  (REAL*8 ARRAY OF DIMENSIONS .GE. N*N) UPPER TRIANGULAR
C             TRANSFORMATION MATRIX.                               
C      R   -  (REAL*8 LINEAR ARRAY OF DIMENSION .GE. N) ROW OF    
C             JACOBIAN TO BE ADDED.                             
C      F   -  (REAL*8 LINEAR ARRAY .GE. TO THE ROW DIMENSION OF THE
C             JACOBIAN) TRANSFORMED DEPENDENT VARIABLE MATRIX.    
C      B   -  (REAL*8) VALUE OF F THAT CORRESPONDS TO THE ADDED  
C             JACOBIAN ROW.                                     
C     GC   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS COSINE TRANSFORMATIONS.
C     GS   -  (REAL*8 LINEAR ARRAY .GE. N) GIVENS SINE TRANSFORMATIONS. 
C--------------------------------------------------------------------
C  AUTHOR : MICHAEL DULICK, Department of Chemistry,
C           UNIVERSITY OF WATERLOO, WATERLOO, ONTARIO N2L 3G1
C--------------------------------------------------------------------
      INTEGER  I,J,K,N,NC,NR
      REAL*8 A(NR,NC), R(N), F(NR), GC(N), GS(N), B, Z(2)
      DO 10 I = 1,N
          Z(1) = R(I)
          J = I - 1
          DO   K = 1,J
            Z(2) = GC(K) * A(K,I) + GS(K) * Z(1)
            Z(1) = GC(K) * Z(1) - GS(K) * A(K,I)
            A(K,I) = Z(2)
            ENDDO
          GC(I) = 1.D0
          GS(I) = 0.D0
          IF(Z(1) .EQ. 0.D0) GOTO 10
          IF(DABS(A(I,I)) .GE. DABS(Z(1))) GOTO 3
          Z(2) = A(I,I) / Z(1)
          GS(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
          GC(I) = Z(2) * GS(I)
          GOTO 4
    3     Z(2) = Z(1) / A(I,I)
          GC(I) = 1.D0 / DSQRT(1.D0 + Z(2) * Z(2))
          GS(I) = Z(2) * GC(I)
    4     A(I,I) = GC(I) * A(I,I) + GS(I) * Z(1)
          Z(2) = GC(I) * F(I) + GS(I) * B
          B = GC(I) * B - GS(I) * F(I)
          F(I) = Z(2)
   10     CONTINUE
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
