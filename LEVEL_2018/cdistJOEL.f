c**************************************************************************
c** This is Joel Tellinghuisen's code for solving the radial eigenvalue 
c  problem and calculating centrifugal distortion constants, whose routine
c  CDNEW  was used as the basis of LeRoy's  CDJOEL  routine for distortion
c  constants.  This code was received  ca. 29 August 1994.
c**************************************************************************
C    Program EIGVCD
C
C
C       This is a special version of EIGVALS designed for accurate computation
c    of eigenvalues, rotational, and centrifugal distortion constants.  The
C    calculations are for selected V of an input potential, but for J = 0 only.
C
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8  MASS1,MASS2,MU
      DIMENSION  EGY(1000),R(1000),PSIB(8000),VB(8000)
      DIMENSION  FORT(10),ROTC(25),VIBC(25),NQ(50),ETRY(50),BVIN(50)
      COMMON /CONE/ AC,BC,CC,DC
      EQUIVALENCE  (PSIB(1),EGY(1)),(PSIB(1001),R(1))
      EXTERNAL  POTSPE,POTMOR
    1 FORMAT(80H                                                     
     1                              )
    2 FORMAT(25I)
    3 FORMAT(8F)
    4 FORMAT(//5X,'Isotopic RHO factor = ',F9.6//5X,'Y00(/cm) = ',F10.3)
    5 FORMAT(/5X,'FUNDAMENTAL CONSTANT = ',1PE17.9//)
    6 FORMAT(//5X,'MASSES OF NUCLEI (BASED ON C12 = 12) ARE ',F12.8,5X,
     1    F12.8//)
    7 FORMAT(///3X,1HV, 6X,'E(Trial)',5X,'E(Found)',5X,'Conv.',3X,
     *    'Iter.',2X,'Bv(in)',5X,'Bv(out)',9X,'Dv',11X,'Hv',10X,'Lv',
     $    10X,'Mv',8X,'<R>'//)
    8 FORMAT(//5X,'INPUT POTENTIAL CURVE FOR BOUND STATE.  R IN ANGSTROM
     1S AND V IN /CM'//)
    9 FORMAT(0PI4,2X,2F13.5,1PE11.2,0PI4,2F12.8,1PE14.6,E13.5,E12.4,
     1    E11.3,0PF10.6)
   10 FORMAT(3(F9.5,2X,F12.3,10X),F9.5,2X,F12.3)
   11 FORMAT(///5X'DISSOCIATION ENERGY (/CM) = ',F12.2/5X,
     1    'CONVERGENCE CRITERION(/CM) = ',1PE9.2 //)
   12 FORMAT(///5X,'POTENTIAL IS GENERATED FROM FUNCTION OF FORM   ',
     *    10A5/10X,'PARAMETERS A,B,C,D ARE  ',4E18.9/)
   15 FORMAT(/////5X,'MAXIMUM NUMBER OF ALLOWED ITERATIONS = ',I5//)
   17 FORMAT(1P8E15.8)
   19 FORMAT(5X,'***  NO CONVERGENCE IN ',I4,'  ITERATIONS.  ON EXIT, E
     1= ',F13.5,5X,'DE = ',F9.5,7X,'NODE COUNT = ',I4)
   20 FORMAT(2I5,F15.6,E20.9)
  204 FORMAT(  5E15.7)
  208 FORMAT(//5X,I5,' VIBRATIONAL CONSTANTS'//)
  209 FORMAT(//5X,I5,' ROTATIONAL CONSTANTS'//)
  210 FORMAT(1H1,//5X,'Eigenvalues, rotational and centrifugal distortio
     1n constants')
  211 FORMAT(10A5)
          PI = 3.141592654
C         H = 6.626196E-27
          H = 6.6260755E-27
          C = 2.99792458E10
C         WT = 1.6605310E-24
          WT = 1.6605402E-24
   28 READ 1
      READ 3,  MASS1,MASS2
      READ 2,  LEVS
      PRINT 1
      PRINT 6,  MASS1,MASS2
          MU = MASS1*MASS2/(MASS1+MASS2)*WT
          CONST = 8.*PI**2*MU*C/H*1.E-16
          CONIN = 1./CONST
      PRINT 5,  CONST
C
C    INPUT INFO FOR STATE 1
C
      READ 2,  NVIBC,NROTC
      READ 211,  FORT
      READ FORT,  (VIBC(I), I=1,NVIBC)
      PRINT 208, NVIBC
      PRINT 17,  (VIBC(I), I=1,NVIBC)
      READ FORT,  (ROTC(I), I=1,NROTC)
      PRINT 209,  NROTC
      PRINT 17,  (ROTC(I), I=1,NROTC)
      READ 3,  RHO,Y00
      PRINT 4,  RHO,Y00
      READ 2,  (NQ(I), I=1,LEVS)
      READ 2,  ITYPE
        IF (ITYPE.NE.0)  GO TO 35
      READ 2,  NEPO
      READ 211,  FORT
      READ FORT ,  (EGY(I),R(I), I=1,NEPO)
      PRINT 8
      PRINT 10,  (R(I),EGY(I), I=1,NEPO)
      READ 3,   DE
      READ 2,  NLEFT
        GO TO 38
   35 READ 3,  AC,BC,CC,DC
      READ 211,  FORT
      READ 3,  DE
      PRINT 12,  FORT,AC,BC,CC,DC
   38 READ 3,  RMIN1,RMAX1,DR,EPS
      READ 2,  MAXIT,NTERP
      PRINT 11,  DE,EPS
      PRINT 15,  MAXIT
          KDIM = 8000
          JROT1 = 0
        IF (ITYPE.NE.0)  GO TO 45
      CALL POTH(RMIN1,RMAX1,DR,DE,CONST,NEPO,R,EGY,NTERP,NLEFT,VB,NBOU,
     1    KDIM,JROT1)
        GO TO 50
   45   IF (ITYPE.EQ.1)  CALL POTFUN(RMIN1,RMAX1,DR,CONST,VB,NBOU,KDIM,
     *    JROT1,POTSPE)
        IF (ITYPE.EQ.2)  CALL POTFUN(RMIN1,RMAX1,DR,CONST,VB,NBOU,KDIM,
     *    JROT1,POTMOR)
   50     IGI = 0
          IBI = 0
      DO 65  I = 1,LEVS
          VVAL = NQ(I)
          ETRY(I) = GVF(NVIBC,VIBC,DE,IGI,VVAL,RHO) + Y00
   65     BVIN(I) = BVF(NROTC,ROTC,IBI,VVAL,RHO)
   72 READ 2,  IT,NPRS,IPS
      PRINT 210
        IF (NPRS.EQ.0)  PRINT 7
      DO 90  I = 1,LEVS
          ETMP = (ETRY(I) - DE)*CONST
          KNQ = NQ(I)
          DELE = EPS*CONST
          MITER = MAXIT
          ISUC = SCHRAD(IT,NPRS,MITER,DELE,IPS,VB,PSIB,NBOU,RMIN1,
     *    RMAX1,DR,KNQ,ETMP)
          ESTO = ETMP
          ETMP = ETMP*CONIN + DE
          DELE = DELE*CONIN
        IF (ISUC.NE.0)  PRINT 19,  MAXIT,ETMP,DELE,KNQ
   82     NQ(I) = KNQ
      CALL PECY(PSIB,RMIN1,DR,NBOU,RX,R2INX)
          BVCOMP = CONIN*R2INX
        IF (NPRS.NE.0)  PRINT 7
      CALL CDNEW(NBOU,VB,PSIB,RMIN1,RMAX1,DR,ESTO,R2INX,
     1    DVV,HVV,VLV,VMV)
          DVO = -DVV*CONIN
          HVO = HVV*CONIN
          VLO = VLV*CONIN
          VMO = VMV*CONIN
      PRINT 9,  KNQ,ETRY(I),ETMP,DELE,MITER,BVIN(I),BVCOMP,DVO,HVO,
     *    VLO,VMO,RX
c     PRINT 101,  DVO,HVO,VLO,VMO
  101 FORMAT(1P4E15.5)
   90 CONTINUE
      READ 2,  IRPT
        IF (IRPT.NE.0)  GO TO 28
      STOP
      END
      SUBROUTINE  PECY(S,RMIN,DR,NP,REX,R2IN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  S(NP)
          A = 0.0
          D = 0.0
      DO 30  I = 1,NP
          R = RMIN + DR*(I-1)
          T = S(I)**2
          P = T*R
          A = A + P
          Q = T/R
          D = D + Q/R
   30 CONTINUE
          REX = A*DR
          R2IN = D*DR
      RETURN
      END
      FUNCTION SCHRAD(NI,NS,MAXIT,EPS,IPQ,V,P,NP,RMIN,RMAX,DR,NQ,E0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       FUNCTION SCHRAD OBTAINS SOLUTIONS TO RADIAL SCHROEDINGER EQUATION.
C
C    NI = 1   PRINT ITERATIONS.
C
C    NS = 1   PRINT WAVE FUNCTION AT EVERY IPQ-TH POINT.
C    NS = 2   PRINT AND PUNCH  **
C    NS = 3   PUNCH  **
C
      DIMENSION  V(1),P(1)
          PMAX = 1.E+16
          PMIN = 1.E-16
          PSTRT = PMIN
          DRSQ = DR**2
          DV = DRSQ/12.
          E = E0
          IT = 0
        IF (NI.NE.1)  GO TO 10
      PRINT 201,  NQ,E0
      PRINT 202
   10     IT = IT + 1
          P(NP) = PSTRT
          V1 = V(NP) - E
          V2 = V(NP-1) - E
        IF (V1.GT.0.0)  GO TO 15
      PRINT 203,  E,V(NP)
          SCHRAD = 2.0
      RETURN
   15     P(NP-1) = P(NP)*DEXP(DSQRT(V1)*DR)
          Y1 = (1.-DV*V1)*P(NP)
          Y2 = (1.-DV*V2)*P(NP-1)
      DO 35  I = 3,NP
          M = NP - I + 1
   20     Y3 = Y2 + Y2 - Y1 + DRSQ*V2*P(M+1)
          V1 = V2
          V2 = V(M) - E
          P(M) = Y3/(1.-DV*V2)
        IF (P(M).LT.PMAX)  GO TO 30
          M1 = M + 1
      DO 25  J = M1,NP
        IF (P(J).LT.1.) GO TO 22
          P(J) = P(J)*PSTRT
        GO TO 25
   22     P(J) = 0.0
   25 CONTINUE
          Y1 = Y1*PSTRT
          Y2 = Y2*PSTRT
          V2 = V1
        GO TO 20
   30   IF (P(M).LT.P(M+1))  GO TO 40
          Y1 = Y2
          Y2 = Y3
   35 CONTINUE
   36 PRINT 204
          SCHRAD = 2.0
      RETURN
   40     MTURN = M
        IF (MTURN.LE.5)  GO TO 36
          PNORM = 1./P(MTURN)
          YIN = Y2*PNORM
      DO 45  J = M,NP
          P(J) = P(J)*PNORM
        IF (P(J).LT.PMIN)  P(J) = 0.0
   45 CONTINUE
          P(1) = PSTRT
          V1 = V(1) - E
        IF (V1.GT.0.0)  GO TO 46
      PRINT 203,  E,V(1)
          SCHRAD = 2.0
      RETURN
   46     P(2) = P(1)*DEXP(DSQRT(V1)*DR)
          Y1 = (1.-DV*V1)*P(1)
          V2 = V(2) - E
          Y2 = (1.-DV*V2)*P(2)
      DO 60  I = 3,MTURN
          M = I - 1
   48     Y3 = Y2 + Y2 - Y1 + DRSQ*V2*P(M)
          V1 = V2
          V2 = V(I) - E
          P(I) = Y3/(1.-DV*V2)
        IF (DABS(P(I)).LT.PMAX)  GO TO 55
      DO 50  J = 1,M
        IF (DABS(P(J)).LT.1.)  GO TO 49
          P(J) = P(J)*PSTRT
        GO TO 50
   49     P(J) = 0.0
   50 CONTINUE
          Y1 = Y1*PSTRT
          Y2 = Y2*PSTRT
          V2 = V1
        GO TO 48
   55     Y1 = Y2
          Y2 = Y3
   60 CONTINUE
          PNORM = 1./P(MTURN)
          YOUT = Y1*PNORM
          YM = Y3*PNORM
      DO 65  J = 1,MTURN
          P(J) = P(J)*PNORM
        IF (DABS(P(J)).LT.PMIN)  P(J) = 0.0
   65 CONTINUE
          DF = 0.0
      DO 70 J = 1,NP
   70     DF = DF + P(J)**2
          F = (2.*YM-YOUT-YIN)/DRSQ + V(MTURN) - E
          DE = F/DF
          DCOMP = DABS(DE/E)
        IF (DCOMP.GT..05)  DE = DE*.05/DCOMP
        IF (NI.NE.1)  GO TO 75
      PRINT 205,  IT,E,F,DF,DE,MTURN
   75     E = E + DE
        IF (IT.LT.2)  GO TO 10
        IF (DABS(DE).LE.EPS)  GO TO 78
        IF (IT.LT.MAXIT)  GO TO 10
          SCHRAD = 1.0
        GO TO 80
   78     SCHRAD = 0.0
   80     NODE = 0
          NM = NP - 1
      DO 90  I = 2,NM
        IF (P(I))  82,84,86
   82   IF (P(I-1).GT.0.0)  NODE = NODE + 1
        GO TO 90
   84   IF (P(I+1))  82,90,86
   86   IF (P(I-1).LT.0.0)  NODE = NODE + 1
   90 CONTINUE
          PNORM = 1./DSQRT(DR*DF)
      DO 95  J = 1,NP
   95     P(J) = P(J)*PNORM
        IF (NS.NE.1.AND.NS.NE.2)  GO TO 105
          NPNT = 1 + (NP-1)/IPQ
          NPCOL = 1 + (NPNT-1)/6
          NPLUS = NPCOL*IPQ
      PRINT 206,  NODE,E
      DO 100  IJK = 1,NPCOL
          JS = 1 + (IJK-1)*IPQ
          JF = MIN0(JS+5*NPLUS,NP)
  100 PRINT 207,  (I,P(I), I=JS,JF,NPLUS)
  105     E0 = E
          NQ = NODE
          EPS = DE
          MAXIT = IT
          ATST = DABS(P(1)/P(MTURN))
          BTST = DABS(P(NP)/P(MTURN))
        IF (ATST.GT.1.E-08)  PRINT 210, ATST
        IF (BTST.GT.1.E-08)  PRINT 211, BTST
        IF (NS.NE.2.AND.NS.NE.3)  RETURN
      WRITE (32,208)  NQ,RMIN,RMAX,DR
      WRITE (32,209)  (I,P(I), I=1,NP,IPQ)
      RETURN
  201 FORMAT(//5X,'RADIAL SCHr.&SOLUTION FOR V = ',I3,'   WITH E(TRIAL)
     1= ',  E15.7/)
  202 FORMAT(4X,'ITER',10X,'E',14X,'F(E)',12X,'DF(E)',11X,'D(E)'//)
  203 FORMAT(/5X,'SCHRAD FAILS BECAUSE E & POTENTIAL AT INTEGRATION END
     1POINT'/10X,'E = ',  E15.6,5X,'POT. = ',  E15.6/)
  204 FORMAT(/5X,'FUNCTION FAILS TO TURN ON INWARD INTEGRATION.  ADJUST
     1INPUT.')
  205 FORMAT(6X,I4,2X,4E16.7,5X,'THE TURNING POINT OCCURS AT',I5/)
  206 FORMAT(//'        SOLUTION OF RADIAL SCHR. EQUATION FOR V = ',I3,
     1  5X,'E = ',  E15.7//6(20H    I       S(I)      )//)
  207 FORMAT(6(I5,  E13.4,2X))
  208 FORMAT(I5,5X,3F10.4)
  209 FORMAT(5(0PI4,1PE11.4))
  210 FORMAT(/5X,'WARNING  --  P(1)/PMAX.GT.E-8, = ',1PE12.4)
  211 FORMAT(/5X,'WARNING  --  P(NP)/PMAX.GT.E-8, = ',1PE12.4)
      END
      SUBROUTINE  POTH(RMIN,RMAX,DR,DE,CONST,N,R,EGY,NTRP,KL,V,NP,KD,
     1    JROT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C    POTH (NOVEMBER,1973) GENERATES WORKING POTENTIALS.
C
      DIMENSION  V(1),R(1),EGY(1),XX(2),YY(2)
    1 FORMAT(  5X,'POTENTIAL OF FORM A + C(1/R)** ',I3,' IS ATTACHED TO
     1 LEFT BRANCH, WITH A = ',1PE15.7,'  AND C = ',1PE15.7/)
    2 FORMAT(  5X,'POTENTIAL OF FORM  C(1/R)**P  IS ATTACHED TO RIGHT BRAN
     1ANCH, WITH P = ',0PF9.5,'  AND  C = ',1PE15.5//)
    3 FORMAT(//5X'TOO MANY POINTS ON FUNCTIONS.  MODIFY  RMIN,RMAX,DR.'/
     1    5X,'PRESENT VALUES ARE ',3F12.5//)
    4 FORMAT(    //5X,'ROUTINE POTH PRODUCES WORKING POTENTIAL CURVE AT'
     $,I5,'  POINTS, WITH'/ 5X,'RMIN(ANGSTROMS) = ',F10.3,10X,'RMAX = ',F
     *  F10.3,10X,'DR = ',F10.5/)
    5 FORMAT( 5X,I5,'  POINTS USED IN DIVIDED DIFFERENCES INTERPOLATION'
     */)
          NP = (RMAX-RMIN)/DR + 1.1
        IF (NP.GT.KD)  GO TO 50
          RL = R(1)
          RR = R(N)
          YY(1) = EGY(1) - DE
          YY(2) = EGY(2) - DE
          XX(1) = 1./R(1)**KL
          XX(2) = 1./R(2)**KL
          CL = (YY(1) - YY(2))/(XX(1)-XX(2))
          AL = YY(1) - CL*XX(1)
          YY(1) = DLOG(DABS(DE-EGY(N)))
          YY(2) = DLOG(DABS(DE-EGY(N-1)))
          XX(1) = DLOG(R(N))
          XX(2) = DLOG(R(N-1))
          PR = (YY(1) - YY(2))/(XX(2) - XX(1))
          CR = -DEXP(YY(1) + PR*XX(1))
        IF (DE.LT.EGY(N))  CR = - CR
          IND = 1
          XAD = 0.0
          ROTM = JROT*(JROT+1)
      DO 20  I = 1,NP
          X = RMIN + DR*(I-1)
        IF (X.LT.RL)  GO TO 12
        IF (X.GT.RR)  GO TO 14
          XT = TERPD(IND,X,N,R,EGY,NTRP) - DE
        GO TO 16
   12     XT = CL/X**KL + AL
        GO TO 16
   14     XT = CR/X**PR
   16   IF (JROT.NE.0)  XAD = ROTM/X**2
          V(I) = XT*CONST + XAD
   20 CONTINUE
      PRINT 4,  NP,RMIN,RMAX,DR
      PRINT 5,  NTRP
      PRINT 1,  KL,AL,CL
      PRINT 2,  PR,CR
      RETURN
   50 PRINT 3,  RMIN,RMAX,DR
      STOP
      END
      SUBROUTINE  POTFUN(RMIN,RMAX,DR,CONST,V,NP,KD,JROT,FPOT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  V(1)
    3 FORMAT(//5X'TOO MANY POINTS ON FUNCTIONS.  MODIFY  RMIN,RMAX,DR.'/
     1    5X,'PRESENT VALUES ARE ',3F12.5//)
    4 FORMAT(    //5X,'ROUTINE POTH PRODUCES WORKING POTENTIAL CURVE AT'
     $,I5,'  POINTS, WITH'/ 5X,'RMIN(ANGSTROMS) = ',F10.3,10X,'RMAX = ',F
     *  F10.3,10X,'DR = ',F10.5/)
          NP = (RMAX-RMIN)/DR + 1.1
        IF (NP.GT.KD)  GO TO 50
          XAD = 0.0
          ROTM = JROT*(JROT+1)
      DO 30  I = 1,NP
          X = RMIN + DR*(I-1)
        IF (JROT.NE.0)  XAD = ROTM/X**2
          V(I) = FPOT(X)*CONST + XAD
   30 CONTINUE
      PRINT 4,  NP,RMIN,RMAX,DR
      RETURN
   50 PRINT 3,  RMIN,RMAX,DR
      STOP
      END
      FUNCTION  TERPD(IND,X1,NN,FX,FY,NTER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       TERPD IS A GENERAL INTERPOLATION ROUTINE BASED ON NEWTONS DIVIDED
C    DIFFERENCES SCHEME.
C
C       IND IS AN INDICATOR USED TO SHORTEN THE SEARCH ROUTINE WHNN A NUMBER
C    OF POINTS ARE TAKEN SUCCESSIVELY BY INTERPOLATING ON A GIVEN ARRAY.
C       X1  IS THE ABSCISSUM OF THE UNKNOWN.
C       NN IS THE NUMBER OF POINTS IN THE ARRAY (FX,FY), WHICH MUST BE ORDERED
C    ACCORDING TO INCREASING VALUES OF FX.
C       NTER IS THE NUMBER OF POINTS USED IN THE INTERPOLATION (VARIABLE,
C    MAXIMUM 10).
C
C       TERPD EMPLOYS THE AUXILIARY FUNCTION TERPX.
C
      DIMENSION  FX(1),FY(1),X(10),Y(10)
          S = X1
        IF (IND.LE.0)  IND = 1
        IF (IND.GT.NN)  IND = NN
          K = IND
          NM1 = NN-1
          NHF = NTER/2
        IF (S-FX(IND))  1,50,2
    1 DO 10  I=1,IND
          K = IND + 1 - I
        IF (S.EQ.FX(K))  GO TO 50
        IF (FX(K).LT.S.AND.FX(K+1).GT.S)  GO TO 15
   10 CONTINUE
        GO TO 22
    2 DO 11  I=IND,NM1
          K = I
        IF (S.EQ.FX(K))  GO TO 50
        IF (FX(K).LT.S.AND.FX(K+1).GT.S)  GO TO 15
   11 CONTINUE
          K = NN
        IF (S.EQ.FX(K))  GO TO 50
        GO TO 24
   15   IF (K-NHF.LT.0)  GO TO 22
        IF (K-NHF+NTER.GT.NN)  GO TO 24
      DO 20  N = 1,NTER
          X(N) = FX(K-NHF+N)
   20     Y(N) = FY(K-NHF+N)
        GO TO 28
   22 DO 23  N = 1,NTER
          X(N) = FX(N)
   23     Y(N) = FY(N)
        GO TO 28
   24 DO 25  N = 1,NTER
          X(N) = FX(N+NN-NTER)
   25     Y(N) = FY(N+NN-NTER)
   28     IND = K
      TERPD = TERPX(NTER,S,X,Y)
      RETURN
   50     IND = K
      TERPD = FY(K)
      RETURN
      END
      FUNCTION TERPX(N,A,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  X(10),Y(10),DD(10,10)
      DO 20  I = 1,N
   20     DD(I,1) = Y(I)
      DO 40  J = 2,N
          K = N + 1 - J
      DO 30  I = 1,K
   30     DD(I,J) = (DD(I,J-1)-DD(I+1,J-1))/(X(I)-X(I+J-1))
   40 CONTINUE
          T = DD(1,N)
      DO 50  I = 2,N
          K = N - I + 1
   50     T = T*(A-X(K)) + DD(1,K)
      TERPX = T
      RETURN
      END
      FUNCTION  POTSPE(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONE/ B,BETA,C4,C3
      POTSPE = B*DEXP(-BETA*X) - 1.16141E+5/X - C3/X**3 - C4/X**4
      POTSPE = POTSPE - 1.4E+6/X**6
C     z = 3./x
C     potspe = 1000.*(z**12 - 2.*z**6)
      RETURN
      END
      FUNCTION  POTMOR(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONE/ DE,RE,BETA
          POTMOR = DE*(1. - DEXP(-BETA*(X-RE)))**2 - DE
      RETURN
      END
      DOUBLE PRECISION FUNCTION GVF(NC,VC,DE,INIT,V,RHO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION VC(25)
          NPOLY = VC(NC) + .1
        IF (NC.EQ.NPOLY+1) GO TO 40
        IF (V.LE.VC(NC-1)) GO TO 40
        IF (INIT.EQ.1) GO TO 20
          POW = VC(NC-2)
          VD = VC(NC-3)
          X0 = VC(NC-4)
          NLR = NC - NPOLY - 5
          INIT = 1
   20     MI = NC-5
          VEF = RHO*V + .5*(RHO-1.)
          A = VD - VEF
          T = VC(MI)
      DO 25 I=2,NLR
          MI = MI-1
   25     T = T*A + VC(MI)
          FG = 1. + T*A
      GVF = DE - X0*A**POW*FG
      RETURN
   40     MI = NPOLY
          A = (V + .5)*RHO
          T = VC(MI)
        IF (NPOLY.EQ.1) GO TO 50
      DO 45 I=2,NPOLY
          MI = MI - 1
   45     T = T*A + VC(MI)
   50 GVF = T*A
      RETURN
      END
      DOUBLE PRECISION FUNCTION BVF(NC,RC,INIT,V,RHO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RC(25)
          NPOLY = RC(NC) + .1
        IF (NC.EQ.NPOLY+1) GO TO 40
        IF (V.LE.RC(NC-1)) GO TO 40
        IF (INIT.EQ.1) GO TO 20
          POW = RC(NC-2) - 2.0
          VD = RC(NC-3)
          X1 = RC(NC-4)
          NLR = NC - NPOLY - 5
          INIT = 1
   20     MI = NC-5
          VEF = RHO*V + .5*(RHO-1.)
          A = VD - VEF
          T = RC(MI)
      DO 25 I=2,NLR
          MI = MI-1
   25     T = T*A + RC(MI)
          FB = DEXP(T*A)
      BVF = X1*A**POW*FB*RHO**2
      RETURN
   40     MI = NPOLY
          T = RC(MI)
        IF (NPOLY.EQ.1) GO TO 50
          A = (V + .5)*RHO
      DO 45 I=2,NPOLY
          MI = MI-1
   45     T = T*A + RC(MI)
   50 BVF = T*RHO**2
      RETURN
      END
      SUBROUTINE CDNEW(NP,V,P0,RMIN,RMAX,DR,E,R2IN,DVV,HVV,VLV,VMV)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       SUBROUTINE CDNEW SOLVES 2ND-ORDER INHOMOGENEOUS EQUATION FOR
C    1ST- AND 2ND-ORDER CORRECTIONS TO THE WAVEFUNCTION, THEN OBTAINS DV, HV,
C    LV, AND MV.
C
C
      DIMENSION V(1),P0(1),P(8000),Q(8000)
          DRSQ = DR*DR
          DV = DRSQ/12.
          IPA = 0
          OV11 = 0.0
          OV12 = 0.0
          OV22 = 0.0
          PER01 = 0.0
          PER02 = 0.0
          PER11 = 0.0
          PER12 = 0.0
          PER22 = 0.0
    3     P(NP) = 0.0
          P(NP-1) = 0.0
          IOT = 0
          V1 = V(NP) - E
        IF (V1.GT.0.0)  GO TO 5
      PRINT 203,  E,V(NP)
      RETURN
    5     V2 = V(NP-1) - E
          R = RMAX
        GO TO (8) IPA
          Y1 = P0(NP)*DV*(R2IN - 1./R**2)
          R = R - DR
          HP2 = P0(NP-1)*(1./R**2 - R2IN)
        GO TO 9
    8     PERT = 1./(R*R) - R2IN
          Y1 = DV*(DVV*P0(NP) - PERT*Q(NP))
          R = R - DR
          PERT = 1./(R*R) - R2IN
          HP2 = PERT*Q(NP-1) - DVV*P0(NP-1)
    9     Y2 = -DV*HP2
          P2 = 0.0
      DO 15  I = 3,NP
          M = NP-I+1
          P22 = HP2 + V2*P2
          Y3 = Y2 + Y2 - Y1 + DRSQ*P22
          R = RMIN + DR*(M-1)
          P03 = P0(M)
          PERT = 1./(R*R) - R2IN
        GO TO (12) IPA
          HP3 = P03*PERT
        GO TO 13
   12     HP3 = PERT*Q(M) - DVV*P03
   13     V3 = V(M) - E
          P3 = (Y3 + DV*HP3)/(1. - DV*V3)
        IF (V3.LT.0.0)  GO TO 20
          P(M) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          HP2 = HP3
c         IOT = IOT + 1
c       IF (IOT.NE.10)  GO TO 15
c         IOT = 0
c     WRITE (22,101)  R,V3,P03,P3
   15 CONTINUE
   20     PRS = P3
          PRT = P(M+1)
          P(1) = 0.0
          P(2) = 0.0
          P2 = 0.0
          IOT = 0.0
          V1 = V(1) - E
        IF (V1.GT.0.0)  GO TO 25
      PRINT 203,  E,V(1)
      RETURN
   25     V2 = V(2) - E
          R = RMIN
        GO TO (26) IPA
          Y1 = P0(1)*DV*(R2IN - 1./R**2)
          R = R + DR
          HP2 = P0(2)*(1./R**2 - R2IN)
        GO TO 27
   26     PERT = 1./(R*R) - R2IN
          Y1 = DV*(DVV*P0(1) - PERT*Q(1))
          R = R + DR
          PERT = 1./(R*R) - R2IN
          HP2 = PERT*Q(2) - DVV*P0(2)
   27     Y2 = -DV*HP2
          AR = 0.0
          M1 = M+1
      DO 30  I = 3,M1
          P22 = HP2 + V2*P2
          Y3 = Y2 + Y2 - Y1 + DRSQ*P22
          R = RMIN + DR*(I-1)
          P03 = P0(I)
        GO TO (28) IPA
          HP3 = P03*(1./R**2 - R2IN)
        GO TO 29
   28     HP3 = (1./(R*R)-R2IN)*Q(I) - DVV*P03
   29     V3 = V(I) - E
          P3 = (Y3 + DV*HP3)/(1. - DV*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          HP2 = HP3
c         IOT = IOT + 1
c       IF (IOT.NE.10)  GO TO 30
c         IOT = 0
c     WRITE (22,101)  R,V3,P03,P3
   30     AR = AR + P03*P3
          AMB2 = (P3-PRT)/P03
          AMB1 = (P(M)-PRS)/P0(M)
c     TYPE 103, AMB1,AMB2,PRS,PRT,P3,P03
C     PRINT 103, AMB1,AMB2,PRS,PRT,P3,P03
          AMB = (AMB1+AMB2)/2.
          M2 = M+2
      DO 35  I = M2,NP
          P03 = P0(I)
          P3 = P(I) + AMB*P03
          P(I) = P3
   35     AR = AR + P3*P03
          OV = AR*DR
C     type 101,  ov
      DO 40  I = 1,NP
          R = RMIN + DR*(I-1)
          P03 = P0(I)
          P3 = P(I) - OV*P03
          P3SQ = P3*P3
          RI2 = 1./(R*R)
        GO TO (37) IPA
          Q(I) = P3
          OV11 = OV11 + P3SQ
          PER01 = PER01 + P3*P03*RI2
          PER11 = PER11 + P3SQ*RI2
        GO TO 40
   37     P(I) = P3
          Q3 = Q(I)
          OV12 = OV12 + P3*Q3
          OV22 = OV22 + P3SQ
          PER02 = PER02 + P3*P03*RI2
          PER12 = PER12 + P3*Q3*RI2
          PER22 = PER22 + P3SQ*RI2
   40 CONTINUE
        GO TO (45) IPA
          DVV = PER01*DR
          HVV = DR*(PER11 - R2IN*OV11)
          IPA = 1
        GO TO 3
   45     HVZ = PER02*DR
          VLV = DR*(PER12 - R2IN*OV12 - DVV*OV11)
          VMV = DR*(PER22 - R2IN*OV22 - 2.*DVV*OV12 - HVV*OV11)
      type 101,  dvv,hvv,hvz,vlv,vmv
C     PRINT 101,  DVV,HVV,HVZ,VLV,VMV
c     DO 60 I=1,NP,10
c         R = RMIN + DR*(I-1)
c  60 WRITE (22,102) R,V(I),P0(I),Q(I),P(I)
c     WRITE (22,101)  OV
c     WRITE (22,101)  (P(I),I=1,NP,10)
  101 FORMAT(2X,'Raw Values of Dv, Hv(1), Hv(2), Lv, Mv'/1P5E16.8)
  102 FORMAT(0PF7.3,1P4E12.4)
  103 format(1p6E20.11)
  203 FORMAT(/5X,'CDNEW FAILS BECAUSE E.GT.V AT END POINT.'/
     1    5X,'E, V = ',1P2E15.6/)
      RETURN
      END

