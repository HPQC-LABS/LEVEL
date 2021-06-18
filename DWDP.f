      SUBROUTINE DWDP(DATN,ISTATE,NPPS,ZMU,vb,Jr,EO,DEDPK,width,dWdPk)
c***********************************************************************
c** This subroutine calculates dW/dP using Pajunen's quadrature method 
c   Ref: JCP. 71(6), 2618, 1979.
c
c** The values of dW/dP(k) for each parameter P(k) are stored in dWdPk.
c
c** On entry:
c    DATN    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopomer being considered.
c    NPPS(i) are the number of parameters varied in each state.
c    ZMU     is the reduced mass of the diatom in atomic units.
c    vb      is the vibrational quantum number.
c    Jr      is the rotational quantum number.
c    EO      is the energy for this level (DE in cm-1).
c    DEDPK(i) are the values of the partial derivative dE/dP(k).
c    
c
c** On exit:
c    width   is the calculated width (DE in cm-1).
c    dWdPk(i) are the values of the partial derivative dW/dP(k).
c
c***********************************************************************
      INCLUDE 'arrsizes.h'

*** Variables declaration
*----------------------------------------------------------------------*
      INTEGER ISTATE,DATN,NPPS(NSTATEMX),vb,Jr,IPV,I,l,j,NOWIDTHS

      COMMON /WIDTHBLK/NOWIDTHS

c ** Width parameters
      REAL*8 width,dWdPk(NPARMX),dvdGv,dvdGv2

c ** Turning points
      REAL*8 R1, R2, R3

c ** Phase integrals
      REAL*8 I0m(NPARMX),I00(NPARMX),I10(NPARMX),I11(NPARMX),I0
 
      REAL*8 ZMU,EO,FCTOR,GAMA,VMAX,DEDPV,DEDPK(NPARMX),Pi,deltaP
      DATA Pi/3.141592653589793D0/
c
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)

      INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)

      CHARACTER*2 NAME(2)

      CHARACTER*2 SLABL(-2:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL

c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      REAL*8 R(NPNTMX,NSTATEMX)
      COMMON /RBLK/R
*----------------------------------------------------------------------*

      FCTOR = SQRT(ZMU/16.857629206d0)


*** Begins by calling subroutine locateTP to find the three turning 
*   points of the quasi-bound level 
*      WRITE (33,150) vb,Jr,EO
* 150  FORMAT ('Level v =',I3,' J =',I3,' Ep =',G16.8)
      CALL locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
c *** if some turning points are not found, ignore this datum
      IF (R1.LT.0 .OR. R2.LT.0 .OR. R3.LT.0) GO TO 195

*** Calls subroutine phaseIntegral to calculate the phase integrals 
*   using Pajunen's quadrature method
      CALL PhaseIntegral(DATN,ISTATE,NPPS,R2,R3,0,-1,EO,DEDPK,I0m)
      CALL PhaseIntegral(DATN,ISTATE,NPPS,R1,R2,0,0,EO,DEDPK,I00)
      CALL PhaseIntegral(DATN,ISTATE,NPPS,R2,R3,1,0,EO,DEDPK,I10)
      CALL PhaseIntegral(DATN,ISTATE,NPPS,R1,R2,1,1,EO,DEDPK,I11)

      l = 0
      DO j = 1, ISTATE-1
        l = l + NPPS(j)
      END DO

*** Calculates the width
      dvdGv = (FCTOR*I00(l+1))/(2*Pi)
      width = 1.0d0/(2*Pi * dvdGv) * EXP(-2.0d0*FCTOR*I0m(l+1))
***      width = 1.0d0/(FCTOR*I00(l+1)) * EXP(-2.0d0*FCTOR*I0m(l+1))

*** Calculates the partial derivatives dWdPk(i)
      DO IPV = l+1, l+NPPS(ISTATE)

        WRITE (33,170) vb,Jr,EO,IPV-l-1,DEDPK(IPV)
 170    FORMAT (' v =',I3,' J =',I3,' E =',G16.8,I2,' DEDPK =',G16.8)
        WRITE (34,180) vb,Jr,EO,I0m(IPV),I00(IPV),I10(IPV),I11(IPV)
 180    FORMAT ('v =',I3,' J =',I3,' E =',G16.8,' I0m =',G14.8,' I00 =',
     &         G14.8,' I10 =',G14.8,' I11 =',G14.8)

******7************* Calculate I11 by difference *********************72
cc        I0 = I00(l+1) 
cc        CALL DeltaBeta(ISTATE,IPV-1,deltaP,.TRUE.)
cc        CALL locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
cc        CALL PhaseIntegral(DATN,ISTATE,NPPS,R1,R2,0,0,EO,DEDPK,I00)
cc        dvdGv2 = (FCTOR*I00(l+1))/(2*Pi)
cc        I11(IPV) = (I00(l+1)-I0)/deltaP
cc        WRITE (35,185) vb,Jr,EO,IPV-1,dvdGv2,dvdGv,(dvdGv2-dvdGv)/deltaP
cc 185    FORMAT ('v =',I3,' J =',I3,' E =',G16.8,I3,' dv/dGv2 =',F12.8,
cc     &         ' dv/dGv =',F12.8,' d(dv/dGv)/dpi =',G16.8)
cc        CALL DeltaBeta(ISTATE,IPV-1,deltaP,.FALSE.)
******7********************** BLOCK END ******************************72

        IF (NOWIDTHS .EQ. 0) THEN
          dWdPk(IPV) = -I10(IPV)/I00(IPV) * EXP(-2*FCTOR*I0m(IPV))
     &    -1.0/(2*FCTOR) * I11(IPV)/I00(IPV)**2 * EXP(-2*FCTOR*I0m(IPV))
        ELSE
          dWdPk(IPV) = -I10(IPV)/I00(IPV) * EXP(-2*FCTOR*I0m(IPV))
        END IF
      END DO
      IF (PSEL(ISTATE) .EQ. 0) THEN
        DO IPV = l+1, l+NPPS(ISTATE)
          dWdPk(IPV) = 0.0d0
        END DO
      END IF
 195  IF (R1.LT.0 .OR. R2.LT.0 .OR. R3.LT.0) THEN
        width = FREQ(DATN)
        DO IPV = l+1, l+NPPS(ISTATE)
          dWdPk(IPV) = 0.0d0
        END DO 
        WRITE (6,197) DATN
 197    FORMAT(' <<Energy of Datum ',I5,' is above barrier maximum, henc
     &e eliminated>>')
      END IF

      WRITE (35,200) vb,Jr,EO,width,1.0d0/dvdGv
 200  FORMAT ('v =',I3,' J =',I3,' E =',G14.8,' width =',G12.6,
     &       ' dGv/dv =',G16.8) 

******7***************************************************************72
* Calculate the partial derivative of dG/dv by difference for testing 
* the Pajunen method
c      DO IPV = l+1, l+NPPS(ISTATE)
c        CALL DeltaBeta(ISTATE,IPV-1,deltaP,.TRUE.)
c        CALL locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
c        CALL PhaseIntegral(DATN,ISTATE,NPPS,R1,R2,0,0,EO,DEDPK,I00)
c        dvdGv2 = (FCTOR*I00(l+1))/(2*Pi)
c        WRITE (35,205) vb,Jr,EO,IPV-1,dvdGv2,dvdGv,(dvdGv2-dvdGv)/deltaP
c 205    FORMAT ('v =',I3,' J =',I3,' E =',G16.8,I3,' dv/dGv2 =',F12.8,
c     &         ' dv/dGv =',F12.8,' d(dv/dGv)/dpi =',G16.8) 
c        WRITE (35,*) 'deltaP = ',deltaP
c        CALL DeltaBeta(ISTATE,IPV-1,deltaP,.FALSE.)
c        CALL locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
c        CALL PhaseIntegral(DATN,ISTATE,NPPS,R1,R2,0,0,EO,DEDPK,I00)
c      END DO
******7********************** BLOCK END ******************************72
      RETURN
      END

******7***************************************************************72
      SUBROUTINE locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
******7***************************************************************72
*  Subroutine locateTP locates the turning points R1, R2 and R3 for    
*  a given energy E.                                                   
******7***************************************************************72
*  To find the R1, R2 and R3 for a given energy E, first obtain trial R
*  values by searching the trial energy on the potential array passed 
*  from the common block from RMIN to RMAX, then compare it with the 
*  given energy. This process is repeated until convergence is achieved.
******7***************************************************************72 
*** On entry
c    DATN    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    vb      is the vibrational quantum number.
c    Jr      is the rotational quantum number.
c    EO      is the energy for this level (DE in cm-1).
*** On exit
*  R1          is the array of inner turning points.                 
*  R2          is the array of second turning points.                
*  R3          is the array of third turning points.                 
**********************************************************************72
*  Re          is the Equilibrium Distance for each state.          
*  R1Up        is the upper trial R1.                              
*  R1Down      is the lower trial R1.                                
*  R2Up        is the upper trial R2.                                
*  R2Down      is the lower trial R2.                                
*  R3Up        is the upper trial R3.                                
*  R3Down      is the lower trial R3.                                
*  RH          is the mesh distance used to locate the trial R.   
******7***************************************************************72
      INCLUDE 'arrsizes.h'
c
c** Type statements & common block for data
c
      REAL*8 FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)

      INTEGER COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)

      CHARACTER*2 NAME(2)

      CHARACTER*2 SLABL(-2:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL

**************************
      REAL*8 RR(1),RM2(1),VV(1),VLIMT
      INTEGER NCNN

c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      REAL*8 R(NPNTMX,NSTATEMX),RMAX(NSTATEMX)
      COMMON /RBLK/R
      COMMON /RMAXBLK/RMAX

      INTEGER I,ISTATE,DATN,vb,Jr,index

      REAL*8 R1,R2,R3,R1Up,R1Down,R2Up,R2Down,R3Up,R3Down,EO,
     &  V1Up,V1Down,V2Up,V2Down,V3Up,V3Down,bMi,bMo
      REAL*8 prevEMinusV,currEMinusV,Rmid,Vmid,Vtemp,Vtotal(NPNTMX)

      COMMON /VBLIK/Vtotal

******7************** Locating inner turning point R1 ******************
      index = 1
  10  IF (EO .LT. Vtotal(index)) THEN
        index = index + 1
        IF (index .GT. NDATPT(ISTATE)) THEN
          WRITE(6,80) vb,Jr
  80      FORMAT('The inner turning point R1 of level v =',I3,
     &           ' J =',I3,' is greater than RMAX')
          R1 = -1.0d0
          GO TO 600
        END IF
        GO TO 10
      END IF
      IF (EO .EQ. Vtotal(index)) THEN
        R1 = R(index,ISTATE)
        GO TO 22
      END IF
      R1Down = R(index,ISTATE)
      R1Up = R(index-1,ISTATE)
      prevEMinusV = Vtotal(index-1) - Vtotal(index) 
  20  Rmid = (R1Up+R1Down)/2.0d0
      IF (PSEL(ISTATE) .NE. 0) THEN
        CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.FALSE.)
      ELSE
        RR(1) = Rmid
        RM2(1) = 1.0d0 / (RR(1) ** 2)
       CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,VV,NCNN)
        Vmid = VV(1) 
        CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.TRUE.)
      END IF
      currEMinusV = EO - Vmid
*      WRITE (30,21) R1Up,Rmid,R1Down,EO,Vmid
*  21  FORMAT ('R1U=',F8.6,' Rmid=',F8.6,' R1Down=',F8.6,' E=',G14.7,
*     &' Vmid=',G14.7)
      IF (DABS(currEMinusV) .GE. DABS(prevEMinusV) .AND. 
     &   ((currEMinusV.GE.0 .AND. prevEMinusV.GE.0) .OR. 
     &   (currEMinusV.LE.0 .AND. prevEMinusV.LE.0))) THEN
        R1 = Rmid
      ELSE IF (currEMinusV .GT. 0) THEN
        R1Down = Rmid
        prevEMinusV = currEMinusV
        GO TO 20 
      ELSE IF (currEMinusV .EQ. 0) THEN
        R1 = Rmid
      ELSE
        R1Up = Rmid
        prevEMinusV = currEMinusV
        GO TO 20
  22  END IF

      WRITE(30,100) vb,Jr,R1,EO,currEMinusV,prevEMinusV
 100  FORMAT('Level v =',I3,' J =',I3,' R1 =',F8.5,' E =',G16.8,
     &  ' EmVc =', G16.8,' EmVp =', G16.8)
******7************** Locating second turning point R2 *****************
      index = IDNINT((RE(ISTATE) - RMIN(ISTATE))/RH(ISTATE))
  25  IF (EO .GT. Vtotal(index)) THEN
        index = index + 1
        IF (index .GT. NDATPT(ISTATE)) THEN    
          WRITE(6,200) vb,Jr
 200      FORMAT('The second turning point R2 of level v =',I3,
     &           ' J =',I3,' is greater than RMAX')
          R2 = -1.0d0
          GO TO 600
        END IF
        GO TO 25
      END IF
      IF (EO .EQ. Vtotal(index)) THEN
        R2 = R(index,ISTATE)
        GO TO 35
      END IF
      R2Up = R(index,ISTATE)
      R2Down = R(index-1,ISTATE)
      prevEMinusV = Vtotal(index) - Vtotal(index-1)
  30  Rmid = (R2Up+R2Down)/2.0d0
      IF (PSEL(ISTATE) .NE. 0) THEN
        CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.FALSE.)
      ELSE
        RR(1) = Rmid
        RM2(1) = 1.0d0 / (RR(1) ** 2)
       CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,VV,NCNN)
        Vmid = VV(1)
        CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.TRUE.)
      END IF
      currEMinusV = EO - Vmid
      IF (DABS(currEMinusV) .GE. DABS(prevEMinusV) .AND. 
     &   ((currEMinusV.GE.0 .AND. prevEMinusV.GE.0) .OR. 
     &   (currEMinusV.LE.0 .AND. prevEMinusV.LE.0))) THEN
        R2 = Rmid
      ELSE IF (currEMinusV .GT. 0) THEN
        R2Down = Rmid
        prevEMinusV = currEMinusV
        GO TO 30
      ELSE IF (currEMinusV .EQ. 0) THEN
        R2 = Rmid
      ELSE 
        R2Up = Rmid
        prevEMinusV = currEMinusV
        GO TO 30
  35  END IF

      WRITE(30,250) vb,Jr,R2,EO,currEMinusV,prevEMinusV
 250  FORMAT('Level v =',I3,' J =',I3,' R2 =',F8.5,' E =',G16.8,
     &  ' EmVc =', G16.8,' EmVp =', G16.8)
******7************** Locating third turning point R3 ******************
  40  IF (EO .LT. Vtotal(index)) THEN
        index = index + 1
        IF (index .GT. NDATPT(ISTATE)) THEN
          WRITE(6,300) vb,Jr
 300      FORMAT('The third turning point R3 of level v =',I3,
     &           ' J =',I3,' is greater than RMAX')
          R3 = -1.0d0
          GO TO 600
        END IF
        GO TO 40
      END IF
      IF (EO .EQ. Vtotal(index)) THEN
        R3 = R(index,ISTATE)
        GO TO 600
      END IF
      R3Down = R(index,ISTATE)
      R3Up = R(index-1,ISTATE)
      prevEMinusV = Vtotal(index-1) - Vtotal(index)   
  50  Rmid = (R3Up+R3Down)/2.0d0
      IF (PSEL(ISTATE) .NE. 0) THEN
        CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.FALSE.)
      ELSE
        RR(1) = Rmid
        RM2(1) = 1.0d0 / (RR(1) ** 2)
       CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,VV,NCNN)
        Vmid = VV(1)
        CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.TRUE.)
      END IF
      currEMinusV = EO - Vmid
      IF (DABS(currEMinusV) .GE. DABS(prevEMinusV) .AND. 
     &   ((currEMinusV.GE.0 .AND. prevEMinusV.GE.0) .OR. 
     &   (currEMinusV.LE.0 .AND. prevEMinusV.LE.0))) THEN
        R3 = Rmid
      ELSE IF (currEMinusV .GT. 0) THEN
        R3Down = Rmid
        prevEMinusV = currEMinusV
        GO TO 50
      ELSE IF (currEMinusV .EQ. 0) THEN
        R3 = Rmid
      ELSE 
        R3Up = Rmid
        prevEMinusV = currEMinusV
        GO TO 50
      END IF
  
      WRITE(30,350) vb,Jr,R3,EO,currEMinusV,prevEMinusV
 350  FORMAT('Level v =',I3,' J =',I3,' R3 =',F8.5,' E =',G16.8,
     &  ' EmVc =', G16.8,' EmVp =', G16.8/)
 600  RETURN
******7***************************************************************72
      END 

******7***************************************************************72
      SUBROUTINE PhaseIntegral(DATN,ISTATE,NPPS,RI,RO,n,k,EO,DEDPK,Ink)
******7***************************************************************72
*  Subroutine PhaseIntegral performs the calculations of the phase     
*  integrals in the case of two adjacent turning points, when integrand
*  has the form f(R)/(E-V(R))**(k+1/2), for k=-1,0 or k>0 (high order).
*  f(R) is some linear combination of powers of R times powers of      
*  (E-V(R)) and its derivatives, where E is energy and V(R) any smooth 
*  potential function.                                                 
*  The formula used to calculate the phase integral Ink is
*                              N
*           Ink = ½ (R2 - R1) Sum(Wi·F(Zi))
*                             i=1
*  in which Wi are the weights.
*  Ref. JCP, 71(6), 2618, 1979
******7***************************************************************72
*** On entry
c    DATN    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    NPPS(i) are the number of parameters varied in each state.
c    RI      is the inner turning point
c    RO      is the outer turning point
c    EO      is the energy for this level (DE in cm-1).
c    n, k    are the powers in the phase integral Ink
c    DEDPK(i) are the values of the partial derivative dE/dP(k).
*** On exit
*   Ink     is the phase integral
******7***************************************************************72
      INCLUDE 'arrsizes.h'

      INTEGER n,k,i,j,l,M,p
 
      REAL*8 Vdi,deltaP

      PARAMETER (M = 32)

      INTEGER ISTATE, DATN, NPPS(NSTATEMX), IPV
 
      REAL*8 RI,RO,EO,Ink(NPARMX),FZi(M,NPARMX),Rd(M),Vd(M),EMinusV,
     &       In1(NPARMX),FZi1(M,NPARMX),FCTOR,DEDPK(NPARMX),bMi,bMo

      REAL*8 Zi(M),Wi(M),Zp11(M),Wp11k1(M),Wp11k3(M),Zp15(M),Wp15k1(M),
     &       Wp15k3(M),Zp21(M),Wp21k1(M),Wp21k3(M),Zgi(M),Wgi(M)        

      DATA Zp11/0.9659258262890683, 0.8660254037844386,
     1          0.7071067811865475, 0.5000000000000000,
     2          0.2588190451025208, 
     3          0.0000000000000000,
     4         -0.2588190451025208,-0.5000000000000000,
     5         -0.7071067811865475,-0.8660254037844386,
     6         -0.9659258262890683,
     c          0.0000000000000000,
     c          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c          0.0,0.0,0.0,0.0/

      DATA Wp11k1/22.06633397656787,-37.69911184307752,
     1            35.60471674068432,-37.69911184307752,
     2            36.57672889044161,
     3           -37.69911184307752,
     4            36.57672889044161,-37.69911184307752,
     5            35.60471674068432,-37.69911184307752,
     6            22.06633397656787,
     c            0.000000000000000,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0,0.0,0.0/
 
      DATA Wp11k3/228.3214545745810, -720.4719152232592,
     1           1139.3509357018984,-1357.1680263507907,
     2           1447.1946273399755,
     3          -1474.4541520848097,
     4           1447.1946273399755,-1357.1680263507907,
     5           1139.3509357018984, -720.4719152232592,
     6            228.3214545745810,
     c            0.000000000000000,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0,0.0,0.0/           
      
      DATA Zp15/0.9807852804032304, 0.9238795325112868, 
     1          0.8314696123025452, 0.7071067811865475,
     2          0.5555702330196022, 0.3826834323650898,
     3          0.1950903220161283, 
     4          0.0000000000000000,
     5         -0.1950903220161283,-0.3826834323650898,
     6         -0.5555702330196022,-0.7071067811865475,
     7         -0.8314696123025452,-0.9238795325112868,
     8         -0.9807852804032304,
     c          0.0000000000000000,
     c          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

      DATA Wp15k1/29.62981929591175,-50.26548245743668,
     1            47.72092686124880,-50.26548245743664,
     2            49.12943331558201,-50.26548245743658,
     3            49.44900912828539,
     4           -50.26548245743656,
     5            49.44900912828539,-50.26548245743658,
     6            49.12943331558201,-50.26548245743664,
     7            47.72092686124880,-50.26548245743668,
     8            29.62981929591175, 
     c            0.000000000000000,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

      DATA Wp15k3/1164.639428963841,-3580.432803129281,
     1            5525.620597073791,-6534.512719466765,
     2            7020.542275282852,-7276.911407677031,
     3            7400.700330802901,
     4           -7439.291403700618,
     5            7400.700330802901,-7276.911407677031,
     6            7020.542275282852,-6534.512719466765,
     7            5525.620597073791,-3580.432803129281,
     8            1164.639428963841,
     c            0.000000000000000,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

      DATA Zp21/0.9898214418809327, 0.9594929736144974,
     1          0.9096319953545184, 0.8412535328311812,
     2          0.7557495743542583, 0.6548607339452851,
     3          0.5406408174555976, 0.4154150130018864,
     4          0.2817325568414297, 0.1423148382732851,
     5          0.0000000000000000,
     6         -0.1423148382732851,-0.2817325568414297,
     7         -0.4154150130018864,-0.5406408174555976,
     8         -0.6548607339452851,-0.7557495743542583, 
     9         -0.8412535328311812,-0.9096319953545184,
     9         -0.9594929736144974,-0.9898214418809327,  
     c          0.0000000000000000,
     c          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c          0.0,0.0/

      DATA Wp21k1/40.91258980361040,-69.11503837897816,
     1            65.80507790523560,-69.11503837898373,
     2            67.78308420797106,-69.11503837899778,
     3            68.30792716563759,-69.11503837900795,
     4            68.49459295516724,-69.11503837901213,
     5            68.54383971474920,
     6           -69.11503837901213, 68.49459295516724,
     7           -69.11503837900795, 68.30792716563759,
     8           -69.11503837899778, 67.78308420797106,
     9           -69.11503837898373, 65.80507790523560,
     9           -69.11503837897816, 40.91258980361040,
     c            0.000000000000000,
     c            0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c            0.0,0.0/ 
     
      DATA Wp21k3/6364.91821744174,-19267.83229098791,
     1           29317.26172550868,-34478.42376106038,
     2           37049.70046340271,-38499.30026119479,
     3           39357.74970048209,-39887.54536375314,
     4           40205.64256060925,-40378.03411697539,
     5           40431.72625305376,
     6          -40378.03411697539, 40205.64256060925,
     7          -39887.54536375314, 39357.74970048209,
     8          -38499.30026119479, 37049.70046340271,
     9          -34478.42376106038, 29317.26172550868,
     9          -19267.83229098791,  6364.91821744174,
     c           0.000000000000000,
     c           0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c           0.0,0.0/

      DATA Zgi/0.095012509837637440185, 0.281603550779258913230,
     1         0.458016777657227386342, 0.617876244402643748447,
     2         0.755404408355003033895, 0.865631202387831743880,
     3         0.944575023073232576078, 0.989400934991649932596,
     4        -0.095012509837637440185,-0.281603550779258913230,
     5        -0.458016777657227386342,-0.617876244402643748447,
     6        -0.755404408355003033895,-0.865631202387831743880,
     7        -0.944575023073232576078,-0.989400934991649932596,
     c         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

      DATA Wgi/0.189450610455068496285, 0.182603415044923588867,
     1         0.169156519395002538189, 0.149595988816576732081,
     2         0.124628971255533872052, 0.095158511682492784810,
     3         0.062253523938647892863, 0.027152459411754094852,
     4         0.189450610455068496285, 0.182603415044923588867,
     5         0.169156519395002538189, 0.149595988816576732081,
     6         0.124628971255533872052, 0.095158511682492784810,
     7         0.062253523938647892863, 0.027152459411754094852,
     c         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     c         0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/

c
c** Type statements & common block for data
c
      REAL*8 FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)

      INTEGER COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)

      CHARACTER*2 NAME(2)

      CHARACTER*2 SLABL(-2:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL
**************************
      REAL*8 RR(1),RM2(1),VV(1),VLIMT
      INTEGER NCNN

      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      REAL*8 R(NPNTMX,NSTATEMX),dVdPk(NPOTMX),Vpt(NPNTMX,NSTATEMX)
      COMMON /RBLK/R
      COMMON /dVdPkBLK/dVdPk
c *** Partial derivatives
      REAL*8 DVBODP(NPNTMX,NPOTMX),DUADP(NPNTMX,NPOTMX),
     1  DUBDP(NPNTMX,NPOTMX),DQADP(NPNTMX,NPOTMX),DQBDP(NPNTMX,NPOTMX)
     2  ,dLdP(NPNTMX,NPOTMX),DBDP(NPNTMX,NPOTMX)

      COMMON /PDVBLK/DVBODP,DUADP,DUBDP,DQADP,DQBDP,DBDP
      COMMON /PdVLBLK/dLdP
      
c      REAL*8 VdT
c      LOGICAL out
c      COMMON /VdTBLK/VdT,out

      DATA Pi/3.141592653589793D0/

      LOGICAL StGsQd,VrGsQd,PajunenK1
******7***************************************************************72
* Using different combinations of logical variable values of StGsQd 
* and VrGsQd can choose different method for integration. 
* e.g.
* If StGsQd = .FALSE. and VrGsQd = .FALSE., 
*    Pajunen method is used for all integrals 
* If StGsQd = .FALSE. and VrGsQd = .TRUE., 
*    Variant Gaussian Quadrature method is used for k = -1 and k = 0
*    Pajunen method is used for k = 1
* If StGsQd = .TRUE. and VrGsQd = .FALSE.,
*    Standard Gaussian Quadrature is used for k = -1
*    Pajunen method is used for k = 0 and k = 1
* If StGsQd = .TRUE. and VrGsQd = .TRUE.,
*    Standard Gaussian Quadrature is used for k = -1
*    Variant Gaussian Quadrature method for k = 0
*    Pajunen method is used for k = 1
***
* In Pajunen method, there is a option to choose different points of 
* weight. To do this, just change the variable names for the Zp* and 
* Wp*, e.g.
* Instead of using 15 points
*                   Zi(i) = Zp15(i) and Wi(i) = -Wp15k1(i), 
* replace them with
*                   Zi(i) = Zp11(i) and Wi(i) = -Wp11k1(i)
* for 11 points.
* Note: when using Pajunen method, the signs of the weights have to be 
*       flipped to get the right values
* Another option in Pajunen method is use the k = 3 points to fake
* k = 1 points; to do this, set PajunenK1 = .FALSE.
******7***************************************************************72  
      PajunenK1 = .TRUE.
      StGsQd = .TRUE.
      VrGsQd = .TRUE.

      FCTOR = SQRT(ZMU/16.857629206d0)

******7************ Choose quadrature method *************************72
* For k = -1, if using Standard Gaussian Quadrature (25.4.29, N.A.)
      IF (k .EQ. -1) THEN
        IF (StGsQd) THEN
          DO i = 1, M
            Zi(i) = Zgi(i)
            Wi(i) = Wgi(i)
          END DO
c *** If Variant Gaussian Quadrature method are used (25.4.38, N.A.)
        ELSE IF (VrGsQd) THEN
          DO i = 1, M
            Zi(i) = DCOS((2.0d0*i-1.0d0)*Pi/(2.0d0*M))
            Wi(i) = Pi/M
          END DO
c *** If Pajunen quadrature points are used (25.4.38, N.A.)
        ELSE
          DO i = 1, M
            Zi(i) = Zp21(i)
            Wi(i) = -Wp21k1(i)
          END DO
        END IF

* For k = 0, if Variant Gaussian Quadrature are used (25.4.38, N.A.)
      ELSE IF (k .EQ. 0) THEN 
        IF (VrGsQd) THEN 
          DO i = 1, M
            Zi(i) = DCOS((2.0d0*i-1.0d0)*Pi/(2.0d0*M))
            Wi(i) = Pi/M
          END DO
c *** If Pajunen quadrature points are used (25.4.38, N.A.)
        ELSE
          DO i = 1, M
            Zi(i) = Zp21(i)
            Wi(i) = -Wp21k1(i)
          END DO
        END IF  

* For k = 1, use Pajunen quadrature points (JCP, 71(6), 2618, 1979)
      ELSE 
******* IF k = 1 points are used
        IF (PajunenK1) THEN
          DO i = 1, M
            Zi(i) = Zp21(i)
            Wi(i) = -Wp21k1(i)
          END DO
******* IF k = 3 points are used (to fake)
        ELSE
          DO i = 1, M
            Zi(i) = Zp21(i)
            Wi(i) = -Wp21k3(i)
          END DO
        END IF
      END IF
******7********************** BLOCK END ******************************72

c *** Calculating distance array
 366  DO i = 1, M
        Rd(i) = (0.5)*(RO+RI) + (0.5)*(RO-RI)*Zi(i)
        Vd(i) = 0
*        WRITE(31,370) RI,RO,Rd(i)
* 370    FORMAT ('RI =',F8.5,' RO =',F8.5,' Rd(i) =',F8.5)
      END DO
*      WRITE (31,*)

      l = 0
      DO j = 1, ISTATE-1
        l = l + NPPS(j)
      END DO
      IF (PSEL(ISTATE) .EQ. 0 .OR. NPPS(ISTATE) .EQ. 0) THEN
        NPPS(ISTATE) = NPPS(ISTATE) + 1
      END IF 
      DO IPV = l+1, l+NPPS(ISTATE)
******7********************* TESTING *********************************72
* Calculate the partial derivative of dV/dPk by difference for testing
* 
c        CALL VGEN(ISTATE,-1.0d0,Vdi,DATN)
c        DO p = 1, NDATPT(ISTATE)
c          Vpt(p,ISTATE) = VPOT(p,ISTATE)
c          WRITE (39,371) R(p,ISTATE),VPOT(p,ISTATE),DVBODP(p,IPV)
c  371     FORMAT ('R(i) =',F10.4,' V =',F16.8,' DVBODP =',G20.10)
c        END DO
c        WRITE (39,*) NDATPT(ISTATE)
c        WRITE (39,*) 
c        CALL DeltaBeta(ISTATE,IPV-1,deltaP,.TRUE.)
c        CALL VGEN(ISTATE,-1.0d0,Vdi,DATN)
c        DO p = 1, NDATPT(ISTATE)
c          WRITE (39,372) R(p,ISTATE),VPOT(p,ISTATE),
c     &                   (VPOT(p,ISTATE)-Vpt(p,ISTATE))/deltaP
c  372     FORMAT ('R(i) =',F10.4,' V =',F16.8,' dVdPi ='G20.10)
c        END DO
c        CALL DeltaBeta(ISTATE,IPV-1,deltaP,.FALSE.)
******7********************* END TESTING *****************************72

        Ink(IPV) = 0.0d0
        In1(IPV) = 0.0d0
        DO i = 1, M
 
******7********************* TESTING *********************************72
* Calculate the partial derivative of dV/dPk by difference for testing
* 
*          out = .TRUE.
*          CALL VGEN(ISTATE,Rd(i),Vd(i),DATN)
*          Vdi = VdT
*          WRITE (40,373) Rd(i),Vd(i),VdT,dVdPk(IPV)
*  373     FORMAT('Rdi = ',F10.4,' Vdi   =',F20.10,19x,' VdT =',F20.10,
*     &          ' dVdPk =',F20.10)
*          CALL DeltaBeta(ISTATE,IPV-1,deltaP,.TRUE.)
*          CALL VGEN(ISTATE,Rd(i),Vd(i),DATN) 
*          WRITE (40,374) Rd(i),Vd(i),VdT,deltaP,(VdT-Vdi)/deltaP
*  374     FORMAT ('Rdi = ',F10.4,' Vd(i) =',F20.10,' VdT =',F20.10,
*     &           ' deltaP =',F9.6,' dV/dPk =',F20.10)
*          CALL DeltaBeta(ISTATE,IPV-1,deltaP,.FALSE.)
*          out = .FALSE.
******7********************* END TESTING *****************************72

c *** Locates the Rd(i) in the potential array
          IF (PSEL(ISTATE) .NE. 0) THEN
            CALL VGEN(ISTATE,Rd(i),Vd(i),DATN)
          ELSE
            RR(1) = Rd(i)
            RM2(1) = 1.0d0 / (RR(1) ** 2)
            CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,VV,
     &                                                             NCNN)
            Vd(i) = VV(1)
            CALL vValue(ISTATE,Rd(i),Vd(i),bMi,bMo,DATN,.TRUE.)
          END IF
          EMinusV = EO-Vd(i)
c *** For classical forbidden regions
          IF (EMinusV .LT. 0) EMinusV = -EMinusV

******7****************** Calculate FZi function *********************72
* For k = -1, if using Standard Gaussian Quadrature (25.4.29, N.A.)
          IF (k .EQ. -1) THEN
            IF (StGsQd) THEN
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n * DSQRT(EminusV)
c *** If Variant Gaussian Quadrature method are used (25.4.38, N.A.)
            ELSE IF (VrGsQd) THEN
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n *
     &                                DSQRT(EminusV)*DSQRT(1-Zi(i)**2) 
c *** If Pajunen quadrature method are used (25.4.38, N.A.)
            ELSE
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n *
     &                            DSQRT(EminusV)*(1-Zi(i)**2)**(1.5)
            END IF 
          ELSE IF (k .EQ. 0) THEN 
* For k = 0, if Variant Gaussian Quadrature are used (25.4.38, N.A.)
            IF (VrGsQd) THEN
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n *
     &                               ((1-Zi(i)**2)/EminusV)**(k+0.5)
c *** If Pajunen quadrature method are used (25.4.38, N.A.)
            ELSE
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n *
     &                        DSQRT((1-Zi(i)**2)/EminusV)*(1-Zi(i)**2)
            END IF 
* For k = 1, use Pajunen quadrature method (JCP, 71(6), 2618, 1979)
          ELSE
* If  k = 1 points are used
            IF (PajunenK1) THEN
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n *
     &                               ((1-Zi(i)**2)/EminusV)**(k+0.5)
* If  k = 3 points are used (to fake)
            ELSE
              FZi(i,IPV) = (dVdPk(IPV)-DEDPK(IPV))**n * (1-Zi(i)**2)**2
     &                                * ((1-Zi(i)**2)/EminusV)**(k+0.5)
            END IF
          END IF
******7********************** BLOCK END ******************************72

c *** Calculate the phase integral
          Ink(IPV) = Ink(IPV) + Wi(i)*FZi(i,IPV)
*          WRITE (32,400) EO,Zi(i),Wi(i),FZi(i,IPV),Ink(IPV),n,k
* 400      FORMAT('EO =',G12.6,' Zi(i) =',G18.12,' Wi(i) =',G18.12,
*     &          ' FZi =',G10.4,' Ink =',G10.4,' n =',I2,' k =',I2)

******7********************** Testing ********************************72
* To test if the Pajunen's method is valid, calculate the partial
* derivatives of dG/dv with respect to parameter pi, and compare its 
* values with the ones calculated by difference from other methods
          IF (k .EQ. 1) THEN 
            IF (PajunenK1) THEN
              FZi1(i,IPV)=dVdPk(IPV)**n*((1-Zi(i)**2)/EminusV)**(k+0.5)      
            ELSE
              FZi1(i,IPV)=dVdPk(IPV)**n * (1-Zi(i)**2)**2
     &                   *((1-Zi(i)**2)/EminusV)**(k+0.5)
            END IF
            In1(IPV) = In1(IPV) + Wi(i)*FZi1(i,IPV)
          END IF

        END DO
        IF (k .EQ. 1) THEN
          In1(IPV) = (0.25)*(RO-RI)*In1(IPV) * 1.0/(4*Pi) * FCTOR
          WRITE (37,401) EO,IPV-l-1,In1(IPV)
 401      FORMAT('EO =',G16.8,I3,' d(dv/dGv)/dpi =',G16.8)      
        END IF
******7********************* END TESTING *****************************72

        Ink(IPV) = (0.5)*(RO-RI)*Ink(IPV)
c *** In Pajunen's method, a factor of (1/2) is incorporated because 
c     the phase integral is a round integration
        IF (k .EQ. 1 .OR. ((.NOT.StGsQd) .AND. (.NOT.VrGsQd))) THEN 
          Ink(IPV) = (0.5)*Ink(IPV)
        END IF

        WRITE (32,403) EO,IPV-l-1,Ink(IPV),n,k
 403    FORMAT('EO =',G16.8,I3,' Ink =',G16.8,' n =',I2,' k =',I2)

      END DO
      IF (PSEL(ISTATE) .EQ. 0 .OR. NPPS(ISTATE) .EQ. 0) THEN
        NPPS(ISTATE) = NPPS(ISTATE) - 1
      END IF
      WRITE (32,*)

      RETURN
******7***************************************************************72
      END      

c**********************************************************************
      SUBROUTINE vValue(ISTATE,Rd,Vd,bMi,bMo,DATN,Forward)
c**********************************************************************
c** This subroutine will generate potential value Vd given a known 
c   distance variable Rd
c   Only the following types of potential form are included in current
c   version
c     - The Long-Range Expanded Morse Oscillator
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 2001 by Yiye Huang
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c
c    ISTATE  is the electronic state being fitted to.
c
c    Rd      is the distance, if Rd > 0, calculate that point potential
c            value Vd
c    DATN    is the date number
c    Forward if Forward = .TRUE., do forward calculation, Vd is an input  
c** On entry via common blocks:
c    PSEL   is used to choose the type of potential used.
c           |PSEL| = 0 : Do forward calculation.
c           |PSEL| = 1 : Use the Generalized Morse Potential.
c           |PSEL| = 2 : Use the Expanded Morse Potential.
c           |PSEL| = 3 : Use the Modified Morse Potential.
c           |PSEL| = 4 : Use the Modified Lennard-Jones Potential.
c           |PSEL| = 5 : Use the Generalized Potential Energy Function.
c           |PSEL| = 6 : Use the Expanded Morse Long-Range Potential. 
c           |PSEL| = 7 : Use the Modified Lennard-Jones-Gaussian Potential. 
c           |PSEL| = 8 : Use the Total Expanded Morse Potential.
c           |PSEL| = 9 : Total Expanded Morse Long-Range Potential(trail ver.).
c           |PSEL| = 10: Use the Total Expanded Morse Long-Range Potential.
c           < 0 : Read in parmaters are in fact the molecular
c                 constants for that molecule and the trial exponential
c                 paramters must be calculated from these values.
c    NBETA  is the number of beta parameters that will be used.
c    NCN    is the power of the dampening term.
c    nPw    is the power of the new variable z.
c    YPARM  is the array of parameters used for the potential.
c    DE     is the Dissociation Energy for each state.
c    RE     is the Equilibrium Distance for each state.
c    BETA   is the array of potential parameters for each state
c           or (if PSEL < 0) selected molecular constants.
c    NDATPT is the number of meshpoints used for the array.
c
c** On exit via common blocks:
c    Vd     is the potential value that is generated.
c
******7*****************************************************************
      INCLUDE 'arrsizes.h'
