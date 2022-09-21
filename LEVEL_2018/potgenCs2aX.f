c***********************************************************************
c     SUBROUTINE POTGENCs2Xa(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)
c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c-----------------------------------------------------------------------
c         **********    Version of  4 February 2013    **********
c      Potentials fitted to synthetic Hutson data from February 4
c
c** Subroutine to generate the effective adiabatic potential energy 
c  function for either of two electronic states of Cs2 
c  using the MLR potential functions:
c  of Baldwin & Le Roy [J.XXXxx (2012 )
c      for the   X(^1\Sigma_g^+) state  [case  ISTATE= 1]
c  and for the   a(^3\Sigma_u^+) state  [case  ISTATE= 2]
c-----------------------------------------------------------------------
c** Note that this is a special 'standalone' version of subroutine 
c  POTGEN from the general-purpose bound-state Franck-Condon code LEVEL.  
c  Replacing the calling statement in line 2 (above) with the (currently)
c  commented-out lines 4 & 5 allows this routine to be compiled and used
c  in the normal way with R.J. Le Roy's program LEVEL.  However, in that
c  case it is necessary to define parameter ISTATE internally at at
c  line 80           as ISTATE=1  for the X(^1\Sigma_g^+) state 
c  or at line 81     as ISTATE=2 for the a(^3\Sigma_u^+}) state.
c==============
c*** On INPUT:
c==============
c*  integer  ISTATE  specifies the choice of PEF/electronic state (above)
c*  integers IMN1 and IMN2  are the mass numbers for the 2 isotopologues
c*  integers IAN1 and IAN2  are the atomic numbers for the isotopologue
c*  integer  NPP  is the dimension of the  REAL*8 radial distance and 
c                 potential energy function arrays
c*  VLIM  is the (externally specified) asymptote energy (in cm-1) for
c         the {79,79}Br_2 isotopologue potential of the chosen state
c*  XO(i)  is the array of  NPP  radial distances (in Angst.) at which
c          the potential enrgy function is to be calculated
c==============
c*** On OUTPUT:
c==============
c*  VV(i)  the array of potential function values (in cm-1) 
c*  RM2(i) is the array of BOB-corrected radial centrifugal strength
c          factors  [1 + g_{na}^{\alpha}(r)]/XO(i)**2  defining the
c   effective centrifugal potential energy function for this isotopologue
c=======================================================================
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,IBOB,IAN1,IAN2,IMN1,IMN2,MN1R,MN2R,IORD,IORDD,IPOTL,
     1  PAD,QAD,PNA,NU1,NU2,NT1,NT2,NCMAX,PPAR,QPAR,NCN,NSR,NLR,NVARB,
     2  NPP,LNPT,GNS,GEL, NCMM,IVSR,LVSR,IDSTT,KDER,MM1, MMLR(9)
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,Asw,Rsw,BETA,BINF,B1,B2,BT,CSAV,U1INF,
     1 U2INF,T1INF,T2INF,YPAD,YQAD,YQADSM,YPNA,YPNASM,ABUND,CNN,
     2 DSCM,DX,DX1,FCT,FC1,FC2,FG1,FG2,MASS1,MASS2,RMASS1,RMASS2,REQ,
     3 Rref,Rinn,Rout,SC1,SC2,SG1,SG2,VLIM,VMIN,XDF,X1,XS,XL,XP1,ZZ,ZP,
     4 ZQ,ZME,
     5 ULR,ULRe,rhoAB,REQP,DM(9),DMP(9),DMPP(9),CMM(9),T0,C6adj,C9adj,
     6 RM3,RET,RETsig,RETpi,RETp,RETm,BFCT2,PPOW,PVSR,
     7 U1(0:NBOB),U2(0:NBOB),T1(0:NBOB),T2(0:NBOB),PARM(50),
     8 XO(NPP),VV(NPP),RM2(NPP), bTT(-1:2),cDS(-2:0),bDS(-2:0)
      SAVE IBOB,IPOTL,IORD,IORDD,PPAR,QPAR,PAD,QAD,PNA,NSR,
     1 NLR,MMLR,NVARB,NCMM
      SAVE DSCM,REQ,Rref,PARM,U1,U2,T1,T2,CSAV,BINF,ALFA,Rsw,ZME,
     2 Rinn,Rout,ULR,ULRe,CMM
c** Define arrays for storing parameter for one or multiple states
      INTEGER ISTATE,PPARA(2),QPARA(2),NLRA(2),MMLRA(4,2),NU1A(2),
     1    NT1A(2)
      REAL*8 DSCMA(2),REQA(2),RrefA(2),CMMA(4,2),PARMA(25,2),
     1    U1A(10,2),T1A(10,2)
c** Damping function parameters for printout .....
      DATA bTT/2.44d0,2.78d0,3.126d0,3.471d0/
      DATA bDS/3.3d0,3.69d0,3.95d0/
      DATA cDS/0.423d0,0.40d0,0.39d0/
      SAVE bTT, bDS, cDS
c** Electron mass, as per 2006 physical constants
      DATA ZME/5.4857990943d-4/
c** when using this with the standard LEVEL program, select ISTATE case here
c*** ISTATE =1 for X(^1\Sigma_g^+),
c           =2 for A(^3\Pi_{1u})
c---------------------------------------------------------------------------
c     ISTATE= 1 
      ISTATE= 2 
ccc Specify parameter values for the two Cs2 states of interest
c*** ISTATE =1 for X(^1\Sigma_g^+),  =2 for a(^3\Sigma_u^+) Cs2,
      DATA PPARA/5,7/, QPARA/5,4/,MMLRA/6,8,10,0, 6,8,10,0/,
     1  NLRA/22,3/, NU1A/-1,-1/, NT1A/ -1,-1/ 
      DATA DSCMA/3.650029360426D+03, 2.785411387409D+02/,
     1  REQA/4.647967805457D+00, 6.277487013349D+00/,
     2  RrefA/6.2d0,13.21124d0/,
c** specify  up to 4 C_m long-range coefficients for up to 2 states
     4  CMMA/3.3188D+07,1.38d09,6.01d10,0.d0,3.3188D+07,
     5  1.38d9,6.01d10,0.d0/,
c** now specify up to 25 exponent parameters \beta_i for 2 states 
     6 PARMA/9.501531033220D-02,-3.722606328132D-01,-4.041904740197D-02,
     7  1.301668046083D-01,1.493847745066D-01,1.736528799714D-01,
     8  2.946191125909D-01,5.203014030305D-01,-1.128813164895D+00,
     9  -2.526260985324D+00,9.933897713444D+00,1.547620902776D+01,
     a  -4.187473012360D+01,-5.275307242344D+01,1.101103556968D+02,
     b  1.136284480558D+02,-1.757316156746D+02,-1.481380993045D+02,
     c  1.647118824973D+02,1.078196559724D+02,-8.087592626435D+01,
     d  -3.370965829120D+01,1.504982030319D+01,
     e  4*0.d0, -3.317815967582D-01,-1.602744130372D-02,
     f  -1.924060381445D-01,-2.596671685451D-01,19*0.d0/

c** Specify up to 10 'adiabatic' BOB parameters for up to 2 states
c     DATA U1A/0.371D0,9*0.D0, 0.050D0,9*0.D0/
c** Specify up to 10 centrifugal BOB parameters for up to 2 states
c     DATA T1A/10*0.D0, 0.D0,1.13D-03,2.945D-02,-1.809D-01,-5.71D-01,
c    1        8.22D+00,-2.934D+01,4.99D+01,-4.19D+01,1.4D+01/
      LNPT= 1
      IAN1= 55
      IAN2= 55
      IPOTL= 4
      IBOB= -1
      IVSR= -2
      IDSTT= 1
      NCMM= 3
c*** Note that rhoAB ia the same for the X and A states
      rhoAB= 0.434d0
      DSCM= DSCMA(ISTATE)
      REQ= REQA(ISTATE)
      Rref= RrefA(ISTATE)
c     IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
cc        READ(5,*) IPOTL, PPAR, QPAR, NSR, NLR, IBOB
cc        READ(5,*) DSCM, REQ, Rref
cc        IF(IPOTL.GE.4) THEN
c** For MLR, DELR or Tiemann-polynomial potentials .....
cc            READ(5,*) NCMM, IVSR, IDSTT, rhoAB
cc            READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
cc            ENDIF
c-----------------------------------------------------------------------
          PPAR= PPARA(ISTATE)
          QPAR= QPARA(ISTATE)
          NLR= NLRA(ISTATE)
          DO  I= 1,NCMM
              MMLR(I)= MMLRA(I,ISTATE)
              CMM(I)= CMMA(I,ISTATE)
              ENDDO
          NSR= NLR
          IORD= NLR
          NVARB= IORD+ 1
          DO  I= 1,IORD+1
              PARM(I)= PARMA(I,ISTATE)
              ENDDO
c-----------------------------------------------------------------------
cc        IF(NVARB.GT.0) READ(5,*) (PARM(I), I=1,IORD+1)
c-----------------------------------------------------------------------
          PNA= 3
c*** Set mass numbers for reference isotopologue
          MN1R= 133
          MN2R= 133
          NU1= NU1A(ISTATE)
          NU2= NU1
          NT1= NT1A(ISTATE)
          NT2= NT1
          PAD= 6
          IF(ISTATE.GE.2) PAD=5
          QAD= PAD
c-----------------------------------------------------------------------
cc            READ(5,*) MN1R, MN2R, PAD, QAD, NU1, NU2, PNA, NT1, NT2
c-----------------------------------------------------------------------
          NCMAX= MAX0(NU1,NU2,NT1,NT2)
c** If appropriate, read parameters & prepare to add mass-dep. BOB corrn
c         CALL MASSES(IAN1,IMN1,NAME1,GEL,GNS,MASS1,ABUND)
c         CALL MASSES(IAN1,MN1R,NAME1,GEL,GNS,RMASS1,ABUND)
c         CALL MASSES(IAN2,IMN2,NAME2,GEL,GNS,MASS2,ABUND)
c         CALL MASSES(IAN2,MN2R,NAME2,GEL,GNS,RMASS2,ABUND)
c  For simplicity, first zero out all correction function coefficients
          DO  I=0,NCMAX
              U1(I)= 0.d0
              U2(I)= 0.d0
              T1(I)= 0.d0
              T2(I)= 0.d0
              ENDDO
          FC1= 0.d0
          FC2= 0.d0
          FG1= 0.d0
          FG2= 0.d0
          U1INF= 0.d0
          U2INF= U1INF
          T1INF= 0.d0
          T2INF= T1INF
          IF(NU1.GE.0) THEN
              DO  I= 0, NU1
                  U1(I)= U1A(I+1,ISTATE)
c ....   for homonuclear molecule
                  U2(I)= U1(I)
                  ENDDO
              ENDIF
          IF(NT1.GE.0) THEN
              DO  I= 0, NT1
                  T1(I)= T1A(I+1,ISTATE)
c ....   for homonuclear molecule
                  T2(I)= T1(I)
                  ENDDO
               ENDIF
c=======================================================================
c** Read actual BOB polynomial expansion coefficients
c=======================================================================
          IF(NU1.GE.0) THEN
c... use Huang/Le Roy form for atom-1 adiabatic potential BOB radial fx.
c-----------------------------------------------------------------------
ccc                   READ(5,*) U1INF,(U1(I), I=0,NU1)
c-----------------------------------------------------------------------
              WRITE(6,630) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1        1,U1INF,PAD,PAD,PAD,PAD,PAD,PAD,NU1,QAD,QAD,QAD,QAD,QAD,
     2                                         NU1+1,(U1(I),I= 0,NU1)
              FC1= 1.d0 - RMASS1/MASS1
              ENDIF
          IF(NU2.GE.0) THEN
c... use Huang/Le Roy form for atom-2 adiabatic potential BOB radial fx.
c-----------------------------------------------------------------------
cc                    READ(5,*) U2INF,(U2(I), I=0,NU2)
c-----------------------------------------------------------------------
              WRITE(6,630) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1        1,U2INF,PAD,PAD,PAD,PAD,PAD,PAD,NU2,QAD,QAD,QAD,QAD,QAD,
     2                                         NU2+1,(U2(I),I= 0,NU2)
              FC2= 1.d0 - RMASS2/MASS2
              ENDIF
          IF(NT1.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-1 ...
c-----------------------------------------------------------------------
cc                    READ(5,*) T1INF,(T1(I), I=0,NT1)
c-----------------------------------------------------------------------
              WRITE(6,634) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1 1,T1INF,PNA,PNA,PNA,PNA,PNA,PNA,NT1,PNA,NT1+1,(T1(I),I= 0,NT1)
              FG1= RMASS1/MASS1
              ENDIF
          IF(NT2.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-2 ...
c-----------------------------------------------------------------------
cc            READ(5,*) T2INF,(T2(I), I=0,NT2)
c-----------------------------------------------------------------------
              WRITE(6,634) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1 2,T2INF,PNA,PNA,PNA,PNA,PNA,PNA,NT2,PNA,NT2+1,(T2(I),I= 0,NT2)
              FG2= RMASS2/MASS2
              ENDIF
          U1INF= U1INF*FC1
          U2INF= U2INF*FC2
          T1INF= T1INF*FG1
          T2INF= T2INF*FG2
c... Now generates scaled expansion parameters all at the same time!
          DO  I=0,NCMAX
              U1(I)= U1(I)*FC1
              U2(I)= U2(I)*FC2
              T1(I)= T1(I)*FG1
              T2(I)= T2(I)*FG2
              ENDDO
cc        ENDIF
c*** If OMEGA
c=======================================================================
c** Generate an MLR potential [as per J.Chem.Phys. 131, 204309 (2009)
c                                 or J.Mol.Spectrosc. XX, xxxx (2011)]
c=======================================================================
      IF(IPOTL.EQ.4) THEN
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= CMM(1)
              ULRe= 0.d0
              IF((NCMM.GE.3).AND.(MMLR(2).LE.0)) THEN
c** For special Aubert-Frecon Li2(^2P + ^2S) {3,0,6,6} type cases
                  RM3= 1.d0/REQ**3
c** Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor here is  2\pi/\lambda  for this case
                  RET= 9.36423830d-4*REQ
                  RETSig= DCOS(RET) + (RET)*DSIN(RET)
                  RETPi= RETSig - RET**2 *DCOS(RET)
                  RETp= RETSig + 0.5d0*RETPi
                  RETm= RETSig - 0.5d0*RETPi
                  C6adj= CMM(3) + CMM(1)**2/(4.d0*DSCM)
                  C9adj = CMM(1)*C6adj/(2.d0*DSCM)
                  IF(MMLR(2).EQ.0) THEN
c ... for Aubert-Frecon 2x2 treatment of {C3,C6,C8} for Alkali A-state
                      T0= CMM(1)*RETm*RM3
                      ULRe= 0.5d0*(-CMM(2) + CMM(1)*RETp*RM3)
                      IF(NCMM.GE.3) THEN
                          T0= T0 + C6adj*RM3**2
                          ULRe= ULRe + 0.5d0*C6adj*RM3**2
                          IF(NCMM.GE.4) THEN
                              T0= T0+ CMM(4)*(RM3/REQ)**2
                              ULRe= ULRe + 0.5d0*CMM(4)*(RM3/REQ)**2
     1                                                  + C9adj*RM3**3
                              ENDIF
                          ENDIF
                      T0= T0/3.d0
                      ULRe= ULRe+0.5d0*DSQRT((T0-CMM(2))**2+ 8.d0*T0**2)
                      ENDIF
                  IF(MMLR(2).EQ.-1) THEN
c ... extension for Li2(c) {3,0,6,6,8,8} 3x3 Aubert-Frecon case ...
                      CALL AF3X3LEV(REQ,CMM(2),CMM(1),C6adj,
     1                                               CMM(4),DSCM,ULRe)
                      ULRe= ULRe + C9adj*RM3**3
                      ENDIF
                ELSE
c*** for 'ordinary' NCMM-term MLR uLR(r) ...  with damping [if rhoAB > 0]
                  IF(rhoAB.GT.0.d0) THEN
                      KDER= 0
                      CALL dampF(REQ,rhoAB,NCMM,MMLR,IVSR,IDSTT,KDER,
     1                                                    DM,DMP,DMPP)
                      ENDIF
                  DO  J= 1,NCMM
                      IF(rhoAB.LE.0.d0) THEN
                          ULRe= ULRe + CMM(J)/REQ**MMLR(J)
                        ELSE
                          ULRe= ULRe + DM(J)*CMM(J)/REQ**MMLR(J)
                        ENDIF
                      ENDDO
                ENDIF
              BINF= DLOG(2.d0*DSCM/ULRe)
              WRITE(6,602) NCN,PPAR,QPAR,DSCM,REQ
c... use THEOCHEM/Huang form:  \beta(yp)= Binf*yp + [1-yp]*{power series in yq}
              WRITE(6,607) PPAR,PPAR,QPAR,NSR,NLR,IORD+1,
     1                                       (PARM(J),J= 1,IORD+1)
              IF(Rref.GT.0) THEN
                  WRITE(6,613) Rref
                ELSE
                  WRITE(6,615) REQ
                  Rref= REQ
                ENDIF
              IF(rhoAB.GT.0.d0) THEN
                  PVSR= 0.5d0*IVSR
                  IF(IDSTT.GT.0) THEN
                      PVSR= 0.5d0*IVSR
                      WRITE(6,664) rhoAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
                    ELSE
                      LVSR= IVSR/2
                      WRITE(6,666) rhoAB,LVSR,bTT(LVSR)
                    ENDIF
                ELSE
                  WRITE(6,668)
                ENDIF
              WRITE(6,617) BINF,MMLR(1),CMM(1),MMLR(1)
              IF(NCMM.GT.1) THEN
                  MM1= 2
                  IF(MMLR(2).LE.0) THEN
                      MM1= 3
                      IF(MMLR(2).EQ.0)
     1                         WRITE(6,623) MMLR(2),CMM(2),MMLR(2)
                      IF(MMLR(2).LT.0) WRITE(6,625) MMLR(2),CMM(2),0
                      ENDIF
                  DO  I= MM1,NCMM
                      IF(MMLR(I).LE.9) WRITE(6,619) MMLR(I),CMM(I)
     1                                                    ,MMLR(I)
                      IF(MMLR(I).GT.9) WRITE(6,621) MMLR(I),CMM(I)
     1                                                    ,MMLR(I)
                      ENDDO
                  ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
              ZZ= (XO(i)**PPAR- REQ**PPAR)/(XO(i)**PPAR+ REQ**PPAR)
              ZP= (XO(i)**PPAR-Rref**PPAR)/(XO(i)**PPAR+Rref**PPAR)
              ZQ= (XO(i)**QPAR-Rref**QPAR)/(XO(i)**QPAR+Rref**QPAR)
              BETA= 0.d0
              IF(ZZ.GT.0) IORDD= NLR
              IF(ZZ.LE.0) IORDD= NSR
              DO  J= IORDD,0,-1
                  BETA= BETA*ZQ+ PARM(J+1)
                  ENDDO
c  Calculate MLR exponent coefficient
              BETA= BINF*ZP + (1.d0- ZP)*BETA
              ULR= 0.d0
c** Calculate local value of uLR(r)
              IF((NCMM.GE.3).AND.(MMLR(2).LE.0)) THEN
c... For special Aubert-Frecon Li2(^2P + ^2S) {3,0,6,6} type case 
                  RM3= 1.d0/XO(I)**3
c... Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor here is  2\pi/\lambda  for Li2(A)
                  RET= 9.36423830d-4*XO(I)
                  RETSig= DCOS(RET) + (RET)*DSIN(RET)
                  RETPi= RETSig - RET**2 *DCOS(RET)
                  RETp= RETSig + 0.5d0*RETPi
                  RETm= RETSig - 0.5d0*RETPi
                  IF(MMLR(2).EQ.0) THEN
c... Aubert-Frecon 2x2 case - for A-state Li2
                      T0= CMM(1)*RETm*RM3
                      ULR= 0.5d0*(-CMM(2) + CMM(1)*RETp*RM3)
                      IF(NCMM.GE.3) THEN
                          T0= T0 + C6adj*RM3**2
                          ULR= ULR + 0.5d0*C6adj*RM3**2
                          IF(NCMM.GE.4) THEN
                              T0= T0+ CMM(4)*(RM3/XO(I))**2
                              ULR= ULR + 0.5d0*CMM(4)*(RM3/XO(I))**2
     1                                                  + C9adj*RM3**3
                              ENDIF
                          ENDIF
                      T0= T0/3.d0
                      ULR= ULR+ 0.5d0*DSQRT((T0- CMM(2))**2+ 8.d0*T0**2)
                      ENDIF
                  IF(MMLR(2).EQ.-1) THEN
c... for Aubert-Frecon 3x3 case yielding lowest (c state) energy
                      CALL AF3X3LEV(XO(I),CMM(2),CMM(1),C6adj,
     1                                                CMM(4),DSCM,ULR)
                      ULR= ULR + C9adj*RM3**3
                      ENDIF
                ELSE
c** For the 'regular' simple inverse-power sum case.
                  IF(rhoAB.GT.0.d0) CALL dampF(XO(I),rhoAB,NCMM,MMLR,
     1                                    IVSR,IDSTT,KDER,DM,DMP,DMPP)
                  DO  J= 1,NCMM
                      IF(rhoAB.LE.0.d0) THEN
                          ULR= ULR + CMM(J)/XO(I)**MMLR(J)
                        ELSE
                          ULR= ULR + DM(J)*CMM(J)/XO(I)**MMLR(J)
                        ENDIF
                      ENDDO
                ENDIF
              BETA= (ULR/ULRe)*DEXP(-BETA*ZZ)
              VV(I)= DSCM*(1.d0 - BETA)**2 - DSCM + VLIM
              ENDDO
          ENDIF
      IF(IBOB.GT.0) THEN
c=======================================================================
c** If appropriate, generate Born-Oppenheimer breakdown correction 
c      functions to rotationless and/or centrifugal potential(s) using
c      LeRoy/Huang radial functions ...
c=======================================================================
          DO  I=1,NPP
              YPAD= (XO(I)**PAD- REQ**PAD)/(XO(I)**PAD+ REQ**PAD)
              YQAD= (XO(I)**QAD- REQ**QAD)/(XO(I)**QAD+ REQ**QAD)
              YPNA= (XO(I)**PNA- REQ**PNA)/(XO(I)**PNA+ REQ**PNA)
              SC1= U1INF*YPAD
              SC2= U2INF*YPAD
              SG1= T1INF*YPNA
              SG2= T2INF*YPNA
              YQADSM= (1.d0- YPAD)
              YPNASM= (1.d0- YPNA)
c ... finally, accumulate overall BOB terms ... all at the same time!
              DO  J= 0,NCMAX
                  SC1= SC1+ YQADSM*U1(J)
                  SC2= SC2+ YQADSM*U2(J)
                  SG1= SG1+ YPNASM*T1(J)
                  SG2= SG2+ YPNASM*T2(J)
                  YQADSM= YQADSM*YQAD
                  YPNASM= YPNASM*YPNA
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          IF((IAN1.EQ.3).AND.(IAN2.EQ.3).AND.(NCMM.GE.3).AND.
     1                                            (MMLR(2).LE.0)) THEN
c!! For mixed isotopopogue {6,7}Li_2(A) state, shift asymptote!
              IF((IMN1.NE.IMN2).AND.(MMLR(2).EQ.0)) THEN
                  DO  I= 1,NPP
                      RM3= (2.d0/3.d0)*CMM(1)/XO(I)**3
                      VV(I)= VV(I)+ RM3- DSQRT(RM3**2+ 3.085959756d-02)
                      ENDDO
                  ENDIF
c** For special case of A and c states of Li2, add BOB centrifugal term
              BFCT2= 2.d0*16.857629206d0*(MASS1+MASS2)/(MASS1*MASS2)
              DO  I= 1, NPP
                  VV(I)= VV(I) + BFCT2/XO(I)**2
                  ENDDO
              ENDIF
          ENDIF
      RETURN
  602 FORMAT(/' MLR(n=',i1,'; p=',I1,', q=',I1,') Potential with:   De='
     1 ,F10.3,'[cm-1]    Re=',F12.8,'[A]')
  607 FORMAT('   with exponent coefficient   beta(r)= beta{INF}*y',I1,
     1  ' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x,'exponent coefft. powe
     2r series orders',I4,' for  R < Re  and',I4,' for  R > Re'/6x,
     3  'and',i3,' coefficients:',1PD16.8,2D16.8:/(10x,4D16.8:))
  613 FORMAT(6x,'with radial variables  y_p & y_q  defined w.r.t.',
     1  '  Rref=',F10.7)
  615 FORMAT(6x,'radial variables  y_p & y_q  defined w.r.t.  Rref= Re='
     1   F10.7)
  617 FORMAT('      betaINF=',f16.12,'  & uLR defined by  C',i1,' =',
     1  1PD13.6,'[cm-1 Ang','^',0PI1,']')
  619 FORMAT(50x,'C',I1,' =',1PD13.6,'[cm-1 Ang','^',0PI1,']')
  621 FORMAT(50x,'C',I2,'=',1PD13.6,'[cm-1 Ang','^',0PI2,']')
  623 FORMAT(3x,'Use Aubert-Frecon 2x2 model for uLR(r) with',
     1  4x,'C',I1,' =',1PD13.6,'[cm-1 Ang^',i1,']'/8x,'including retarda
     2tion & B(r) within V_{ad}, plus Delta(V)_{gu}^{(6,7)})')
  625 FORMAT(3x,'Use Aubert-Frecon 3x3 model for uLR(r) with',
     1  3x,'C',I2,' =',1PD13.6,'[cm-1 Ang^',i1,']'/6x,'including retarda
     2tion & B(r) within V_{ad}')
  630 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1 '  of mass ',f15.11/'   consists of mass factor  [1- MASS(',I3,
     2 A2,')/MASS(',I3,A2,')]  multiplying all of:'/5x,'u',i1,'INF=',
     3 f11.6,'  times  y',i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,
     4 ' + Re**',i1,')]'/5x,'plus  [1 - y',i1,']  times an order',I3,
     5 ' polynomial in'/7x,'y',i1,'=[(r**',i1,' - Re**',i1,')/(r**',i1,
     6 ' + Re**',i1,')]  with the ',i3,' coefficients:'/(3x,4G17.9:))
  634 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor  [MASS(',I3,A2,')/MASS(',I3,
     2 A2,')]  multiplying all of:'/5x,'q',i1,'INF=',F11.6,' times  y',
     3 i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',i1,')]'/
     4 5x,'plus [1 - y',i1,'] times an order',I3,' polynomial in y',i1,
     5 ' with the',i3,' coefficients:'/(3x,4G17.9:))
  664 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2,
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
  666 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     3)^k/k!}]   where   bTT=',f6.3,'*rhoAB')
  668 FORMAT(4x,'uLR inverse-power terms incorporate NO damping function
     1s')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c---------------------- RJL Version of 06 July 2010 --------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IDF requires damping to be defined s.th.  Dm(r)/r^m --> r^{IDF/2}
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c* KDER:  if KDER.GT.0  the first derivative is also calculated 
c*        if KDER.GT.1  the second derivative is also calculated 
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m) - The first derivative of the damping function  DM(m)
c  DMPP(m) - The second derivative of the damping function  DM(m)
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IDF,IDSTT,KDER,IDFF,FIRST,
     1  Lsr,m,MM,MMAX
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-2:0), cpm(20,-2:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.126d0,3.471d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.40d0,0.39d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
      IF(RHOab.LE.0) THEN
          WRITE(6,602) RHOab
          STOP
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          IF((IDF.LT.-4).OR.(IDF.GT.4)) THEN
                WRITE(6,600) IDSTT,IDF
                STOP
                ENDIF
          Lsr= IDF/2
          MMAX= MMLR(NCMM) + Lsr - 1
          aTT= RHOab*bTT(Lsr)
          br= aTT*r
          XP= DEXP(-br)
          SM(-3)= 0.d0
          SM(-2)= 0.d0
          SM(-1)= 0.d0
          SM(0)=  1.d0
          TK= 1.d0
          IF(br.GT.0.5d0) THEN
              DO  m= 1,MMAX
                  TK= TK*br/DFLOAT(m)
                  SM(m)= SM(m-1)+ TK
                  ENDDO
              DO m= 1, NCMM
                  MM= MMLR(m) - 1 + Lsr
                  DM(m)= 1.d0 - XP*SM(MM)
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                      IF(KDER.GT.1) DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                      ENDIF
                  ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r.  However, at very small r that algorithm
c  becomes unstable due to numerical noise.  To avoid this, if the 
c  argument is very small it is re-evaluated as a finite sum ...
c-----------------------------------------------------------------------
            ELSE
              MMAX= MMAX+5
              DO  m= 1, MMAX
c... NOTE that here SM(m) is the m'th term  (b*r)^m/m!  [not a sum]
                  SM(m)= SM(m-1)*br/DFLOAT(m)
                  ENDDO
              DO  m= 1, NCMM
                  MM= MMLR(m) + Lsr
                  DM(m)= XP*(SM(MM)+ SM(MM+1)+ SM(MM+2)+ SM(MM+3) 
     1                                                     + SM(MM+4))
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*SM(m-1)
                      IF(KDER.GT.1)DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                      ENDIF
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IDF.LT.-4).OR.(IDF.GT.0)) THEN
              WRITE(6,600) IDSTT,IDF
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -2,0
                      bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IDF) + cpm(MM,IDF)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
c... Actually ...  DM(m)= YP**(MM + IDF/2)  :  set it up this way to 
c   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IDF.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IDF.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
              IF(KDER.GT.0) THEN
                  TK= bpm(MM,IDF) + 2.d0*cpm(MM,IDF)*br
                  DMP(m) = ZK*XP*rhoAB*TK*DM(m)/YP
                  IF(KDER.GT.1) THEN
c ... if desired ... calculate second derivative [for DELR case] {check this!}
                      DMPP(m)= (ZK-1.d0)*XP*TK*DMP(m)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IDF)*rhoAB**2/TK
                      ENDIF
                  ENDIF
              ENDDO   
          ENDIF  
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IDF=',i3,'  no dampi
     1ng function is defined')
  602 FORMAT( /,' ***ERROR ***  rhoAB=', F7.4,'  yields an invalid Dampi
     1ng Function definition')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c=======================================================================
      SUBROUTINE AF3X3LEV(RDIST,DELTAE,C3val,C6val,C8val,De,ULR)
c=======================================================================
c*** Simplified version of subroutiune AF3x3potRet which locates the
c  lowest eigenvalue of the 3x3 long-range Li2 interaction matrix of 
c  Eq.(25) of J.Mol.Spectrosc. 268, 199 (2011), and does not return 
c  derivatives w.r.t. parameters
c*****               Version of  10 Sept. 2011
c-----------------------------------------------------------------------
      REAL*8  H(3,3),DM1(3,3),DM3(3,3),DM5(3,3),Q(3,3), EIGVEC(3,1),
     1 RESID(3,1), W(3) 
      REAL*8  RDIST,RDIST2,RDIST3,DELTAE,C3val,C6val,C8val,De,ULR,
     1   RET,RETSig,RETPi,Modulus,M1,M3,M5,Z
      INTEGER  I,J,L,K
      M1= C3val
      M3= C6val
      M5= C8val
      RET= 9.36423830d-4*RDIST
      RETSig= DCOS(RET) + (RET)*DSIN(RET)
      RETPi= RETSig - RET**2 *DCOS(RET)
      RDIST2= RDIST**2
      RDIST3= RDIST*RDIST2
c      WRITE(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)" ' 
c      WRITE(25,*) 'zone T = "U(r)"'
c  Initialize interaction matrix to 0.d0
      DO  I= 1,3
          H(I,I)=0.0D0
          ENDDO
ccccc Prepare interation matrix  H 
      H(1,1)= -(M1*RETSig+ M3/(RDIST3)+M5/(RDIST3*RDIST2))/(3.d0*RDIST3)
      H(1,2)= -(DSQRT(2.D0))*H(1,1)
      H(2,1)= H(1,2)
      H(1,3)= M1*RETPi/(DSQRT(6.D0)*RDIST3)
      H(3,1)= H(1,3)
      H(2,2)= 2*H(1,1) + DELTAE
      H(2,3)= H(1,3)/DSQRT(2.d0)
      H(3,2)= H(2,3)
      H(3,3)= DELTAE
cccccc Call subroutine to prepare and invert interaction matrix  H
      CALL ZHEEVJ3(H,Q,W)
      L=1
ccc Nor - identify the lowest eigenvalue of  H  and label it  L
      DO J=2,3
          IF (W(J) .LT. W(L)) THEN
              L=J
              ENDIF
          ENDDO  
      ULR= -W(L)
c     WRITE(25,600) RDIST ,ULR 
c 600 FORMAT(2D16.7)
      Modulus = (Z)**2 
      RETURN
      CONTAINS
c=======================================================================
      SUBROUTINE ZHEEVJ3(H,Q,W)
c=======================================================================
c** Subroutine to setup and invert the matrix  H  and return 
c   eigenvalues W and eigenvector matric  Q
      INTEGER   N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8    H(3,3),Q(3,3), W(3)
      REAL*8    SD,SO,S,T,C,G,B,Z,THRESH
c Initialize Q to the identitity matrix
c --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO
c Initialize W to diag(A)
      DO  X = 1, N
          W(X) = H(X, X)
          ENDDO
c Calculate SQR(tr(A))
      SD= 0.0D0
      DO  X = 1, N
          SD= SD + ABS(W(X))
          ENDDO
      SD = SD**2
c Main iteration loop
      DO  I = 1, 30
c Test for convergence
          SO = 0.0D0
          DO  X = 1, N
              DO  Y = X+1, N
                  SO = SO + ABS(H(X, Y))
                  ENDDO
              ENDDO
          IF(SO.EQ.0.0D0) RETURN
          IF (I .LT. 4) THEN
              THRESH = 0.2D0 * SO / N**2
            ELSE
              THRESH = 0.0D0
            END IF
c Do sweep
          DO  X= 1, N
              DO  Y= X+1, N
                  G= 100.0D0*(ABS(H(X, Y)))
                  IF((I.GT.4).AND.((ABS(W(X))+G).EQ.ABS(W(X)))
     $                         .AND.((ABS(W(Y))+G).EQ.ABS(W(Y)))) THEN
                      H(X, Y)= 0.0D0
                    ELSEIF(ABS(H(X, Y)).GT.THRESH) THEN
c Calculate Jacobi transformation
                      B= W(Y) - W(X)
                      IF((ABS(B)+G).EQ.ABS(B)) THEN
                          T= H(X, Y) / B
                        ELSE
                          IF(B .LE. 0.0D0) THEN
                              T= -2.0D0 * H(X, Y)
     $                       /(SQRT(B**2 + 4.0D0*H(X, Y)**2) - B)
                            ELSE IF (B .EQ. 0.0D0) THEN
                              T= H(X, Y) * (1.0D0 / ABS(H(X, Y)))
                            ELSE
                              T= 2.0D0 * H(X, Y)
     $                       /(SQRT(B**2 + 4.0D0*H(X, Y)**2) + B)
                            ENDIF
                        ENDIF
                      C= 1.0D0 / SQRT( 1.0D0 + T**2 )
                      S= T * C
                      Z= T * (H(X, Y))
c Apply Jacobi transformation
                      H(X, Y) = 0.0D0
                      W(X)    = W(X) - Z
                      W(Y)    = W(Y) + Z
                      DO  R = 1, X-1
                          T       = H(R, X)
                          H(R, X) = C * T - (S) * H(R, Y)
                          H(R, Y) = S * T + C * H(R, Y)
                          ENDDO
                      DO  R = X+1, Y-1
                          T       = H(X, R)
                          H(X, R) = C * T - S * (H(R, Y))
                          H(R, Y) = S * (T) + C * H(R, Y)
                          ENDDO
                      DO  R = Y+1, N
                          T       = H(X, R)
                          H(X, R) = C * T - S * H(Y, R)
                          H(Y, R) = (S) * T + C * H(Y, R)
                          ENDDO
c eigenvectors
c This loop can be omitted if only the eigenvalues are desired ---
                      DO  R = 1, N
                          T       = Q(R, X)
                          Q(R, X) = C * T - (S) * Q(R, Y)
                          Q(R, Y) = S * T + C * Q(R, Y)
                          ENDDO
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
cc    PRINT *, "ZHEEVJ3: No convergence.   I=", I
      END SUBROUTINE ZHEEVJ3
      END SUBROUTINE AF3X3LEV
c***********************************************************************
c***********************************************************************
