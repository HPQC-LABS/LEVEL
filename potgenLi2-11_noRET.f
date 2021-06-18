c***********************************************************************
c     SUBROUTINE POTGENLI2(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)
c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c***********************************************************************
c                       Version of  11 Septemper 2011
c** Subroutine to generate the effective adiabatic potential energy 
c  function for any one of four electronic states of Li2 and any one 
c  of the three Li2 isotopologus using the MLR potential functions:
c  of Le Roy et al.[J.Chem.Phys. 131, 204309 (2009)]  
c      for the   X(^1\Sigma_g^+) state  [case  ISTATE= 1]
c      for the   A(^1\Sigma_u^+) state  [case  ISTATE= 3]
c  and those of Dattani & Le Roy [J.Mol.Spectrosc. 268, 199 (2011)]
c      for the   a(^3\Sigma_u^+) state    [case  ISTATE= 2]
c      for the   c(1 ^3\Sigma_g^+) state  [case  ISTATE= 4]
c-----------------------------------------------------------------------
c** Note that this is a special 'standalone' version of subroutine 
c  POTGEN from the general-purpose bound-state Franck-Condon code LEVEL.  
c  Replacing the calling statement in line 2 (above) with the (currently)
c  commented-out lines 4 & 5 allows this routine to be compiled and used
c  in the normal way with R.J. Le Roy's program LEVEL, once one removes
c  the 'cc' in columns #1 & 2 of Line #82 and inserts a chosen value for 
c  parameter ISTATE
c==============
c*** On INPUT:
c==============
c*  integer  ISTATE  specifies the choice of PEF/electronic state (above)
c*  integers IMN1 and IMN2  are the mass numbers for the isotopologue
c*  integer  NPP  is the dimension of the  REAL*8 radial distance and 
c                 potential energy function arrays
c*  VLIM  is the (externally specified) asymptote energy (in cm-1) for
c         the {7.7}Li_2 isotopologue potential of the chosen state
c*  XO(i)  is the array of  NPP  radial distances (in Angst.) at which
c          the potential enrgy function is to be calculated
c==============
c*** On OUTPUT:
c==============
c*  VV(i)  the array of potential function values (in cm-1) 
c*  RM2(i) is the array of BOB-corrected radial centrifugal strength
c          factors  [1 + g_{na}^{\alpha}(r)]/XO(i)**2  defining the
c          centrifugal potential energy function for this isotopologue
c** NOTE  that the effective adiabatic potentials for  ISTATE= 3 & 4
c  include both retardation effects and the non-adiabatic  2*B(r)  term
c  of Eq.(31) of JCP 131, 204309 (2009), so the cenrifugal potential to 
c  be used with these potentials is:   [hbar*2/(2*mu)] [J(J+1)]*RM2(i)
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
     4 ZQ,ZME,Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,fad1e,fad2e,FSW,
     5 ULR,ULRe,rhoAB,REQP,DM(9),DMP(9),DMPP(9),CMM(9),T0,C6adj,C9adj,
     6 RM3,RET,RETsig,RETpi,RETp,RETm,BFCT2,PPOW,PVSR,
     7 U1(0:NBOB),U2(0:NBOB),T1(0:NBOB),T2(0:NBOB),PARM(50),
     8 XO(NPP),VV(NPP),RM2(NPP), bTT(-1:2),cDS(-2:0),bDS(-2:0)
      SAVE IBOB,IPOTL,IORD,IORDD,PPAR,QPAR,PAD,QAD,PNA,NSR,
     1 NLR,MMLR,NVARB,NCMM
      SAVE DSCM,REQ,Rref,PARM,U1,U2,T1,T2,CSAV,BINF,ALFA,Rsw,ZME,
     2 Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,Rinn,Rout,ULR,ULRe,CMM
c** Damping function parameters for printout .....
      DATA bTT/2.44d0,2.78d0,3.126d0,3.471d0/
      DATA bDS/3.3d0,3.69d0,3.95d0/
      DATA cDS/0.423d0,0.40d0,0.39d0/
      SAVE bTT, bDS, cDS
c** Electron mass, as per 2006 physical constants
      DATA ZME/5.4857990943d-4/
ccc Define variables for the four Li2 states of interest
      INTEGER ISTATE,PPARA(4),QPAR4(4),NLRA(4),MMLRA(4,4),NU1A(4)
      REAL*8 DSCMA(4),REQA(4),RrefA(4),CMMA(4,4),PARMA(20,4),U1A(6,4)
c** when using this with the standard LEVEL program, select ISTATE case here
c*** ISTATE =1 for X(^1\Sigma_g^+),
c           =2 for a(^3\Sigma_u^+),
c           =3 for A(^1Pi_u),  
c           =4 for c(1 ^3\Sigma_g^+)
c-----------------------------------------------------------------------
cc    ISTATE= 4     !! a choice of ISTATE required when using with LEVEL
      ISTATE= 4     !! a choice of ISTATE required when using with LEVEL
c-----------------------------------------------------------------------
ccc Specify parameter values for one of the four Li2 states of interest
      DATA PPARA/5,5,6,6/, MMLRA/6,8,10,0, 6,8,10,0, 3,0,6,8,3,-1,6,8/,
     1  NLRA/16,3,16,9/, NU1A/2,0,5,3/ 
      DATA DSCMA/8516.709d0,333.758d0,9353.167d0,7093.44d0/,
     1  REQA/2.672993d0,4.17005d0,3.107856d0,3.06514d0/,
     2  RrefA/3.85d0,8.0d0,4.4d0,3.6d0/, 
c** specify the 4 C_m long-range coefficients for the 4 states
     4  CMMA/6.71527d6,1.12588d8,2.78604d9,0.d0, 6.7185d6,1.12629d8,
     5    2.78683d9,0.d0, 3.57829d5,0.335338d0,1.000045d7,3.7020d8,
     6    3.57557d5,0.335338d0,1.00054d7,3.69953d8/,
c** now specify 20 exponent parameters \beta_i for the 4 states
     7    PARMA/-2.89828701D0,-1.309265D0,-2.018507D0,-1.38162D0,
     8   -1.21933D0,3.463D-1,1.061D-1,-1.886D-1,1.6710D0,1.0943D1,
     9   3.944D0,-2.723D1,-1.149D1,5.67D1,3.89D1,-3.74D1,-3.30D1,3*0.0,
     a   -5.1608D-1,-9.782D-2,1.137D-1,-2.48D-2, 16*0.d0,
     b   -1.75783235D0,-1.034889D0,-1.811924D0,-1.62341D0,-1.44347D0,
     c   1.23978D0,2.9415D0,2.4661D0,1.9912D0,1.6947D1,2.0506D1,
     d  -2.826D1,-5.565D1,9.68D0,4.61D1,2.D0,-1.43D1, 3*0.d0,
     e   -1.6373863D0,2.9197D-1,-5.5544D-1,-2.794D-1,-1.5993D0,
     f   -6.73D-1,-1.23D0,-1.29D0,5.D-1,2.6D+0, 10*0.d0/
c** Specify up to 10 'adiabatic' BOB parameters for up to 4 states
      DATA U1A/0.194d0,-0.01d0,0.39d0, 3*0.d0, 0.059d0, 5*0.d0,
     1    1.066d0,2.98d0,-0.32d0,2.3d0,-7.5d0,3.3d0,
     2    1.367d0,2.7d0,-1.3d0,-1.8d0, 2*0.d0/
      LNPT= 1
      IAN1= 3
      IAN2= 3
      IPOTL= 4
      IBOB= 1
      IVSR= -2
      IDSTT= 1
      QPAR= 3
      NCMM= 3
      IF(ISTATE.GE.3) NCMM= 4
      rhoAB= -1.d0
      IF(ISTATE.EQ.2) rhoAB= 0.54d0
      DSCM= DSCMA(ISTATE)
      REQ= REQA(ISTATE)
      Rref= RrefA(ISTATE)
c
      IF(LNPT.GT.0) THEN
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
          NLR= NLRA(ISTATE)
          DO  I= 1,NCMM
              MMLR(I)= MMLRA(I,ISTATE)
              CMM(I)= CMMA(I,ISTATE)
              ENDDO
          NCN= MMLR(1)
          CNN= CMM(1)
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
          NT1= -1
          IF(ISTATE.EQ.3) NT1= 1
c*** Set mass numbers for reference isotopologue
          NT2= NT1
          MN1R= 7
          MN2R= 7
          NU1= NU1A(ISTATE)
          NU2= NU1
          PAD= 6
          IF(ISTATE.GE.3) PAD=3
          QAD= PAD
c-----------------------------------------------------------------------
cc            READ(5,*) MN1R, MN2R, PAD, QAD, NU1, NU2, PNA, NT1, NT2
c-----------------------------------------------------------------------
          NCMAX= MAX0(NU1,NU2,NT1,NT2)
c** If appropriate, read parameters & prepare to add mass-dep. BOB corrn
          CALL MASSES(IAN1,IMN1,NAME1,GEL,GNS,MASS1,ABUND)
          CALL MASSES(IAN1,MN1R,NAME1,GEL,GNS,RMASS1,ABUND)
          CALL MASSES(IAN2,IMN2,NAME2,GEL,GNS,MASS2,ABUND)
          CALL MASSES(IAN2,MN2R,NAME2,GEL,GNS,RMASS2,ABUND)
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
          U2INF= 0.d0
          T1INF= 0.d0
          T2INF= 0.d0
          Aad1= 0.d0
          Aad2= 0.d0
          Ana1= 0.d0
          Ana1= 0.d0
          Rad1= 0.d0
          Rad2= 0.d0
          Rna1= 0.d0
          Rna2= 0.d0
          IF(ISTATE.EQ.3) THEN
              T1inf= 0.d0
              T1(0)= 0.d0
              T1(1)= 9.2d-05
              T2inf= 0.d0
              T2(0)= 0.d0
              T2(1)= 9.2d-05
              ENDIF
          DO  I= 0, NU1
              U1(I)= U1A(I+1,ISTATE)
              U2(I)= U1(I)
              ENDDO
          U1INF= 0.d0
          IF(ISTATE.GE.3) U1INF= 1.05574d0
          U2INF= U1INF
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
c
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
c
          IF(NT1.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-1 ...
c-----------------------------------------------------------------------
cc                    READ(5,*) T1INF,(T1(I), I=0,NT1)
c-----------------------------------------------------------------------
              WRITE(6,634) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1 1,T1INF,PNA,PNA,PNA,PNA,PNA,PNA,NT1,PNA,NT1+1,(T1(I),I= 0,NT1)
              FG1= RMASS1/MASS1
              ENDIF
c
          IF(NT2.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-2 ...
c-----------------------------------------------------------------------
cc                    READ(5,*) T2INF,(T2(I), I=0,NT2)
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
ccc           ENDIF
          ENDIF
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
                  RET= 0.d0                            !! testing no RET
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
                  RET= 0.d0                               !! testing no RET
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
              BFCT2= 2.d0*16.85762920d0*(MASS1+MASS2)/(MASS1*MASS2)
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
     1 '  Rref=',F10.7)
  615 FORMAT(6x,'radial variables  y_p & y_q  defined w.r.t.',
     1 '  Rref= Re=' F10.7)
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
      RET= 0.d0                                        !! testing no RET
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
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
 
c=======================================================================
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,GNS,MASS,ABUND)
c=======================================================================
c** For isotope with (input) atomic number IAN and mass number IMN,
c  return (output):  (i) as the right-adjusted 2-character variable NAME
c  the alphabetic symbol for that element,  (ii) the ground state
c  electronic degeneracy GELGS, (iii) the nuclear spin degeneracy GNS,
c  (iv) the atomic mass MASS [amu], and  (v) the natural isotopic
c  abundance ABUND [in percent].   GELGS values based on atomic states
c  in Moore's "Atomic Energy Level" tables, the isotope masses are taken
c  from the 2003 mass table [Audi, Wapstra & Thibault, Nucl.Phys. A729,
c  337-676 (2003)] and other quantities from Tables 6.2 and 6.3 of 
c  "Quantities, Units and Symbols in Physical Chemistry", by Mills et 
c  al. (Blackwell, 2'nd Edition, Oxford, 1993).
c** If the input value of IMN does not equal one of the tabulated values
c  for atomic species IAN, return the abundance-averaged standard atomic
c  weight of that atom and set GNS=-1 and ABUND=-1.
c                          COPYRIGHT 2005-2010
c** By R.J. Le Roy (with assistance from G.T. Kraemer & J.Y. Seto).
c           Last modified  25 July 2010 (added proton, d, t)
c***********************************************************************
      REAL*8 zm(123,0:10),mass,ab(123,10),abund
      INTEGER i,ian,imn,gel(123),nmn(123),mn(123,10),ns2(123,10),
     1        gelgs, gns
      CHARACTER*2 NAME,AT(123)
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,6)/' H',2,6,1,2,3,4,5,6/
      DATA  (zm(1,i),i=0,6)/1.00794d0, 1.00782503207d0, 2.0141017778d0,
     1                 3.0160492777d0,1.00727646677d0,2.013553212724d0,
     2                 3.0155007134d0/
      DATA  (ns2(1,i),i=1,3)/1,2,1/
      DATA  (ab(1,i),i=1,3)/99.985d0,0.015d0,0.d0/
c
      DATA  at(2),gel(2),nmn(2),(mn(2,i),i=1,2)/'He',1,2,3,4/
      DATA  (zm(2,i),i=0,2)/4.002602d0, 3.0160293191d0, 4.00260325415d0/
      DATA  (ns2(2,i),i=1,2)/1,0/
      DATA  (ab(2,i),i=1,2)/0.000137d0,99.999863d0/
c
      DATA  at(3),gel(3),nmn(3),(mn(3,i),i=1,2)/'Li',2,2,6,7/
      DATA  (zm(3,i),i=0,2)/6.941d0, 6.015122795d0, 7.01600455d0/
      DATA  (ns2(3,i),i=1,2)/2,3/
      DATA  (ab(3,i),i=1,2)/7.5d0,92.5d0/
c
      DATA  at(4),gel(4),nmn(4),(mn(4,i),i=1,1)/'Be',1,1,9/
      DATA  (zm(4,i),i=0,1)/9.012182d0, 9.0121822d0/
      DATA  (ns2(4,i),i=1,1)/3/
      DATA  (ab(4,i),i=1,1)/100.d0/
c
      DATA at(5),gel(5),nmn(5),(mn(5,i),i=1,2)/' B',2,2,10,11/
      DATA (zm(5,i),i=0,2)/10.811d0, 10.0129370d0, 11.0093054d0/
      DATA  (ns2(5,i),i=1,2)/6,3/
      DATA  (ab(5,i),i=1,2)/19.9d0,80.1d0/
c
      DATA at(6),gel(6),nmn(6),(mn(6,i),i=1,3)/' C',1,3,12,13,14/
      DATA (zm(6,i),i=0,3)/12.011d0, 12.d0, 13.0033548378d0, 
     1                      14.003241989d0/
      DATA  (ns2(6,i),i=1,3)/0,1,0/
      DATA  (ab(6,i),i=1,3)/98.90d0,1.10d0, 0.d0/
c
      DATA at(7),gel(7),nmn(7),(mn(7,i),i=1,2)/' N',4,2,14,15/
      DATA (zm(7,i),i=0,2)/14.00674d0, 14.0030740048d0, 15.0001088982d0/
      DATA (ns2(7,i),i=1,2)/2,1/
      DATA (ab(7,i),i=1,2)/99.634d0,0.366d0/
c
      DATA at(8),gel(8),nmn(8),(mn(8,i),i=1,3)/' O',5,3,16,17,18/
      DATA (zm(8,i),i=0,3)/15.9994d0, 15.99491461956d0, 16.99913170d0,
     1                      17.9991610d0/
      DATA (ns2(8,i),i=1,3)/0,5,0/
      DATA (ab(8,i),i=1,3)/99.762d0, 0.038d0, 0.200d0/
c
      DATA at(9),gel(9),nmn(9),(mn(9,i),i=1,1)/' F',4,1,19/
      DATA (zm(9,i),i=0,1)/18.9984032d0, 18.99840322d0/
      DATA (ns2(9,i),i=1,1)/1/
      DATA (ab(9,i),i=1,1)/100.d0/
c
      DATA at(10),gel(10),nmn(10),(mn(10,i),i=1,3)/'Ne',1,3,20,21,22/
      DATA (zm(10,i),i=0,3)/20.1797d0, 19.9924401754d0, 20.99384668d0,
     1                       21.991385114d0/
      DATA (ns2(10,i),i=1,3)/0,3,0/
      DATA (ab(10,i),i=1,3)/90.48d0, 0.27d0, 9.25d0/
c
      DATA at(11),gel(11),nmn(11),(mn(11,i),i=1,1)/'Na',2,1,23/
      DATA (zm(11,i),i=0,1)/22.989768d0, 22.9897692809d0/
      DATA (ns2(11,i),i=1,1)/3/
      DATA (ab(11,i),i=1,1)/100.d0/
c
      DATA at(12),gel(12),nmn(12),(mn(12,i),i=1,3)/'Mg',1,3,24,25,26/
      DATA (zm(12,i),i=0,3)/24.3050d0, 23.985041700d0, 24.98583692d0,
     1                       25.982592929d0/
      DATA (ns2(12,i),i=1,3)/0,5,0/
      DATA (ab(12,i),i=1,3)/78.99d0, 10.00d0, 11.01d0/
c
      DATA at(13),gel(13),nmn(13),(mn(13,i),i=1,1)/'Al',2,1,27/
      DATA (zm(13,i),i=0,1)/26.981539d0, 26.98153863d0/
      DATA (ns2(13,i),i=1,1)/5/
      DATA (ab(13,i),i=1,1)/100.d0/
c
      DATA at(14),gel(14),nmn(14),(mn(14,i),i=1,3)/'Si',1,3,28,29,30/
      DATA (zm(14,i),i=0,3)/28.0855d0, 27.9769265325d0, 28.976494700d0,
     1                       29.97377017d0/
      DATA (ns2(14,i),i=1,3)/0,1,0/
      DATA (ab(14,i),i=1,3)/92.23d0, 4.67d0, 3.10d0/
 
      DATA at(15),gel(15),nmn(15),(mn(15,i),i=1,1)/' P',4,1,31/
      DATA (zm(15,i),i=0,1)/30.973762d0, 30.97376163d0/
      DATA (ns2(15,i),i=1,1)/1/
      DATA (ab(15,i),i=1,1)/100.d0/
c
      DATA at(16),gel(16),nmn(16),(mn(16,i),i=1,4)/' S',5,4,32,33,34,36/
      DATA (zm(16,i),i=0,4)/32.066d0, 31.97207100d0, 32.97145876d0,
     1                       33.96786690d0, 35.96708076d0/
      DATA (ns2(16,i),i=1,4)/0,3,0,0/
      DATA (ab(16,i),i=1,4)/95.02d0, 0.75d0, 4.21d0, 0.02d0/
c
      DATA at(17),gel(17),nmn(17),(mn(17,i),i=1,2)/'Cl',4,2,35,37/
      DATA (zm(17,i),i=0,2)/35.4527d0, 34.96885268d0, 36.96590259d0/
      DATA (ns2(17,i),i=1,2)/3,3/
      DATA (ab(17,i),i=1,2)/75.77d0, 24.23d0/
c
      DATA at(18),gel(18),nmn(18),(mn(18,i),i=1,3)/'Ar',1,3,36,38,40/
      DATA (zm(18,i),i=0,3)/39.948d0, 35.967545106d0, 37.9627324d0,
     1                       39.9623831225d0/
      DATA (ns2(18,i),i=1,3)/0,0,0/
      DATA (ab(18,i),i=1,3)/0.337d0, 0.063d0, 99.600d0/
c
      DATA at(19),gel(19),nmn(19),(mn(19,i),i=1,3)/' K',2,3,39,40,41/
      DATA (zm(19,i),i=0,3)/39.0983d0, 38.96370668d0, 39.96399848d0,
     1                       40.96182576d0/
      DATA (ns2(19,i),i=1,3)/3,8,3/
      DATA (ab(19,i),i=1,3)/93.2581d0, 0.0117d0, 6.7302d0/
 
      DATA at(20),gel(20),nmn(20),(mn(20,i),i=1,6)/'Ca',1,6,40,42,43,44,
     1                                              46,48/
      DATA (zm(20,i),i=0,6)/40.078d0, 39.96259098d0, 41.95861801d0,
     1         42.9587666d0, 43.9554818d0, 45.9536926d0, 47.952534d0/
      DATA (ns2(20,i),i=1,6)/0,0,7,0,0,0/
      DATA (ab(20,i),i=1,6)/96.941d0, 0.647d0, 0.135d0, 2.086d0,
     1                      0.004d0, 0.187d0/
c
      DATA at(21),gel(21),nmn(21),(mn(21,i),i=1,1)/'Sc',4,1,45/
      DATA (zm(21,i),i=0,1)/44.955910d0, 44.9559119d0/
      DATA (ns2(21,i),i=1,1)/7/
      DATA (ab(21,i),i=1,1)/100.d0/
c
      DATA at(22),gel(22),nmn(22),(mn(22,i),i=1,5)/'Ti',5,5,46,47,48,49,
     1                                              50/
      DATA (zm(22,i),i=0,5)/47.88d0, 45.9526316d0, 46.9517631d0,
     1         47.9479463d0, 48.9478700d0, 49.9447912d0/
      DATA (ns2(22,i),i=1,5)/0,5,0,7,0/
      DATA (ab(22,i),i=1,5)/8.0d0, 7.3d0, 73.8d0, 5.5d0, 5.4d0/
c
      DATA at(23),gel(23),nmn(23),(mn(23,i),i=1,2)/' V',4,2,50,51/
      DATA (zm(23,i),i=0,2)/50.9415d0, 49.9471585d0, 50.9439595d0/
      DATA (ns2(23,i),i=1,2)/12,7/
      DATA (ab(23,i),i=1,2)/0.250d0, 99.750d0/
c
      DATA at(24),gel(24),nmn(24),(mn(24,i),i=1,4)/'Cr',7,4,50,52,53,54/
      DATA (zm(24,i),i=0,4)/51.9961d0, 49.9460442d0, 51.9405075d0,
     1                       52.9406494d0, 53.9388804d0/
      DATA (ns2(24,i),i=1,4)/0,0,3,0/
      DATA (ab(24,i),i=1,4)/4.345d0, 83.789d0, 9.501d0, 2.365d0/
c
      DATA at(25),gel(25),nmn(25),(mn(25,i),i=1,1)/'Mn',6,1,55/
      DATA (zm(25,i),i=0,1)/54.93805d0, 54.9380451d0/
      DATA (ns2(25,i),i=1,1)/5/
      DATA (ab(25,i),i=1,1)/100.d0/
c
      DATA at(26),gel(26),nmn(26),(mn(26,i),i=1,4)/'Fe',9,4,54,56,57,58/
      DATA (zm(26,i),i=0,4)/55.847d0, 53.9396105d0, 55.9349375d0,
     1                       56.9353940d0, 57.9332756d0/
      DATA (ns2(26,i),i=1,4)/0,0,1,0/
      DATA (ab(26,i),i=1,4)/5.8d0, 91.72d0, 2.2d0, 0.28d0/
c
      DATA at(27),gel(27),nmn(27),(mn(27,i),i=1,1)/'Co',10,1,59/
      DATA (zm(27,i),i=0,1)/58.93320d0, 58.9331950d0/
      DATA (ns2(27,i),i=1,1)/7/
      DATA (ab(27,i),i=1,1)/100.d0/
c
      DATA at(28),gel(28),nmn(28),(mn(28,i),i=1,5)/'Ni',9,5,58,60,61,62,
     1                                              64/
      DATA (zm(28,i),i=0,5)/58.69d0, 57.9353429d0, 59.9307864d0,
     1         60.9310560d0, 61.9283451d0, 63.9279660d0/
      DATA (ns2(28,i),i=1,5)/0,0,3,0,0/
      DATA (ab(28,i),i=1,5)/68.077d0,26.223d0,1.140d0,3.634d0,0.926d0/
c
      DATA at(29),gel(29),nmn(29),(mn(29,i),i=1,2)/'Cu',2,2,63,65/
      DATA (zm(29,i),i=0,2)/63.546d0, 62.9295975d0,64.9277895d0/
      DATA (ns2(29,i),i=1,2)/3,3/
      DATA (ab(29,i),i=1,2)/69.17d0, 30.83d0/
c
      DATA at(30),gel(30),nmn(30),(mn(30,i),i=1,5)/'Zn',1,5,64,66,67,68,
     1                                              70/
      DATA (zm(30,i),i=0,5)/65.40d0, 63.9291422d0, 65.9260334d0,
     1         66.9271273d0, 67.9248442d0, 69.9253193d0/
      DATA (ns2(30,i),i=1,5)/0,0,5,0,0/
      DATA (ab(30,i),i=1,5)/48.6d0, 27.9d0, 4.1d0, 18.8d0, 0.6d0/
c
      DATA at(31),gel(31),nmn(31),(mn(31,i),i=1,2)/'Ga',2,2,69,71/
Coxon   DATA (zm(31,i),i=0,2)/69.723d0, 68.925581d0, 70.9247073d0/
      DATA (zm(31,i),i=0,2)/69.723d0, 68.9255736d0, 70.9247013d0/
      DATA (ns2(31,i),i=1,2)/3,3/
      DATA (ab(31,i),i=1,2)/60.108d0, 39.892d0/
c
      DATA at(32),gel(32),nmn(32),(mn(32,i),i=1,5)/'Ge',1,5,70,72,73,74,
     1                                              76/
      DATA (zm(32,i),i=0,5)/72.61d0, 69.9242474d0, 71.9220758d0,
     1         72.9234589d0, 73.9211778d0, 75.9214026d0/
      DATA (ns2(32,i),i=1,5)/0,0,9,0,0/
      DATA (ab(32,i),i=1,5)/21.23d0, 27.66d0, 7.73d0, 35.94d0, 7.44d0/
c
      DATA at(33),gel(33),nmn(33),(mn(33,i),i=1,1)/'As',4,1,75/
      DATA (zm(33,i),i=0,1)/74.92159d0, 74.9215965d0/
      DATA (ns2(33,i),i=1,1)/3/
      DATA (ab(33,i),i=1,1)/100.d0/
c
      DATA at(34),gel(34),nmn(34),(mn(34,i),i=1,6)/'Se',5,6,74,76,77,78,
     1                                              80,82/
      DATA (zm(34,i),i=0,6)/78.96d0, 73.9224764d0, 75.9192136d0,
     1         76.9199140d0, 77.9173091d0, 79.9165213d0, 81.9166994d0/
      DATA (ns2(34,i),i=1,6)/0,0,1,0,0,0/
      DATA (ab(34,i),i=1,6)/0.89d0, 9.36d0, 7.63d0, 23.78d0, 49.61d0,
     1                      8.73d0/
c
      DATA at(35),gel(35),nmn(35),(mn(35,i),i=1,2)/'Br',4,2,79,81/
      DATA (zm(35,i),i=0,2)/79.904d0, 78.9183371d0, 80.9162906d0/
      DATA (ns2(35,i),i=1,2)/3,3/
      DATA (ab(35,i),i=1,2)/50.69d0, 49.31d0/
c
      DATA at(36),gel(36),nmn(36),(mn(36,i),i=1,6)/'Kr',1,6,78,80,82,83,
     1                                              84,86/
      DATA (zm(36,i),i=0,6)/83.80d0, 77.9203648d0, 79.9163790d0,
     1          81.9134836d0, 82.914136d0, 83.911507d0, 85.91061073d0/
      DATA (ns2(36,i),i=1,6)/0,0,0,9,0,0/
      DATA (ab(36,i),i=1,6)/0.35d0, 2.25d0, 11.6d0, 11.5d0, 57.0d0,
     1                      17.3d0/
c
      DATA at(37),gel(37),nmn(37),(mn(37,i),i=1,2)/'Rb',2,2,85,87/
      DATA (zm(37,i),i=0,2)/85.4678d0, 84.911789738d0, 86.909180527d0/
      DATA (ns2(37,i),i=1,2)/5,3/
      DATA (ab(37,i),i=1,2)/72.165d0, 27.835d0/
c
      DATA at(38),gel(38),nmn(38),(mn(38,i),i=1,4)/'Sr',1,4,84,86,87,88/
      DATA (zm(38,i),i=0,4)/87.62d0, 83.913425d0, 85.9092602d0,
     1                      86.9088771d0, 87.9056121d0/
      DATA (ns2(38,i),i=1,4)/0,0,9,0/
      DATA (ab(38,i),i=1,4)/0.56d0, 9.86d0, 7.00d0, 82.58d0/
c
      DATA at(39),gel(39),nmn(39),(mn(39,i),i=1,1)/' Y',4,1,89/
      DATA (zm(39,i),i=0,1)/88.90585d0, 88.9058483d0/
      DATA (ns2(39,i),i=1,1)/1/
      DATA (ab(39,i),i=1,1)/100.d0/
c
      DATA at(40),gel(40),nmn(40),(mn(40,i),i=1,5)/'Zr',5,5,90,91,92,94,
     1                                              96/
      DATA (zm(40,i),i=0,5)/91.224d0, 89.9047044d0, 90.9056458d0,
     1                      91.9050408d0, 93.9063152d0, 95.9082734d0/
      DATA (ns2(40,i),i=1,5)/0,5,0,0,0/
      DATA (ab(40,i),i=1,5)/51.45d0, 11.22d0, 17.15d0, 17.38d0, 2.80d0/
c
      DATA at(41),gel(41),nmn(41),(mn(41,i),i=1,1)/'Nb',2,1,93/
      DATA (zm(41,i),i=0,1)/92.90638d0, 92.9063781d0/
      DATA (ns2(41,i),i=1,1)/9/
      DATA (ab(41,i),i=1,1)/100.d0/
c
      DATA at(42),gel(42),nmn(42),(mn(42,i),i=1,7)/'Mo',7,7,92,94,95,96,
     1                                              97,98,100/
      DATA (zm(42,i),i=0,7)/95.94d0, 91.906811d0, 93.9050883d0,
     1        94.9058421d0, 95.9046795d0, 96.9060215d0, 97.9054082d0,
     2        99.907477d0/
      DATA (ns2(42,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(42,i),i=1,7)/14.84d0, 9.25d0, 15.92d0, 16.68d0, 9.55d0,
     1                      24.13d0, 9.63d0/
c
      DATA at(43),gel(43),nmn(43),(mn(43,i),i=1,1)/'Tc',6,1,98/
      DATA (zm(43,i),i=0,1)/97.907215d0, 97.907216d0/
      DATA (ns2(43,i),i=1,1)/12/
      DATA (ab(43,i),i=1,1)/100.d0/
c
      DATA at(44),gel(44),nmn(44),(mn(44,i),i=1,7)/'Ru',11,7,96,98,99,
     1                                              100,101,102,104/
      DATA (zm(44,i),i=0,7)/101.07d0, 95.907598d0, 97.905287d0,
     1     98.9059393d0, 99.9042195d0, 100.9055821d0, 101.9043493d0,
     2     103.905433d0/
      DATA (ns2(44,i),i=1,7)/0,0,5,0,5,0,0/
      DATA (ab(44,i),i=1,7)/5.52d0, 1.88d0, 12.7d0, 12.6d0, 17.0d0,
     1                      31.6d0, 18.7d0/
c
      DATA at(45),gel(45),nmn(45),(mn(45,i),i=1,1)/'Rh',10,1,103/
      DATA (zm(45,i),i=0,1)/102.90550d0, 102.905504d0/
      DATA (ns2(45,i),i=1,1)/1/
      DATA (ab(45,i),i=1,1)/100.d0/
c
      DATA at(46),gel(46),nmn(46),(mn(46,i),i=1,6)/'Pd',1,6,102,104,105,
     1                                              106,108,110/
      DATA (zm(46,i),i=0,6)/106.42d0, 101.905609d0, 103.904036d0,
     1       104.905085d0, 105.903486d0, 107.903892d0, 109.905153d0/
      DATA (ns2(46,i),i=1,6)/0,0,5,0,0,0/
      DATA (ab(46,i),i=1,6)/1.02d0, 11.14d0, 22.33d0, 27.33d0, 26.46d0,
     1                      11.72d0/
c
      DATA at(47),gel(47),nmn(47),(mn(47,i),i=1,2)/'Ag',2,2,107,109/
      DATA (zm(47,i),i=0,2)/107.8682d0, 106.905097d0, 108.904752d0/
      DATA (ns2(47,i),i=1,2)/1,1/
      DATA (ab(47,i),i=1,2)/51.839d0, 48.161d0/
c
      DATA at(48),gel(48),nmn(48),(mn(48,i),i=1,8)/'Cd',1,8,106,108,110,
     1                                             111,112,113,114,116/ 
      DATA (zm(48,i),i=0,8)/112.411d0, 105.906459d0, 107.904184d0, 
     1       109.9030021d0, 110.9041781d0, 111.9027578d0, 112.9044017d0,
     2       113.9033585d0, 115.904756d0/
      DATA (ns2(48,i),i=1,8)/0,0,0,1,0,1,0,0/
      DATA (ab(48,i),i=1,8)/1.25d0, 0.89d0, 12.49d0, 12.80d0, 24.13d0,
     1                      12.22d0, 28.73d0, 7.49d0/
c
      DATA at(49),gel(49),nmn(49),(mn(49,i),i=1,2)/'In',2,2,113,115/
      DATA (zm(49,i),i=0,2)/114.818d0, 112.904058d0, 114.903878d0/
      DATA  (ns2(49,i),i=1,2)/9,9/
      DATA (ab(49,i),i=1,2)/4.3d0, 95.7d0/
c
      DATA at(50),gel(50),nmn(50),(mn(50,i),i=1,10)/'Sn',1,10,112,114,
     1                                 115,116,117,118,119,120,122,124/
      DATA (zm(50,i),i=0,10)/118.710d0, 111.904818d0, 113.902779d0,
     1     114.903342d0, 115.901741d0, 116.902952d0, 117.901603d0,
     2     118.903308d0, 119.9021947d0, 121.9034390d0, 123.9052739d0/
      DATA (ns2(50,i),i=1,10)/0,0,1,0,1,0,1,0,0,0/
      DATA (ab(50,i),i=1,10)/0.97d0, 0.65d0, 0.34d0, 14.53d0, 7.68d0,
     1                       24.23d0, 8.59d0, 32.59d0, 4.63d0, 5.79d0/
c
      DATA at(51),gel(51),nmn(51),(mn(51,i),i=1,2)/'Sb',4,2,121,123/
      DATA (zm(51,i),i=0,2)/121.757d0, 120.9038157d0, 122.9042140d0/
      DATA (ns2(51,i),i=1,2)/5,7/
      DATA (ab(51,i),i=1,2)/57.36d0, 42.64d0/
c
      DATA at(52),gel(52),nmn(52),(mn(52,i),i=1,8)/'Te',5,8,120,122,123,
     1                                             124,125,126,128,130/
      DATA (zm(52,i),i=0,8)/127.60d0, 119.904020d0, 121.9030439d0,
     1    122.9042700d0, 123.9028179d0, 124.9044307d0, 125.9033117d0,
     2    127.9044631d0, 129.9062244d0/
      DATA (ns2(52,i),i=1,8)/0,0,1,0,1,0,0,0/
      DATA (ab(52,i),i=1,8)/0.096d0, 2.603d0, 0.908d0, 4.816d0,
     1                      7.139d0, 18.95d0, 31.69d0, 33.80d0/
c
      DATA at(53),gel(53),nmn(53),(mn(53,i),i=1,2)/' I',4,2,127,129/
      DATA (zm(53,i),i=0,2)/126.90447d0, 126.904473d0, 128.904988d0/
      DATA (ns2(53,i),i=1,2)/5,7/
      DATA (ab(53,i),i=1,2)/100.d0,0.d0/
c
      DATA at(54),gel(54),nmn(54),(mn(54,i),i=1,9)/'Xe',1,9,124,126,128,
     1                                          129,130,131,132,134,136/
      DATA (zm(54,i),i=0,9)/131.29d0, 123.9058930d0, 125.904274d0,
     1    127.9035313d0, 128.9047794d0, 129.9035080d0, 130.9050824d0,
     2    131.9041535d0, 133.9053945d0, 135.907219d0/
      DATA (ns2(54,i),i=1,9)/0,0,0,1,0,3,0,0,0/
      DATA (ab(54,i),i=1,9)/0.10d0, 0.09d0, 1.91d0, 26.4d0, 4.1d0,
     1                      21.2d0, 26.9d0, 10.4d0, 8.9d0/
c
      DATA at(55),gel(55),nmn(55),(mn(55,i),i=1,1)/'Cs',2,1,133/
      DATA (zm(55,i),i=0,1)/132.90543d0, 132.905451933d0/
      DATA (ns2(55,i),i=1,1)/7/
      DATA (ab(55,i),i=1,1)/100.d0/
c
      DATA at(56),gel(56),nmn(56),(mn(56,i),i=1,7)/'Ba',1,7,130,132,134,
     1                                             135,136,137,138/
      DATA (zm(56,i),i=0,7)/137.327d0, 129.9063208d0, 131.9050613d0,
     1    133.9045084d0, 134.9056886d0, 135.9045759d0, 136.9058274d0,
     2    137.9052472d0/
      DATA (ns2(56,i),i=1,7)/0,0,0,3,0,3,0/
      DATA (ab(56,i),i=1,7)/0.106d0, 0.101d0, 2.417d0, 6.592d0, 
     1                      7.854d0, 11.23d0, 71.70d0/
c
      DATA at(57),gel(57),nmn(57),(mn(57,i),i=1,2)/'La',4,2,138,139/
      DATA (zm(57,i),i=0,2)/138.9055d0, 137.907112d0, 138.9063533d0/
      DATA (ns2(57,i),i=1,2)/10,7/ 
      DATA (ab(57,i),i=1,2)/0.0902d0, 99.9098d0/
c
      DATA at(58),gel(58),nmn(58),(mn(58,i),i=1,4)/'Ce',9,4,136,138,140,
     1                                             142/
      DATA (zm(58,i),i=0,4)/140.115d0, 135.907172d0, 137.905991d0,
     1    139.9054387d0, 141.909244d0/
      DATA (ns2(58,i),i=1,4)/0,0,0,0/
      DATA (ab(58,i),i=1,4)/0.19d0, 0.25d0, 88.48d0, 11.08d0/
c
      DATA at(59),gel(59),nmn(59),(mn(59,i),i=1,1)/'Pr',10,1,141/
      DATA (zm(59,i),i=0,1)/140.90765d0, 140.9076528d0/
      DATA (ns2(59,i),i=1,1)/5/
      DATA (ab(59,i),i=1,1)/100.d0/
c
      DATA at(60),gel(60),nmn(60),(mn(60,i),i=1,7)/'Nd',9,7,142,143,144,
     1                                             145,146,148,150/
      DATA (zm(60,i),i=0,7)/144.24d0, 141.9077233d0, 142.9098143d0,
     1    143.9100873d0, 144.9125736d0, 145.9131169d0, 147.916893d0,
     2    149.920891d0/
      DATA (ns2(60,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(60,i),i=1,7)/27.13d0, 12.18d0, 23.80d0, 8.30d0, 17.19d0,
     1                       5.76d0, 5.64d0/
c
      DATA at(61),gel(61),nmn(61),(mn(61,i),i=1,1)/'Pm',6,1,145/
      DATA (zm(61,i),i=0,1)/144.912743d0, 144.912749d0/
      DATA (ns2(61,i),i=1,1)/5/
      DATA (ab(61,i),i=1,1)/100.d0/
c
      DATA at(62),gel(62),nmn(62),(mn(62,i),i=1,7)/'Sm',1,7,144,147,148,
     1                                             149,150,152,154/
      DATA (zm(62,i),i=0,7)/150.36d0, 143.911999d0, 146.9148979d0,
     1    147.9148227d0, 148.9171847d0, 149.9172755d0, 151.9197324d0,
     2    153.9222093d0/
      DATA (ns2(62,i),i=1,7)/0,7,0,7,0,0,0/
      DATA (ab(62,i),i=1,7)/3.1d0, 15.0d0, 11.3d0, 13.8d0, 7.4d0,
     1                      26.7d0, 22.7d0/
c
      DATA at(63),gel(63),nmn(63),(mn(63,i),i=1,2)/'Eu',8,2,151,153/
      DATA (zm(63,i),i=0,2)/151.965d0, 150.9198502d0, 152.9212303d0/
      DATA (ns2(63,i),i=1,2)/5,5/
      DATA (ab(63,i),i=1,2)/47.8d0, 52.2d0/
c
      DATA at(64),gel(64),nmn(64),(mn(64,i),i=1,7)/'Gd',5,7,152,154,155,
     1                                              156,157,158,160/
      DATA (zm(64,i),i=0,7)/157.25d0, 151.9197910d0, 153.92086560d0,
     1    154.9226220d0, 155.9221227d0, 156.9239601d0, 157.9241039d0,
     2    159.9270541d0/
      DATA (ns2(64,i),i=1,7)/0,0,3,0,3,0,0/
      DATA (ab(64,i),i=1,7)/0.20d0, 2.18d0, 14.80d0, 20.47d0, 15.65d0,
     1                      24.84d0, 21.86d0/
c
      DATA at(65),gel(65),nmn(65),(mn(65,i),i=1,1)/'Tb',16,1,159/
      DATA (zm(65,i),i=0,1)/158.92534d0, 158.9253468d0/
      DATA (ns2(65,i),i=1,1)/3/
      DATA (ab(65,i),i=1,1)/100.d0/
c
      DATA at(66),gel(66),nmn(66),(mn(66,i),i=1,7)/'Dy',17,7,156,158,
     1                                           160,161,162,163,164/
      DATA (zm(66,i),i=0,7)/162.50d0, 155.924283d0, 157.924409d0,
     1    159.9251975d0, 160.9269334d0, 161.9267984d0, 162.9287312d0,
     2    163.9291748d0/
      DATA (ns2(66,i),i=1,7)/0,0,0,5,0,5,0/
      DATA (ab(66,i),i=1,7)/0.06d0, 0.10d0, 2.34d0, 18.9d0, 25.5d0,
     1                      24.9d0, 28.2d0/
c
      DATA at(67),gel(67),nmn(67),(mn(67,i),i=1,1)/'Ho',16,1,165/
      DATA (zm(67,i),i=0,1)/164.93032d0, 164.9303221d0/
      DATA (ns2(67,i),i=1,1)/7/
      DATA (ab(67,i),i=1,1)/100.d0/
     
      DATA at(68),gel(68),nmn(68),(mn(68,i),i=1,6)/'Er',13,6,162,164,
     1                                            166,167,168,170/
      DATA (zm(68,i),i=0,6)/167.26d0, 161.928778d0, 163.929200d0,
     1    165.9302931d0, 166.9320482d0, 167.9323702d0, 169.9354643d0/
      DATA (ns2(68,i),i=1,6)/0,0,0,7,0,0/
      DATA (ab(68,i),i=1,6)/0.14d0, 1.61d0, 33.6d0, 22.95d0, 26.8d0,
     1                      14.9d0/
c
      DATA at(69),gel(69),nmn(69),(mn(69,i),i=1,1)/'Tm',8,1,169/  
      DATA (zm(69,i),i=0,1)/168.93421d0, 168.9342133d0/
      DATA (ns2(69,i),i=1,1)/1/
      DATA (ab(69,i),i=1,1)/100.d0/
c
      DATA at(70),gel(70),nmn(70),(mn(70,i),i=1,7)/'Yb',1,7,168,170,171,
     1                                            172,173,174,176/
      DATA (zm(70,i),i=0,7)/173.04d0, 167.933897d0, 169.9347618d0,
     1    170.936323580d0, 171.9363815d0, 172.9382108d0, 173.9388621d0,
     2    175.9425717d0/
      DATA (ns2(70,i),i=1,7)/0,0,1,0,5,0,0/
      DATA (ab(70,i),i=1,7)/0.13d0, 3.05d0, 14.3d0, 21.9d0, 16.12d0,
     1                      31.8d0, 12.7d0/
c
      DATA at(71),gel(71),nmn(71),(mn(71,i),i=1,2)/'Lu',4,2,175,176/
      DATA (zm(71,i),i=0,2)/174.967d0, 174.9407718d0, 175.9426863d0/
      DATA (ns2(71,i),i=1,2)/7,14/
      DATA (ab(71,i),i=1,2)/97.41d0, 2.59d0/
c
      DATA at(72),gel(72),nmn(72),(mn(72,i),i=1,6)/'Hf',5,6,174,176,177,
     1                                             178,179,180/
      DATA (zm(72,i),i=0,6)/178.49d0, 173.940046d0, 175.9414086d0,
     1    176.9432207d0, 177.9436988d0, 178.9458161d0, 179.9465500d0/
      DATA (ns2(72,i),i=1,6)/0,0,7,0,9,0/
      DATA (ab(72,i),i=1,6)/0.162d0, 5.206d0, 18.606d0, 27.297d0,
     1                      13.629d0, 35.100d0/
c
      DATA at(73),gel(73),nmn(73),(mn(73,i),i=1,2)/'Ta',4,2,180,181/
      DATA (zm(73,i),i=0,2)/180.9479d0, 179.9474648d0, 180.9479958d0/
      DATA (ns2(73,i),i=1,2)/16,7/
      DATA (ab(73,i),i=1,2)/0.012d0, 99.988d0/
c
      DATA at(74),gel(74),nmn(74),(mn(74,i),i=1,5)/' W',1,5,180,182,183,
     1                                             184,186/
      DATA (zm(74,i),i=0,5)/183.84d0, 179.946704d0, 181.9482042d0,
     1    182.9502230d0, 183.9509312d0, 185.9543641d0/
      DATA (ns2(74,i),i=1,5)/0,0,1,0,0/
      DATA (ab(74,i),i=1,5)/0.13d0, 26.3d0, 14.3d0, 30.67d0, 28.6d0/
c
      DATA at(75),gel(75),nmn(75),(mn(75,i),i=1,2)/'Re',6,2,185,187/
      DATA (zm(75,i),i=0,2)/186.207d0, 184.9529550d0, 186.9557531d0/
      DATA (ns2(75,i),i=1,2)/5,5/
      DATA (ab(75,i),i=1,2)/37.40d0, 62.60d0/
c
      DATA at(76),gel(76),nmn(76),(mn(76,i),i=1,7)/'Os',9,7,184,186,187,
     1                                             188,189,190,192/
      DATA (zm(76,i),i=0,7)/190.23d0, 183.9524891d0, 185.9538382d0,
     1    186.9557505d0, 187.9558382d0, 188.9581475d0, 189.9584470d0,
     2    191.9614807d0/
      DATA (ns2(76,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(76,i),i=1,7)/0.02d0, 1.58d0, 1.6d0, 13.3d0, 16.1d0,
     1                      26.4d0, 41.0d0/
c
      DATA at(77),gel(77),nmn(77),(mn(77,i),i=1,2)/'Ir',10,2,191,193/
      DATA (zm(77,i),i=0,2)/192.22d0, 190.9605940d0, 192.9629264d0/
      DATA (ns2(77,i),i=1,2)/3,3/
      DATA (ab(77,i),i=1,2)/37.3d0, 62.7d0/
c
c
      DATA at(78),gel(78),nmn(78),(mn(78,i),i=1,6)/'Pt',7,6,190,192,194,
     1                                            195,196,198/
      DATA (zm(78,i),i=0,6)/195.08d0, 189.959932d0, 191.9610380d0,
     1    193.9626803d0, 194.9647911d0, 195.9649515d0, 197.967893d0/
      DATA (ns2(78,i),i=1,6)/0,0,0,1,0,0/
      DATA (ab(78,i),i=1,6)/0.01d0,0.79d0,32.9d0,33.8d0,25.3d0,7.2d0/
c
      DATA at(79),gel(79),nmn(79),(mn(79,i),i=1,1)/'Au',2,1,197/
      DATA (zm(79,i),i=0,1)/196.96654d0, 196.9665687d0/
      DATA (ns2(79,i),i=1,1)/3/
      DATA (ab(79,i),i=1,1)/100.d0/
c
      DATA at(80),gel(80),nmn(80),(mn(80,i),i=1,7)/'Hg',1,7,196,198,199,
     1                                            200,201,202,204/
      DATA (zm(80,i),i=0,7)/200.59d0, 195.965833d0, 197.9667690d0,
     1    198.9682799d0, 199.9683260d0, 200.9703023d0, 201.9706430d0,
     2    203.9734939d0/
      DATA (ns2(80,i),i=1,7)/0,0,1,0,3,0,0/
      DATA (ab(80,i),i=1,7)/0.15d0, 9.97d0, 16.87d0, 23.10d0, 13.18d0,
     1                      29.86d0, 6.87d0/
c
      DATA at(81),gel(81),nmn(81),(mn(81,i),i=1,2)/'Tl',2,2,203,205/
      DATA (zm(81,i),i=0,2)/204.3833d0, 202.9723442d0, 204.9744275d0/
      DATA (ns2(81,i),i=1,2)/1,1/
      DATA (ab(81,i),i=1,2)/29.524d0, 70.476d0/
c
      DATA at(82),gel(82),nmn(82),(mn(82,i),i=1,4)/'Pb',1,4,204,206,207,
     1                                             208/
      DATA (zm(82,i),i=0,4)/207.2d0, 203.9730436d0, 205.9744653d0,
     1    206.9758969d0, 207.9766521d0/
      DATA (ns2(82,i),i=1,4)/0,0,1,0/
      DATA (ab(82,i),i=1,4)/1.4d0, 24.1d0, 22.1d0, 52.4d0/
c
      DATA at(83),gel(83),nmn(83),(mn(83,i),i=1,1)/'Bi',4,1,209/
      DATA (zm(83,i),i=0,1)/208.98037d0, 208.9803987d0/
      DATA (ns2(83,i),i=1,1)/9/
      DATA (ab(83,i),i=1,1)/100.d0/
c
      DATA at(84),gel(84),nmn(84),(mn(84,i),i=1,1)/'Po',5,1,209/
      DATA (zm(84,i),i=0,1)/208.982404d0, 208.9824304d0/
      DATA (ns2(84,i),i=1,1)/1/
      DATA (ab(84,i),i=1,1)/100.d0/
c
      DATA at(85),gel(85),nmn(85),(mn(85,i),i=1,1)/'At',-1,1,210/
      DATA (zm(85,i),i=0,1)/209.987126d0, 209.987148d0/
      DATA (ns2(85,i),i=1,1)/10/
      DATA (ab(85,i),i=1,1)/100.d0/
c
      DATA at(86),gel(86),nmn(86),(mn(86,i),i=1,1)/'Rn',1,1,222/
      DATA (zm(86,i),i=0,1)/222.017571d0, 222.0175777d0/
      DATA (ns2(86,i),i=1,1)/0/
      DATA (ab(86,i),i=1,1)/100.d0/
c
      DATA at(87),gel(87),nmn(87),(mn(87,i),i=1,1)/'Fr',-1,1,223/
      DATA (zm(87,i),i=0,1)/223.019733d0, 223.0197359d0/
      DATA (ns2(87,i),i=1,1)/3/
      DATA (ab(87,i),i=1,1)/100.d0/
c
      DATA at(88),gel(88),nmn(88),(mn(88,i),i=1,1)/'Ra',1,1,226/
      DATA (zm(88,i),i=0,1)/226.025403d0, 226.0254098d0/
      DATA (ns2(88,i),i=1,1)/0/
      DATA (ab(88,i),i=1,1)/100.d0/
c
      DATA at(89),gel(89),nmn(89),(mn(89,i),i=1,1)/'Ac',4,1,227/
      DATA (zm(89,i),i=0,1)/227.027750d0, 227.0277521d0/
      DATA (ns2(89,i),i=1,1)/3/
      DATA (ab(89,i),i=1,1)/100.d0/
c
      DATA at(90),gel(90),nmn(90),(mn(90,i),i=1,1)/'Th',-1,1,232/
      DATA (zm(90,i),i=0,1)/232.038d0, 232.0380553d0/
      DATA (ns2(90,i),i=1,1)/0/
      DATA (ab(90,i),i=1,1)/100.d0/
c
      DATA at(91),gel(91),nmn(91),(mn(91,i),i=1,1)/'Pa',-1,1,231/
      DATA (zm(91,i),i=0,1)/231.03588d0, 231.0358840d0/
      DATA (ns2(91,i),i=1,1)/3/
      DATA (ab(91,i),i=1,1)/100.d0/
c
      DATA at(92),gel(92),nmn(92),(mn(92,i),i=1,4)/' U',-1,4,233,234,
     1                                             235,238/
      DATA (zm(92,i),i=0,4)/238.0289d0, 233.0396352d0, 234.0409521d0,
     1    235.0439299d0, 238.0507882d0/
      DATA (ns2(92,i),i=1,4)/5,0,7,0/
      DATA (ab(92,i),i=1,4)/0.d0, 0.0055d0, 0.7200d0, 99.2745d0/
c
      DATA at(93),gel(93),nmn(93),(mn(93,i),i=1,1)/'Np',-1,1,237/
      DATA (zm(93,i),i=0,1)/237.0481678d0, 237.0481734d0/
      DATA (ns2(93,i),i=1,1)/5/
      DATA (ab(93,i),i=1,1)/100.d0/
c
      DATA at(94),gel(94),nmn(94),(mn(94,i),i=1,1)/'Pu',-1,1,244/
      DATA (zm(94,i),i=0,1)/244.064199d0, 244.064204d0/
      DATA (ns2(94,i),i=1,1)/0/
      DATA (ab(94,i),i=1,1)/100.d0/
c
      DATA at(95),gel(95),nmn(95),(mn(95,i),i=1,1)/'Am',-1,1,243/
      DATA (zm(95,i),i=0,1)/243.061375d0, 243.0613811d0/
      DATA (ns2(95,i),i=1,1)/5/
      DATA (ab(95,i),i=1,1)/100.d0/
c
      DATA at(96),gel(96),nmn(96),(mn(96,i),i=1,1)/'Cm',-1,1,247/
      DATA (zm(96,i),i=0,1)/247.070347d0, 247.070354d0/
      DATA (ns2(96,i),i=1,1)/9/
      DATA (ab(96,i),i=1,1)/100.d0/
c
      DATA at(97),gel(97),nmn(97),(mn(97,i),i=1,1)/'Bk',-1,1,247/
      DATA (zm(97,i),i=0,1)/247.070300d0, 247.070307d0/
      DATA (ns2(97,i),i=1,1)/3/
      DATA (ab(97,i),i=1,1)/100.d0/
c
      DATA at(98),gel(98),nmn(98),(mn(98,i),i=1,1)/'Cf',-1,1,251/
      DATA (zm(98,i),i=0,1)/251.079580d0, 251.079587d0/
      DATA (ns2(98,i),i=1,1)/1/
      DATA (ab(98,i),i=1,1)/100.d0/
c
      DATA at(99),gel(99),nmn(99),(mn(99,i),i=1,1)/'Es',-1,1,252/
      DATA (zm(99,i),i=0,1)/252.082944d0, 252.082980d0/
      DATA (ns2(99,i),i=1,1)/10/
      DATA (ab(99,i),i=1,1)/100.d0/
c
      DATA at(100),gel(100),nmn(100),(mn(100,i),i=1,1)/'Fm',-1,1,257/
      DATA (zm(100,i),i=0,1)/257.095099d0, 257.095105d0/
      DATA (ns2(100,i),i=1,1)/9/
      DATA (ab(100,i),i=1,1)/100.d0/
c
      DATA at(101),gel(101),nmn(101),(mn(101,i),i=1,1)/'Md',-1,1,258/
      DATA (zm(101,i),i=0,1)/258.09857d0, 258.098431d0/
      DATA (ns2(101,i),i=1,1)/16/
      DATA (ab(101,i),i=1,1)/100.d0/
c
      DATA at(102),gel(102),nmn(102),(mn(102,i),i=1,1)/'No',-1,1,259/
      DATA (zm(102,i),i=0,1)/259.100931d0, 259.101030d0/
      DATA (ns2(102,i),i=1,1)/9/
      DATA (ab(102,i),i=1,1)/100.d0/
c
      DATA at(103),gel(103),nmn(103),(mn(103,i),i=1,1)/'Lr',-1,1,260/
      DATA (zm(103,i),i=0,1)/260.105320d0, 260.105500d0/
      DATA (ns2(103,i),i=1,1)/-1/
      DATA (ab(103,i),i=1,1)/100.d0/
c
      DATA at(104),gel(104),nmn(104),(mn(104,i),i=1,1)/'Rf',-1,1,261/
      DATA (zm(104,i),i=0,1)/261.10869d0, 261.108770d0/
      DATA (ns2(104,i),i=1,1)/-1/
      DATA (ab(104,i),i=1,1)/100.d0/
c
      DATA at(105),gel(105),nmn(105),(mn(105,i),i=1,1)/'Db',-1,1,262/
      DATA (zm(105,i),i=0,1)/262.11376d0, 262.114080d0/
      DATA (ns2(105,i),i=1,1)/-1/
      DATA (ab(105,i),i=1,1)/100.d0/
c
      DATA at(106),gel(106),nmn(106),(mn(106,i),i=1,1)/'Sg',-1,1,263/
      DATA (zm(106,i),i=0,1)/263.11822d0, 263.118320d0/
      DATA (ns2(106,i),i=1,1)/-1/
      DATA (ab(106,i),i=1,1)/100.d0/
c
      DATA at(107),gel(107),nmn(107),(mn(107,i),i=1,1)/'Bh',-1,1,262/
      DATA (zm(107,i),i=0,1)/262.12293d0, 262.122890d0/
      DATA (ns2(107,i),i=1,1)/-1/
      DATA (ab(107,i),i=1,1)/100.d0/
c
      DATA at(108),gel(108),nmn(108),(mn(108,i),i=1,1)/'Hs',-1,1,265/
      DATA (zm(108,i),i=0,1)/265.13016d0, 265.130090d0/
      DATA (ns2(108,i),i=1,1)/-1/
      DATA (ab(108,i),i=1,1)/100.d0/
c
      DATA at(109),gel(109),nmn(109),(mn(109,i),i=1,1)/'Mt',-1,1,266/
      DATA (zm(109,i),i=0,1)/266.13764d0, 266.137300d0/
      DATA (ns2(109,i),i=1,1)/-1/
      DATA (ab(109,i),i=1,1)/100.d0/
c
      IF((IAN.LE.0).OR.(IAN.GT.109)) THEN
          MASS= 0.d0
          NAME= 'XX'
          IMN= 0
          WRITE(6,601) IAN
          RETURN
        ELSE
          NAME= AT(IAN)
        ENDIF
      IF((IAN.EQ.1).AND.(IMN.NE.1)) THEN
c** Special case: insert common name for deuterium or tritium
          IF(IMN.EQ.2) NAME=' D'
          IF(IMN.EQ.3) NAME=' T'
          IF(IMN.EQ.4) NAME=' p'
          IF(IMN.EQ.5) NAME=' d'
          IF(IMN.EQ.6) NAME=' t'
          ENDIF
      GELGS= GEL(IAN)
      MASS= -1.d0
      GNS= -1
      ABUND = -1.d0
      DO  I= 1,NMN(IAN)
          if(i.gt.10)  write(6,606) ian,imn,nmn(ian)
  606  format(3i9)
          IF(IMN.EQ.MN(IAN,I)) THEN
              MASS= ZM(IAN,I)
              GNS= NS2(IAN,I)+1
              ABUND = AB(IAN,I)
              ENDIF
          ENDDO
      IF(MASS.LT.0.d0) THEN
          MASS= ZM(IAN,0)
          IF(IMN.NE.0) WRITE(6,602) AT(IAN),IMN
          IMN= 0
          ENDIF
      RETURN
  601 FORMAT(' *** MASSES Data base does not include Atomic Number=',i4)
  602 FORMAT(' *** MASSES Data base does not include ',A2,'(',i3,
     1 '), so use average atomic mass.')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
c***********************************************************************

