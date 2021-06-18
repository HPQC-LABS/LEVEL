c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c** Generate analytic potential  VV(i)  as specified by the choice
c  of parameter IPOTL (see comments in PREPOT (& in main program))
c** All potentials generated in units cm-1 with absolute asymptote at
c  (input) energy VLIM for distance array  X0(i) Angstroms.
c** Return with NCN equal to power of asymptotically dominant inverse
c  power term in long range part of potential
c** Born-Oppenheimer correction functions in IPOTL=3 option may have up
c  to NBOBmx+1 terms.  ||    ****** last updated  9 Sept 2016  *********
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'arrsizes.h'       !! import array dimension parameters
      INTEGER  I,J,M,IBOB,IAN1,IAN2,IMN1,IMN2,MN1R,MN2R,IORD,IORDD,
     1 IPOTL,IMIN,PAD,QAD,QNA,NU1,NU2,NT1,NT2,NCMAX,PPAR,QPAR,NCN,Nbeta,
     2 APSE,NVARB,NPP,LNPT,GNS,GEL,NCMM,MCMM,sVSR2,IDSTT,MM1,
     3 MMLR(NCMMAX)
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,AT,BT,BETA,BINF,B1,B2,CSAV,U1INF,U2INF,
     1 T1INF,T2INF,YPAD,YQAD,YQADSM,YQNA,YQNASM,ABUND,CNN,DSCM,DX,DX1,
     2 FCT,FC1,FC2,FG1,FG2,MASS1,MASS2,RMASS1,RMASS2,REQ,Rref,Rinn,
     3 Rout,SC1,SC2,SG1,SG2,VLIM,DVLIM,VMIN,XDF,X1,XS,XL,XP1,ZZ,ZP,ZQ,
     4 ZME,ULR,ULRe,rhoAB,rhoINT,nDSCM,nREQ,Scalc,XXQ,REQp,REQq,RREFp,
     5 RREFq,DSUM,DSUMP,bohr,Rbohr,T0,dULRdR,RM3,BFCT2,RH,f2,f2p,GAMMA,
     6 dULRdCm(NCMMAX),
     7 DM(NCMMAX),DMP(NCMMAX),DMPP(NCMMAX),CmVAL(NCMMAX),CmEFF(NCMMAX),
     6 U1(0:NBOBmx),U2(0:NBOBmx),T1(0:NBOBmx),T2(0:NBOBmx),
     9 PARM(NbetaMX),XPARM(NbetaMX),rKL(NbetaMX,NbetaMX),XO(NDIMR),
     a VV(NDIMR),RM2(NDIMR),bTT(-1:2),cDS(-2:0),bDS(-2:0)
      SAVE IBOB,IPOTL,PPAR,QPAR,PAD,QAD,QNA,Nbeta,MMLR,NVARB,NCMM
      SAVE DSCM,REQ,Rref,PARM,U1,U2,T1,T2,CSAV,BINF,ALFA,ZME,
     1 Rinn,Rout,ULR,ULRe,CmVAL,XPARM
c** Damping function parameters for use and printout .....
      DATA bTT/2.44d0,2.78d0,3.13d0,3.47d0/
      DATA bDS/3.3d0,3.69d0,3.95d0/
      DATA cDS/0.423d0,0.40d0,0.39d0/
      SAVE bTT, bDS, cDS
c** Electron mass, as per 2010 physical constants
      DATA ZME/5.4857990946d-4/,bohr/0.52917721092d0/
c
      IF(LNPT.GT.0) THEN
c** Most parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
          READ(5,*) IPOTL, QPAR, PPAR, Nbeta, APSE, IBOB
          READ(5,*) DSCM, REQ, Rref
          IF(IPOTL.GE.4) THEN
c** For MLR, DELR, HFD, Tang-Toennies or Tiemann-polynomial potentials .....
c   For each long-range term read power  MMLR(i)  & coefficient CmVAL(i)
c** For special Aubert-Frecon 2x2 cases,  NCMM= 7,  MMLR= {x,3,3,6,6,8,8},
c   with x= 0 for the A state, x= -1 for the b state, and  CmVAL= {Aso, 
c   C3Sig, C3Pi, C6Sig, C6Pi, C8Sig, C8Pi}, 
c* while for the 3x3 diagonalization cases, NCMM=10, MMLR= {x,3,3,3,6,6,6,
c   8,8,8}  with x= -2 for the (lowest eigenvalue) c(1\,^3\Sigma_g^+ state,  
c   x= -3 for the (middle root) B^1\Pi_u state, and x=-4  for the 
c   highest-root state, while CmVal= {Aso, C3Sig, C3Pi1, C3Pi3, C6Sig, C6Pi1,
c    C6Pi3, C8Sig, C8Pi1, C8Pi3}
c=======================================================================
              READ(5,*) NCMM, rhoAB, sVSR2, IDSTT
              DO m=1, NCMM
                  READ(5,*) MMLR(m), CmVAL(m)
                  CmEFF(m)= CmVAL(m)
                  ENDDO
              MCMM= NCMM
              ENDIF
c-----------------------------------------------------------------------
          IF(IPOTL.EQ.1) NVARB= 0
          IF(IPOTL.EQ.2) THEN
              NVARB= Nbeta+2
              ENDIF
          IF(IPOTL.EQ.3) THEN
              NVARB= Nbeta+1
              IF(QPAR.LE.0) NVARB=2
              ENDIF
          IF(IPOTL.EQ.4) THEN
              NVARB= Nbeta+ 1
              IF(APSE.GT.0) NVARB= Nbeta
              ENDIF
          IF(IPOTL.EQ.5) THEN
              IORD= Nbeta
              NVARB= IORD+ 1
              ENDIF
          IF(IPOTL.EQ.6) NVARB= Nbeta
          IF(IPOTL.EQ.7) NVARB= 9
          IF(IPOTL.EQ.8) NVARB= Nbeta+ 4
c-----------------------------------------------------------------------
          IF(NVARB.GT.0) THEN
              IF((IPOTL.EQ.4).AND.(APSE.GT.0)) THEN
                  DO I=1, NVARB
                      READ(5,*) XPARM(I),PARM(I)
                      ENDDO
                ELSE
                  READ(5,*) (PARM(I),I=1,NVARB)
                ENDIF  
              ENDIF 
c-----------------------------------------------------------------------
          IF(IBOB.GT.0) THEN
c-----------------------------------------------------------------------
              READ(5,*) MN1R, MN2R, qAD, pAD, NU1, NU2, qNA, NT1, NT2
c-----------------------------------------------------------------------
              NCMAX= MAX0(NU1,NU2,NT1,NT2)
              IF(NCMAX.LT.0) THEN
                  IBOB= 0
                ELSE
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
                  DVLIM= 0.d0
c=======================================================================
c** Read actual BOB polynomial expansion coefficients
c=======================================================================
                  IF(NU1.GE.0) THEN
c... use Huang/Le Roy form for atom-1 adiabatic potential BOB radial fx.
c-----------------------------------------------------------------------
                      READ(5,*) (U1(I), I=0,NU1)
                      READ(5,*) U1INF
c-----------------------------------------------------------------------
                      FC1= 1.d0 - RMASS1/MASS1
                      DVLIM= DVLIM + FC1*U1INF
                      WRITE(6,630) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1        1,U1INF,PAD,PAD,PAD,PAD,PAD,PAD,NU1,QAD,QAD,QAD,QAD,QAD,
     2                                         NU1+1,(U1(I),I= 0,NU1)
                      ENDIF
c
                  IF(NU2.GE.0) THEN
c... use Huang/Le Roy form for atom-2 adiabatic potential BOB radial fx.
c-----------------------------------------------------------------------
                      READ(5,*) (U2(I), I=0,NU2)
                      READ(5,*) U2INF
c-----------------------------------------------------------------------
                      FC2= 1.d0 - RMASS2/MASS2
                      DVLIM= DVLIM + FC2*U2INF
                      WRITE(6,630) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1        1,U2INF,PAD,PAD,PAD,PAD,PAD,PAD,NU2,QAD,QAD,QAD,QAD,QAD,
     2                                         NU2+1,(U2(I),I= 0,NU2)
                      ENDIF
c
                  IF(NT1.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-1 ...
c-----------------------------------------------------------------------
                      READ(5,*) (T1(I), I=0,NT1)
                      READ(5,*) T1INF
c-----------------------------------------------------------------------
                      WRITE(6,634) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1 1,T1INF,QNA,QNA,QNA,QNA,QNA,QNA,NT1,QNA,NT1+1,(T1(I),I= 0,NT1)
                      FG1= RMASS1/MASS1
                      ENDIF
c
                  IF(NT2.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-2 ...
c-----------------------------------------------------------------------
                      READ(5,*) (T2(I), I=0,NT2)
                      READ(5,*) T2INF
c-----------------------------------------------------------------------
                      WRITE(6,634) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1 2,T2INF,QNA,QNA,QNA,QNA,QNA,QNA,NT2,QNA,NT2+1,(T2(I),I= 0,NT2)
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
                ENDIF
              ENDIF
          IF(IPOTL.GE.4) THEN           !! now describe long-range tail
              IF(rhoAB.GT.0.d0) THEN
                  IF(IDSTT.GT.0) WRITE(6,660) rhoAB,sVSR2,bDS(sVSR2),
     1                                                cDS(sVSR2),sVSR2
                  IF(IDSTT.LE.0) THEN
                      IF(IPOTL.NE.7) WRITE(6,662) rhoAB,sVSR2/2,
     1                                                    bTT(sVSR2/2)
                      IF(IPOTL.EQ.7) WRITE(6,663) rhoAB,sVSR2/2
                      ENDIF
                ELSE
                  WRITE(6,664) 
                ENDIF
              IF(MMLR(1).LE.0) THEN
c** uLR printout for Lyon 2x2 or 3x3 treatment of 2S + 2p alkali dimers ...
                  IF((NCMM.NE.7).AND.(NCMM.NE.10)) THEN
                      WRITE(6,666) MMLR(1),NCMM
                      STOP
                      ENDIF
                  IF(MMLR(1).EQ.0) WRITE(6,668) 'A-state',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
                  IF(MMLR(1).EQ.-1) WRITE(6,668) 'b-state',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
c... For Lyon treatment of b-state alkali dimers ...
                  IF(MMLR(1).EQ.-2) WRITE(6,670) 'c-state',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
                  IF(MMLR(1).EQ.-3) WRITE(6,670) 'B-state',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
                  IF(MMLR(1).EQ.-4) WRITE(6,670) '1 ^3Pi',CmVAL(1),
     1                                    CmVAL(2),(CmVAL(m),m=3,NCMM)
                ELSE
c... uLR printout for 'conventional' (damped or non-damped) inverse-power sum
                  WRITE(6,672) NCMM,(MMLR(m),CmEFF(m),m= 1,NCMM)
                ENDIF
              ENDIF
          ENDIF
c
c=======================================================================
c** Generate a  Lennard-Jones(QPAR,PPAR)  potential here.
c=======================================================================
      IF(IPOTL.EQ.1) THEN 
          XS= QPAR
          XL= PPAR
          XDF= DSCM/(XS-XL)
          IF(LNPT.GT.0) WRITE(6,600) QPAR, PPAR, DSCM, REQ
          CNN= XS*XDF*REQ**QPAR
          NCN= PPAR
          DO  I= 1,NPP
              VV(I)= (XL*(REQ/XO(I))**QPAR - XS*(REQ/XO(I))**PPAR)*XDF
     1                  +VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.2) THEN
c=======================================================================
c** Generate Seto-modified form of Surkus' GPEF function which includes
c  Dunham, SPF and OT forms as special cases.
c=======================================================================
          VMIN= VLIM
          X1= 1.d0
          A0= DSCM
          IF((PPAR.NE.0).AND.(DABS(PARM(Nbeta+1)).GT.0.d0)) THEN
              FCT= 1.d0/PARM(Nbeta+1)
              DO  J=1,IORD
                  X1= X1+ PARM(J)*FCT**J
                  ENDDO
c... Actual Dissoc. limit for this GPEF power series potential
              DSCM= A0*X1*FCT**2 + VMIN
              ENDIF
          IF(PPAR.EQ.1) THEN
c  Cases with power =1 (including Dunham, SPF & O-T expansions).
              IF(DABS(PARM(Nbeta+1)).LE.0.d0) THEN
c ... print for Dunham expansion ...
                  WRITE(6,612) PARM(Nbeta+2),REQ,VMIN,A0,Nbeta,
     1                                              (PARM(I),I= 1,Nbeta)
                  NCN= -99
                  CNN= 0.d0
                  ENDIF
              IF(DABS(PARM(Nbeta+2)).LE.0.d0) THEN
c ... print for Simons-Parr-Finlan expansion ...
                  WRITE(6,614) PARM(Nbeta+1),REQ,DSCM,A0,Nbeta,
     1                                              (PARM(I),I= 1,Nbeta)
                  NCN= 1
                  ENDIF
              IF(DABS(PARM(Nbeta+2)-PARM(Nbeta+1)).LE.0.d0) THEN
c ... print for Ogilvie-Tipping expansion ...
                  WRITE(6,616) PARM(Nbeta+2),REQ,DSCM,A0,Nbeta,
     1                                              (PARM(I),I= 1,Nbeta)
                  NCN= 1
                  ENDIF
              ENDIF
          IF((PPAR.NE.0).AND.((PPAR.NE.1).OR.
     1               ((DABS(PARM(Nbeta+2)-PARM(Nbeta+1)).GT.0.d0).AND.
     2             (DABS(PARM(Nbeta+2)*PARM(Nbeta+1)).GT.0.d0)))) THEN
c ... print for general GPEF expansion variable ...
              IF(PPAR.LT.0) THEN
c ... for negative PPAR, convert to equivalent positive PPAR case
                  PPAR= -PPAR
                  A1= PARM(Nbeta+2)
                  PARM(Nbeta+2)= -PARM(Nbeta+1)
                  PARM(Nbeta+1)= -A1
                  ENDIF
              WRITE(6,618) PPAR,PPAR,PARM(Nbeta+1),PPAR,PARM(Nbeta+2),
     1                     PPAR,REQ,DSCM,A0,Nbeta,(PARM(I),I= 1,Nbeta)
              NCN= PPAR
              ENDIF
          IF(PPAR.EQ.0) THEN
c** For case of simple power series in  R  itself
              NCN= -1
              WRITE(6,620) Nbeta,VMIN,(PARM(I),I= 1,Nbeta)
              DO  I= 1, NPP
                  ZP= 1.d0
                  A1= VMIN
                  DO  J= 1,Nbeta
                      ZP= ZP*XO(I)
                      A1= A1+ PARM(J)*ZP
                      ENDDO
                  VV(I)= A1
                  ENDDO
c ...Reset asymptote to avoid spurious  E > VLIM  warnings (e.g. for HO)
cc            VLIM= VV(NPP)
              RETURN
              ENDIF
c ... otherwise - generate potential as a GPEF-type expansion
          DO  I= 1, NPP
              ZZ= (XO(I)**PPAR - REQ**PPAR)/(PARM(Nbeta+1)*XO(I)**PPAR
     1                                        + PARM(Nbeta+2)*REQ**PPAR)
              A1= 1.d0
              ZP= 1.d0
              DO  J=1, Nbeta
                  ZP= ZP*ZZ
                  A1= A1+ PARM(J)*ZP
                  ENDDO
              VV(I)= A0*ZZ*ZZ*A1 + VMIN
              ENDDO
c ...Reset asymptote to avoid spurious  E > VLIM  warnings (e.g. for HO)
          IF(DABS(PARM(Nbeta+1)).LE.0) VLIM= VMIN + MIN(VV(NPP),VV(1))
          ENDIF
c
c=======================================================================
c** Generate a simple Morse, or Extended (EMOp) Morse potential, or as
c  a special cases, Wei Hua's 4-parameter generalized Morse
c=======================================================================
      IF(IPOTL.EQ.3) THEN
          IF(Rref.LE.0.d0) Rref= REQ
          BETA= PARM(1)
          NCN= 99
          IF(LNPT.GE.0) THEN
              IF(QPAR.GT.0) THEN
c... Normal case is Morse or EMO
                  IF(Nbeta.EQ.0) THEN
                      WRITE(6,606) DSCM,REQ,BETA
                    ELSE
                      WRITE(6,608) QPAR,DSCM,REQ,Rref,Nbeta,QPAR,QPAR,
     1                            QPAR,QPAR,NVARB,(PARM(i),i= 1,NVARB)
                    ENDIF
                ELSE
c... Option to generate Wei Hua's extended 4-parameter Morse-type potl.
                  CSAV= PARM(2)
                  WRITE(6,605) DSCM,REQ,CSAV,BETA
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
c... for Wei Hua's extended Morse function ...
              IF(QPAR.LE.0) THEN
                  VV(I)= DSCM*((1.d0 - DEXP(-BETA*(XO(I)-REQ)))/(1.d0 
     1                - CSAV*DEXP(-BETA*(XO(I)-REQ))))**2 - DSCM+ VLIM
                ELSE 
c... for proper Morse or EMO function ...
                  IF(Nbeta.GE.1) THEN
                      ZZ= (XO(I)- Rref)/(XO(I)+ Rref)
c... for proper LeRoy-Huang yp(r) expansion ...
                      IF(QPAR.GT.1) ZZ= (XO(i)**QPAR - Rref**QPAR)/
     1                                  (XO(i)**QPAR + Rref**QPAR)
                      BETA= 0.d0
                      DO  J= Nbeta,0,-1
                          BETA= BETA*ZZ+ PARM(J+1)
                          ENDDO
                      ENDIF
                  VV(I)=  DSCM*(1.d0 - DEXP(-BETA*(XO(I)-REQ)))**2 
     1                                                    - DSCM+ VLIM
                ENDIF
              ENDDO
          ENDIF
c=======================================================================
c** Generate an MLR potential [as per J.Chem.Phys. 131, 204309 (2009)]
c=======================================================================
      IF(IPOTL.EQ.4) THEN
          IF(LNPT.GT.0) THEN
c** for a new case ... define ULRE and print potential description
              NCN= MMLR(1)
              IF(NCN.LE.0) NCN= MMLR(2)
              CNN= CmVAL(1)
              ULRe= 0.d0
c*** print for MLR form
              WRITE(6,602) QPAR,PPAR,DSCM,REQ
c... for Huang form: \beta(yp)= Binf*yp + [1-yp]*{power series in yq}
              IF(APSE.LE.0) WRITE(6,607) PPAR,PPAR,QPAR,Nbeta,RREF,
     1                                  Nbeta+1,(PARM(J),J= 1,Nbeta+1)
c... print for Asen Pashov Spline Exponent (APSE > 0) MLR form
              IF(APSE.GT.0) THEN 
                  WRITE(6,604) PPAR,Nbeta,(PARM(J),J= 1,Nbeta) 
                  WRITE(6,610) QPAR,Rref,(XPARM(J),J= 1,Nbeta)
c** Prepare Asen's Rlk array for later use in generating Spline fx. 
                  CALL Lkoef(Nbeta,XPARM,rKL)
                  ENDIF
c=======================================================================
              CALL quadCORR(NCMM,MCMM,NCMMAX,MMLR,DSCM,CmVAL,CmEFF)
c=======================================================================
c** Now - initialize at r= REQ 
              IF(MMLR(1).LE.0) THEN
c..... for AF 2x2 or 3x3 case ...
                  CALL AFdiag(REQ,VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                sVSR2,IDSTT,ULRe,dULRdCm,dULRdR)
                ELSE
c..... or for 'simple' (damped) inverse-power sum 
                  CALL dampF(REQ,rhoAB,MCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
                  DO  J= 1,MCMM
                      ULRe= ULRe + DM(J)*CmEFF(J)/REQ**MMLR(J)
                      ENDDO
                ENDIF
              BINF= DLOG(2.d0*DSCM/ULRe)
              WRITE(6,674) BINF,ULRe
              REQp= REQ**PPAR
              RREFp= Rref**PPAR
              RREFq= Rref**QPAR
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
              ZZ= (XO(i)**PPAR- REQp)/(XO(i)**PPAR+ REQp)
              ZP= (XO(i)**PPAR- RREFp)/(XO(i)**PPAR+ RREFp)
              ZQ= (XO(i)**QPAR- RREFq)/(XO(i)**QPAR+ RREFq)
              IF(APSE.LE.0) THEN
c... for Huang/THEOCHEM constrained polynomial for \beta(r) ...
                  BETA= 0.d0
                  DO  J= Nbeta,0,-1
                      BETA= BETA*ZQ+ PARM(J+1)
                      ENDDO
c...  calculate constrained polynomial MLR exponent coefficient
                  BETA= BINF*ZP + (1.d0- ZP)*BETA
                ELSE
c... calculate Pashov cubic spline exponent coefficient  ...
                  BETA= 0.d0
                  DO m= 1,Nbeta
                      BETA= BETA + Scalc(ZQ,m,Nbeta,XPARM,rKL,NbetaMX)
     1                                                        *PARM(m)
                      ENDDO
                ENDIF
c** Now Calculate local value of uLR(r)
              ULR= 0.d0
              IF(MMLR(1).LE.0) THEN
c..... for AF 2x2 or 3x3 case ...
                  CALL AFdiag(XO(i),VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                 sVSR2,IDSTT,ULR,dULRdCm,dULRdR)
                ELSE
c..... or for 'simple' (damped) inverse-power sum 
                  CALL dampF(XO(i),rhoAB,MCMM,NCMMAX,MMLR,sVSR2,IDSTT,
     1                                                    DM,DMP,DMPP)
                  ULR= 0.d0
                  DO  J= 1,MCMM
                      ULR= ULR + DM(J)*CmEFF(J)/XO(I)**MMLR(J)
                      ENDDO
                ENDIF
              BETA= (ULR/ULRe)*DEXP(-BETA*ZZ)
              VV(I)= DSCM*(1.d0 - BETA)**2 - DSCM + VLIM
c???????? Print for testing !!!!!!!
cc            write(8,777)xo(i),ulr,ulre,beta,(ULR/ULRe)*DEXP(-BETA*ZZ)
cc   1                                                          ,VV(I)
cc777 format('    r=',f9.4,'   ulr=',1pd12.5,'   uLRe=',d12.5,'   beta='
cc   1   ,d12.5,'   XP=',D14.7, '   V=',d14.7)         
c???????? Print for testing !!!!!!!
              ENDDO
          ENDIF
c=======================================================================
c** Generate a DELR potential [as per JCP 119, 7398 (2003) {revised}] 
c=======================================================================
      IF(IPOTL.EQ.5) THEN
          IF(LNPT.GT.0) THEN
              REQq= REQ**QPAR
              RREFq= Rref**QPAR 
              ZZ= (REQq - RREFq)/(REQq + RREFq)
              BETA= 0.d0
                  DO  J= Nbeta,0,-1
                      BETA= BETA*ZZ+ PARM(J+1)
                      ENDDO
              ULRe= 0.0d0
              B1= 0.0d0
c... First, calculations @ Re to get  AA & BB 
              CALL dampF(REQ,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
              IF(MMLR(1).LE.0) THEN  !! for A-F 2x2 or 3x3 uLR fx,
                  CALL AFdiag(REQ,VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                sVSR2,IDSTT,ULRe,dULRdCm,dULRdR)
                  B1= dULRdR
                ELSE             !! For conventional inverse-power sum
                  DO  J= 1,MCMM
                      T0= CmEFF(J)/REQ**MMLR(J)
                      ULRe= ULRe+ T0*DM(J)
                      B1= B1+ T0*(DMP(J) - DM(J)*MMLR(J)/REQ)
                      ENDDO
                ENDIF
              A1= DSCM - ULRe - B1/BETA
              B1= 2.d0*A1 + B1/BETA
              WRITE(6,650) QPAR,DSCM,REQ,Nbeta,(PARM(I),I= 1,IORD+1)
              WRITE(6,652) QPAR,QPAR,QPAR,QPAR,QPAR
              IF(Rref.GT.0.d0) WRITE(6,654) Rref
              IF(Rref.LE.0.d0) WRITE(6,656) REQ
              WRITE(6,658) A1,B1,NCMM
              ENDIF
c** Now ... generate potential function array for DELR form
          DO  I= 1, NPP
              XXQ= XO(I)**QPAR
              ZZ= (XXQ - RREFq)/(XXQ + RREFq)
              BETA= 0.d0
c ... calculate the exponent
              DO  J= Nbeta,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
               BETA= DEXP(-BETA*(XO(I)-REQ))
c ... calculate the (damped) long-range tail
              CALL dampF(XO(I),rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
              IF(MMLR(1).LE.0) THEN  !! for A-F 2x2 or 3x3 uLR fx,
                  CALL AFdiag(XO(i),VLIM,NCMM,NCMMax,MMLR,CmEFF,rhoAB,
     1                                 sVSR2,IDSTT,ULR,dULRdCm,dULRdR)
                ELSE
                  ULR= 0.0d0
                  DO  J= 1, MCMM
                      ULR= ULR+ DM(J)*CmEFF(J)/XO(I)**MMLR(J)
                      ENDDO
                ENDIF
              VV(I)=  (A1*BETA - B1)*BETA - ULR + VLIM
              ENDDO
          ENDIF
c
      IF((IPOTL.EQ.6).AND.(Nbeta.EQ.5)) THEN
c=======================================================================
c** For generalized  HFDB(m= MMLR(j), j=1,NCMM) potential 
c          V(r) = ALFA*(r/R_e)**PARM(5) * exp[-BETR*r - PARM(4)*r**2] 
c - D(r)*[CmEFF(1)/r**MMLR(1)+ CmEFF(2)/r**sMMLR(2)+ CmEFF(3)/r**MMLR(3)+ ...
c     and    D(r) = 1 for r > PARM(2)   and
c            D(x)= exp[-PARM(1)*(PARM(2)/r - 1)**PARM(3)] for  r < PARM(2)
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= CmEFF(1)
              A1= PARM(1)
              A2= PARM(2)
              A3= PARM(3)
              B2= PARM(4)
              DX= 1.d0
              DX1= 0.d0
              IF(A2.GT.1.d0) THEN     !!!!!!!!!!!!!!!!!!! GT.REQ) THEN
                  DX= DEXP(-A1*(A2/REQ - 1.d0)**A3)
                  DX1= A1*A2*A3*DX*(A2/REQ - 1.d0)**(A3- 1.d0)/REQ**2
                  ENDIF
              DSUM= 0.d0
              DSUMP= 0.d0
              DO  J= 1, NCMM
                  B1= CmEFF(J)/REQ**MMLR(j)
                  DSUM= DSUM + B1
                  DSUMP= DSUMP + B1*(DX1 - DX*MMLR(j)/REQ)
                  ENDDO
              ALFA= DSUM*DX -DSCM
              IF(ALFA.LE.0.d0) THEN
                  WRITE(6,622) ALFA,(MMLR(J),CmEFF(J),J= 1, NCMM)
                  STOP
                  ENDIF
              B1= PARM(5)/REQ - 2.d0*B2*REQ - DSUMP/ALFA
              ALFA= ALFA*DEXP(REQ*(B1 + B2*REQ))
              WRITE(6,624) A1,A2,A3 
              WRITE(6,626) 'ABC',PARM(5),DSCM,REQ,B1,B2,ALFA 
              ENDIF
          DO  I= 1,NPP
              X1= XO(I)
              XP1= 0.0D0
              IF((X1*(B1+ B2*X1)).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*(X1/REQ)**PARM(5)
              FC1= 0.d0
              DO  J= 1, NCMM
                  FC1= FC1 + CmEFF(J)/X1**MMLR(J)
                  ENDDO
              IF(X1.LT.A2) FC1= FC1*DEXP(-A1*(A2/X1- 1.d0)**A3)
              VV(I)= ALFA*XP1- FC1 + VLIM
              ENDDO
          ENDIF
c
      IF((IPOTL.EQ.6).AND.(Nbeta.EQ.2)) THEN
c=======================================================================
c** For generalized  HFD-ID(m= MMLR(j), j=1,NCMM) potential 
c   V(r) = ALFA*x**PARM(5) * exp[-BETR*r - PARM(4)*r**2]  - f2(\rho*r) *
c    \sum_m{ D_m^{ds}(\rho*r)*CmEFF(m)/r**MMLR(m)} with  x=r/R_e  and
c    f2(\rho*Rbohr)= (\rho*Rbohr)^{1.68} * exp{-0.78*\rho*Rbohr}
c=======================================================================
          IF(LNPT.GT.0) THEN
              B2= PARM(1)
              GAMMA= PARM(2)
              NCN= MMLR(1)
              CNN= CmEFF(1)
              sVSR2= 0
              Rbohr= REQ/bohr
              f2= 1.d0 - (rhoAB*Rbohr)**1.68d0 *EXP(-0.78d0*rhoAB*Rbohr)
              f2p= (f2 - 1.d0)*(1.68d0/REQ - 0.78d0*rhoAB/bohr)
              CALL dampF(REQ,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
              DSUM= 0.d0
              DSUMP= 0.d0
              DO  m= 1, NCMM
                  B1= CmEFF(m)/REQ**MMLR(m)
                  DSUM= DSUM + DM(m)*B1
                  DSUMP= DSUMP + B1*(f2p*DM(m) + f2*(DMP(m) 
     1                                           - DM(m)*MMLR(m)/REQ))
                  ENDDO
              ALFA= f2*DSUM - DSCM
              IF(ALFA.LE.0.d0) THEN
                  WRITE(6,622) ALFA,(MMLR(J),CmEFF(J),J= 1, NCMM)
                  STOP
                  ENDIF
              B1= GAMMA/REQ - 2.d0*B2*REQ - DSUMP/ALFA
              ALFA= ALFA*DEXP(REQ*(B1 + B2*REQ))
              WRITE(6,625) 
              WRITE(6,626) 'ID ',GAMMA,DSCM,REQ,B1,B2,ALFA
              ENDIF
          DO  I= 1,NPP
              X1= XO(I)
              Rbohr= X1/bohr
              f2= 1.d0 - (rhoAB*Rbohr)**1.68d0 *EXP(-0.78d0*rhoAB*Rbohr)
              CALL dampF(X1,rhoAB,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
              DSUM= 0.d0
              DO  m= 1, NCMM
                  DSUM = DSUM+ DM(m)*CmEFF(m)/X1**MMLR(m)
                  ENDDO
              XP1= 0.0D0
              IF((X1*(B1+ B2*X1)).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*(X1/REQ)**GAMMA
              VV(I)= ALFA*XP1- F2*DSUM + VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.7) THEN
c=======================================================================
c** Generate Generalized Tang-Toennies (TT) type potential as desribed
c   in the LEVEL manual: JQSRT(submitted Feb. 2016)
c   NCMM = number of inverse-power long-range terms and NVARB = 9.
c   DSCM and Re are the reported PEC minimum parameters.  The powers and
c   coefficients of the NCMM inverse-power long-range terms are MMCM(j) 
c   and CmEFF(j), with damping fx defined by rhoAB, IDSTT & sVSR2
c=======================================================================
          NCN= MMLR(1)
          CNN= CmEFF(1)
          IDSTT= 0
          sVSR2= 2
c** Define  rhoINT for consistency with conventional TT(sVSR2=+2) damping fx.
          rhoINT= rhoAB/3.13d0       
          VMIN= VLIM
          IMIN= 1
          DO I= 1, NPP
c....generate potential function array
              CALL dampF(XO(I),rhoINT,NCMM,NCMMAX,MMLR,sVSR2,IDSTT,DM,
     1                                                       DMP,DMPP)
c....calculate the (damped) long range tail
              A3= 0.d0
              DO J= 1, NCMM
                  A3= A3+ DM(J)*CmEFF(J)/XO(I)**MMLR(J)
                  ENDDO
c....For Generalized TT model
              XP1= PARM(1)*XO(I)+ PARM(2)*XO(I)**2+ PARM(3)/XO(I)
     1                                              + PARM(4)/XO(I)**2
              VV(I)= (PARM(5) + PARM(6)*XO(I) + PARM(7)/XO(I) + 
     1       PARM(8)*XO(I)**2 +PARM(9)*XO(I)**3)*DEXP(-XP1) - A3 + VLIM
              IF(VV(I).LE.VMIN) THEN
c... search for potential minimum ...
                  VMIN= VV(I)
                  IMIN= I
                  ENDIF
              ENDDO
          IF(LNPT.GT.0) THEN
              WRITE(6,628) (PARM(i),i=1,9)
c*** Use quadratic approximation to determine REQ and DSCM
              IF(IMIN.GT.1) THEN
                  A1= VV(IMIN-1)
                  A2= VV(IMIN)
                  A3= VV(IMIN+1)
                  RH= XO(IMIN) - XO(IMIN-1)
                  B1= (A3- 2.d0*A2 + A1)/(2.d0*RH**2)             !! curvature
                  nREQ= XO(IMIN) + 0.5d0*RH - (A3-A2)/(2.d0*RH*B1)
                  A2= A2- B1*(XO(IMIN)-nREQ)**2
                  nDSCM= VLIM - A2
                  IF(LNPT.GT.0) WRITE(6,629) DSCM,REQ, nDSCM,nREQ
                ELSE  
                        WRITE(6,629) DSCM,REQ
                ENDIF
              ENDIF
          ENDIF
c
      IF(IPOTL.EQ.8) THEN
c=======================================================================
c** Generate Tiemann-type polynomial potential attached to inverse-power
c  tail and 1/R^{12} (or exponential) inner wall [PRA 63, 012710 (2000)]
c  Polynomial expansion variable is  z= [R - Rm]/[R + b*Rm] where 
c  expansion has constant and linear terms.  The read-in DSCM= De (well
c  depth), but  Rm (read in as REQ) is not precisely Re (for a1 .neq. 0)
c  NCMM= number of inverse-power long-range terms;  
c  NVARB= (polynomial order) + 4.  [PPAR and APSE are dummy parameters]
c** Read-in parameters PARM(i) are in order: the  (Nbeta+1)  polynomial
c  coefficients  a(0) - a(Nbeta), the expansion variable denominator
c  factor b=PARM(Nbeta+2), and the the inner and outer bounds on the 
c  polynomial domain, Tiemann's Rinn=PARM(Nbeta+3) & Rout=PARM(Nbeta+4),
c  respectively.  The powers and coefficients (-ve if attractive) of the
c  NCMM inverse-power long-range terms are MMCM(j) and CmEFF(j).
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= -CmEFF(1)
              A0= VLIM- DSCM
              BT= PARM(Nbeta+2)
              Rinn= PARM(Nbeta+3)
              Rout= PARM(Nbeta+4)
c** Determine analytic function attaching smoothly to inner wall of 
c  polynomial expansion at  R= Rinn < Rm
              ZZ= (Rinn - REQ)/(Rinn+ BT*REQ)
              ZP= 1.d0
              A1= PARM(1)
              A2= 0.d0
              DO  J= 1,Nbeta
                  A2= A2+ J*ZP*PARM(J+1)
                  ZP= ZP*ZZ
                  A1= A1+ ZP*PARM(J+1)
                  ENDDO
              A2= A2*(REQ+ BT*REQ)/(Rinn + BT*REQ)**2
c* If inward extrapolation is exponential:   A1*exp(-A2*(R-Rinn))
              A2= -A2/A1
c* If inward extrapolation is inverse-power:   A1 + A2/R**12
c*** To invoke this version, comment out precious line and UNcomment 
c     the next 2 lines
c             A2= -A2*Rinn**13/12.d0
c             A1= A1 - A2/Rinn**12 + VLIM - DSCM
c** With long-range tail an NCMM-term inverse-power sum, add 1 additional
c   higher-power term to ensure continuity (not smoothness) at  Rout
c** NOTE attractive long-range terms have negative (-) coefficients!
              ZZ= (Rout - REQ)/(Rout+ BT*REQ)
              ZP= 1.d0
              B1= PARM(1)
              DO  J= 1,Nbeta
                  ZP= ZP*ZZ
                  B1= B1+ ZP*PARM(J+1)
                  ENDDO
              A3= DSCM
              DO  J= 1,NCMM
                  A3= A3+ CmEFF(J)/Rout**MMLR(J)
                  ENDDO
              PPAR= NCMM+ 1
              MMLR(PPAR)= MMLR(NCMM)+ 2
              CmEFF(PPAR)= (B1-A3)*Rout**MMLR(PPAR)
c*** Print for Tiemann-type potential
              IF(LNPT.GE.0) THEN
                  WRITE(6,640) DSCM,REQ,PARM(Nbeta+2),Nbeta,Nbeta+1, 
     1                                            (PARM(J),J= 1,Nbeta+1)
ccc               IF(XO(1).LT.Rinn) WRITE(6,642) PARM(Nbeta+3),A1,A2,A0
                  IF(XO(1).LT.Rinn) WRITE(6,642) PARM(Nbeta+3),A1,A2
                  IF(XO(NPP).GT.Rout) WRITE(6,644) PARM(Nbeta+4),
     1                                     (CmEFF(J),MMLR(J),J= 1, PPAR)
                  ENDIF
              ENDIF
c ... now generate potential as a Tiemann-type expansion
          DO  I= 1, NPP
              IF(XO(I).LE.Rinn) THEN
c ... for exponential inward extrapolation ... for consistency with manual
                  VV(I)= A1*DEXP(-A2*(XO(I)- Rinn)) + A0
c ... for   A + B/R**12  inward extrapolation ... possible alternative
c                 VV(I)= A1 + A2/XO(I)**12
                ELSEIF(XO(I).LE.Rout) THEN
                  ZZ= (XO(I) - REQ)/(XO(I) + BT*REQ)
                  A3= A0 + PARM(1)
                  ZP= 1.d0
                  DO  J= 1,Nbeta
                      ZP= ZP*ZZ
                      A3= A3+ PARM(J+1)*ZP
                      ENDDO
                  VV(I)= A3
                ELSEIF(XO(I).GT.Rout) THEN
                  A3= VLIM
                  DO  J= 1, PPAR
                      A3= A3+ CmEFF(J)/XO(I)**MMLR(J)
                      ENDDO
                  VV(I)= A3
                ENDIF
              ENDDO
          ENDIF
c
      IF(IBOB.GT.0) THEN
c=======================================================================
c** If appropriate, generate Born-Oppenheimer breakdown correction 
c      functions to rotationless and/or centrifugal potential(s) using
c      LeRoy/Huang radial functions ...
c=======================================================================
          DO  I=1,NPP
              YPAD= (XO(I)**PAD- REQ**PAD)/(XO(I)**PAD+ REQ**PAD)
              YQAD= (XO(I)**QAD- REQ**QAD)/(XO(I)**QAD+ REQ**QAD)
              YQNA= (XO(I)**QNA- REQ**QNA)/(XO(I)**QNA+ REQ**QNA)
              SC1= U1INF*YPAD
              SC2= U2INF*YPAD
              SG1= T1INF*YQNA
              SG2= T2INF*YQNA
              YQADSM= (1.d0- YPAD)
              YQNASM= (1.d0- YQNA)
c ... finally, accumulate overall BOB terms ... all at the same time!
              DO  J= 0,NCMAX
                  SC1= SC1+ YQADSM*U1(J)
                  SC2= SC2+ YQADSM*U2(J)
                  SG1= SG1+ YQNASM*T1(J)
                  SG2= SG2+ YQNASM*T2(J)
                  YQADSM= YQADSM*YQAD
                  YQNASM= YQNASM*YQNA
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          VLIM= VLIM+ DVLIM
          IF((IPOTL.EQ.4).AND.(MMLR(1).LE.0)) THEN
c!! For mixed isotopopogue {6,7}Li_2(A) state, shift asymptote! ??? HUH ???
              IF(IMN1.NE.IMN2) THEN
                  DO  I= 1,NPP
                      RM3= (2.d0/3.d0)*CmEFF(1)/XO(I)**3
                      VV(I)= VV(I)+ RM3- DSQRT(RM3**2+ 3.085959756d-02)
                      ENDDO
                  VLIM= VLIM + DSQRT(3.085959756d-02)
                  ENDIF
c** For special case of A and c states of Li2, add BOB centrifugal term
              IF((MMLR(1).EQ.0).OR.(MMLR(1).EQ.-2)) THEN
                  BFCT2= 2.d0*16.857629206d0*(MASS1+MASS2)/(MASS1*MASS2)
                  DO  I= 1, NPP
                      VV(I)= VV(I) + BFCT2/XO(I)**2   !!!  ??? HUH ???
                      ENDDO
                  ENDIF
                  WRITE(6,646)
              ENDIF
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  602 FORMAT(/' MLR(q=',I1,', p=',I1,') Potential with:   De='
     1 ,F10.4,'[cm-1]    Re=',F12.8,'[A]')
  604 FORMAT('   with SE-MLR exponent coefft   beta(r)='/22x,'y',I1,
     1  '^{eq} *{Spline through the',I3,' function values} beta_i ='/
     2  (10x,4D16.8:))
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F12.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  607 FORMAT('  with PE-MLR exponent coefft:  beta(r)= beta{INF}*y',I1,
     1  ' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x,'exponent power series
     2 of order',I3,' in a variable in which   Rref=',f8.5/
     3   6x,'with',i3,' coefficients:',1PD17.9,2D17.9:/(10x,4D17.9:))
  608 FORMAT(/' EMO_',i1,' Potential with   De=',F11.4,'    Re=',F11.8,
     1 '   Rref=',F11.8/3x,'Exponent coeft: order-',i2,
     2 ' power series in  y=(r**',i1,' - Rref**',i1,')/(r**',i1,
     3 ' + Rref**',i1,')'/'   with',I3,' coefficients:',1x,1PD17.9,
     4 2D17.9:/(7X,4D17.9:))
  610 FORMAT(5x,'at distances defined by y_',I1,'(r; RREF) ='/
     1  (10x,4D16.8:))
  612 FORMAT(/' Potential is a Dunham expansion in  (r-Re)/(',f5.2,
     1  ' * Re)  with   Re=',f12.9/'  V(Re)=',f12.4,'    a0=',1PD16.9,
     2  '   and',i3,'  a_i coefficients:'/(5D16.8))
  614 FORMAT(/' Potential is an SPF expansion in  (r-Re)/(',F5.2,
     1  '* r)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',
     2  1PD16.9,'   and',i3,'  b_i  coefficients:'/(5D16.8))
  616 FORMAT(/' Potential is an O-T expansion in  (r-Re)/[',f5.2,
     1  '*(r+Re)]  with   Re=',f12.9/5x,'De=',G18.10,
     2  '   c0=',1PD16.9,'   and',i3,'  c_i coefficients:'/(5D16.8))
  618 FORMAT(/' Potential is a general GPEF expansion in  (r**',i1,
     1  ' - Re**',i1,')/(',SP,F5.2,'*r**',SS,i1,SP,F6.2,'*Re**',SS,i1,
     2  ')'/5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/
     3  5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
  620 FORMAT(/' Potential is a power series in  r  of  order',i3,
     1 ' with   V(r=0)=',f11.4/3x,'& coefficients (from linear term):',
     2 1P2d16.8:/(5x,4D16.8:))
  617 FORMAT('      betaINF=',f16.12,'  & uLR defined by  C',i1,' =',
     1  1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  622 FORMAT(/' *** ERROR in generating HFD potential *** generate   ALF
     1A=',1PD15.7,'  from reduced  Cm  coefficients:'/(3x,3('   C',I2,
     2 '=',D15.7:)) )
  624 FORMAT(15x,'and overall damping function:'/20x,'D(r)= exp[ -',
     1  0P,f8.6,'*(',f11.8,'/r -1.0)**',f5.2,']')
  625 FORMAT(15x,'and overall damping function:'/20x,'f2(r)= 1 - rhoAB*r
     3[bohr]^1.68 *exp{0.78*rhoAB*r[bohr]}')
  626 FORMAT(/' Potential is Generalized HFD-',a3,'  with   radial power
     1   gamma=',F9.6/ '   De=',f10.4,'[cm-1]   Re=',f9.6,'[Ang.],   wit
     2h  exponential-term factors:'
     3 5x,'beta1=',f11.8,'   beta=',f11.8,'   and A(pre-exp)=',1PD16.9)
  628 FORMAT(/' Generalized Tang-Tonnies Potential function with exponen
     1t function'/' - {{',SP,F15.11,'*r',F15.11,'*r^2',F15.11,'/r',
     2  F15.11,'/r^2}}'/' and pre-exp factor:'/3x,'{{',SP,1PD15.8,D16.8,
     3  '*r',d16.8,'/r',d16.8,'*r^2'/21x,D16.8,'*r^3}}',S)
  629 FORMAT(/10x,'Input    DSCM=',F10.4,'   REQ=',f9.6:/ 10x,
     1 'Actual   DSCM=',F10.4,'   REQ=',f9.6)
  630 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1 '  of mass ',f15.11/'   consists of mass factor  [1- MASS(',I3,
     2 A2,')/MASS(',I3,A2,')]  multiplying all of:'/5x,'u',I1,'INF=',
     3 f11.6,'  times  y',i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1, 
     4 ' + Re**',i1,')]  plus'/7x,'[1 - y',i1,']  times an order',I3, 
     5 ' polynomial in'/7x,'y',i1,'=[(r**',i1,' - Re**',i1,')/(r**',i1, 
     6 ' + Re**',i1,')]  with the ',i3,' coefficients:'/1P,(3x,4D17.9:))
  634 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor  [MASS(',I3,A2,')/MASS(',I3,
     2 A2,')]  multiplying all of:'/5x,'q',i1,'INF=',1PD17.9,
     3 ' times  y',i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',
     4 i1,')]'/ 3x,'plus [1 - y',i1,'] times an order',I3,' polynomial i
     6n y',i1,  '(r) with the',i3,' coefficients:'/(3x,1P,4D17.9:))
  636 FORMAT(3x,'where   fsw(r) = 1/[1 - exp{',f7.4,'*(r -',f7.4,')}]')
  638 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor   [mass{electron}/MASS(',I3,
     2 A2,')]'/'   multiplying   q',i1,'INF=',1PD17.9,'  times [1 - fsw(
     3r)/fsw(Re)]'/ '   plus  fsw(r)  times an order',0P,i3,' polynomial
     4 in z{O-T} with coefficients:'/ 1P,(3x,4D17.9:))
  640 FORMAT(/' Tiemann-type potential with   De=',F11.4,'   Rm=',f9.6,
     1 '   is a power series'/10x,'in  (r - Re)/(r ',SP,F9.5, 
     2 '*Re) of order',SS,I3,'  with the',I3,' coefficients:'/(5D16.8))
c 642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',1PD13.6,'*exp[-',
c    1 0P,F9.6,'*(r - Rinn)] ',SP,F10.3)
  642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',SP,F12.4,1x,1PD13.6,
     1  '/R**12' )
  644 FORMAT('  and  for  r > Rout=',F7.3,'   V= VLIM ',
     1 (SP,1PD14.6,'/r**',SS,I2):/(39x,SP,1PD14.6,'/r**',SS,I2))
  646 FORMAT(2x,'Potential fx. includes centrifugal BOB term  +2*B(r)')
  650 FORMAT(/' DELR(q=',i2,') Potential with   De=', F11.4,'[cm-1]   Re
     1=',F11.8,'[A]   where'/3x,'exponent coefft. has power series order
     2',I4/6x,'with polynomial coefficients',8x,1PD17.8,D17.8/ 
     3 (8x,4D17.8))
  652 FORMAT(6x,'where the radial variable   y_',I1,'= (r**',I1,' - Rref
     4**',i1,')/(r**',I1,' + Rref**',i1, ')')
  654 FORMAT(10x,'is defined w.r.t.   Rref=',F11.8) 
  656 FORMAT(10x,'is defined w.r.t.   Rref= Re= ',F11.8) 
  658 FORMAT(3x,'Generate A(DELR)=',1Pd17.9,'   B(DELR)=',D17.9/
     1 6x,'from uLR defined by',I2,' inverse-power terms')
  660 FORMAT(/' uLR inverse-power terms incorporate DS-type damping with
     1   rhoAB=',f9.6/8x,'defined to give very short-range  Dm(r)*Cm/r^m
     2  behaviour   r^{',SS,I2,'/2}'/8x,'Dm(r)= [1 - exp(-',f5.2,
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,I3,'/2}')
  662 FORMAT(/' uLR inverse-power terms incorporate TT-type damping with
     1   rhoAB=',f9.6/8x,'defined to give very short-range  Dm(r)*Cm/r^m
     2  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r)^
     3k/k!}]   where   bTT=',f6.3,'*rhoAB')
  663 FORMAT(/' uLR inverse-power terms incorporate TT-type damping with
     1   rhoAB=',f13.10/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     3)^k/k!}]   where   bTT= rhoAB')
  664 FORMAT(' uLR(r) inverse-power terms inlude NO individual-term damp
     1ing')
  666 FORMAT(4x,'*** ERROR ***  MMLR(1)=',I3,' A-F diagonalization not d
     1efined for  NCMM=', I3)
  668 FORMAT(5x,'Use Lyon 2x2  ',A7,'  uLR(r)  with   Aso=',F11.6/
     1  47x,'C_3(^1Sig)=',1P,D15.7:/47x,'C_3(^3Pi) =',D15.7:/
     1  47x,'C_6(^1Sig)=',1PD15.7:/47x,'C_6(^3Pi) =',D15.7:/
     1  47x,'C_8(^1Sig)=',1PD15.7:/47x,'C_8(^3Pi) =',D15.7)
  670 FORMAT(' Use Lyon 3x3 ',A7,'  uLR(r)  with   Aso=',F11.6 /
     1  47x,'C_3(^3Sig)=',1P,D15.7:/47x,'C_3(^1Pi) =',D15.7:/
     2  47x,'C_3(^3Pi) =',D15.7:/
     3  47x,'C_6(^3Sig)=',D15.7:/47x,'C_6(^1Pi) =',D15.7:/
     4  47x,'C_6(^3Pi) =',D15.7:/
     5  47x,'C_8(^3Sig)=',D15.7:/47x,'C_8(^1Pi) =',D15.7:/
     6  47x,'C_8(^3Pi) =',D15.7)
  672 FORMAT(' uLR(r) has ',I3,' inverse-power terms:',4x,'C',I1,
     1  ' =',1PD16.8:/40x,'C',i1,' =',D16.8:/(40x,'C',i2,'=',D16.8:))
  674 FORMAT(5x,'Generate   betaINF=',f16.12,'  from uLR(Re)=',1PD17.10)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

