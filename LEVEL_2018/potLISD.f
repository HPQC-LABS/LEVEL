c***********************************************************************
c  Special version for Lisdat potentials for node-defined levl
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
c  to NBOB+1 terms.  ||    ****** last updated  06 June 2010***********
c-----------------------------------------------------------------------
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
     6 RM3,RET,RETsig,RETpi,RETp,RETm,BFCT,PPOW,PVSR,
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
co



      DO  I= 1, NPP

          VDISP= 0.d0
          DO  n= 6,12,2
              FN= (1.d0 - dexp(-0.91286d0*(XO(I) - n*0.37595d0)))**n
              VDISP= VDISP + FN*ZZ


