c***********************************************************************
      SUBROUTINE POTGEN_CH+XA(ISTATE,IMN1,IMN2,NPP,VLIM,XO,RM2,VV)
c***********************************************************************
c     SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
c    1                                                        NCN,CNN)
c***********************************************************************
c                       Version of  10 December 2015
c** Subroutine to generate the effective adiabatic potential energy 
c  and centrifugal potential functions for the X and A states of CH+
c  reported by Cho and Le Roy [JCP 144, xxxx (2016)]. 
c----------------------------------------------------------------------
c** Note that this is a special 'standalone' version of subroutine 
c  POTGEN from the general-purpose bound-state Franck-Condon code LEVEL.  
c  Replacing the calling statement in line 2 (above) with the (currently)
c  commented-out lines 4 & 5 allows this routine to be compiled and used
c  in the normal way with R.J. Le Roy's program LEVEL, AFTER (!!) one 
c  inserts a value for parameters ISTATE @ line 90 **** (!!)
c==============
c*** On INPUT:
c==============
c*  integer  ISTATE  specifies the choice of PEF/electronic state (above)
c*  integers IMN1 and IMN2  are the mass numbers for the isotopologue
c*  integer  NPP  is the dimension of the  REAL*8 radial distance and 
c                 potential energy function arrays
c*  VLIM  is the (externally specified) asymptote energy (in cm-1) for
c         the {6,6}Li_2 isotopologue potential of the chosen state
c*  XO(i)  is the array of  NPP  radial distances (in Angst.) at which
c          the potential energy function is to be calculated
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
c** NCN is the power of asymptotically dominant inverse power term in 
c   the long-range potential tail
c** Born-Oppenheimer correction functions in IPOTL=3 option may have up
c  to NBOB+1 terms.  ||  Allow for up to  NbetaMX  EMO/MLR coeffocients 
c-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NBOB,NbetaMX,ISMX
      PARAMETER (NBOB=10,NbetaMX=20,ISMX=2)
      INTEGER  I,J,M,IBOB,IAN1,IAN2,IMN1,IMN2,MN1R,MN2R,IORD,IORDD,
     1 IPOTL,IMIN,PAD,QAD,QNA,NU1,NU2,NT1,NT2,NCMAX,PPAR,QPAR,NCN,Nbeta,
     2 NSR,NVARB,NPP,LNPT,GNS,GEL,NCMM,MCMM,sVSR2,LVSR,IDSTT,MM1,MMLR(9)
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,AT,BT,BETA,BINF,B1,B2,CSAV,U1INF,U2INF,
     1  T1INF,T2INF,YPAD,YQAD,YQADSM,YQNA,YQNASM,ABUND,CNN,DSCM,DX,DX1,
     2  FCT,FC1,FC2,FG1,FG2,MASS1,MASS2,RMASS1,RMASS2,REQ,Rref,Rinn,
     3  Rout,SC1,SC2,SG1,SG2,VLIM,DVLIM,VMIN,XDF,X1,XS,XL,XP1,ZZ,ZP,ZQ,
     4  ZME,ULR,ULRe,rhoAB,REQP,DM(9),DMP(9),DMPP(9),CMM(9),CmEFF(9),T0,
     5  C6adj,C9adj,C11adj,RM3,BFCT2,PVSR,RH,Scalc,XXQ,REQq,RREFq,
     7  U1(0:NBOB),U2(0:NBOB),T1(0:NBOB),T2(0:NBOB),
     8  PARM(NbetaMX),XPARM(NbetaMX),rKL(NbetaMX,NbetaMX), XO(NPP),
     9  VV(NPP),RM2(NPP),bTT(-1:2),cDS(-2:0),bDS(-2:0)
      SAVE IBOB,IPOTL,PPAR,QPAR,PAD,QAD,QNA,Nbeta,MMLR,NVARB,NCMM
      SAVE DSCM,REQ,Rref,PARM,U1,U2,T1,T2,CSAV,BINF,ALFA,ZME,
     2 Rinn,Rout,ULR,ULRe,CMM,XPARM
c** Damping function parameters for use and printout .....
      DATA bTT/2.44d0,2.78d0,3.126d0,3.471d0/
      DATA bDS/3.3d0,3.69d0,3.95d0/
      DATA cDS/0.423d0,0.40d0,0.39d0/
      SAVE bTT, bDS, cDS
c** Electron mass, as per 2010 physical constants
      DATA ZME/5.4857990946d-4/
c-----------------------------------------------------------------------
ccc Define variables for the ISMX states of interest
      INTEGER ISTATE,PSEL(ISMX),PPARA(ISMX),QPARA(ISMX),NSRA(ISMX),
     1 NbetaA(ISMX),NCMMA(ISMX),sVSR2A(ISMX),IDSTTA(ISMX),MMLRA(9,ISMX),
     2 MN1RA(ISMX),MN2RA(ISMX),PADA(ISMX),QADA(ISMX),QNAA(ISMX),
     3 NU1A(ISMX),NU2A(ISMX),NT1A(ISMX),NT2A(ISMX)
      REAL*8 DSCMA(ISMX),REQA(ISMX),RrefA(ISMX),rhoABA(ISMX),
     1 CMMA(9,ISMX),PARMA(NbetaMX,ISMX),XPARMA(NbetaMX,ISMX),
     2 U1infA(ISMX),U2infA(ISMX),T1INFA(ISMX),T2infA(ISMX),
     3 U1A(0:NBOB,ISMX),U2A(0:NBOB,ISMX),T1A(0:NBOB,ISMX),
     4 T2A(0:NBOB,ISMX)
c** Begin by identifying the chemical species by the IAN1 & IAN2 values
      IAN1= 6
      IAN2= 1
c** when using this with the standard LEVEL program, select ISTATE case here
c*** ISTATE =1 for X(^1\Sigma^+) CH+
c           =2 for A(^1\Sigma^+) CH+ 
c-----------------------------------------------------------------------
      ISTATE= 1                     !! {for example ... for the X-state}
c-----------------------------------------------------------------------
ccc Store parameter values for the specific ISMX states of interest
      DATA PSEL/4,4/, PPARA/5,6/,QPARA/1,1/,NSRA/6,8/,NbetaA/6,8/,
     1  sVSR2A/-2,-2/, IDSTTA/1,1/, NCMMA/3,3/, MMLRA/4,6,8,6*0, 4,6,8,
     2  6*0/, MN1RA/12,12/, MN2RA/1,1/,PADA/9,9/,QADA/1,1/, QNAA/7,4/,
     3  NU1A/3,2/, Nu2A/6,6/, Nt1A/-1,-1/, Nt2A/1,2/
      DATA DSCMA/34362.8d0, 10303.7d0/, REQA/1.1284625d0, 1.235896d0/,
     3  RrefA/2.55d0,2.96d0/, rhoABA/1.19d0,1.19d0/
c** specify up to 9 C_m long-range coefficients for the ISMX states
     4  CMMA/3.87847236D+04, 4.3D+03, 1.6D+08, 6*0.d0, 3.87847236D+04,
     5      1.07D+05,3.3D+08,6*0.d0/, 
c** now specify up to 20 exponent parameters \beta_i for the ISMX states
     7 PARMA/-1.50240186D-01, 1.37297788D+01, 7.1225228D+01, 
     8 2.133232D+02, 3.8073D+02, 3.7529D+02, 1.579D+02, 13*0.d0, 
c... now for state-2 ...
     a  -7.2644835D-01, 1.9172913D+01, 1.1120663D+02, 3.904437D+02,
     b  1.004442D+03, 2.07555D+03, 3.172D+03, 2.9267D+03, 1.167D+03,
     b  11*0.D0/,
c... now up to 20 (dummy) values of XPARMX for each of the ISMX states
     e  XPARMA/20*1.d0, 20*1.d0/
c** Next the asympptotic limits for the BOB functions
     f U1infA/0.d0,0.d0/,U2infA/0.d0,0.d0/,T1infA/0.d0,0.d0/,
     g T2infA/0.d0,0.d0/,
c... and FINALLY up to 11 BOB expansion parameters for each of ISMX states
     h U1A/0.0D+0,  7.759D+01, 5.17D+02, 1.7D+04, 7*0.D0, -2.01D-01, 
     i  4.83D+01, 1.42D+02, 8*0.d0/,
     i U2A/0.0D+0, 2.288D+02, 1.7048D+03, 3.723D+04, 9.6D+04, 2.29D+05, 
     j   5.59D+06, 4*0.D0, 1.1989D+01, 4.379D+02, 2.4704D+03, 4.129D+04,
     k   3.44D+04, 5.71D+05, 2.8D+06, 4*0.d0/, 
     l T1A/11*0.d0, 11*0.d0/, T2A/-4.17D-3, -2.613D-3, 9*0.D0,
     m   0.0d+0, -5.6d-05, 1.727d-2, 8*0.d0/
C==========END OF CASE SPECIFIC PARAMETER DEFINITIONS===================
c** and now use 'general' potgen structure - replacing READs with DATA
      LNPT= 1
c
      IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
ccc       READ(5,*) IPOTL, PPAR, QPAR, NSR, Nbeta, IBOB
ccc       READ(5,*) DSCM, REQ, Rref
c-----------------------------------------------------------------------
          IPOTL= PSEL(ISTATE)
          PPAR= PPARA(ISTATE)
          QPAR= QPARA(ISTATE)          
          NSR= NSRA(ISTATE)          
          Nbeta= NbetaA(ISTATE)          
          DSCM= DSCMA(ISTATE)          
          REQ= REQA(ISTATE)          
          Rref= RrefA(ISTATE)          
          IBOB= 0
          IF((NU1A(ISTATE).GE.0).OR.(NU1A(ISTATE).GE.0).OR.
     1              (NT1A(ISTATE).GE.0).OR.(NT2A(ISTATE).GE.0)) IBOB=1
          IF(IPOTL.GE.4) THEN
c-----------------------------------------------------------------------
c** For MLR, DELR, HFD, Tang-Toennies or Tiemann-polynomial potentials .....
ccc           READ(5,*) NCMM, rhoAB, sVSR2, IDSTT
ccc           READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
c-----------------------------------------------------------------------
              rhoAB= rhoABA(ISTATE)
              sVSR2= sVSR2A(ISTATE)
              IDSTT= IDSTTA(ISTATE)
              NCMM= NCMMA(ISTATE)
              DO  I= 1,NCMM
                  MMLR(I)= MMLRA(I,ISTATE)
                  CMM(I)= CMMA(I,ISTATE)
                  ENDDO
              ENDIF
c-----------------------------------------------------------------------
          IF(IPOTL.EQ.1) NVARB= 0
          IF(IPOTL.EQ.2) THEN
              NVARB= Nbeta+2
              ENDIF
          IF(IPOTL.EQ.3) THEN
              NVARB= Nbeta+1
              IF(PPAR.LE.0) NVARB=2
              ENDIF
          IF(IPOTL.EQ.4) THEN
              NVARB= Nbeta+ 1
              IF(NSR.LT.0) NVARB= Nbeta - NSR
              ENDIF
          IF(IPOTL.EQ.5) THEN
              IORD= Nbeta
              NVARB= Nbeta+ 1
              ENDIF
          IF(IPOTL.EQ.6) NVARB= 5
          IF(IPOTL.EQ.7) THEN
              NVARB= 2
              IF(QPAR.GT.2) NVARB= QPAR
              ENDIF
          IF(IPOTL.EQ.8) NVARB= Nbeta+ 4
c-----------------------------------------------------------------------
          IF(NVARB.GT.0) THEN
ccc           IF((IPOTL.EQ.4).AND.(NSR.LT.0)) THEN
ccc               READ(5,*) (XPARM(I),PARM(I), I=1,NVARB)
ccc             ELSE
ccc               READ(5,*) (PARM(I), I=1,NVARB)
ccc             ENDIF
c-----------------------------------------------------------------------
              DO  I= 1,NVARB
                  PARM(I)= PARMA(I,ISTATE)
                  XPARM(I)= XPARMA(I,ISTATE)
                  ENDDO
              ENDIF 
c-----------------------------------------------------------------------
          IF(IBOB.GT.0) THEN
c-----------------------------------------------------------------------
ccc           READ(5,*) MN1R, MN2R, PAD, QAD, NU1, NU2, QNA, NT1, NT2
c-----------------------------------------------------------------------
              MN1R= MN1RA(ISTATE)
              MN2R= MN2RA(ISTATE)
              PAD= PADA(ISTATE)
              QAD= QADA(ISTATE)
              NU1= NU1A(ISTATE)
              NU2= NU2A(ISTATE)
              QNA= QNAA(ISTATE)
              NT1= NT1A(ISTATE)
              NT2= NT2A(ISTATE)
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
ccc                   READ(5,*) U1INF,(U1(I), I=0,NU1)
c-----------------------------------------------------------------------
                      DO I=0,NU1
                        U1(I)= U1A(I,ISTATE)
                        ENDDO
                      U1inf= U1infA(ISTATE)
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
ccc                   READ(5,*) U2INF,(U2(I), I=0,NU2)
c-----------------------------------------------------------------------
                      DO I=0,NU2
                        U2(I)= U2A(I,ISTATE)
                        ENDDO
                      U2inf= U2infA(ISTATE)
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
ccc                   READ(5,*) T1INF,(T1(I), I=0,NT1)
c-----------------------------------------------------------------------
                      DO I=0,NT1
                        T1(I)= U1A(I,ISTATE)
                        ENDDO
                      T1inf= T1infA(ISTATE)
c-----------------------------------------------------------------------
                      WRITE(6,634) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1 1,T1INF,QNA,QNA,QNA,QNA,QNA,QNA,NT1,QNA,NT1+1,(T1(I),I= 0,NT1)
                      FG1= RMASS1/MASS1
                      ENDIF
c
                  IF(NT2.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-2 ...
c-----------------------------------------------------------------------
ccc                   READ(5,*) T2INF,(T2(I), I=0,NT2)
c-----------------------------------------------------------------------
                      DO I=0,NT2
                        T2(I)= T2A(I,ISTATE)
                        ENDDO
                      T2inf= T2infA(ISTATE)
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
          ENDIF
c
c=======================================================================
c** Generate a  Lennard-Jones(PPAR,QPAR)  potential here.
c=======================================================================
      IF(IPOTL.EQ.1) THEN 
          XS= PPAR
          XL= QPAR
          XDF= DSCM/(XS-XL)
          IF(LNPT.GT.0) WRITE(6,600) PPAR,QPAR,DSCM,REQ
          CNN= XS*XDF*REQ**QPAR
          NCN= QPAR
          DO  I= 1,NPP
              VV(I)= (XL*(REQ/XO(I))**PPAR - XS*(REQ/XO(I))**QPAR)*XDF
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
              VLIM= VV(NPP)
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
c** for a new case ... define ULRE an print potential description
              NCN= MMLR(1)
              CNN= CMM(1)
              ULRe= 0.d0
              IF((NCMM.GE.3).AND.(MMLR(1).LE.0)) THEN
c** Get uLR(Re) For special Aubert-Frecon Li2(^2P + ^2S) {3,0,6,6} type cases
                  RM3= 1.d0/REQ**3
c** QED retardation for Li2(A) deleted here  5 lines & RETp, RETm vbles removed
                  C6adj= CMM(3) + CMM(1)**2/(4.d0*DSCM)
                  C9adj = CMM(1)*C6adj/(2.d0*DSCM)
                  IF((MMLR(2).EQ.0).or.(MMLR(2).EQ.-2)) THEN
c ... for Aubert-Frecon 2x2 treatment of {C3,C6,C8} for Alkali A-state
                      T0= CMM(1)*RM3
                      ULRe= 0.5d0*(-CMM(2) + CMM(1)*RM3)
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
c.... ajdustment to get the upper root ...for the b-state
                      IF(MMLR(2).EQ.-2) ULRe= ULRe - T0
                      ENDIF
                  IF(MMLR(2).EQ.-1) THEN
c ... extension for Li2(c) {3,0,6,6,8,8} 3x3 Aubert-Frecon case ...
                      CALL AF3X3LEV(REQ,CMM(2),CMM(1),C6adj,
     1                                               CMM(4),DSCM,ULRe)
                      ULRe= ULRe + C9adj*RM3**3
                      ENDIF
                ELSE
c*** for 'ordinary' NCMM-term MLR uLR(r) ...  with damping [if rhoAB > 0]
c!! Initialize CmEFF on VERY first call from MAIN prog.
                  DO m= 1, NCMM
                      CmEFF(m)= CMM(m)
                      ENDDO
                  MCMM= NCMM
c=======================================================================
cc   ! As appropriate - make (& write) 'Dattani' Cm{adj} correctons here
          IF((MMLR(1).EQ.6).AND.(NCMM.GE.4)) THEN
c... First, consider C6/C12adj(C14adj) for MMLR(m)={6,8,10,(11),12,14} case
              IF(MMLR(4).EQ.12) THEN      ! explicitly MMLR(4)=12
                  CmEFF(4)= CMM(4) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN    ! assuming MMLR(5)=14
                      CmEFF(5)= CMM(5) + 0.5d0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                      ENDIF
                ELSE      !! Assuming explicitly MMLR(2)=11 & MMLR(5)=12
                  CmEFF(5)= CMM(5) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN             ! implicitly MMLR(6)=14
                      CmEFF(6)= CMM(6)+0.5D0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                ENDIF
              ENDIF
          IF((MMLR(1).EQ.5).AND.(NCMM.GE.4)) THEN
c... Then, consider C5/C10adj + C12adj for MMLR(m)={5,6,8,10,12,14} cases
              CmEFF(4)= CMM(4) + 0.25D0*CMM(1)**2/DSCM
              WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
              IF(NCMM.GE.5) THEN                   ! introduce C12^{adj}
                  CmEFF(5)= CMM(5) + 0.25D0*CMM(2)**2/DSCM
                  WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                  IF(NCMM.GE.6) THEN               ! introduce C14^{adj}
                      CmEFF(6)= CMM(6) + 0.5D0*CMM(2)*CMM(3)/DSCM
                      WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                      ENDIF
                  ENDIF
              ENDIF
          IF((MMLR(1).EQ.4).AND.(NCMM.GE.3)) THEN
              IF(MMLR(3).EQ.8) THEN
c... First, consider C4/C8adj + C12adj for MMLR(m)={4,6,8,10,12,14} cases
                  CmEFF(3)= CMM(3) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,712) MMLR(3),MMLR(3),CmEFF(3)
                  IF(NCMM.GE.4) THEN           ! implicitly MMLR(4)=10
                      CmEFF(4)= CMM(4) + 0.5D0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710) MMLR(4),MMLR(4),CmEFF(4)
                      IF(NCMM.GE.5) THEN       ! implicitly MMLR(5)=12
                          CmEFF(5)= CMM(5) + 0.5D0*CMM(1)*CMM(3)/DSCM
     1                                         + 0.25D0*CMM(2)**2/DSCM
                          WRITE(6,710) MMLR(5),MMLR(5),CmEFF(5)
                          IF(NCMM.GE.6) THEN     ! implicitly MMLR(6)=14
                              CmEFF(6)= CMM(6)+0.5D0*CMM(2)*CMM(3)/DSCM
     1                                      + 0.5D0*CMM(1)*CMM(4)/DSCM
                              WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                              ENDIF
                          ENDIF
                      ENDIF
                ELSEIF((MMLR(3).EQ.7).AND.(NCMM.GT.3)) THEN
c... Then, consider C4/C8adj + C12adj for MMLR(m)={4,6,7,8,10,12,14} cases
c...       allowing for a C7:   C8 is m=4, C10 is m=5 ... etc.
                  CmEFF(4)= CMM(4) + 0.25D0*CMM(1)**2/DSCM
                  WRITE(6,712) MMLR(4),MMLR(4),CmEFF(4)
                  IF(NCMM.GE.5) THEN             ! implicitly MMLR(5)=10
                      CmEFF(5)= CMM(5) + 0.5D0*CMM(1)*CMM(2)/DSCM
                      WRITE(6,710)MMLR(5),MMLR(5), CmEFF(5)
                      IF(NCMM.GE.6) THEN         ! implicitly MMLR(6)=12
                          CmEFF(6)= CMM(6) + 0.5D0*CMM(1)*CmEFF(4)/DSCM
     1                                         + 0.25D0*CMM(3)**2/DSCM
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                          IF(NCMM.GE.7) THEN     ! implicitly MMLR(7)=14
                              CmEFF(7)= CMM(7) + 0.25d0*CMM(3)**2/DSCM
     2           + 0.5D0*CMM(2)*CMM(4)/DSCM + 0.5D0*CMM(1)*CMM(5)/DSCM
                              WRITE(6,710) MMLR(7),MMLR(7),CmEFF(7)
                              ENDIF
                          ENDIF
                      ENDIF
                ENDIF
              ENDIF
          IF((MMLR(1).EQ.3).AND.(NCMM.GE.2)) THEN
c... Then, consider C3/C6adj & C9adj for MMLR(m)={3,6,8,(9),10,(11),12,14} cases
              CmEFF(2)= CMM(2) + 0.25D0*CMM(1)**2/DSCM
              WRITE(6,712) MMLR(2),MMLR(2),CmEFF(2)
              IF(NCMM.GE.3) THEN                !! introduce C9adj & MMLR=9
                  C9adj= 0.5d0*CMM(1)*CmEFF(2)/DSCM
                  IF((NCMM.EQ.3).OR.(MMLR(4).NE.9)) THEN
                      MCMM= NCMM+1             !! adding m=9 as last power
                      MMLR(MCMM)= 9
                      CMM(MCMM)= C9adj
                      WRITE(6,714) MMLR(MCMM),CmEFF(MCMM)
                    ELSE                        !! case of MMLR(4)=9) 
                      CmEFF(4)= CmEFF(4) + C9adj
                      WRITE(6,712) MMLR(4),MMLR(4),CmEFF(4)
                    ENDIF      !! finished for MMLR={3,6,8} or {3,6,8 ...}
                  ENDIF         !!  and for MMLR={3,6,8,9, ...}    
              IF(NCMM.GE.4) THEN               !! first - create C11adj
                  C11adj= 0.5d0*CMM(1)*CmEFF(3)/DSCM
                  IF((MMLR(4).EQ.10).OR.
     1           ((NCMM.GT.4).AND.(MMLR(4).EQ.9))) THEN
                      MCMM= NCMM+1          ! adding m=11 as last power
                      MMLR(MCMM)= 11
                      IF(NCMM.GE.6) THEN
                          CmEFF(6)= C11adj+ CMM(6)
                          WRITE(6,710) MMLR(6),MMLR(6),CmEFF(6)
                        ELSE
                          CmEFF(MCMM)= C11adj
                          WRITE(6,716) MMLR(MCMM),CmEFF(MCMM)
                        ENDIF
                      ENDIF
                  ENDIF
              ENDIF
c** End of  CmEFF= Cm + CmADJ  setup ===================================
  710 Format("  'Dattani adjustment' for MLR  C",I2,'  yields   internal
     1  C',  I2,'{eff}=',1PD14.7)
  712 Format("  'Dattani adjustment' for MLR  C",I1, '   yields internal
     1   C',  I1,'{eff}=',1PD15.8)
  714 Format("  'Dattani adjustment' for MLR(m1=3) yields internal   C",
     1  I1,'{eff}=',1PD15.8)
  716 Format("  'Dattani adjustment' for MLR(m1=3) yields internal   C",
     1  I2,'{eff}=',1PD15.8)
c** Now - initialize at r= REQ for 'simple' inverse-power sum case
                  IF(rhoAB.GT.0.d0) CALL dampF(REQ,rhoAB,MCMM,MMLR,
     1                               sVSR2,IDSTT,DM,DMP,DMPP)
                  DO  J= 1,MCMM
                      IF(rhoAB.LE.0.d0) THEN
                          ULRe= ULRe + CmEFF(J)/REQ**MMLR(J)
                        ELSE
                          ULRe= ULRe + DM(J)*CmEFF(J)/REQ**MMLR(J)
                        ENDIF
                      ENDDO
                ENDIF
              BINF= DLOG(2.d0*DSCM/ULRe)
c*** print for MLR form
              WRITE(6,602) PPAR,QPAR,DSCM,REQ
c... for Huang form: \beta(yp)= Binf*yp + [1-yp]*{power series in yq}
              IF(NSR.GT.0) WRITE(6,607) PPAR,PPAR,QPAR,Nbeta,Nbeta+1,
     1                                       (PARM(J),J= 1,Nbeta+1)
c... print for Asen Pashov Spline Exponent (NSR < 0) MLR form
              IF(NSR.LT.0) THEN 
                  WRITE(6,604) PPAR,Nbeta,(PARM(J),J= 1,Nbeta) 
                  WRITE(6,610) QPAR,Rref,(XPARM(J),J= 1,Nbeta)
c** Prepare Asen's Rlk array for later use in generating Spline fx. 
                  CALL Lkoef(Nbeta,XPARM,rKL,NbetaMX)
                  ENDIF
              IF(Rref.GT.0) THEN
                  WRITE(6,613) Rref
                ELSE
                  WRITE(6,615) REQ
                  Rref= REQ
                ENDIF  
              IF(rhoAB.GT.0.d0) THEN      !! describe type of Damping Fx
                  PVSR= 0.5d0*sVSR2
                  IF(IDSTT.GT.0) THEN   !! first option: Douketis-Scoles
                      PVSR= 0.5d0*sVSR2
                      WRITE(6,664) rhoAB,PVSR,bDS(sVSR2),cDS(sVSR2),PVSR
                    ELSE                 !! second option: Tang-Toennies
                      LVSR= sVSR2/2
                      WRITE(6,666) rhoAB,LVSR,bTT(LVSR)
                    ENDIF
                ELSE                     !!  ELSE ... for no damping ...
                  WRITE(6,668)
                ENDIF
              WRITE(6,617) BINF,MMLR(1),CmEFF(1),MMLR(1)
              IF(NCMM.GT.1) THEN
                  MM1= 2
                  IF(MMLR(2).LE.0) THEN
                      MM1= 3
                      IF(MMLR(2).EQ.0)
     1                         WRITE(6,623) MMLR(2),CMM(2),MMLR(2),0
                      IF(MMLR(2).EQ.-2)
     1                         WRITE(6,627) MMLR(2),CMM(2),MMLR(2),0
                      IF(MMLR(2).EQ.-1) WRITE(6,625) MMLR(2),CMM(2),0
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
              IF(NSR.GE.0) THEN
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
              IF((NCMM.GE.3).AND.(MMLR(1).LE.0)) THEN
c... For special Aubert-Frecon Li2(^2P + ^2S) {3,0,6,6} type case 
                  RM3= 1.d0/XO(I)**3
c** QED retardation for Li2(A) deleted here  5 lines & RETp, RETm vbles removed
                  IF((MMLR(1).EQ.0).OR.(MMLR(1).EQ.-2)) THEN
c... Aubert-Frecon 2x2 case - for Li2(A) or Li2(b)
                      T0= CMM(1)*RM3
                      ULR= 0.5d0*(-CMM(2) + CMM(1)*RM3)
                      IF(NCMM.GE.3) THEN
                          T0= T0 + C6adj*RM3**2
                          ULR= ULR + 0.5d0*C6adj*RM3**2
                          IF(NCMM.GE.4) THEN
                              T0= T0+ CMM(4)*(RM3/XO(I))**2
                              ULR= ULR + 0.5d0*CMM(4)*(RM3/XO(I))**2
     1                                                  + C9adj*RM3**3
                              ENDIF
                          ENDIF
c.....  adjustment for the b-state .........
                      IF(MMLR(1).EQ.-2) ULR= ULR- T0
                      T0= T0/3.d0
                      ULR= ULR+ 0.5d0*DSQRT((T0- CMM(2))**2+ 8.d0*T0**2)
                      ENDIF
                  IF(MMLR(1).EQ.-1) THEN
c... for Aubert-Frecon 3x3 case yielding lowest (c state) root
                      CALL AF3X3LEV(XO(I),CMM(2),CMM(1),C6adj,
     1                                                CMM(4),DSCM,ULR)
                      ULR= ULR + C9adj*RM3**3
                      ENDIF
                ELSE
c*** for 'ordinary' NCMM-term MLR uLR(r) ...  with damping [if rhoAB > 0]
                  CALL dampF(XO(I),rhoAB,NCMM,MMLR,
     1                                         sVSR2,IDSTT,DM,DMP,DMPP)
                  DO  J= 1,MCMM
                      IF(rhoAB.LE.0.d0) THEN
                          ULR= ULR + CmEFF(J)/XO(I)**MMLR(J)
                        ELSE
                          ULR= ULR + DM(J)*CmEFF(J)/XO(I)**MMLR(J)
                        ENDIF
                      ENDDO
                ENDIF
              BETA= (ULR/ULRe)*DEXP(-BETA*ZZ)
              VV(I)= DSCM*(1.d0 - BETA)**2 - DSCM + VLIM
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
              CALL dampF(REQ,rhoAB,NCMM,MMLR,sVSR2,IDSTT,DM,DMP,DMPP)
              DO  J= 1,NCMM
                  T0= CMM(J)/REQ**MMLR(J)
                  ULRe= ULRe+ T0*DM(J)
                  B1= B1+ T0*(DMP(J) - DM(J)*MMLR(J)/REQ)
                  ENDDO
              A1= DSCM - ULRe - B1/BETA
              B1= 2.d0*A1 + B1/BETA
              WRITE(6,650) QPAR,DSCM,REQ,Nbeta,(PARM(I),I= 1,IORD+1)
              WRITE(6,652) QPAR,QPAR,QPAR,QPAR,QPAR
              IF(Rref.GT.0.d0) WRITE(6,654) Rref
              IF(Rref.LE.0.d0) WRITE(6,656) REQ
              WRITE(6,658) A1,B1,NCMM,(MMLR(J),CMM(J),J= 1,NCMM)
              IF(IDSTT.GT.0) THEN
                  PVSR= 0.5d0*sVSR2
                  WRITE(6,664) rhoAB,PVSR,bDS(sVSR2),cDS(sVSR2),PVSR
                ELSE
                  LVSR= sVSR2/2
                  WRITE(6,666) rhoAB,LVSR,bTT(LVSR)
                ENDIF
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
              ULR= 0.0d0
              CALL dampF(XO(I),rhoAB,NCMM,MMLR,sVSR2,IDSTT,DM,DMP,DMPP)
              DO  J= 1, NCMM
                  ULR= ULR+ DM(J)*CMM(J)/XO(I)**MMLR(J)
                  ENDDO
              VV(I)=  (A1*BETA - B1)*BETA - ULR + VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.6) THEN
c=======================================================================
c** For generalized  HFDB(m= MMLR(j), j=1,NCMM) potential with reduced 
c  form   VBAR = ALFA*x**PARM(5) * exp[-BETR*x - PARM(4)*x**2] - D(x)*
c      [CMM(1)/x**MMLR(1) + CMM(2)/x**sMMLR(2) + CMM(3)/x**MMLR(3) + ...
c  x=r/R_e ,  VBAR= V/De   and    D(x) = 1 for x > PARM(2)   and
c      D(x)= exp[-PARM(1)*(PARM(2)/x - 1)**PARM(3)] for  x < PARM(2)
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= CMM(1)*DSCM*REQ**MMLR(1)
              A1= PARM(1)
              A2= PARM(2)
              A3= PARM(3)
              B2= PARM(4)
              DX= 1.d0
              DX1= 0.d0
              IF(A2.GT.1.d0) THEN
                  DX= DEXP(-A1*(A2- 1.d0)**A3)
                  DX1= A1*A2*A3*DX*(A2- 1.d0)**(A3- 1.d0)
                  ENDIF
              ALFA= 0.d0
              DO  J= 1, NCMM
                  ALFA= ALFA+ CMM(J)
                  ENDDO
              ALFA= ALFA*DX -1.D0
              IF(ALFA.LE.0.d0) THEN
                  WRITE(6,622) ALFA,(MMLR(J),CMM(J),J= 1, NCMM)
                  STOP
                  ENDIF
              B1= 0.d0
              DO  J= 1, NCMM
                  B1= B1+ (MMLR(J)*DX - DX1)*CMM(J)
                  ENDDO
              B1= B1/ALFA + PARM(5) - 2.d0*B2
              ALFA= ALFA*DEXP(B1+B2)
              WRITE(6,624) PARM(5),B1,B2,ALFA*DSCM,
     1                                     (MMLR(J),CMM(J),J= 1, NCMM)
              WRITE(6,626) DSCM,REQ,A1,A2,A3
              ENDIF
          DO  I= 1,NPP
              X1= XO(I)/REQ
              XP1= 0.0D0
              IF((B1*X1+ B2*X1**2).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*X1**PARM(5)
              FC1= 0.d0
              DO  J= 1, NCMM
                  FC1= FC1 + CMM(J)/X1**MMLR(J)
                  ENDDO
              IF(X1.LT.A2) FC1= FC1*DEXP(-A1*(A2/X1- 1.d0)**A3)
              VV(I)= DSCM*(ALFA*XP1- FC1) + VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.7) THEN
c=======================================================================
c** Generate Tang-Toennies (TT) type potential as per JCP 118, 11 (2003) 
c   or generalized Bich-type Tang-Toennies fx.
c NCMM = number of inverse-power long-range terms and NVARB = 2.
c Here DSCM and Re are dummy parameters [equal to 1.d0]. The powers and
c coefficients of the NCMM inverse-power long-range terms are
c MMCM(j) and CMM(j).
c=======================================================================
          NCN= MMLR(1)
          CNN= CMM(1)
          AT= PARM(1)
          BT= PARM(2)
          IDSTT= 0
          sVSR2= 2
c** Define  rhoAB for consistency with conventional TT(sVSR2=+2) damping fx.
          rhoAB= BT/3.126D0
          IF((QPAR.LE.2).AND.(AT.LT.0.d0)) THEN
c** For conventonal TT potential, if input  AT  not specified, generate
c   AT & BT from REQ & DSCM 
              CALL dampF(REQ,rhoAB,NCMM,MMLR,sVSR2,IDSTT,DM,DMP,DMPP)
              A2= 0.d0
              A3= 0.d0
              DO J=1, NCMM
                 A2= A2+ DM(J)*CMM(J)/REQ**MMLR(J)
                 A3= A3+(CMM(J)/REQ**MMLR(J))*(DMP(J)-MMLR(J)*DM(J)/REQ)
                 ENDDO
              BT= A3/(DSCM-A2)
              AT= DEXP(BT*REQ)*A3/BT
              PARM(1)= AT
              PARM(2)= BT
              ENDIF
          VMIN= 9.d99
          IMIN=1
          DO I= 1, NPP
c....generate potential function array
              CALL dampF(XO(I),rhoAB,NCMM,MMLR,sVSR2,IDSTT,DM,DMP,DMPP)
c....calculate the (damped) long range tail
              A3= 0.d0
              DO J= 1, NCMM
                  A3= A3+ DM(J)*CMM(J)/XO(I)**MMLR(J)
                  ENDDO
              IF(QPAR.LE.2) THEN
c....For 'conventional' TT Born-Meyer repulsion
                  XP1= -BT*XO(i)
                ELSE
c....For Bich/Vogel modified TT model
                  XP1= PARM(3)*XO(I)+ PARM(4)*XO(I)**2+ PARM(5)/XO(I)
     1                                              + PARM(6)/XO(I)**2
                ENDIF
              VV(I)= AT*DEXP(XP1) - A3 + VLIM
              IF(VV(I).LE.VMIN) THEN
c... search for potential minimum ...
                  VMIN= VV(I)
                  IMIN= I
                  ENDIF
              ENDDO
c*** Use quadratic approximation to determine REQ and DSCM
          IF(IMIN.EQ.1) IMIN=2
          A1= VV(IMIN-1)
          A2= VV(IMIN)
          A3= VV(IMIN+1)
          RH= XO(IMIN) - XO(IMIN-1)
          B1= (A3- 2.d0*A2 + A1)/(2.d0*RH**2)
          REQ= XO(IMIN) + 0.5d0*RH - (A3-A2)/(2.d0*RH*B1)
c...... If ?????
cc        REQ= XO(IMIN) + (RH/2.d0)*(A1-A3)/(A1-2.d0*A2+A3)
cc        B1= (A1-A2)/((2.d0*REQ- 2.d0*XO(IMIN)+ RH)*RH)
c................
          A2= A2- B1*(XO(IMIN)-REQ)**2
          DSCM= VLIM - A2
          IF(QPAR.LE.2) WRITE(6,628) PARM(1),PARM(2),REQ, DSCM,
     1                                     (MMLR(J),CMM(J),J= 1, NCMM)
          IF(QPAR.GT.2) WRITE(6,629) PARM(1),PARM(2),REQ, DSCM,PARM(3),
     1             PARM(4),PARM(5),PARM(6),(MMLR(J),CMM(J),J= 1, NCMM)
          ENDIF
c
      IF(IPOTL.EQ.8) THEN
c=======================================================================
c** Generate Tiemann-type polynomial potential attached to inverse-power
c  tail and 1/R^{12} (or exponential) inner wall [PRA 63, 012710 (2000)].
c  Polynomial expansion variable is  z= [R - Rm]/[R + b*Rm] where 
c  expansion has constant and linear terms.  The read-in DSCM= De (well
c  depth), but  Rm (read in as REQ) is not precisely Re (for a1 .neq. 0).
c  NCMM= number of inverse-power long-range terms;  
c  NVARB= (polynomial order) + 4.  [PPAR and NSR are dummy parameters]
c** Read-in parameters PARM(i) are in order: the  (Nbeta+1)  polynomial
c  coefficients  a(0) - a(Nbeta), the expansion variable denominator
c  factor b=PARM(Nbeta+2), and the the inner and outer bounds on the 
c  polynomial domain, Tiemann's Rinn= PARM(Nbeta+3) & Rout= PARM(Nbeta+4), 
c  respectively.  The powers and coefficients (-ve if attractive) of the
c  NCMM inverse-power long-range terms are MMCM(j) and CMM(j).
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= -CMM(1)
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
                  A3= A3+ CMM(J)/Rout**MMLR(J)
                  ENDDO
              PPAR= NCMM+ 1
              MMLR(PPAR)= MMLR(NCMM)+ 2
              CMM(PPAR)= (B1-A3)*Rout**MMLR(PPAR)
c*** Print for Tiemann-type potential
              IF(LNPT.GE.0) THEN
                  WRITE(6,640) DSCM,REQ,PARM(Nbeta+2),Nbeta,Nbeta+1, 
     1                                            (PARM(J),J= 1,Nbeta+1)
ccc               IF(XO(1).LT.Rinn) WRITE(6,642) PARM(Nbeta+3),A1,A2,A0
                  IF(XO(1).LT.Rinn) WRITE(6,642) PARM(Nbeta+3),A1,A2
                  IF(XO(NPP).GT.Rout) WRITE(6,644) PARM(Nbeta+4),
     1                                     (CMM(J),MMLR(J),J= 1, PPAR)
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
                      A3= A3+ CMM(J)/XO(I)**MMLR(J)
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
          IF((IAN1.EQ.3).AND.(IAN2.EQ.3).AND.(NCMM.GE.3).AND.
     1                                            (MMLR(2).LE.0)) THEN
c!! For mixed isotopopogue {6,7}Li_2(A) state, shift asymptote!
              IF((IMN1.NE.IMN2).AND.(MMLR(2).EQ.0)) THEN
                  DO  I= 1,NPP
                      RM3= (2.d0/3.d0)*CMM(1)/XO(I)**3
                      VV(I)= VV(I)+ RM3- DSQRT(RM3**2+ 3.085959756d-02)
                      ENDDO
                  VLIM= VLIM + DSQRT(3.085959756d-02)
                  ENDIF
c** For special case of A and c states of Li2, add BOB centrifugal term
              BFCT2= 2.d0*16.857629206d0*(MASS1+MASS2)/(MASS1*MASS2)
              DO  I= 1, NPP
                  VV(I)= VV(I) + BFCT2/XO(I)**2
                  ENDDO
              ENDIF
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  602 FORMAT(/' MLR(p=',I1,', q=',I1,') Potential with:   De='
     1 ,F10.4,'[cm-1]    Re=',F12.8,'[A]')
  604 FORMAT('   with exponent coefficient   beta(r)= y',I1,'^{eq} *{Spl
     1ine through the',I3,' function values: beta_i ='/(10x,4D16.8:))
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F11.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  607 FORMAT('   with exponent coefficient   beta(r)= beta{INF}*y',I1,
     1  ' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x,'exponent coefft. powe
     2r series order',I3/6x,'and',i3,' coefficients:',1PD17.9,2D17.9:/
     3  (10x,4D17.9:))
  608 FORMAT(/' EMO_',i1,' Potential with   De=',F11.4,'    Re=',F11.8,
     1 '   Rref=',F11.8/3x,'Exponent coeft: order-',i2,
     2 ' power series in  y=(r**',i1,' - Rref**',i1,')/(r**',i1,
     3 ' + Rref**',i1,')'/'   with',I3,' coefficients:',1x,1PD17.9,
     4 2D18.9:/(7X,4D17.9:))
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
  613 FORMAT(6x,'with radial variables  y_p and/or y_q  defined w.r.t.',
     1 '  Rref=',F10.7)
  615 FORMAT(6x,'radial variables  y_p and/or y_q  defined w.r.t.',
     1 '  Rref= Re=' F10.7)
  617 FORMAT('      betaINF=',f16.12,'  & uLR defined by  C',i1,' =',
     1  1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  619 FORMAT(50x,'C',I1,' =',1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  621 FORMAT(50x,'C',I2,'=',1PD13.6,'[cm-1 Ang','^',0P,I2,']')
  623 FORMAT(3x,'Use LOWER root of Aubert-Frecon 2x2 model for uLR(r) wi
     1th',4x,'C',I1,' =',1PD13.6,'[cm-1 Ang^',0P,I1,']'/8x,'including re
     2tardation & B(r) within V_{ad}, plus Delta(V)_{gu}^{(6,7)}')
  627 FORMAT(3x,'Use UPPER root of Aubert-Frecon 2x2 model for uLR(r) wi
     1th',4x,'C',I1,' =',1PD13.6,'[cm-1 Ang^',0P,I1,']'/8x,'including re
     2tardation & B(r) within V_{ad}, plus Delta(V)_{gu}^{(6,7)}')
  625 FORMAT(3x,'Use Aubert-Frecon 3x3 model for uLR(r) with',
     1  3x,'C',I2,' =',1PD13.6,'[cm-1 Ang^',0P,I1,']'/6x,'including reta
     2rtion & B(r) within V_{ad}')
  622 FORMAT(/' *** ERROR in generating HFD potential *** generate   ALF
     1A=',1PD15.7,'  from reduced  Cm  coefficients:'/(3x,3('   C',I2,
     2 '=',D15.7:)) )
  624 FORMAT(/' Potential is Generalized HFD with exponent factors   gam
     1ma=',f9.6/'   beta1=',f12.8,'   beta2=',f9.6,'   A=',1PD16.9, 
     2 " & reduced Cm's:"/(3x,3('   C',I2,' =',D15.7:)) )
  626 FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]   and'/
     1  '     Damping function  D(r)= exp[ -',0P,f6.4,'*(',f7.4, 
     2  '/X -1.0)**',f5.2,']')
  628 FORMAT(/' Potential is conventional Tang-Tonnies potential with
     1A=',1PD16.8/'   b=',0Pf11.8,'   Re=',f12.8,'   DSCM=',f12.4/ 
     2 " Cm's[cm-1 Ang^m]:"/(3x,3('   C',I2,'=',1PD15.8)))
  629 FORMAT(/' Potential is Bich modified Tang-Tonnies potential with',
     1 ' A=',1PD16.8/'   b=',0P,f11.8,'   Re=',f12.8,'   DSCM=',f12.4/
     2  '   and exponent:',SP1PD17.9,'*R ',D17.9,'*R^2 ',D17.9,'/R '/
     3  58x,D17.9,S,'/R^2' /
     4 "   Cm's[cm-1 Ang^m]:"/(3x,3('   C',0P,i2,'=',1PD15.8)))
  630 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1 '  of mass ',f15.11/'   consists of mass factor  [1- MASS(',I3,
     2 A2,')/MASS(',I3,A2,')]  multiplying all of:'/5x,'u',I1,'INF=',
     3 f11.6,'  times  y',i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1, 
     4 ' + Re**',i1,')]'/4x,'plus  [1 - y',i1,']  times an order',I3, 
     5 ' polynomial in'/7x,'y',i1,'=[(r**',i1,' - Re**',i1,')/(r**',i1, 
     6 ' + Re**',i1,')]  with the ',i3,' coefficients:'/(3x,4D17.9:))
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
     4 in z{O-T} with coefficients:'/ (3x,4D17.9:))
  640 FORMAT(/' Tiemann-type potential with   De=',F11.4,'   Rm=',f9.6,
     1 '   is a power series'/10x,'in  (r - Re)/(r ',SP,F9.5, 
     2 '*Re) of order',SS,I3,'  with the',I3,' coefficients:'/(5D16.8))
c 642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',1PD13.6,'*exp[-',
c    1 0P,F9.6,'*(r - Rinn)] ',SP,F10.3)
  642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',SP,F12.4,1x,1PD13.6,
     1  '/R**12' )
  644 FORMAT('  and  for  r > Rout=',F7.3,'   V= VLIM ',
     1 (SP,1PD14.6,'/r**',SS,I2):/(39x,SP,1PD14.6,'/r**',SS,I2))
  650 FORMAT(/' DELR(p=',i2,') Potential with   De=', F11.4,'[cm-1]   Re
     1=',F11.8,'[A]   where'/3x,'exponent coefft. has power series order
     2',I4/6x,'with polynomial coefficients',8x,1PD17.8,D17.8/ 
     3 (8x,4D17.8))
  652 FORMAT(6x,'where the radial variable   y_',I1,'= (r**',I1,' - Rref
     4**',i1,')/(r**',I1,' + Rref**',i1, ')')
  654 FORMAT(10x,'is defined w.r.t.   Rref=',F11.8) 
  656 FORMAT(10x,'is defined w.r.t.   Rref= Re= ',F11.8) 
  658 FORMAT(3x,'Generate A(DELR)=',1Pd17.9,'   B(DELR)=',D17.9/
     1 6x,'from uLR defined by',I2," inverse-power terms with coeffts (+
     2 've repulsive):"/(5x,3(5x,'C',0P,i2,' =',1Pd14.6:)))
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
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,sVSR2,IDSTT,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. r of the chosen form of the damping
c    function, for  m= 1 to MMAX.
c---------------------- RJL Version of 30 January 2014 ----------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!) 
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* sVSR2 defines damping s.th.  Dm(r)/r^m --> r^{sVSR2/2} as r --> 0
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form 
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term 
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m): first derivative of the damping function  DM(m) w.r.t. r
c  DMPP(m): second derivative of the damping function  DM(m) w.r.t. r
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),sVSR2,IDSTT,sVSR2,FIRST, Lsr,m,
     1  MM,MMAX,MMTEMP
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:4),bDS(-4:4),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-4:4), cpm(20,-4:4),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.13d0,3.47d0/
       DATA bDS/2.50d0,2.90d0,3.30d0,3.69d0,3.95d0,0.d0,4.53d0,0.d0,
     1          4.99d0/         
       DATA cDS/0.468d0,0.446d0,0.423d0,0.405d0,0.390d0,0.d0,0.360d0,
     1            0.d0,0.340d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
      MMTEMP = MMLR(1)
      IF(MMLR(1).LE.0) MMLR(1) = 1

      IF(RHOab.LE.0) THEN
          DO  m=1,NCMMax
              DM(m)=1.d0
              DMP(m)= 0.d0
              DMPP(m)= 0.d0
              ENDDO 
          WRITE(6,602) RHOab
          RETURN
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          Lsr= sVSR2/2
          IF((sVSR2.LT.-4).OR.(sVSR2.GT.4).OR.((2*LSR).NE.sVSR2)) THEN
                WRITE(6,600) 'TT',sVSR2
                STOP
                ENDIF
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
                  DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                  DMPP(m)= -aTT*aTT*XP*(SM(MM) 
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
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
                  DMP(m)= aTT*XP*SM(m-1)
                  DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((sVSR2.LT.-4).OR.(sVSR2.GT.0).OR.(sVSR2.EQ.1).OR.
     1                                              (sVSR2.EQ.3)) THEN
              WRITE(6,600) 'DS',sVSR2
              STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  sVSR2= -4,0
                      bpm(m,sVSR2)= bDS(sVSR2)/DFLOAT(m)
                      cpm(m,sVSR2)= cDS(sVSR2)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0 
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,sVSR2) + cpm(MM,sVSR2)*br)*br)
              YP= 1.d0 - XP
              ZK= MM + 0.5d0*sVSR2
              DM(m)= YP**ZK
              TK= (bpm(MM,sVSR2) + 2.d0*cpm(MM,sVSR2)*br)*rhoAB
              DMP(m) = ZK*XP*TK*DM(m)/YP
c ... calculate second derivative [for DELR case] {check this!}
              DMPP(m)= (ZK-1.d0)*XP*TK*DMP(m)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,sVSR2)*rhoAB**2/TK
              ENDDO   
          ENDIF  
      MMLR(1) = MMTEMP
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  ',A2,'-damping functions not yet de
     1fined for   sVSR2=',i3)
  602 FORMAT( /,' ***ERROR ***  should not call dampF when rhoAB=',
     1  F7.4)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
