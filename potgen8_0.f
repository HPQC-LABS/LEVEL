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
c  to NBOB+1 terms.  ||    ****** updated  01 November 2010 ******
c-----------------------------------------------------------------------
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,IBOB,IAN1,IAN2,IMN1,IMN2,MN1R,MN2R,IORD,IPOTL,
     1  PAD,MAD,PNA,NU1,NU2,NT1,NT2,NCMAX,MPAR,pPOW,NCN,NSR,NLR,
     2  NVARB,NPP,LNPT,GNS,GEL, NCMM,IDF, Jinn, MMLR(9)
      CHARACTER*2 NAME1,NAME2
      CHARACTER*3 MLname
      REAL*8  A0,A1,A2,A3,ALFA,Asw,Rsw,BETA,BINF,B1,B2,BT,CSAV,U1INF,
     1 U2INF,T1INF,T2INF,YPAD,YMAD1,YMAD2,YMNA1,YMNA2,YPNA,ABUND,CNN,
     2 DSCM,DX,DX1,FCT,FC1,FC2,FG1,FG2,MASS1,MASS2,RMASS1,RMASS2,REQ,
     3 Rinn,Rout,SC1,SC2,SG1,SG2,VLIM,VMIN,XDF,X1,XS,XL,XP1,ZZ,ZP,ZME,
     4 Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,fad1e,fad2e,
     5 ULR,ULRe,RHOd,REQP,Dm,Dmp,Dmpp,CMM(9),
     6 U1(0:NBOB),U2(0:NBOB),T1(0:NBOB),T2(0:NBOB),PARM(50),
     7 XO(NPP),VV(NPP),RM2(NPP)
      SAVE IBOB,IPOTL,IORD,MPAR,pPOW,PAD,MAD,PNA,NSR,NLR,MMLR,
     1 NVARB,NCMM
      SAVE DSCM,REQ,PARM,U1,U2,T1,T2,CSAV,BINF,ALFA,Rsw,ZME,
     2 Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,Rinn,Rout,ULR,ULRe,CMM
c** Electron mass, as per 2002 physical constants
      DATA ZME/5.4857990945d-4/
c
      IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
          READ(5,*) IPOTL, MPAR, NSR, NCMM, NVARB, IBOB, DSCM, REQ
          IF(IPOTL.GE.4) READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
          IF(IPOTL.EQ.1) NVARB= 0
          IF(IPOTL.EQ.6) NVARB= 5
          IF((IPOTL.EQ.3).AND.(MPAR.EQ.-1)) NVARB= 2
          IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
          IF(IBOB.GT.0) THEN
              READ(5,*) MN1R, MN2R, PAD, MAD, NU1, NU2, PNA, NT1, NT2
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
                  Aad1= 0.d0
                  Aad2= 0.d0
                  Ana1= 0.d0
                  Ana1= 0.d0
                  Rad1= 0.d0
                  Rad2= 0.d0
                  Rna1= 0.d0
                  Rna2= 0.d0
c=======================================================================
c** Read actual BOB polynomial expansion coefficients
c=======================================================================
                  IF(NU1.GE.0) THEN
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
                          READ(5,*) U1INF,(U1(I), I=0,NU1)
                        ELSE
                          READ(5,*) U1INF,(U1(I),I=0,NU1),Aad1,Rad1
                        ENDIF
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
c ... for Huang/Le Roy form for the adiabatic potential BOB radial fx.
                          WRITE(6,630) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1        1,U1INF,MAD,MAD,MAD,MAD,MAD,MAD,NU1,PAD,PAD,PAD,PAD,PAD,
     2                                         NU1+1,(U1(I),I= 0,NU1)
                          FC1= 1.d0 - RMASS1/MASS1
                        ELSE
c ... for Coxon/Hajigeorgiou switching fx. & mass scaling fx. for atom-1
                          WRITE(6,632) 1,MASS1,IMN1,NAME1,MN1R,NAME1,
     1                                   1,U1INF,NU1,(U1(I),I= 0,NU1)
                          WRITE(6,636) Aad1, Rad1
                          FC1= ZME*(RMASS1 - MASS1)/(RMASS1*MASS1)
                        ENDIF
                      ENDIF
c
                  IF(NU2.GE.0) THEN
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
                          READ(5,*) U2INF,(U2(I), I=0,NU2)
                        ELSE
                          READ(5,*) U2INF,(U2(I),I=0,NU2),Aad2,Rad2
                        ENDIF
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
c ... for Huang/Le Roy adiabatic radial strength function for atom-2
                          WRITE(6,630) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1        2,U2INF,MAD,MAD,MAD,MAD,MAD,MAD,NU2,PAD,PAD,PAD,PAD,PAD,
     2                                         NU2+1,(U2(I),I= 0,NU2)
                          FC2= 1.d0 - RMASS2/MASS2
                        ELSE
c ... for Coxon/Hajigeorgiou switching fx. & mass scaling fx. for atom-1
                          WRITE(6,632) 2,MASS2,IMN2,NAME2,MN2R,NAME2,
     1                                    2,U2INF,NU2,(U2(I),I= 0,NU2)
                          WRITE(6,636) Aad2, Rad2
                          FC2= ZME*(RMASS2 - MASS2)/(RMASS2*MASS2)
                        ENDIF
                      ENDIF
c
                  IF(NT1.GE.0) THEN
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
                          READ(5,*) T1INF,(T1(I), I=0,NT1)
                        ELSE
                          READ(5,*) T1INF,(T1(I),I=0,NT1),Ana1,Rna1
                        ENDIF
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
c ... for Huang/Le Roy centrifugal BOB radial functionfor atom-1 ...
                          WRITE(6,634) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1 1,T1INF,PNA,PNA,PNA,PNA,PNA,PNA,NT1,PNA,NT1+1,(T1(I),I= 0,NT1)
                          FG1= RMASS1/MASS1
                        ELSE
c ... for Coxon O-T polynomial multiplied by  m_e/M_i
                          WRITE(6,638) 1,MASS1,IMN1,NAME1,1,T1INF,NT1,
     1                                               (T1(I),I= 0,NT1)
                          WRITE(6,636) Ana1, Rna1
                          FG1= ZME/MASS1
                        ENDIF                     
                      ENDIF
c
                  IF(NT2.GE.0) THEN
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
                          READ(5,*) T2INF,(T2(I), I=0,NT2)
                        ELSE
                          READ(5,*) T2INF,(T2(I),I=0,NT2),Ana2,Rna2
                        ENDIF
c-----------------------------------------------------------------------
                      IF(PAD.GT.0) THEN
c ... for Huang/Le Roy centrifugal BOB radial function for atom-2 ...
                          WRITE(6,634) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1 2,T2INF,PNA,PNA,PNA,PNA,PNA,PNA,NT2,PNA,NT2+1,(T2(I),I= 0,NT2)
                          FG2= RMASS2/MASS2
                        ELSE
c ... for Coxon O-T polynomial multiplied by  m_e/M_i
                          WRITE(6,638) 2,MASS2,IMN2,NAME2,2,T2INF,NT2, 
     1                                               (T2(I),I= 0,NT2)
                          WRITE(6,636) Ana2, Rna2
                          FG2= ZME/MASS2
                        ENDIF
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
c** Generate a  Lennard-Jones(NSR,MPAR)  potential here.
c=======================================================================
      IF(IPOTL.EQ.1) THEN 
          NCN= MPAR
          XS= NSR
          XL= MPAR
          XDF= DSCM/(XS-XL)
          IF(LNPT.GE.0) WRITE(6,600) NSR,MPAR,DSCM,REQ
          CNN= XS*XDF*REQ**NCN
          DO  I= 1,NPP
              VV(I)= (XL*(REQ/XO(I))**NSR - XS*(REQ/XO(I))**MPAR)*XDF
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
          A0= DSCM
          IORD= NVARB-2
          X1= 1.d0
          FCT= PARM(NVARB-1)
          IF((MPAR.NE.0).AND.(DABS(FCT).GT.0.d0)) THEN
              FCT= 1.d0/PARM(NVARB-1)
              DO  J=1,IORD
                  X1= X1+ PARM(J)*FCT**J
                  ENDDO
              DSCM= DSCM*X1*FCT**2 + VMIN
              ENDIF
          IF(MPAR.EQ.1) THEN
c  Cases with power =1 (including Dunham, SPF & O-T expansions).
              IF(DABS(PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Dunham expansion ...
                  WRITE(6,612) PARM(NVARB),REQ,VMIN,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= -99
                  CNN= 0.d0
                  ENDIF
              IF(DABS(PARM(NVARB)).LE.0.d0) THEN
c ... print for Simons-Parr-Finlan expansion ...
                  WRITE(6,614) PARM(NVARB-1),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              IF(DABS(PARM(NVARB)-PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Ogilvie-Tipping expansion ...
                  WRITE(6,616) PARM(NVARB),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              ENDIF
          IF((MPAR.NE.0).AND.((MPAR.NE.1).OR.
     1                ((DABS(PARM(NVARB)-PARM(NVARB-1)).GT.0.d0).AND.
     2               (DABS(PARM(NVARB)*PARM(NVARB-1)).GT.0.d0)))) THEN
c ... print for general GPEF expansion variable ...
              IF(MPAR.LT.0) THEN
c ... for negative MPAR, convert to equivalent positive MPAR case
                  MPAR= -MPAR
                  A1= PARM(NVARB)
                  PARM(NVARB)= -PARM(NVARB-1)
                  PARM(NVARB-1)= -A1
                  ENDIF
              WRITE(6,618) MPAR,MPAR,PARM(NVARB-1),MPAR,PARM(NVARB),
     1                 MPAR,REQ,DSCM,A0,NVARB-2,(PARM(I),I= 1,NVARB-2)
              NCN= MPAR
              ENDIF
          IF(MPAR.EQ.0) THEN
c** For case of simple power series in  R  itself
              NCN= -1
              WRITE(6,620) NVARB,VMIN,(PARM(I),I= 1,NVARB)
              DO  I= 1, NPP
                  ZP= 1.d0
                  A1= VMIN
                  DO  J= 1,NVARB
                      ZP= ZP*XO(I)
                      A1= A1+ PARM(J)*ZP
                      ENDDO
                  VV(I)= A1
                  ENDDO
              VLIM= VV(NPP)
              RETURN
              ENDIF
c ... otherwise - generate potential as a GPEF-type expansion
          DO  I= 1, NPP
              ZZ= (XO(I)**MPAR - REQ**MPAR)/(PARM(NVARB-1)*XO(I)**MPAR
     1                                       + PARM(NVARB)*REQ**MPAR)
              A1= 1.d0
              ZP= 1.d0
              DO  J=1, NVARB-2
                  ZP= ZP*ZZ
                  A1= A1+ PARM(J)*ZP
                  ENDDO
              VV(I)= A0*ZZ*ZZ*A1 + VMIN
              ENDDO
          IF(DABS(PARM(NVARB-1)).GT.0) THEN
c ...Reset asymptote to avoid spurious  E > VLIM  warnings (e.g. for HO)
              VLIM= VV(NPP)
            ELSE
              VLIM= VV(1)
            ENDIF
          ENDIF
c
c=======================================================================
c** Generate a simple Morse, or Extended (EMOp) Morse potential, or as
c  special cases, Coxon's GMO or Wei Hua's generalized Morse
c=======================================================================
      IF(IPOTL.EQ.3) THEN
          BETA= PARM(1)
          NCN= 99
          IF(LNPT.GE.0) THEN
              IF(MPAR.EQ.-1) THEN
c** Option to generate Wei Hua's extended 4-parameter Morse-type potl.
                  CSAV= PARM(2)
                  WRITE(6,605) DSCM,REQ,CSAV,BETA
                ELSE 
                  IF(NVARB.LE.1) WRITE(6,606) DSCM,REQ,BETA
                  IF(NVARB.GT.1) THEN
                      IF(MPAR.GT.0) THEN
                          WRITE(6,608) MPAR,DSCM,REQ,NVARB-1,MPAR,MPAR,
     1                            MPAR,MPAR,NVARB,(PARM(i),i= 1,NVARB)
                          IF((NSR.GE.NVARB).OR.(NSR.LE.0)) NSR= NVARB
                          IF(NSR.LT.NVARB-1) WRITE(6,611) NSR
                          ENDIF
c ... Coxon's original GMO ...
                      IF(MPAR.EQ.-2) WRITE(6,610) DSCM,REQ,NVARB-1,
     1                                            (PARM(i),i= 1,NVARB)
                      ENDIF
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          NLR= NVARB-1
          DO  I= 1,NPP
c ... for Wei Hua's extended Morse function ...
              IF(MPAR.EQ.-1) THEN
                  VV(I)= DSCM*((1.d0 - DEXP(-BETA*(XO(I)-REQ)))/(1.d0 
     1                - CSAV*DEXP(-BETA*(XO(I)-REQ))))**2 - DSCM+ VLIM
                ELSE 
                  IF(NVARB.GT.1) THEN
                      ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
                      IF(MPAR.GT.1) ZZ= (XO(i)**MPAR - REQ**MPAR)/
     1                                  (XO(i)**MPAR + REQ**MPAR)
c ... for Coxon-Hajigeorgiou "GMO" potential
                      IF(MPAR.EQ.-2) ZZ= (XO(I)- REQ)
                      BETA= 0.d0
                      IORD= NLR+1
                      IF((ZZ.LT.0.d0).AND.(NSR.LT.NLR)) IORD= NSR+1
                      DO  J= IORD,1,-1             
                          BETA= BETA*ZZ+ PARM(J)
                          ENDDO
                      ENDIF
                  VV(I)=  DSCM*(1.d0 - DEXP(-BETA*(XO(I)-REQ)))**2 
     1                                                    - DSCM+ VLIM
                ENDIF
              ENDDO
          ENDIF
c
c=======================================================================
c** Generate an MLR (or MLJ) potential [as per Mol.Phys. 105, 691 (2007)
c   and JCP 112, 3949 (2000)]] 
c=======================================================================
      IF(IPOTL.EQ.4) THEN
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= CMM(1)
              pPOW= MPAR
              IF(MPAR.EQ.0) pPOW= 1
              IF(MPAR.LT.0) pPOW= -MPAR
              MLname= 'MLJ'
              IF(NCMM.GT.1) MLname= 'MLR'
              ULRe= 0.d0
              DO  I= 1,NCMM
                  ULRe= ULRe + CMM(I)/REQ**MMLR(I)
                  ENDDO
              BINF= DLOG(2.d0*DSCM/ULRe)
              WRITE(6,602) MLname,pPOW,NCN,DSCM,REQ
              IF(MPAR.LE.0) THEN
c  If appropriate, prepare to calculate switching function
                  NLR= NVARB- 4
                  Asw= PARM(NVARB-2)
                  Rsw= PARM(NVARB-1)
                  Rinn= PARM(NVARB)
                  IF(MPAR.NE.0) WRITE(6,603) NLR,1,1,1,1,1,NLR+1,
     1                                            (PARM(J),J= 1,NLR+1)
                  IF(MPAR.EQ.0) THEN
                      WRITE(6,598) NLR,1,1,1,1,1,NLR+1,
     1                                            (PARM(J),J= 1,NLR+1)
                      BINF= 0.5d0*BINF
                      ENDIF
                  WRITE(6,604) NCN,NCN,NCN,CMM(1),Asw,Rsw,BINF
                  IF(Rinn.GT.XO(1)) WRITE(6,596) Rinn
                ELSE
c  ... and for the two non-switching function cases ...
                  IF(NSR.LE.0) THEN
c  For case of simple exponent power series of order (NVARB-1) ...
                      NLR= NVARB-1
                      BINF= 0.D0
                      DO  I= 1,NVARB
                          BINF= BINF+ PARM(I)
                          ENDDO
                      CMM(1)= 2.d0*DSCM*REQ**NCN *DEXP(-BINF)
                      WRITE(6,603) NLR,MPAR,MPAR,MPAR,MPAR,MPAR,
     1                                     NLR+1,(PARM(J),J= 1,NLR+1)
                      WRITE(6,601) BINF,NCN,CMM(1)
                    ELSE
c  For THEOCHEM/Huang form:  \beta(yp)= Binf*yp + [1-yp]*{power series}
                      NLR= NVARB-1
                      IF(NSR.GT.NLR) NSR= NLR
                      WRITE(6,607) MPAR,MPAR,MPAR,NLR,MPAR,NLR+1,
     1                                           (PARM(J),J= 1,NLR+1)
                      IF(NSR.LT.NLR) WRITE(6,609) NSR
                      WRITE(6,613) MPAR,MPAR,MPAR,MPAR,MPAR,BINF,NCMM,
     1                                      (MMLR(I),CMM(I),I= 1,NCMM)
                    ENDIF
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          IORD= NLR
          DO  I= 1,NPP
              ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
c** MPAR=0  signals use of original O-T version of p=1 variable
              IF(MPAR.EQ.0) ZZ= 2.d0*ZZ
              IF(pPOW.GE.2) ZZ= (XO(i)**pPOW - REQ**pPOW)/
     1                          (XO(i)**pPOW  + REQ**pPOW)
              IORD= NLR
              IF((ZZ.LT.0).AND.(NSR.GT.0)) IORD= NSR
              BETA= 0.d0
              DO  J= IORD,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
c  Calculate and apply switching function to MLJ exponent coefficient
              IF(MPAR.LE.0) BETA= BINF+ (BETA- BINF)/
     1                                 (1.d0+ DEXP(Asw*(XO(I) - Rsw)))
              IF((MPAR.GT.0).AND.(NSR.GT.0)) BETA= BINF*ZZ + 
     1                                                 (1.d0- ZZ)*BETA
              ULR= 0.d0
              DO  J= 1,NCMM
                  ULR= ULR+ CMM(J)/XO(I)**MMLR(J)
                  ENDDO
              BETA= (ULR/ULRe) *DEXP(-BETA*ZZ)
              VV(I)= DSCM*(1.d0 - BETA)**2 - DSCM + VLIM
              ENDDO
          IF((MPAR.LE.0).AND.(XO(1).LT.Rinn)) THEN
c.... If necessary, prevent potential turnover in Coxon form by extrapolating 
c     exponent coeficient linearly inward from  R{inn}= PARM(NVARB)
              Jinn= INT((Rinn-XO(1))/(XO(2)-XO(1)) + 1.01)
              ZZ= 2.d0*(XO(Jinn+20) - REQ)/(XO(Jinn+20) + REQ)
              B2= 0.d0
              DO  J= IORD,0,-1
                  B2= B2*ZZ+ PARM(J+1)
                  ENDDO
              B2= BINF+(B2-BINF)/(1.d0+ DEXP(Asw*(XO(Jinn+20)- Rsw)))
              ZZ= 2.d0*(XO(Jinn) - REQ)/(XO(Jinn) + REQ)
              BETA= 0.d0
              DO  J= IORD,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
              BETA= BINF+ (BETA- BINF)/(1.d0+ DEXP(Asw*(XO(Jinn)- Rsw)))
              B1= (BETA-B2)/20.d0
              DO  I= Jinn-1,1,-1
                  BETA= BETA + B1
                  ZZ= 2.d0*(XO(I) - REQ)/(XO(I) + REQ)
                  ULR= 0.d0
                  DO  J= 1,NCMM
                      ULR= ULR+ CMM(J)/XO(I)**MMLR(J)
                      ENDDO
                  B2= (ULR/ULRe) *DEXP(-BETA*ZZ)
                  VV(I)= DSCM*(1.d0 - B2)**2 - DSCM + VLIM
                  ENDDO
              ENDIF
          ENDIF
c
c=======================================================================
c** Generate a DELR potential [as per JCP 119, 7398 (2003)] 
c=======================================================================
      IF(IPOTL.EQ.5) THEN
          IF(LNPT.GT.0) THEN
              NLR= NVARB - 3
              RHOd= PARM(NVARB-1)
              IDF= 1
              IF(PARM(NVARB).LE.0.d0) IDF= 2
              PPOW= MPAR
              REQP= REQ**PPOW
              A1= 0.0d0
              B1= 0.0d0
c... first, get  AA & BB and their derivatives!
              DO  J= 1,NCMM
                  CALL dampF(REQ,RHOd,MMLR(J),Dm,Dmp,Dmpp,IDF)
                  A1= A1+ CMM(J)*Dm/REQ**MMLR(J)*(1.d0+ Dmp/(PARM(1)*Dm)
     1                                        - MMLR(J)/(PARM(1)*REQ))
                  B1= B1+ CMM(J)*Dm/REQ**MMLR(J)*(2.d0+ Dmp/(PARM(1)*Dm)
     1                                        - MMLR(J)/(PARM(1)*REQ))
                  ENDDO
              A1 = A1 + DSCM
              B1 = B1 + 2.0d0*DSCM
              WRITE(6,650)  PPOW,DSCM,REQ,NSR,NLR,A1,B1,NCMM,
     1                                      (MMLR(J),CMM(J),J= 1,NCMM)
              WRITE(6,652)  NLR+1,(I-1,PARM(I),I= 1,NLR+1)
              IF(IDF.EQ.1) WRITE(6,654) RHOd
              IF(IDF.EQ.2) WRITE(6,656) RHOd
              ENDIF
c** Now ... generate potential function array for DELR form
          DO  I= 1, NPP
              ZZ= (XO(I)**PPOW -REQP)/(XO(I)**PPOW +REQP)
              IORD= NLR+1
              IF((ZZ.LT.0.d0).AND.(NSR.LT.NLR)) IORD= NSR+1
              BETA= 0.d0
c ... calculate the exponent
              DO  J= IORD,1,-1
                  BETA= BETA*ZZ+ PARM(J)
                  ENDDO
              BETA= DEXP(-BETA*(XO(I)-REQ))
c ... calculate the (damped) long-range tail
              A3= 0.0d0
              DO  J= 1, NCMM
                  CALL dampF(XO(I),RHOd,MMLR(J),Dm,Dmp,Dmpp,IDF)
                  A3= A3+ Dm*CMM(J)/XO(I)**MMLR(J)
                  ENDDO
              VV(I)=  (A1*BETA - B1)*BETA + A3 + VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.6) THEN
c=======================================================================
c** For generalized  HFD(m= MMLR(j), j=1,NCMM) potential with reduced 
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
c** Generate Tiemann-type polynomial potential attached to inverse-power
c  tail and 1/R^{12} (or exponential) inner wall [PRA 63, 012710 (2000)].
c  Polynomial expansion variable is  z= [R - Rm]/[R + b*Rm] where 
c  expansion has constant and linear terms.  The read-in DSCM= De (well
c  depth), but  Rm (read in as REQ) is not precisely Re (for a1 .neq. 0).
c  NCMM= number of inverse-power long-range terms;  
c  NVARB= (polynomial order) + 4.  [MPAR and NSR are dummy parameters]
c** Read-in parameters PARM(i) are in order: the  (IORD+1)=(NVARB-3) 
c  polynomial coefficients  a(0) - a(IORD), the expansion variable 
c  denominator factor b, and the the inner and outer bounds on the 
c  polynomial domain, Tiemann's Rinn & Rout, resp.  The powers and 
c  coefficients (-ve if attractive) of the NCMM inverse-power long-range
c  terms are MMCM(j) and CMM(j).
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= -CMM(1)
              A0= VLIM- DSCM
              IORD= NVARB - 4
              BT= PARM(NVARB-2)
              Rinn= PARM(NVARB-1)
              Rout= PARM(NVARB)
c** Determine analytic function attaching smoothly to inner wall of 
c  polynomial expansion at  R= Rinn < Rm
              ZZ= (Rinn - REQ)/(Rinn+ BT*REQ)
              ZP= 1.d0
              A1= PARM(1)
              A2= 0.d0
              DO  J= 1,IORD
                  A2= A2+ J*ZP*PARM(J+1)
                  ZP= ZP*ZZ
                  A1= A1+ ZP*PARM(J+1)
                  ENDDO
              A2= A2*(REQ+ BT*REQ)/(Rinn + BT*REQ)**2
c* If inward extrapolation is exponential:   A1*exp(-A2*(R-Rinn))
c             A2= -A2/A1
c* If inward extrapolation is inverse-power:   A1 + A2/R**12
              A2= -A2*Rinn**13/12.d0
              A1= A1 - A2/Rinn**12 + VLIM - DSCM
c** With long-range tail an NCMM-term inverse-power sum, add 1 additional
c   higher-power term to ensure continuity (not smoothness) at  Rout
c** NOTE attractive long-range terms have negative (-) coefficients!
              ZZ= (Rout - REQ)/(Rout+ BT*REQ)
              ZP= 1.d0
              B1= PARM(1)
              DO  J= 1,IORD
                  ZP= ZP*ZZ
                  B1= B1+ ZP*PARM(J+1)
                  ENDDO
              A3= DSCM
              DO  J= 1,NCMM
                  A3= A3+ CMM(J)/Rout**MMLR(J)
                  ENDDO
              MPAR= NCMM+ 1
              MMLR(MPAR)= MMLR(NCMM)+ 2
              CMM(MPAR)= (B1-A3)*Rout**MMLR(MPAR)
c*** Print for Tiemann-type potential
              IF(LNPT.GE.0) THEN
                  WRITE(6,640) DSCM,REQ,PARM(IORD+2),IORD,IORD+1, 
     1                                           (PARM(J),J= 1,IORD+1)
ccc               IF(XO(1).LT.Rinn) WRITE(6,642) PARM(NVARB-1),A1,A2,A0
                  IF(XO(1).LT.Rinn) WRITE(6,642) PARM(NVARB-1),A1,A2
                  IF(XO(NPP).GT.Rout) WRITE(6,644) PARM(NVARB),
     1                                   (CMM(J),MMLR(J),J= 1, MPAR)
                  ENDIF
              ENDIF
c ... now generate potential as a Tiemann-type expansion
          DO  I= 1, NPP
              IF(XO(I).LE.Rinn) THEN
c ... for exponential inward extrapolation ...
c                 VV(I)= A1*DEXP(-A2*(XO(I)- Rinn)) + A0
c ... for   A + B/R**12  inward extrapolation ...
                  VV(I)= A1 + A2/XO(I)**12
                ELSEIF(XO(I).LE.Rout) THEN
                  ZZ= (XO(I) - REQ)/(XO(I) + BT*REQ)
                  A3= A0 + PARM(1)
                  ZP= 1.d0
                  DO  J= 1,IORD
                      ZP= ZP*ZZ
                      A3= A3+ PARM(J+1)*ZP
                      ENDDO
                  VV(I)= A3
                ELSEIF(XO(I).GT.Rout) THEN
                  A3= VLIM
                  DO  J= 1, MPAR
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
c      functions to rotationless and/or centrifugal potential(s).
c  If  PAD > 0, use LeRoy-Huang radial functions; else Coxon-Hajig.
c=======================================================================
          IF(PAD.LE.0) THEN
              fad1e= 1.d0/(1.d0 + DEXP(Aad1*(REQ- Rad1)))
              fad2e= 1.d0/(1.d0 + DEXP(Aad2*(REQ- Rad2)))
              ENDIF
          DO  I=1,NPP
              IF(PAD.GT.0) THEN
c ... for LeRoy/Huang radial functions ...
                  YPAD= (XO(I)**PAD- REQ**PAD)/(XO(I)**PAD+ REQ**PAD)
                  YMAD1= (XO(I)**MAD- REQ**MAD)/(XO(I)**MAD+ REQ**MAD)
                  YMAD2= YMAD1
                  YPNA= (XO(I)**PNA- REQ**PNA)/(XO(I)**PNA+ REQ**PNA)
                  SC1= U1INF*YMAD1
                  SC2= U2INF*YMAD2
                  SG1= T1INF*YPNA
                  SG2= T2INF*YPNA
                  YMAD1= (1.d0- YMAD1)
                  YMAD2= (1.d0- YMAD2)
                  YMNA1= (1.d0- YPNA)
                  YMNA2= (1.d0- YPNA)
                ELSE
c** For [Coxon] BOB radial functions ....
                  YPAD= 2.d0*(XO(I) - REQ)/(XO(I) + REQ)
                  YPNA= YPAD
ccc               YPNA= (XO(i) - REQ)/REQ
                  YMAD1= 1.d0/(1.d0 + DEXP(Aad1*(XO(i)- Rad1)))
                  YMAD2= 1.d0/(1.d0 + DEXP(Aad2*(XO(i)- Rad2)))
                  IF(Aad1.LE.0.d0) YMAD1= 1.d0
                  IF(Aad2.LE.0.d0) YMAD2= 1.d0
                  YMNA1= 1.d0/(1.d0 + DEXP(Ana1*(XO(i)- Rna1)))
                  YMNA2= 1.d0/(1.d0 + DEXP(Ana2*(XO(i)- Rna2)))
                  IF(Ana1.LE.0.d0) YMNA1= 1.d0
                  IF(Ana2.LE.0.d0) YMNA2= 1.d0
                  SC1= U1INF*(1.d0 - YMAD1/fad1e)
                  SC2= U2INF*(1.d0 - YMAD2/fad2e)
                  SG1= 0.d0
                  SG2= 0.d0
                ENDIF
c ... finally, accumulate overall BOB terms ... all at the same time!
              DO  J= 0,NCMAX
                  SC1= SC1+ YMAD1*U1(J)
                  SC2= SC2+ YMAD2*U2(J)
                  SG1= SG1+ YMNA1*T1(J)
                  SG2= SG2+ YMNA2*T2(J)
                  YMAD1= YMAD1*YPAD
                  YMAD2= YMAD2*YPAD
                  YMNA1= YMNA1*YPNA
                  YMNA2= YMNA2*YPNA
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  601 FORMAT(3x,'which sum to   beta_{inf}=',f12.8,'   so that   C',i1,
     1  '=',1Pd14.7)
  602 FORMAT(/' Use an  ',A3,'_',i1,'(n=',i2,')  potential with   De=',
     1  F10.3,'[cm-1]    Re=',F12.8,'[A]')
  598 FORMAT(3x,'with exponent an order-',i2,' polynomial in   y',I1,
     1 '= 2*(r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',i1,')'/
     2 '   with',i3,' coefficients:',1PD16.8,2D16.8:/(8x,4D16.8:))
  596 FORMAT(5x,'and for   r <',f8.4,'   extrapolate inward linearly')
  603 FORMAT(3x,'with exponent an order-',i2,' polynomial in   y',I1,
     1 ' = (r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',i1,')'/
     2 '   with',i3,' coefficients:',1PD16.8,2D16.8:/(8x,4D16.8:))
  604 FORMAT(' exponent switching function yields limiting  C',i1,
     1 '/r**',i1,'  with  C_',i1,'=',1PD14.7/8x,'defined by   A_sw=',
     2  0Pf10.6,'   R_sw=',f10.6,'   and BINF=',1PD14.7)
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F11.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  607 FORMAT('   with  beta(y',i1,')= beta{INF}*y',i1,' + [1-y',i1,
     1 ']x{order-',i2,' polynomial in y',i1,'}'/'   with',i3,
     2 ' coefficients:',1PD16.8,2D16.8:/(8x,4D16.8:))
  608 FORMAT(/' Potential is an  EMO_',i1,'  with   De=',F11.4,
     1 '    Re=',F12.9/3x,'Exponent factor is order-',i2, ' power series
     2 in  y=(r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',i1,')'/
     3 '   with',I3,' coefficients:',1x,1PD18.9,2D18.9:/(7X,4D18.9:))
  609 FORMAT('   where for  r < Re  polynomial order is truncated to ord
     1er   NSR=',i2,'   with')
  610 FORMAT(/' Potential is Generalized Morse Oscillator with   De=',
     1 F10.3,'   Re=',F11.8/4x,'Exponent factor "beta" is',i3,' order po
     2wer series in (r-Re) with coefficients:'/4x,1PD18.9,3D18.9:/
     3 (4X,4D18.9:))
  611 FORMAT('   where for  r < Re  polynomial order is truncated to ord
     1er   NSR=',i2) 
  612 FORMAT(/' Potential is a Dunham expansion in  (r-Re)/(',f5.2,
     1  ' * Re)  with   Re=',f12.9/'  V(Re)=',f12.4,'    a0=',1PD16.9,
     2  '   and',i3,'  a_i coefficients:'/(5D16.8))
  613 FORMAT(6x,'y',i1,'= (r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',
     1 i1,')     where   phiINF=',F12.8,'   and'/  '   uLR(r) is a sum o
     2f',I3," term(s) with coefficients (+'ve if attractive):"/
     3 (3x,3(5x,'C',i2,' =',1PD14.6:)))
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
  622 FORMAT(/' *** ERROR in generating HFD potential *** generate   ALF
     1A=',G15.7,'  from reduced  Cm  coefficients:'/
     2   (3x,3('   C',i2,' =',1PD15.7:)) )
  624 FORMAT(/' Potential is Generalized HFD with exponent factors   gam
     1ma=',f9.6/'   beta1=',f12.8,'   beta2=',f9.6,'   A=',1PD16.9,
     2 "   & reduced Cm's:"/(3x,3('   C',i2,' =',D15.7:)) )
  626 FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]   and'/
     1  '     Damping function  D(r)= exp[ -',0Pf6.4,'*(',f7.4,
     1  '/X -1.0)**',f5.2,']')
  630 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1 '  of mass ',f15.11/'   consists of mass factor  [1- MASS(',I3,
     2 A2,')/MASS(',I3,A2,')]  multiplying all of:'/5x,'u',i1,'INF=',
     3 f11.6,'  times  y',i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,
     4 ' + Re**',i1,')]'/5x,'plus  [1 - y',i1,']  times an order',I3,
     5 ' polynomial in'/7x,'y',i1,'=[(r**',i1,' - Re**',i1,')/(r**',i1,
     6 ' + Re**',i1,')]  with the ',i3,' coefficients:'/(3x,4G17.9:))
  632 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1  '  of mass ',f15.11/'   consists of mass factor  m{electron}*[1/
     2MASS(',I3,A2,') - 1/MASS(',I3,A2,')]'/5x,'multiplying   u',i1,
     3 'INF=',1PD17.9,'  times [1 - fsw(r)/fsw(Re)]'/'   plus  fsw(r)  t
     4imes an order',0P,i3,' polynomial in z{Dun} with coefficients:' / 
     5  (3x,1P4D17.9:)) 
  634 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor  [MASS(',I3,A2,')/MASS(',I3,
     2 A2,')]  multiplying all of:'/5x,'q',i1,'INF=',F11.6,' times  y',
     3 i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',i1,')]'/
     4 5x,'plus [1 - y',i1,'] times an order',I3,' polynomial in y',i1,
     5 ' with the',i3,' coefficients:'/(3x,4G17.9:))
  636 FORMAT(3x,'where   fsw(r) = 1/[1 + exp{',f7.4,'*(r -',f7.4,')}]')
  638 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor   [mass{electron}/MASS(',I3,
     2 A2,')]'/'   multiplying   q',i1,'INF=',1PD17.9,'  times [1 - fsw(
     3r)/fsw(Re)]'/ '   plus  fsw(r)  times an order',0P,i3,' polynomial
     4 in z{O-T} with coefficients:'/ (3x,4G17.9:)) 
  640 FORMAT(/' Tiemann-type potential with   De=',F11.4,'   Rm=',f9.6,
     1 '   is a power series'/10x,'in  (r - Re)/(r ',SP,F9.5,
     2 '*Re) of order',SS,I3,'  with the',I3,' coefficients:'/(5D16.8))
c 642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',1PD13.6,'*exp[-',
c    1  0PF9.6,'*(r - Rinn)] ',SP,F10.3)
  642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',SP,F12.4,1x,1PD13.6,
     1  '/R**12' )
  644 FORMAT('  and  for  r > Rout=',F7.3,'   V= VLIM ',
     1 (SP,1PD14.6,'/r**',SS,I2):/(39x,SP,1PD14.6,'/r**',SS,I2))
  650 FORMAT(/' DELR(p=',i2,') potential with   DSCM=', F11.4,'   REQ=',
     1 F11.8,'   and exponent'/8x,'coefficient power series of order',
     2 i3,' for  r < Re','  and',i3,' for  r >= Re'/3x,'Generate   A(DEL
     3R)=',1Pd17.9,'   B(DELR)=',D17.9/6x,'where ULR defined by',I2,
     4 " inverse-power terms with coeffts (+'ve repulsive):"/
     5 (5x,3(5x,'C',0P,i2,' =',1Pd14.6:)))
  652 FORMAT(3x,'Exponent defined by',I3,' expansion coefficients:',6x,
     1 'phi(',i2,')=',1PD17.10/(2x,3(3x,'phi(',0P,i2,')=',1PD15.8:)))
  654 FORMAT(3x,'U_{LR} term uses Tang-Toennies damping function [JCP 80
     1,3726(1984)]'/10x,'[1 - exp(-RHOd*r)*SUM{(RHOd*r)^k/k!}]',
     2  '  with   RHOd =',F9.6)
  656 FORMAT(3x,'U_{LR} term uses Douketis dispersion damping function [
     1Mol.P. 52,763(1984)]'/5x,'[1 - exp(-3.97*RHOd*r/m - 0.39*(RHOd*r)^
     22/sqrt{m})]^m  with  RHOd=',F9.6)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

