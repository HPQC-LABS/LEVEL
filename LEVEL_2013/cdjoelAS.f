c***********************************************************************
      SUBROUTINE CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V,WF0,RM2,RCNST)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Subroutine solving the linear inhomogeneous differential equations
c  formulated by J.M. Hutson [J.Phys.B14, 851 (1982)] for treating 
c  centrifugal distortion as a perturbation, to determine centrifugal 
c  distortion constants of a diatomic molecule.  Uses the algorithm of
c  J. Tellinghuisen [J.Mol.Spectrosc. 122, 455 (1987)].  The current
c  version calculates Bv, Dv, Hv, Lv, Mv, Nv and Ov and writes them out, 
c  but does not return values to the calling program.
c
c** On entry:   EO    is the eigenvalue (in units [cm-1])
c               NBEG & NEND  the mesh point range over which the input
c wavefunction  WF0  (in units 1/sqrt(Ang))  has non-negligible values
c               BvWn  is the numerical factor (hbar^2/2mu) [cm-1 Ang^2]
c               YH    is the integration stepsize 
c               WARN  is an integer flag: > 0 print internal warnings,
c               V(i)  is the effective potential (including centrifugal
c                     term if calculation performed at  J > 0) in 
c                     'internal' units, including the factor  YH**2/BvWN
c               RM2(i) is the array  (r')^2/(distance**2) 
c** On exit:    RCNST(i)  is the set of 7 rotational constants: Bv, -Dv,
c                       Hv, Lv, Mv, Nv & Ov
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 1994  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Authors: R.J. Le Roy & J. Tellinghuisen         Version of 30/09/1999
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Dimension:  potential arrays  and  vib. level arrays.
      INTEGER NDMINT
      PARAMETER (NDMINT= 200001)
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
      REAL*8 V(NEND),WF0(NEND),RM2(NEND),P(NDMINT),WF1(NDMINT),
     1                                            WF2(NDMINT),RCNST(7)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,YH,YH2,
     1  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,
     2  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,
     3  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,
     4  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33
c
c!!
      INTEGER NDIMR
      PARAMETER (NDIMR=200001)
      REAL*8 pRV,aRV,RFN(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
     1                                                    SDRDY(NDIMR)
      COMMON /BLKAS/pRV,aRV,RFN,YVB,DRDY2,FAS,SDRDY
c!!
      IF(NEND.GT.NDMINT) THEN
          WRITE(6,602) NEND,NDMINT
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      YH2 = YH*YH
      DV = YH2*ZTW
      E= EO*YH2/BvWN
      IPASS = 1
      OV01 = 0.D0
      OV02 = 0.D0
      OV03 = 0.D0
      OV11 = 0.D0
      OV22 = 0.D0
      OV12 = 0.D0
      OV33 = 0.D0
      OV23 = 0.D0
      OV13 = 0.D0
      PER01 = 0.D0
      PER02 = 0.D0
      PER03 = 0.D0
      PER11 = 0.D0
      PER12 = 0.D0
      PER13 = 0.D0
      PER22 = 0.D0
      PER23 = 0.D0
      PER33 = 0.D0
c** First, calculate the expectation value of  1/r**2  and hence Bv
      R2IN= 0.5D0*(RM2(NBEG)*WF0(NBEG)**2 + RM2(NEND)*WF0(NEND)**2)
      DO   I= NBEG+1, NEND-1
         R2IN= R2IN+ RM2(I)*WF0(I)**2
         ENDDO
      R2IN = R2IN*YH
      RCNST(1)= R2IN*BvWN
c
c** On First pass  IPASS=1  and calculate first-order wavefx., Dv & Hv
c  On second pass  IPASS=2  and calculate second-order wavefx., Lv & Mv
c  On third pass   IPASS=3  and calculate third-order wavefx., Nv & Ov
c
   10 P1= 0.D0
      P2= 0.D0
c
c     P1= WF0(NEND)
c     P2= WF0(NEND-1)
c
      P(NEND) = P1
      P(NEND-1) = P2
      V1 = V(NEND) - E*DRDY2(NEND)
      V2 = V(NEND-1) - E*DRDY2(NEND-1)
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) 
     1                   - DV*(RM2(NEND) - R2IN*DRDY2(NEND))*WF0(NEND)
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF0(NEND-1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF0(NEND)
     1                     - (RM2(NEND) - R2IN*DRDY2(NEND))*WF1(NEND))
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF1(NEND-1) 
     1                                 - DVV*WF0(NEND-1)*DRDY2(NEND-1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF1(NEND) - HVV*WF0(NEND)
     1                     - (RM2(NEND) - R2IN*DRDY2(NEND))*WF2(NEND))
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF2(NEND-1)
     1             - (DVV*WF1(NEND-1) + HVV*WF0(NEND-1))*DRDY2(NEND-1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
c** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + YH2*G2 + V2*P2
          R2XX= R2IN*DRDY2(M)
          IF(IPASS.EQ.1) G3 = (RM2(M)- R2XX)*WF0(M)
          IF(IPASS.EQ.2) G3 = (RM2(M)-R2XX)*WF1(M) - DVV*WF0(M)*DRDY2(M)
          IF(IPASS.EQ.3) G3 = (RM2(M)- R2XX)*WF2(M) 
     1                            - (DVV*WF1(M) + HVV*WF0(M))*DRDY2(M)
          V3 = V(M) - E*DRDY2(M)
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          IF(V3.LT.0.D0)  GO TO 32
          P(M) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          ENDDO
      GO TO 90
c** Escaped loop at outer turning point:  initialize outward integration
   32 PRS = P3
      PRT = P(M+1)
      P1 = 0.D0
      P2 = 0.D0
c
c     P1 = WF0(NBEG)
c     P2 = WF0(NBEG+1)
c
      P(NBEG) = P1
      P(NBEG+1) = P2
      V1 = V(NBEG) - E*DRDY2(NBEG)
      V2 = V(NBEG+1) - E*DRDY2(NBEG+1)
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) 
     1                   - DV*(RM2(NBEG) - R2IN*DRDY2(NBEG))*WF0(NBEG)
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF0(NEND)*DRDY2(NBEG)
     1                     - (RM2(NBEG) - R2IN*DRDY2(NBEG))*WF1(NBEG))
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF1(NBEG+1)
     1                                 - DVV*WF0(NBEG+1)*DRDY2(NBEG+1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*((DVV*WF1(NEND) + HVV*WF0(NEND))
     1       *DRDY2(NBEG)  - (RM2(NBEG) - R2IN*DRDY2(NBEG))*WF2(NBEG))
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF2(NBEG+1)
     1             - (DVV*WF1(NBEG+1) + HVV*WF0(NBEG+1))*DRDY2(NBEG+1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
c** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + YH2*G2 + V2*P2
          P0 = WF0(I)
          R2XX= R2IN*DRDY2(I)
          IF(IPASS.EQ.1) G3 = (RM2(I)-R2XX)*P0
          IF(IPASS.EQ.2) G3 = (RM2(I)-R2XX)*WF1(I) - DVV*P0*DRDY2(I)
          IF(IPASS.EQ.3) G3 = (RM2(I)-R2XX)*WF2(I) 
     1                                - (DVV*WF1(I) + HVV*P0)*DRDY2(I)
          V3 = V(I) - E*DRDY2(I)
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          AR = AR + P0*P3*DRDY2(I)
          ENDDO
c** Average for 2 adjacent mesh points to get Joel's "(a-b)"
      AMB2 = (P3-PRT)/P0
      AMB1 = (P(M)-PRS)/WF0(M)
      AMB = (AMB1+AMB2)*0.5D0
      M2 = M+2
c** Find the rest of the overlap with zero-th order solution ...
      DO  I = M2,NEND
          P0 = WF0(I)
          PI = P(I) + AMB*P0
          P(I) = PI
          AR = AR + PI*P0*DRDY2(I)
          ENDDO
      OV = AR*YH
      DO  I = NBEG,NEND
          P0 = WF0(I)
c ... and project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*RM2(I)
          IF(IPASS.EQ.1) THEN
c** Now - on first pass accumulate integrals for Dv and Hv
              WF1(I) = PI
              OV01 = OV01 + PI*P0 * drdy2(i)
              OV11 = OV11 + PI*PI * drdy2(i)
              PER01 = PER01 + PIF*P0
              PER11 = PER11 + PI*PIF
            ELSEIF(IPASS.EQ.2) THEN
c ... and on next pass, accumulate integrals for Lv and Mv
              WF2(I) = PI
              P1 = WF1(I)
              OV02 = OV02 + PI*P0 * drdy2(i)
              OV12 = OV12 + PI*P1 * drdy2(i)
              OV22 = OV22 + PI*PI * drdy2(i)
              PER02 = PER02 + PIF*P0
              PER12 = PER12 + PIF*P1
              PER22 = PER22 + PI*PIF
            ELSEIF(IPASS.EQ.3) THEN
c ... and on next pass, accumulate integrals for Nv and Ov
              P1 = WF1(I)
              P2 = WF2(I)
              OV03 = OV03 + PI*P0 * drdy2(i)
              OV13 = OV13 + PI*P1 * drdy2(i)
              OV23 = OV23 + PI*P2 * drdy2(i)
              OV33 = OV33 + PI*PI * drdy2(i)
              PER03 = PER03 + PIF*P0
              PER13 = PER13 + PIF*P1
              PER23 = PER23 + PIF*P2
              PER33 = PER33 + PIF*PI
            ENDIF
          ENDDO
      IF(IPASS.EQ.1) THEN
          DVV = YH*PER01
          HVV = YH*(PER11 - R2IN*OV11)
          IPASS = 2
          RCNST(2) = DVV*BvWN
          RCNST(3) = HVV*BvWn
          GO TO 10
        ELSEIF(IPASS.EQ.2) THEN
          HV2 = YH*PER02*BvWN
          LVV = YH*(PER12 - R2IN*OV12 - DVV*OV11)
          MVV = YH*(PER22 - R2IN*OV22 - 2.D0*DVV*OV12 - HVV*OV11)
          IPASS = 3
          RCNST(4) = LVV*BvWN
          RCNST(5) = MVV*BvWN
          GO TO 10
        ELSEIF(IPASS.EQ.3) THEN
          LV2 = YH*PER03*BvWN
          MV2 = YH*(PER13 - R2IN*OV13 - DVV*OV12 - HVV*OV11)*BvWN
          NVV = YH*(PER23 - R2IN*OV23 - DVV*(OV13 + OV22) 
     1                                     - 2.D0*HVV*OV12 - LVV*OV11)
          OVV = YH*(PER33 - R2IN*OV33 - 2.D0*DVV*OV23 
     1             - HVV*(2.D0*OV13+ OV22) - 2.D0*LVV*OV12 - MVV*OV11)
          RCNST(6) = NVV*BvWN
          RCNST(7) = OVV*BvWN
        ENDIF
      IF(WARN.GT.0) THEN
          IF(DMAX1(DABS(OV01),DABS(OV02),DABS(OV01)).GT.1.D-9)
     1                                     WRITE(6,604) OV01,OV02,OV03
          TSTHV= dabs(RCNST(3)/HV2-1.D0)
          TSTLV= dabs(RCNST(4)/LV2-1.D0)
          TSTMV= dabs(RCNST(5)/MV2-1.D0)
          IF(DMAX1(TSTHV,TSTLV,TSTMV).GT.1.d-5)
     1                                  WRITE(6,603) TSTHV,TSTLV,TSTMV
          ENDIF
      DO  M= 2, 7
c** Kill nonsensical high-order CDCs (which can occur in double-well cases)
          IF(DABS(RCNST(M)).GT.DABS(RCNST(M-1))) THEN
              DO I= M, 7
                  RCNST(I)= 0.d0
                  ENDDO
              EXIT
              ENDIF
          ENDDO
      RETURN
   90 WRITE(6,601) EO
      RETURN
  601 FORMAT(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,
     1   '  never reach outer turning point')
  602 FORMAT(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,
     1  ' > NDMINT=',i6)
  603 FORMAT(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',
     1 3(1Pd9.1))
  604 FORMAT(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03
     1:',3(1Pd9.1))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

