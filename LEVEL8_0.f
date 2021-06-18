c***********************************************************************
c******  Eigenvalue program  LEVEL 8.0 : as of  8 April 2007  ********
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2005-07  by  Robert J. Le Roy               +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Please inform me of any bugs, by phone at: (519)888-4567, ext. 4051 +
c++++++++ by e-mail to: leroy@UWaterloo.ca , or write me at: +++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Program for calculating eigenvalues and eigenfunctions (and if
c  desired, also various expectation values & matrix elements) of a
c  one-dimensional potential, and/or matrix elements (& Franck-Condon
c  factors) between levels of two different potentials.
c** As with most similar codes, the original version of this program was
c  based on the Franck-Condon intensity program of R.N. Zare, report
c  UCRL-10925(1963), but the present version is massively modified.
c** This program is unique in that it can:  (1) automatically locate &
c      calculate the widths of quasibound levels (orbiting resonances);
c  (2) can calculate diatomic molecule centrifugal distortion constants;
c  (3) can find levels in either well of a double minimum potential;
c  (4) starting from a single suitable (almost arbitrary) trial energy,
c      it will also automatically generate the eigenvalues etc. for all
c      vibrational and/or rotational levels of a given well-behaved
c      single-minimum potential.
c***** Main calling and I/O routines.  Last Updated  8 April 2007 *****
c----------------------------------------------------------------------
c** Dimension for  potential arrays  and  vib. level arrays.
      INTEGER NDIMR,VIBMX,MORDRMX,RORDR,NTPMX
      PARAMETER (NDIMR=80001, VIBMX=400,RORDR=7,MORDRMX=20,NTPMX= 1600)
c
      INTEGER I,J,M,III,IJD,ILEV1,ILEV2,IOMEG1,IOMEG2,INNOD1,INNOD2,
     1 INNER,SINNER,IQT,IWR,IRFN,IVD,IVS,IAN1,IAN2,IMN1,IMN2,GEL1,GEL2,
     2 GNS1,GNS2,JDJR,JCT,J2DL,J2DU,J2DD,JROT,JROT2,JLEV,JREF, ICOR,
     3 CHARGE, KV,KV2,KVIN,LCDC,LPRWF,LRPT,LXPCT,MORDR,NUSEF,ILRF,IR2F,
     4 NUMPOT,NBEG,NBEG2,NEND,NEND2,NPP,NCN1,NCN2,NCNF,NLEV,NLEV1,
     5 NLEV2,NJM,NJMM,NFP,NLP,NRFN,NROW,WARN,VMAX,VMAX1,VMAX2,AFLAG,
     6 AUTO1,AUTO2, IV(VIBMX),IJ(VIBMX),IV2(VIBMX),JWR(VIBMX),
     7 INNR1(0:VIBMX),INNR2(0:VIBMX), NoPRWF
c
      REAL*8 ZK1(0:VIBMX,0:RORDR),ZK2(0:VIBMX,0:RORDR),RCNST(RORDR),
     1 V1(NDIMR),V2(NDIMR),VJ(NDIMR),WF1(NDIMR),WF2(NDIMR)
c
      REAL*8  RFN(NDIMR),RM2(NDIMR),RM22(NDIMR),GV(0:VIBMX), 
     2  ESOLN(VIBMX),ESLJ(VIBMX), XIF(NTPMX),YIF(NTPMX),
     4  ABUND1,ABUND2,MASS1,MASS2,DM(0:MORDRMX)
c
      REAL*8 BZ,BvWN,BFCT,BEFF,DEJ,EPS,EO,EO2,EJ,EJ2,EJP,EJREF,GAMA,
     1  MEL,PMAX1,PMAX2,PW,RMIN,RMAX,RH,RMINN,RREF,RREFP,RFLIM,CNNF,
     2  RFACTF,MFACTF,SOMEG1,SOMEG2,VLIM1,VLIM2,VD,VDMV,XX,ZMU, 
     3  GI,GB,GBB,WV
c
      CHARACTER*78 TITL
      CHARACTER*2 NAME1,NAME2
c
      DATA MEL/5.4857990945d-4/
c** Default (Q-branch) defining J-increments for matrix element calcn.
      DATA J2DL,J2DU,J2DD/0,0,1/
      NoPRWF= 0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Begin by reading in the (integer) atomic numbers and mass numbers
c  defining the effective reduced mass of the system considered.
c** IAN1 & IAM2, and IMN1 & IMN2 are, respectively, the atomic numbers
c    and the mass numbers identifying the atoms forming the molecule.
c    Their masses are extracted from data subroutine MASSES and used
c    to generate the the reduced mass ZMU.
c** If  IMN1  or  IMN2  lie outside the range of mass numbers for normal
c  stable isotopes of that species, subroutine MASSES returns the 
c  average atomic mass based on the natural isotope abundance.
c** If the read-in value of IAN1 and/or IAN2 is .LE.0, then instead of
c  using the MASS table, read an actual particle mass for it/them.
c** CHARGE (integer) is the charge on the molecule (=0 for neutral). If
c   (CHARGE.ne.0)  generate & use Watson's charge-adjusted reduced mass.
c** Parameter NUMPOT specifies whether to calculate eigenvalues etc. for
c  a single potential (when NUMPOT.LE.1), or to generate two independent
c  potentials & calculate matrix elements coupling levels of one to
c  levels of the other (for NUMPOT.GE.2).
c----------------------------------------------------------------------
    2 READ(5,*,END=999) IAN1, IMN1, IAN2, IMN2, CHARGE, NUMPOT
c----------------------------------------------------------------------
c** Subroutine MASSES returns the names of the atoms NAMEi,ground
c  electronic state degeneracy GELi, nuclear spin degeneracy GNSi,
c  mass MASSi, and isotopic abundance ABUNDi for a given atomic isotope.
      IF((IAN1.GT.0).AND.(IAN1.LE.109)) THEN
          CALL MASSES(IAN1,IMN1,NAME1,GEL1,GNS1,MASS1,ABUND1)
        ELSE
c** If particle-i is not a normal atomic isotope, read a 2-character
c   name (enclosed between '', as in 'mu') and its actual mass.
c----------------------------------------------------------------------
          READ(5,*) NAME1, MASS1
c----------------------------------------------------------------------
        ENDIF
      IF((IAN2.GT.0).AND.(IAN2.LE.109)) THEN
          CALL MASSES(IAN2,IMN2,NAME2,GEL2,GNS2,MASS2,ABUND2)
        ELSE
c----------------------------------------------------------------------
          READ(5,*) NAME2, MASS2
c----------------------------------------------------------------------
        ENDIF
      ZMU= MASS1*MASS2/(MASS1+MASS2- CHARGE*MEL)
c=======================================================================
c TITL is a title or output header of up to 78 characters, read on a
c   single line enclosed between single quotes: e.g.  'title of problem'
c=======================================================================
      READ(5,*) TITL
c----------------------------------------------------------------------
c** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
c  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920D0
      BvWN= 1.D0/BZ
      WRITE(6,605) TITL,ZMU,BZ,MASS1,MASS2
      IF(CHARGE.NE.0) WRITE(6,624) CHARGE,CHARGE
      EJ= 0.D0
      EJ2= 0.D0
      LRPT= 1
c** Lower limit (RMIN) and increment (RH) of integration in (Angstroms).
c** The upper limit (RMAX) of the integration range is automatically
c  set at the SMALLER of: (i) the read-in value, and  (ii) the largest 
c  value allowed by the array-size dimensions.
c* A hard wall boundary condition may be imposed at a smaller distance
c  using an appropriate choice of the read-in level parameter IV (below)
c** EPS (cm-1) is the desired eigenvalue convergence criterion
c---------------------------------------------------------------------
      READ(5,*) RH, RMIN, RMAX, EPS
c---------------------------------------------------------------------
c** NPP = no. of points in potential and wavefunction array.
      NPP= INT((RMAX-RMIN)/RH+ 1.00001)
      NPP= MIN0(NDIMR,NPP)
      RMINN= RMIN-RH
      RMAX= RMINN+ NPP*RH
      WRITE(6,604) RMIN,RMAX,RH,NAME1,IMN1,NAME2,IMN2
      DO  I= 2,NPP
          RFN(I)= RMINN+I*RH
          WF1(I)= RFN(I)
          RM2(I)= 1.D0/RFN(I)**2
          RM22(I)= RM2(I)
          ENDDO
      RFN(1)= RMIN
      WF1(1)= RMIN
      RM2(1)= RM2(2)
      IF(RMIN.GT.0.D0) RM2(1)= 1.D0/RMIN**2
      RM22(1)= RM2(1)
c
c++ Begin reading appropriate parameters & preparing potential(s)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+  Subroutine "PREPOT" prepares (and if desired, writes) the potential 
c+  array V(i) (cm-1)  at the NPP distances RFN(i) (Angst).
c** NPP = no. of points in potential and wavefunction array.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c* If NTP > 0 :  define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
c   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
c* If LPPOT > 0 : at every |LPPOT|-th point, print potential and
c      derivatives-by-differences. ***  If  LPPOT < 0  write potential
c      at every |LPPOT|-th point to channel-8 in a compact format **
c* OMEGA is the electronic contribution to the angular momentum such
c  that the reduced centrifugal potential is:  (J*(J+1)-OMEGA**2)/R**2
c* Set (OMEGA.GE.99) if wish to use centrifugal factor for rotation
c  in two dimensions:   (J**2 - 1/4)/R**2  .
c* VLIM (cm-1) is the energy associated with the potential asymptote.
c-----------------------------------------------------------------------
c++   READ(5,*) NTP, LPPOT, OMEGA, VLIM
c----------------------------------------------------------------------
c** For pointwise potentials, PREPOT uses subroutine GENINT to read
c  points and conditions and interpolate (using subroutines NTRPSR,
c  SPLINE & SPLINT) and extrapolate to get full potential.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP (read above) is number of turning points (XI,YI) to be read in.
c** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
c    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
c    interpolate with cubic spline instead of local polynomials.
c** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
c   This may help if interpolation has trouble on steep repulsive wall.
c** ILR specifies how to extrapolate beyond largest input distance XI(i)
c  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
c  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
c* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c* If ILR > 3 : successive higher power terms differ by factor  1/R
c
c** RFACT & EFACT are factors required to convert units of input turning
c       points (XI,YI) to Angstroms & cm-1, respectively (often = 1.d0)
c** Turning points (XI,YI) must be ordered with increasing XI(I)
c** Energy VSHIFT (cm-1) is added to the input potential points to
c   make their absolute energy consistent with VLIM (often VSHIFT=Te).
c-----------------------------------------------------------------------
c++   READ(5,*) NUSE, IR2, ILR, NCN, CNN
c++   READ(5,*) RFACT, EFACT, VSHIFT
c++   READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-----------------------------------------------------------------------
c** NCN1 (returned by PREPOT) is the power of the asymptotically-
c  dominant inverse-power long range potential term.   
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RFN,RM2,VLIM1,V1,
     1                                                            NCN1)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If (NTP.le.0) PREPOT uses subroutine POTGEN to generate a fully 
c  analytic potential defined by the following read-in parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
c  and for all cases except IPOTL=2, the potential asymptote energy is
c  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
c  potential minimum and  DSCM  the leading (quadratic) potential coeft.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR, NSR & NCMM are integers cwcharacterizing the chosen potential
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass
c      dependent Born-Oppenheimer breakdown corrections will be included
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NCN) potential.
c** If IPOTL=2  use Seto's modification of Surkus' GPEF expansion in
c       z = [R**MPAR - Re**MPAR]/[a*R**MPAR + b*Re**MPAR] where
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NCN=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NCN=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NCN=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus MPAR < 0 case:  z(MPAR,a,b)= z(|MPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c      with exponent factor "beta" defined as a power series of order
c      (NVARB-1) in  y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR)
c      with NVARB coefficients PARM(i).    [!! MPAR .ge.1 !!]
c    * For conventional "simple" Morse potential,  NVARB=1 & MPAR dummy
c*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  generate an MLR potential [Mol.Phys. 105, 691 (2007)]
c      If MPAR > 0  exponent parameter defined in terms of a polynomial
c           of order (NVARB-2) with the NVARB coefficients  PARM(j).
c           in y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR),
c           and long-range defined by NCN inverse-power terms
c      If MPAR = 0  exponent polynomial variable is  2*y_{1}= y{O-T}
c      If MPAR < 0  exponent polynomial variable is  y_{|MPAR|}
c      If MPAR.le.0  exponent polynomial connected to limiting inverse-
c           power potential exponent by exponential switching function
c           with parameters  Asw= PARM(NVARB-1)  and  RSW= PARM(NVARB).
c** If IPOTL=5  generate a Double-Exponential Long-Range (DELR)
c       potential [JCP 119, 7398 (2003)] with additive long-range part
c       defined by a sum of NCMM damped inverse-power terms, & exponent
c       polynomial radial variable defined by parameter MPAR (=p)
c** If IPOTL=6  generate generalized HFD({m_i},i=1,NCMM) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c** If IPOTL=7  use Tiemann-type polynomial potential attached to an
c     inverse-power long-range tail and an 1/R^{12} (or exponential)
c     inner wall.
c----------------------------------------------------------------------
c++     READ(5,*) IPOTL, MPAR, NSR, NCMM, NVARB, IBOB, DSCM, REQ
c++     IF(IPOTL.GE.4) READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
c++     IF((IPOTL.EQ.4).OR.(IPOT.EQ.7)) READ(5,*) (MMLR(I),CMM(I),I= 1,MPAR)
c++     IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
c++     IF(IBOB.GT.0) THEN
c++         READ(5,*) MN1R, MN2R, PAD, MAD, NU1, NU2, PNA, NT1, NT2
c++         IF(PAD.GT.0) THEN
c++              IF(NU1.GE.0) READ(5,*) U1INF, (U1(I), I=0,NU1)
c++              IF(NU2.GE.0) READ(5,*) U2INF, (U2(I), I=0,NU2)
c++              IF(NT1.GE.0) READ(5,*) T1INF, (T1(I), I=0,NT1)
c++              IF(NT2.GE.0) READ(5,*) T2INF, (T2(I), I=0,NT2)
c++           ELSE
c++              IF(NU1.GE.0) READ(5,*) U1INF,(U1(I), I=0,NU1), Aad1, Rad1
c++              IF(NU2.GE.0) READ(5,*) U2INF,(U2(I), I=0,NU2), Aad2, Rad2
c++              IF(NT1.GE.0) READ(5,*) T1INF,(T1(I), I=0,NT1), Ana1, Rna1
c++              IF(NT2.GE.0) READ(5,*) T2INF,(T2(I), I=0,NT2), Ana2, Rna2
c++           ENDIF
c++         ENDIF
c++     ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      BFCT= BZ*RH*RH
      PW= 2.D0
      IF((NCN1.GT.0).AND.(NCN1.NE.2)) PW= 2.D0*NCN1/(NCN1-2.D0)
      DO  I= 1,NPP
          V1(I)= V1(I)*BFCT
          V2(I)= V1(I)
          ENDDO
      VLIM2= VLIM1
      IF(NUMPOT.LE.1) THEN
          WRITE(6,636)
          IOMEG2= IOMEG1
c** For case in which Potl-1 has centrifugal BOB function, systematize
c   definition of centrifugal potential for 'case-2'
          DO  I= 1,NPP
              RM22(I)= RM2(I)
              ENDDO
        ELSE
          WRITE(6,635)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For 2-potential Franck-Condon factor calculation, get the second
c  potential in this second call to PREPOT (uses the same parameter
c  reading sequence so exhaustively described immediately above).
          CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RFN,RM22,
     1                                                  VLIM2,V2,NCN2)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Convert potential (in (cm-1)) to form appropriate for SCHRQ
          DO  I=1,NPP
              V2(I)= V2(I)*BFCT
              ENDDO
        ENDIF
c
c** NLEV1 is the no. of levels {v=IV(i), J=IJ(i)} of potential-1 which
c   we wish to find.
c* IF(NLEV1=0) calculate (and print?) potential, and then quit.
c* If read-in value of  NLEV1 < 0 , and program (attempts to) find all 
c  vibrational levels of potential-1 up to  v = |NLEV1|.  [This case
c  assumes  AUTO1 > 0.]
c** If (AUTO1.gt.0) read in only (v,J) quantum numbers of NLEV1 desired 
c  level & subroutine ALF tries to locate them (normal preferred case).
c   If (AUTO1.le.0) also read in trial energy for each level.  In this
c      case, the  NLEV.le.0  option does not work.
c   If (AUTO1.le.0) and vib. quant. No. IV < 0, seek level nearest to 
c      given trial energy but for whatever q. number shows up
c** If(LCDC.gt.0) calculate centrifugal distortion constants for each 
c  level via the Tellinghuisen implementation of Hutson's method.
c  If(LCDC.GE.2) also write them compactly to channel-9.
c** IF(LXPCT=0) calculate no expectation values or matrix elements.
c* IF(LXPCT = -1) only calculate and write compactly to channel-7 the
c  eigenvalues and level widths.
c* IF(LXPCT= 1,2 or -2) calculate expectation values, or if |LXPCT| > 2
c  the off-diagonal matrix elements, of powers of the distance 
c  coordinate or radial function defined by parameters IRFN & RREF.  
c* For   LXPCT > 0  write all these results to channel-6;  otherwise
c                  supress most such printing to channel-6.
c* For  |LXPCT| = 2  write eigenvalues and expectation values in 
c          	        compact form on channel-7.
c* For  |LXPCT| > 2  calculate matrix elements coupling each level 
c  to all (up to VIBMX) preceeding levels of the same potential (for 
c  NUMPOT.le.1), or to NLEV2 (see below) vib. levels of potential-2
c  (for NUMPOT.ge.2), and if (|LXPCT| > 3) write the overall 
c  off-diagonal matrix elements on channel-8. 
c* For  |LXPCT| > 4  also write to channel-7 the matrix elements of the
c  individual powers of the chosen distance coordinate (or radial fx.)
c* For  |LXPCT| > 5  WRITE(7,xx) only those matrix element components.
c** IF(NJM > 0), for each vibrational level, calculate all rotational
c  levels up to  J=NJM  or predissociation, whichever comes first.
c  Note that  AUTO1.le.0  forces  NJM= 0
c** When (NJM.GT.0) increase J in increments of JDJR.
c** IF(IWR.NE.0) print error & warning descriptions
c  IF(IWR.GE.1) also print final eigenvalues & node count.
c  IF(IWR.GE.2) also show end-of-range wave function amplitudes
c  IF(IWR.GE.3) print also intermediate trial eigenvalues, etc.
c** IF(LPRWF.GT.0) print wave function every LPRWF-th  point.
c** IF(LPRWF.LT.0) compactly write to channel-7 every |LPRWF|-th 
c  wave function value.  **  The first line of this file identifies the level, 
c  gives the position of 1-st point and radial mesh, & states No. of  points
c=======================================================================
c** INNOD1 specifies wave fx. initiation at RMIN.  Normal case of
c  INNOD1 > 0  gives initiation with wave fx. node @ RMIN.
c  INNOD1.le.0  give initiation with  zero slope @ RMIN.  This determines
c    symmetric eigenfunctions for rare special case when input potential
c    is half of a precisely symmetric potential with mid-point at RMIN.
c-----------------------------------------------------------------------
c ... comment this out if use first version of following READ statement
      INNOD1= 1
      INNOD2= INNOD1
c-----------------------------------------------------------------------
c     READ(5,*) NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF, INNOD1
c-----------------------------------------------------------------------
      READ(5,*) NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF
c-----------------------------------------------------------------------
c** SINNER specifies whether wave function matching occurs at outermost 
c   (SINNER.le.0) or innermost well turning point, to facilitate finding
c   inner vs. outer wells of a double well potential; Normally controlled 
c    automatically, 
      SINNER= 0
      INNER= SINNER
      IF(INNOD1.GT.0) WRITE(6,686) 1
      IF(INNOD1.LE.0) WRITE(6,688) 1
      IF(JDJR.LE.0) JDJR=1
ccc   IF(AUTO1.LE.0) NJM= 0 
      WRITE(6,612) EPS
      NLEV= NLEV1
      IF(NLEV1.LE.0) NLEV= 1
      SOMEG1= IOMEG1**2
      IF(IOMEG1.LT.99) WRITE(6,608) 1,IOMEG1,SOMEG1
      IF(IOMEG1.GE.99) WRITE(6,609)
      VMAX1= 0
c** Read the vibrational & rotational quantum numbers IV(i) & IJ(i) [and
c  if AUTO1.le.0 also trial energy GV(I)] of the NLEV levels to be found
c** For  IV(i)  values < -10,  SCHRQ  imposes a hard wall boundary
c  condition (i.e., a node) at mesh point # |-IV(i)| .
c-----------------------------------------------------------------------
      IF(AUTO1.GT.0) READ(5,*) (IV(I), IJ(I), I= 1,NLEV)
      IF(AUTO1.LE.0) READ(5,*) (IV(I), IJ(I), GV(I), I= 1,NLEV)
c-----------------------------------------------------------------------
      IF(NLEV1.GT.0) THEN
          IF(AUTO1.GT.0) WRITE(6,607) NLEV,(IV(I),IJ(I),I=1,NLEV)
          IF(AUTO1.LE.0) THEN
              WRITE(6,6607) NLEV,(IV(I),IJ(I),GV(I),I=1,NLEV)
              DO  I= 1,NLEV1
                  ZK1(IV(I),0)= GV(I)
                  ENDDO
              ENDIF
          DO I= 1,NLEV
              VMAX1= MAX(VMAX1,IV(I))
              ENDDO
          JREF= 0
        ELSE
          IF(NLEV1.LT.-399) NLEV1= -399
          VMAX1= -NLEV1
          NLEV= VMAX1+ 1
          WRITE(6,625) IJ(1),NLEV
          JREF= IJ(1)
          DO  I= 1,NLEV
              IV(I)= I-1
              IJ(I)= JREF
              ENDDO
        ENDIF
      IF(NJM.GT.IJ(1)) WRITE(6,638) JDJR,NJM
      IF(LCDC.GT.1)  WRITE(9,901) TITL
      IF(LXPCT.EQ.-1) WRITE(7,723) TITL
c** MORDR is the highest power of the radial function (or distance 
c  coordinate whose expectation values or matrix elements are to be 
c  calculated.  Program currently dimensioned for (MORDR.LE.10).  To
c  calculate only F-C Factors (when LXPCT>2), set  MORDR = -1.
c** IRFN & RREF specify the definition of the radial function or 
c  distance coordinate  RFN(R), powers of which are averaged over in 
c  expectation value or matrix element calculations.
c* If(IRFN .le. -10) utilize the USER-CHOSEN and CODED radial function
c                 generated in Lines #500-504 (below)
c* If(IRFN = -4)  the function is a power series in  R  premultiplying a
c          first derivative operator acting on the wavefx of Potential-2
c* If(IRFN = -3)  the function is the inverse power   1/R**3
c* If(IRFN = -2)  the function is the inverse power   1/R**2
c* If(IRFN = -1)  the function is the Dunham coordinate  X=(R-RREF)/RREF
c* If(IRFN =  0)  the function  RFN(R)  is the distance  R  itself.
c* If(IRFN = 1-9)  use the Surkus-type variable  
c                 X=(R^p - RREF^p)/(R^p + RREF^p)  where  p= IRFN
c* For  IRFN = -1 or 1-9,  if  RREF.gt.0  the read-in RREF value is the
c  reference length used to define the distance coordinate, while 
c  if  RREF.le.0  determine the value of this reference length by 
c  requiring that the expectation value  X**1  of the distance 
c  coordinate for the first level considered be identically zero.
c* IF(IRFN > 10) define  RFN(R)   by reading in, interpolating over (and
c  extrapolating beyond) IRFN read-in values of some known radial 
c  (transition moment) function, whose asymptotic value is RREF.  Do 
c  this in using the same read statements and GENINT subroutine calls 
c  used for generating a numerical potential.
      IF((LXPCT.NE.0).AND.(LXPCT.NE.-1)) THEN
c-----------------------------------------------------------------------
          READ(5,*) MORDR, IRFN, RREF
c-----------------------------------------------------------------------
          IF(MORDR.GT.MORDRMX) MORDR= MORDRMX
          IF(IABS(LXPCT).EQ.2) WRITE(7,724) TITL,MORDR
          IF((IABS(LXPCT).EQ.4).OR.(IABS(LXPCT).EQ.5)) WRITE(8,824) TITL
ccc       IF((IABS(LXPCT).EQ.4).OR.(IABS(LXPCT).EQ.5)) WRITE(11,811)TITL
          IF(IABS(LXPCT).GE.5) WRITE(7,725) TITL,MORDR
          IF(IABS(IRFN).GE.10) THEN
              MORDR= 1
              DM(0)= 0.d0
              DM(1)= 1.d0
            ELSE
              IF(MORDR.GE.0) THEN
c** Overall calculated matrix elements are for a power series in the
c  radial function  RFN(i)  (specified by IRFN & RREF), so must read
c  coefficients  DM(J)  of this power series.
c-----------------------------------------------------------------------
                  READ(5,*) (DM(J), J= 0,MORDR)
c-----------------------------------------------------------------------
                ELSE
                  DO  I= 1,NPP
                      RFN(I)= 1.D0
                      ENDDO
                  IF(MORDR.LT.0) WRITE(6,617)
                ENDIF 
            ENDIF
c** Define radial function (distance coordinate) operator  RFN(R)  for 
c  expectation values or matrix elements.
c** First ... for matrix elements of an operator consisting of a power 
c    series in  R  premultiplying the radial derivative of the wavefx.
          IF(IRFN.EQ.-4) WRITE(6,650) MORDR
          IF(MORDR.GT.0) THEN
c** If  RFN(R)  is the distance itself ...	
              IF(IRFN.EQ.0) WRITE(6,614)
              IF((IRFN.EQ.-2).OR.(IRFN.EQ.-3)) THEN
              IF((IRFN.EQ.0).OR.(IRFN.EQ.-2).OR.(IRFN.EQ.-3)) RREF=0.D0
c** If  RFN(R)  is   1/(distance)**|IRFN|  ....
                  J= -IRFN
                  WRITE(6,616) -IRFN
                  DO  I= 1,NPP
                      RFN(I)= 1.d0/RFN(I)**J
                      ENDDO
                  ENDIF
c%% Any other user-defined matrix element argument radial function 
c   may be introduced to the code here, and invoked by:  IRFN= -4 
c  Note that the existing  RFN(i)  array is the radial distances  R .
              IF(IRFN.LE.-10) THEN
c&&
c&& Illustrative user-defined analysis RFN(R) function
c&&               WRITE(6,*) 'Print description of function introduced'
c&&               WRITE(6,*) 'Use Freedman Pade DMF for CO'
c&&               DO  I= 1,NPP
c%%---------------------------------------------------------------------
c%%                  RFN(I)= {calculate users chosen radial function}
c%% Freedman's DMF for CO  ---------------------------------------------
c%%                  RFN(I)= {calculate users chosen radial function}
c&&   data coeff_new /-24.6005858d0,-109.5939637d0,-524.8233323d0,
c&&  +                 4.5194090d0,19.7954955d0,
c&&  +                 6.6011985d0,19.7206690d0/
c&&   dm = -0.122706d0*(1.+coeff(1)*x+coeff(2)*x*x+coeff(3)*x**3)/
c&&  +                   (1.+coeff(4)*x+coeff(5)*x*x+coeff(6)*x**3
c&&  +                    + coeff(7)*x**6)
c&&---------------------------------------------------------------------
c&&                   XX= RFN(I)/1.128322714d0 - 1.d0
c&&                   RFN(I)= -0.122706d0*(1.d0+ XX*(-24.6005858d0 
c&&  1                  + XX*(-109.5939637d0 + XX*(-524.8233323d0))))/
c&&  2                    (1.d0 + XX*(4.5194090d0 + XX*(19.7954955d0 +
c&&  3                        XX*(6.6011985d0 + 19.7206690d0*XX**3))))
c&&                   ENDDO
                  ENDIF
              IF((IRFN.EQ.-1).OR.((IRFN.GE.1).AND.(IRFN.LE.9))) THEN
c** If  RFN(R)  is the Dunham or Surkus-type distance coordinate 
                  IF(IRFN.EQ.-1) WRITE(6,615)
                  IF((IRFN.GE.1).AND.(IRFN.LE.9)) WRITE(6,611) 
     1                                             IRFN,IRFN,IRFN,IRFN
                  IF(RREF.GT.0.D0) THEN
                      RREFP= RREF**IRFN
                      WRITE(6,613) RREF
                      DO  I=1,NPP
                          XX= RMINN+I*RH
                          IF(IRFN.EQ.-1) RFN(I)= (XX- RREF)/RREF
                          IF(IRFN.GE.1) RFN(I)= (XX**IRFN- RREFP)
     1                                        /(XX**IRFN + RREFP)
                          ENDDO
                    ELSE
                      WRITE(6,610)
                    ENDIF
                  ENDIF
c
cQQQQQQQQQQQ new ... the following not fully tested ...
c
c** If  RFN(R)  is defined by interpolating over read-in points, use
c  potential generating routine to do interpolation/extrapolation.
              IF(IRFN.GE.10) THEN
                  MORDR= 1
                  DM(0)= 0.d0
                  DM(1)= 1.d0
                  WRITE(6,603) 
c** If the expectation value/matrix element radial function argument to
c   be defined by interpolating/extrapolating over read-in points, then 
c   read input analogous to that for a pointwise potential, and then call 
c   interpolation/extrapolation routine GENINT (from PREPOT package)
c* NRFN is the number of points [XIF(i),YIF(i)] to be read in
c* RFLIM  is the limiting asymptotic value imposed on the extrapolation
c* Interpolate with NUSEF-point piecewise polynomials (or splines for
c    NUSEF.le.0), which are extrapolated to the asymptote as specified by
c    parameters ILRF, NCNF & CNNF (see read #20).
c* RFACTF - factor converts read-in distances XIF(i) to angstroms
c* MFACTF - factor converts read-in moment values YIF(i) to debye.
c-----------------------------------------------------------------------
                  READ(5,*) NRFN, RFLIM
                  READ(5,*) NUSEF, ILRF, NCNF, CNNF
                  READ(5,*) RFACTF, MFACTF
                  READ(5,*) (XIF(I), YIF(I), I= 1,NRFN)
c-----------------------------------------------------------------------
                  WRITE(6,810) NRFN, RFLIM
                  IF(NUSEF.GT.0) WRITE(6,812) NUSEF, NRFN
                  IF(NUSEF.LE.0) WRITE(6,814) NRFN
                  IF((ILRF.GT.1).AND.(DABS(CNNF).GT.0.D0))
     1                                         WRITE(6,816) CNNF, NCNF
                  WRITE(6,818) RFACTF, MFACTF
                  NROW= (NRFN+ 2)/3
                  DO  J= 1,NROW
                      WRITE(6,820) (XIF(I), YIF(I), I=J, NRFN, NROW)
                      ENDDO
                  DO  I= 1,NRFN
                      XIF(I)= XIF(I)*RFACTF
                      YIF(I)= YIF(I)*MFACTF
                      ENDDO
  810 FORMAT(' Transition moment function defined by interpolating over'
     1  ,I4,' read-in points'/5x,'and approaching the asymptotic value',
     2  f12.6)
  812 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  814 FORMAT(' Perform cubic spline interpolation over the',I5,' input p
     1oints' )
  816 FORMAT('- Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
  818 FORMAT(' Scale input points:  (distance)*',1PD16.9,'   &  (moment)
     1*',D16.9/4x,'to get required units  [Angstroms & debye]'/
     3  3('      R(i)         Y(i)  ')/3(3X,11('--')))
  820 FORMAT((3(F12.6,F13.6)))
                  IR2F= 0
                  CALL GENINT(LRPT,NPP,WF1,RFN,NUSEF,IR2F,NRFN,XIF,YIF,
     1                                           RFLIM,ILRF,NCNF,CNNF)
                  ENDIF
              ENDIF
          IF((MORDR.GE.0).AND.(IABS(IRFN).LE.9)) 
     1                                  WRITE(6,602) (DM(J),J=0,MORDR)
          ENDIF
c** For matrix element calculation, couple each level of potential-1 to
c  up to (see below) NLEV2 other vibrational levels, subject to 
c  rotational selection rules:  DELTA(J)= J2DL to J2DU with increment
c  J2DD (e.g., -1,+1,+2 for P- and R-branches).
c** If (AUTO2.gt.0) read in only (v,J) quantum numbers of desired levels
c  and trust subroutine ALF to locate them (normal preferred case).
c   If (AUTO2.le.0) also read in a trial pure vib energy for each level.
c* For the one-potential case (NUMPOT.LE.1), automatically truncate to
c  avoid redundancy and consider only emission into these NLEV2 levels.
c* Trial level energies are generated internally.
c**  IV2(i) are the vibrational quantum numbers of the Potential-2
c  levels for which matrix elements are desired.
c**  ZK(IV2(i),0) are the associated pure vibrational trial energies 
c  (which are only read in if AUTO2.le.0!)
c=======================================================================
c** INNOD2 specified wave fx. initiation at RMIN.  Normal case of
c  INNOD2 > 0  gives initiation with wave fx. node @ RMIN.
c  INNOD2.le.0  give initiation with  zero slope @ RMIN.  This determines
c    symmetric eigenfunctions for rare special case when input potential
c    is half of a precisely symmetric potential with mid-point at RMIN.
c=======================================================================
      IF(IABS(LXPCT).GE.3) THEN
c-----------------------------------------------------------------------
ccc       READ(5,*) NLEV2, AUTO2, J2DL, J2DU, J2DD, INNOD2
c-----------------------------------------------------------------------
          READ(5,*) NLEV2, AUTO2, J2DL, J2DU, J2DD
c-----------------------------------------------------------------------
          IF(NLEV2.GT.VIBMX) NLEV2= VIBMX
          IF(NLEV2.LE.0) THEN
              WRITE(6,644) NLEV2
              STOP
              ENDIF
c----------------------------------------------------------------------
          IF(AUTO2.GT.0) READ(5,*) (IV2(I), I= 1,NLEV2)
          IF(AUTO2.LE.0) THEN
              READ(5,*) (IV2(I), ZK2(I,1), I= 1,NLEV2)
c----------------------------------------------------------------------
c** Give potential-2 trial energy the correct vibrational label
              DO  I= 1,NLEV2
                  ZK2(IV2(I),0)= ZK2(I,1)
                  ENDDO
              ENDIF
          IF(NUMPOT.GT.1) THEN
              IF(INNOD2.GT.0) WRITE(6,686) 2
              IF(INNOD2.LE.0) WRITE(6,688) 2
              ENDIF
          VMAX2= 0
          DO  ILEV2= 1,NLEV2
              VMAX2= MAX(VMAX2,IV2(ILEV2))
              ENDDO
          IF(MORDR.LT.0) DM(1)= 1.d0
          SOMEG2= IOMEG2**2
          IF(J2DD.EQ.0) J2DD= 1
          IF(AUTO2.GT.0) WRITE(6,634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),
     1                                                     I= 1,NLEV2)
          IF(AUTO2.LE.0) WRITE(6,6634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),
     1                                      ZK2(IV2(I),0), I= 1,NLEV2)
          IF(NUMPOT.GE.2) THEN
              IF(IOMEG2.LT.99) WRITE(6,608) 2,IOMEG2,SOMEG2
              IF(IOMEG2.GE.99) WRITE(6,609)
              ENDIF
          WRITE(6,632)
          ENDIF
c
      IF(AUTO1.GT.0) THEN
c** If using automatic search for desired levels, subroutine ALF gets
c  eigenvalues ZK1(v,0) for desired vibrational levels of Potential-1,
c  centrifugally-distorted to J=JREF.
          EJREF= JREF*(JREF+1)*RH**2
          DO  I= 1,NPP
              VJ(I)= V1(I)+EJREF*RM2(I)
              ENDDO
          IF((NLEV1.EQ.1).AND.(IV(1).gt.998)) THEN
c** Option to search for the very highest level (within 0.001 of Disoc)
              EO= VLIM1- 0.001d0
              KV= IV(1)
              CALL SCHRQ(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,
     1                   RMIN,RH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
              IV(1)= KV
              IF(KV.GE.0) THEN
                  WRITE(6,622) IJ(1),KV,VLIM1-EO
                  GV(KV)= EO
                  VMAX1= KV
                ELSE
                  WRITE(6,626) J, 0.001d0
                  STOP
                ENDIF
            ELSE
              VMAX= VMAX1
              AFLAG= JREF
              IF((IABS(LXPCT).GT.2).AND.(NUMPOT.EQ.1)) VMAX=
     1                                                MAX(VMAX1,VMAX2)
              CALL ALF(NPP,RMIN,RH,VJ,WF1,VLIM1,VMAX,AFLAG,ZMU,EPS,NCN1,
     1                                            GV,INNOD1,INNR1,IWR)
              VMAX1= VMAX
            ENDIF
c** Get band constants for v=0-VMAX1 for generating trial eigenvalues
          WARN= 0
          DO  ILEV1= 0,VMAX
              KV= ILEV1
              EO= GV(KV)
              INNER= INNR1(KV)
              CALL SCHRQ(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,
     1                 RMIN,RH,NPP,NBEG,NEND,INNOD1,INNER,WARN,NoPRWF)
              CALL CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,VJ,WF1,RM2,RCNST)
              IF(NLEV1.LT.0) THEN
                  IV(ILEV1+1)= KV
                  IJ(ILEV1+1)= JREF
                  ENDIF
              ZK1(ILEV1,0)= GV(ILEV1)
              DO  M= 1,7
                  ZK1(ILEV1,M)= RCNST(M)
                  ENDDO
              ENDDO
          ENDIF
      IF(IABS(LXPCT).GT.2) THEN
          IF(AUTO2.GT.0) THEN
c** If using automatic location for levels of potential-2 (AUTO2 > 0)
c  for matrix element calculation, also need Potential-2 band constants
c  (rotational energy derivatives) ... again, calculate them at J=JREF
              IF(NUMPOT.GT.1) THEN
                  AFLAG= JREF
                  DO  I= 1,NPP
                      VJ(I)= V2(I)+EJREF*RM2(I)
                      ENDDO
                  CALL ALF(NPP,RMIN,RH,VJ,WF1,VLIM2,VMAX2,AFLAG,ZMU,
     1                                   EPS,NCN2,GV,INNOD2,INNR2,IWR)
                  ENDIF
              ENDIF 
          DO  ILEV2= 1,NLEV2
              IF(NUMPOT.EQ.1) THEN
c** For matrix elements within a single potl., copy above band constants
                  DO  M= 0,7
                      ZK2(IV2(ILEV2),M)= ZK1(IV2(ILEV2),M)
                      ENDDO
                ELSE
c ... otherwise, generate them (as above) with SCHRQ & CDJOEL
                  KV= IV2(ILEV2)
                  IF(AUTO2.GT.0) EO= GV(KV)
                  IF(AUTO2.LE.0) EO= ZK2(KV,0) 
                  INNER= INNR2(KV)
                  CALL SCHRQ(KV,JREF,EO,GAMA,PMAX2,VLIM2,VJ,WF1,
     1        BFCT,EPS,RMIN,RH,NPP,NBEG,NEND,INNOD2,INNER,WARN,NoPRWF)
                  CALL CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,VJ,WF1,
     1                                                      RM2,RCNST)
                  ZK2(IV2(ILEV2),0)= EO
                  DO  M= 1,7
                      ZK2(IV2(ILEV2),M)= RCNST(M)
                      ENDDO
                ENDIF
              ENDDO
          ENDIF 
      WARN= 1
      EJREF= EJREF/RH**2
      IF(NLEV1.LE.0) NLEV= VMAX1+1
c
c===== Begin Actual Potential-1 Eigenvalue Calculation Loop Here =======
c** Loop to compute eigenvalues ... etc. for NLEV levels of Potential-1
      DO 190 ILEV1= 1,NLEV
          KV= IV(ILEV1)
          NJMM= MAX(NJM,IJ(ILEV1))
          JROT= IJ(ILEV1)- JDJR
          IQT= 0
          JCT= 0
c** If NJM > IJ(ILEV1) loop over range of rotational levels too
          DO  JLEV= IJ(ILEV1),NJMM,JDJR
              JROT= JROT+ JDJR
              EJ= JROT*(JROT+1)-SOMEG1
              IF(IOMEG1.GE.99) EJ= JROT*JROT-0.25D0
c** If appropriate (AUTO1>0) use ALF results to generate trial eigenvalue
              IF(AUTO1.GT.0) THEN
                  EO= ZK1(KV,0)
                  DEJ= EJ- EJREF
                  EJP= 1.d0
                  DO M= 1,7
                      EJP= EJP*DEJ
                      EO= EO+ EJP*ZK1(KV,M)
                      ENDDO
                ELSE
                  EO= ZK1(IV(ILEV1),0)
                ENDIF
c ... or if JLEV > IJ(ILEV1) ... use local Beff to estimate next level
              IF(JLEV.GT.IJ(ILEV1)) THEN
                  BEFF= 0.d0
                  DO  I= NBEG,NEND
                      BEFF= BEFF+ WF1(I)**2*RM2(I)
                      ENDDO
                  BEFF= BEFF*RH*BvWN
                  EO=  ESLJ(JCT)+ (2*JLEV+ 1- JDJR)*JDJR*BEFF
                  ENDIF
c** Now add centrifugal term to get effective (radial) potential
              EJ= EJ*RH**2
              DO  J=1,NPP
                  VJ(J)=V1(J)+EJ*RM2(J)
                  ENDDO
c** Set wall outer boundary condition, if specified by input IV(ILEV1)
              IF(KV.LT.-10) THEN
                  WF1(-IV(ILEV1))= 0.D0
                  WF1(-IV(ILEV1)-1)= -1.D0
                  ENDIF
              KVIN= KV
              INNER= INNR1(KV)
              IF(SINNER.NE.0) INNER= SINNER
c** Call SCHRQ to find Potential-1 eigenvalue EO and eigenfn. WF1(i)
  100         CALL SCHRQ(KV,JROT,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,
     1                   RMIN,RH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
              IF(KV.LT.0) THEN
c** SCHRQ  error condition is  (KV.LT.0) .
                  IF(NJM.GT.IJ(ILEV1)) THEN
c ... in automatic search for ever-higher J levels
                      IF(IQT.LE.0) THEN
c ... try one more time with E(trial) slightly below barrier maximum
                          IQT= 1
                          EO= PMAX1- 0.1d0
                          GO TO 100
                        ELSE
                          KV= KVIN
                          GO TO 130
                        ENDIF
                      ENDIF
                  GO TO 122
                  ENDIF 
              IF((KV.NE.KVIN).AND.
     1                        ((AUTO1.GT.0))) THEN
c             IF(KV.NE.KVIN) THEN
c** If got wrong vib level, do a brute force ALF calculation to find it.
                  KV= KVIN
                  AFLAG= JROT
                  CALL ALF(NPP,RMIN,RH,VJ,WF1,VLIM1,KV,AFLAG,ZMU,
     1                                   EPS,NCN1,GV,INNOD1,INNR1,IWR)
                  IF(KV.EQ.KVIN) THEN 
                      EO= GV(KVIN)
                      GO TO 100
                    ELSE
                      WRITE(6,618) KVIN,JROT,KV
                      KV= KVIN
                      GO TO 130
                    ENDIF
                  ENDIF
              if(kv.ne.iv(ilev1)) iv(ilev1)= KV
c** If desired, calculate rotational & centrifugal distortion constants
              IF(LCDC.GT.0) THEN
                  IF((IOMEG1.GT.0).AND.(JROT.EQ.0)) THEN
c** Calculate 'true' rotational constants for rotationless IOMEG>0 case
                      CALL SCHRQ(KV,0,EO,GAMA,PMAX1,VLIM1,V1,WF1,BFCT,
     1               EPS,RMIN,RH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
                      CALL CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,V1,WF1,
     1                                                      RM2,RCNST)
                    ELSE
c** Calculate rotational constants for actual (v,J) level of interest.
                      CALL CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,VJ,WF1,
     1                                                      RM2,RCNST)
                    ENDIF
                  WRITE(6,606) KV,JROT,EO,(RCNST(M),M=1,7)
                  IF(LCDC.GT.1) WRITE(9,902) KV,JROT,EO,(RCNST(M),M=1,7)
                  ENDIF
              IF(LXPCT.EQ.-1)  WRITE(7,703) KV,JROT,EO,GAMA
              IF(((LXPCT.EQ.1).OR.(IABS(LXPCT).EQ.2)).OR.
     1                  ((IABS(LXPCT).GT.2).AND.((IRFN.EQ.-1).OR.
     2          (IRFN.GE.1).AND.(IRFN.GE.9)).AND.(RREF.LE.0.d0))) THEN
c** Calculate various expectation values in LEVXPC 
                  CALL LEVXPC(KV,JROT,EO,GAMA,NPP,WF1,VJ,VLIM1,RFN,RMIN,
     1                       RH,RREF,NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
                  IF((LXPCT.GT.0).AND.(MORDR.GT.0)) WRITE(6,632)
                  ENDIF
              IF((IABS(LXPCT).LE.2).OR.(NLEV2.LE.0)) GO TO 122
c** If desired, now calculate off-diagonal matrix elements, either
c  between levels of different potentials, IF(NUMPOT.GE.2), or between
c  levels of a single potential, for (NUMPOT.LE.1).
c** First prepare centrifugally distorted potential, trial energy, etc.,
c  and calculate second wave function and matrix element(s)
              DO 120 ILEV2= 1,NLEV2
c** For case of a single potential, avoid redundancy by considering 
c  only emission
                  IF((NUMPOT.LE.1).AND.(IV2(ILEV2).GT.KV)) GO TO 120
c** Loop over J2's allowed by given selection rule.
                  DO 116 IJD= J2DL,J2DU,J2DD
                      KV2= IV2(ILEV2)
                      KVIN= KV2
                      JROT2= JROT+IJD
                      IF(JROT2.LT.0) GO TO 116
                      IF((NUMPOT.LE.1).AND.(IV2(ILEV2).EQ.KV).AND.
     1                                      (JROT2.GT.JROT)) GO TO 116
                      EJ2= JROT2*(JROT2+1)- SOMEG2
                      IF(IOMEG2.GE.99) EJ2=JROT2**2-0.25D0
                      EO2= ZK2(KV2,0)
                      DEJ= EJ2- EJREF
                      EJP= 1.d0
c** Use calculated state-2 CDC's to predict trial eigenvalue
                      DO  M= 1,7
                          EJP= EJP*DEJ
                          EO2= EO2+ EJP*ZK2(KV2,M)
                          ENDDO
c** Now ... update to appropriate centrifugally distorted potential
                      EJ2= EJ2*RH*RH
                      DO  I=1,NPP
                          VJ(I)= V2(I)+ EJ2*RM22(I)
                          ENDDO
                      INNER= INNR2(KV2)
                      IF(SINNER.NE.0) INNER= SINNER
                      ICOR= 0
  110                 CALL SCHRQ(KV2,JROT2,EO2,GAMA,PMAX2,VLIM2,VJ,WF2,
     1        BFCT,EPS,RMIN,RH,NPP,NBEG2,NEND2,INNOD2,INNER,IWR,LPRWF)
                      IF(KV2.NE.KVIN) THEN
                          IF(KV2.LT.0) GO TO 114
c** If Using CDC's to estimate trial eigenvalue failed:
                          ICOR= ICOR+1
                          IF(ICOR.LE.2) THEN
c ... first correction attempt ... use semiclassical dv/dE to improve
                              GB= -1.d0
                              GI= -1.d0
                              WV= 0.d0
                              XX= EO2*BFCT
                              DO  I= NBEG2,NEND2
                                  GBB= GB
                                  GB= GI
                                  GI= XX-VJ(I)
                                  IF((GBB.GT.0.d0).AND.(GI.GT.0.d0))
     1                                          WV= WV+ 1.d0/DSQRT(GB)
                                  ENDDO
                              WV= 6.2832d0/(BFCT*WV)
                              EO2= EO2+ WV*(KVIN- KV2)
                              GO TO 110
                              ENDIF
                          WRITE(6,633) IV2(ILEV2),JROT2,KV2
c ... if that fails, do a brute force ALF calculation to find it.
  114                     KV2= KVIN
                          AFLAG= JROT2
                          CALL ALF(NPP,RMIN,RH,VJ,WF2,VLIM2,KV2,AFLAG,
     1                               ZMU,EPS,NCN2,GV,INNOD2,INNR2,IWR)
                          IF(KV2.EQ.KVIN) THEN
                              EO2= GV(KV2)
                              INNER= INNR2(KV2)
                              GO TO 110
                            ELSE
                              WRITE(6,618) KVIN,JROT,KV2
                              GO TO 116
                            ENDIF
                          ENDIF
                      IF(NBEG.GT.NBEG2) NBEG2= NBEG
                      IF(NEND.LT.NEND2) NEND2= NEND
c                     IF((NUMPOT.LE.1).AND.(EO2.GT.(EO+EPS))) GO TO 120
                      CALL MATXEL(KV,JROT,IOMEG1,EO,KV2,JROT2,IOMEG2,
     1             IRFN,EO2,NBEG2,NEND2,LXPCT,MORDR,DM,RH,WF1,WF2,RFN)
  116                 CONTINUE
c** End of Potential-2 rotational selection level loop
  120             CONTINUE
c++++ End of Potential-2 vibrational level matrix element loop +++++++++
c
  122         CONTINUE
              JCT= JCT+1
c ... check to avoid array overflow
              IF(JCT.GT.VIBMX) THEN
                  WRITE(6,637)  VIBMX
                  STOP
                  ENDIF 
              JWR(JCT)= JROT
              ESLJ(JCT)= EO
              ENDDO
c++ End of Potential-1 loop over NJM-specified J-sublevels
  130     IF(NJM.GT.IJ(ILEV1)) THEN
c** Print rotational sublevels generated for vibrational level  ILEV
              NROW=(JCT+4)/5
              WRITE(6,627) KV
              DO  J=1,NROW
                  WRITE(6,628) (JWR(I),ESLJ(I),I=J,JCT,NROW)
                  ENDDO
              WRITE(6,641)
              ENDIF
          ESOLN(ILEV1)= ESLJ(1)
  190     CONTINUE
c++ End of loop over the NLEV Potential-1 input levels
      IF(NLEV1.LT.0) THEN
          NROW=(NLEV+3)/4
          WRITE(6,623) NLEV,IJ(1)
          DO  J=1,NROW
              WRITE(6,630) (IV(I),ESOLN(I),I=J,NLEV,NROW)
              ENDDO
          IF((NLEV.GT.1).AND.(IJ(1).EQ.0).AND.(NCN1.GT.0)
     1                               .AND.(ESOLN(NLEV).LT.VLIM1)) THEN
c** In (NLEV1 < 0) option, estimate vD using the N-D theory result that:
c     (vD - v) {is proportional to} (binding energy)**((NCN-2)/(2*NCN))
              VDMV=1.D0/(((VLIM1-ESOLN(NLEV-1))/
     1                         (VLIM1-ESOLN(NLEV)))**(1.D0/PW) - 1.D0)
c** Use empirical N-D Expression to predict number and (if there are
c  any) energies of missing levels
              VD= IV(NLEV)+VDMV
              IVD= INT(VD)
              IF(IVD.GE.VIBMX) IVD= VIBMX-1
              IVS= IV(NLEV)+1
              WRITE(6,620) NCN1,VD
              IF((IVD.GE.IVS).AND.(DBLE(IV(NLEV))/VD.GT.0.9d0)) THEN
                  NFP= NLEV+1
                  DO  I= IVS,IVD
                      NLEV= NLEV+1
                      IV(NLEV)= IV(NLEV-1)+1
                      ESOLN(NLEV)= VLIM1 - (VLIM1 - ESOLN(NLEV-1))*
     1                                          (1.D0 - 1.D0/VDMV)**PW
                      VDMV= VDMV-1.D0
                      ENDDO
                  NLP= NLEV-NFP+1
                  NROW= (NLP+3)/4
                  WRITE(6,621) NLP
                  DO  J= 1,NROW
                      III= NFP+J-1
                      WRITE(6,630) (IV(I),ESOLN(I),I= III,NLEV,NROW)
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF
      IF((NJM.LE.0).AND.(NLEV1.GE.0)) THEN
          NROW=  (NLEV+2)/3
          WRITE(6,619) NLEV
          DO  J= 1,NROW
              WRITE(6,631) (IV(I),IJ(I),ESOLN(I),I= J,NLEV,NROW)
              ENDDO
          ENDIF
      WRITE(6,601)
      GO TO 2
  999 STOP
c-------------------------------------------------------------------
  601 FORMAT(1x,79('=')////)
  602 FORMAT( ' Coefficients of expansion for radial matrix element/expe
     1ctation value argument:'/(5X,5(1PD14.6)))
  603 FORMAT(/' Expectation value/matrix element arguments are powers of
     1 a radial function'/5x,'defined by interpolating over read-in poin
     2ts'//' Transition moment function:'/1x,9('==='))
  604 FORMAT(' Integrate from  RMIN=',f7.3,'  to  RMAX=',f7.2,
     1  '  with mesh  RH=',f9.6,'(Angst)'//' Potential #1 for ',A2,'(',
     2  I3,')-',A2,'(',I3,')'/1x,32('='))
  605 FORMAT(/A78/40('==')/' Generate   ZMU=',F15.11,'(u)',
     1  '   &   BZ=',1PD16.9,'((1/cm-1)(1/Ang**2))'/
     2  10x,'from atomic masses:',0Pf16.11,'  & ',F16.11,'(u)')
  606 FORMAT(' E(v=',i3,', J=',i3,')=',f10.3,'   Bv=',F11.7,
     1  '  -Dv=',1PD12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,
     2  '   Mv=',D12.4,'   Nv=',D12.4,'   Ov=',D12.4)
  607 FORMAT(/' Solve for the',i4,' vibration-rotation levels of Potenti
     1al-1:'/'   (v,J) =',6('  (',i3,',',i3,')':)/(10x,6('  (',i3,',',
     2  i3,')':)))
 6607 FORMAT(/' Solve for',i4,' vibration-rotation levels of Potential-1
     1 using Trial energies:'/(2x,3('   E(',I3,',',I3,') =',
     2 F11.2:)))
  608 FORMAT(/' State-',I1,' electronic angular momentum  OMEGA=',I2/
     1  9x,'yields centrifugal potential  [J*(J+1) -',F5.2,']/r**2' )
  609 FORMAT('  Use centrifugal potential for rotation in two dimensions
     1:   (J**2 - 1/4)/r**2')
  610 FORMAT(5X, 'where RREF defined by requiring  <X**1> = 0  for first
     1 level considered')
  611 FORMAT(/' Matrix element argument expansion vble is   X = ((r^',
     1  i1,' - RREF^',i1,')/(r^',i1,' + RREF^',i1,'))')
  612 FORMAT(/' Eigenvalue convergence criterion is   EPS=',1PD8.1,
     1 '(cm-1)'/' Airy function at 3-rd turning point is quasibound oute
     2r boundary condition')
  613 FORMAT(5X,'where reference length is held fixed at   RREF =',
     1 F13.10,'(Angstroms)')
  614 FORMAT(/' Matrix element arguments are powers of the distance  r (
     1in Angstroms)')
  615 FORMAT(/' Matrix element argument expansion variable is:    X = (r
     1 - RREF)/RREF')
  616 FORMAT(/' Matrix element arguments are powers of the squared inver
     1se distance  X = 1/r**',i1)
  617 FORMAT(/' Matrix element argument is fixed as a constant = 1')
  618 FORMAT(' *** PROBLEM *** Searching for  v=',i3,' , J=',i3,
     1 '  ALF only found to  v=',i3)
  619 FORMAT(/' Find the',i4,' vibration-rotation levels:'/
     1  3('     v   J      E(v)   ')/3(2x,7('---')))
  620 FORMAT(/' An  n=',I2,'  N-D theory extrapolation from last 2 level
     1s implies   vD =',F8.3)
  621 FORMAT(5X,'with the',I4,' missing level(s) predicted to be:'/
     1  4('     v     E(v)   ')/4(4x,7('--')))
  622 FORMAT(/' Search for highest bound  J=',i3,'  level finds  E(v=',
     1  i3,') = VLIM -',1PD12.5/)
  623 FORMAT(/' Find',I4,' Potential-1 vibrational levels with  J=',i3/
     1  4('     v     E(v)   ')/4(4x,7('--')))
  624 FORMAT(4x,'Since the molecule is an ion with charge',SP,I3/6x,"use
     1 Watson's charge-adjusted reduced mass   mu = M1*M2/[M1 + M2 - (",
     2  i2,')*me]')
  625 FORMAT(' For  J=',i3,', try to find the first',i4,' vibrational le
     1vels of Potential-1')
  626 FORMAT(/' *** FAIL to find highest bound J=',i3,'  level from tria
     1l   E = VLIM -',1PD11.4)
  627 FORMAT(/' For vibrational level  v =',I3,'   of Potential-1'/
     1 1X,5('  J',6X,'E',7X)/1X,5(7('--'),2X))
  628 FORMAT((1X,5(I3,F11.3,2X)))
  630 FORMAT((4(I6,F12.4:)))
  631 FORMAT((3(I6,I4,F13.5:)))
  632 FORMAT(1X,79('-'))
  633 FORMAT(' **** Caution: Search for   v=',I3,'   J=',i3,
     1  '  on potential-2 actually found   v=',I3)
  634 FORMAT(/' Using the rotational selection rule:  delta(J)=',
     1 i3,' to',i2,' with increment',i2/'   calculate matrix elements fo
     2r coupling to the',I4,' vibrational levels of'/
     3 '   Potential-2:   v =',14I4:/(21x,14i4:))
 6634 FORMAT(/' Using the rotational selection rule:  delta(J)=',
     1 i3,' to',i2,' with increment',i2/'   calculate matrix elements fo
     2r coupling to the',I4,' vibrational levels of'/'   Potential-2 usi
     3ng trial energies:',2('   E(',I3,')=',F9.2:)/4('   E(',I3,')=',
     4 F9.2:))
  635 FORMAT(/' Get matrix elements between levels of Potential-1 (above
     1) & Potential-2 (below)'/1X,39('--')/' For Potential #2:'/
     2  1x,17('='))
  636 FORMAT(/' Calculate properties of the single potential described a
     1bove')
  637 FORMAT(/' *** Array Dimension OVERFLOW ***   (Number of J sublevel
     1s) > VIBMX=',i4)
  638 FORMAT('   and automatically increment  J  in steps of',i3, ' to a
     1 maximum value of',i4)
  641 FORMAT(1X,39('++'))
  644 FORMAT(/' *** Input data ERROR *** matrix element calculation need
     1s  NLEV2=',i3,' > 0')
  650 FORMAT(/' Matrix element argument is radial first derivative opera
     1tor premultiplied by'/5x,'a power series in  r  of order',i3)
  686 FORMAT(' Potential-',i1,' uses inner boundary condition of  zero v
     1alue  at  RMIN')
  688 FORMAT(' Potential-',i1,' uses symmetric-well inner boundary condi
     1tion  of zero slope at RMIN')
  703 FORMAT(1X,I4,I5,F13.4,G13.5)
  723 FORMAT(/A78/1x,'Output values of:  v, J, E & (Level Width)')
  724 FORMAT(//A78//'   v   J    E(v,J)     Width       <KE>',
     1  6x,'<M(r)>  &  <XI**k>  for k=1 to',i3/2x,38('=='))
  725 FORMAT(//A78//"   v'  J'",'  v"  J"     FREQ',"    <v',J'| XI**k",
     1  ' |v",J">  for  k=0  to  MORDR=',i2/2x,37('=='))
  824 FORMAT(//A78/30('==')/" Note that (v',J') &",' (v",J") strictly la
     1bel the upper and lower levels, resp.,'/6x,'and  E(lower)=E"'/
     2 ' but  E(2)-E(1)  is:  (energy of State-2 level) - (energy of Sta
     3te-1 level)'//12x,'Band'/' dJ(J")',4x,7hv'   v",'  E(lower)  E(2)-
     4E(1) A(Einstein)  F-C Factor  ',13h<v'j'|M|v"j"> /
     5 1x,3('--'),('   -------'),'  --------',3x,
     6 3('---'),2x,11('-'),2x,11('-'),2x,11('-') )
cc811 FORMAT(//A78/30('==')//12x,"   v'","  J'",'    v"','  J"',
cc   1 '   position    E(upper)    E(lower)',16h   <v'j'|M|v"j">/
cc   2 1x,68('-') )
  901 FORMAT(//A78/1x,62('==')/'   v    J',7x,'E',10x,'Bv',11x,'-Dv',
     1  13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'/1x,62('=='))
  902 FORMAT(I4,I5,f12.4,f14.10,6(1PD15.7))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE LEVXPC(KV,JR,EPR,GAMA,NPP,WF,V,VLIM,RFN,RMIN,RH,RREF,
     1      NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculates expectation values of the kinetic energy and of X**IP 
c  (IP=1,MORDR), denoted XPTKE and XPCTR(IP), respectively, for level  
c  v=KV, J=JR, E=EPR(cm-1), using wave function WF(i), (i=NBEG,NEND).
c** Assumes units of length are (Angstroms) .
c** Division by BFCT converts potential V(I) to units (cm-1).
c** If (|LXPCT| = 2  or  4), "punch" (WRITE(7,XXX)) results to channel-7
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER I,K,IRFN,IPNCH,ITRY,JR,KV,LXPCT,
     1    NPP,NBEG,NEND,MORDR
      REAL*8  WF(NPP),V(NPP),RFN(NPP),XPCTR(0:11),DM(0:20)
      REAL*8 BFCT,DS,DRT,DMR,DER,EPR,EINN,GAMA,RMIN,RMINN,RH,RREF,
     1  RR,RXPCT,SS2,SF2,VLIM,XPTKE
c
      EINN= BFCT*EPR
      IPNCH=0
      IF((IABS(LXPCT).EQ.2).OR.(IABS(LXPCT).GE.4)) IPNCH=1
c** MORDR is the highest-power expectation value considered.
      IF(MORDR.GT.11) MORDR=11
      ITRY=20
      IF(((IRFN.EQ.-1).OR.((IRFN.GE.1).AND.(IRFN.LE.9)))
     1                                    .AND. (RREF.LE.0.D0)) ITRY=0
c** Start by calculating contributions at end points
    2 SS2=WF(NBEG)**2
      SF2=WF(NEND)**2
      XPTKE= 0.5D0*(SS2*(EINN-V(NBEG)) + SF2*(EINN-V(NEND)))
      IF(MORDR.GT.0) THEN
          XPCTR(0)= 1.d0/RH
          DO  K=1,MORDR
              SS2=SS2*RFN(NBEG)
              SF2=SF2*RFN(NEND)
              XPCTR(K)=0.5D0*(SS2+SF2)
              ENDDO
          ENDIF
      IF(IRFN.GT.-4) THEN
c** For regular expectation values of a radial function ...
          DO  I=NBEG+1,NEND-1
              DS=WF(I)**2
              XPTKE= XPTKE+ DS*(EINN-V(I))
              IF(MORDR.GT.0) THEN
                  RR= RFN(I)
                  DO  K=1,MORDR
                      DS=DS*RR
                      XPCTR(K)=XPCTR(K)+DS
                      ENDDO
                  ENDIF
              ENDDO
        ELSE 
c** For expectation values involving partial derivative operator ...
          DO  K= 0,MORDR
              XPCTR(K)= 0.d0
              ENDDO
          DO  I=NBEG+1,NEND-1
              DS=WF(I)**2
              XPTKE= XPTKE+ DS*(EINN-V(I))
              DS= WF(I)*(WF(I+1)- WF(I-1))
              IF(MORDR.GT.0) THEN
                  RR= RFN(I)
                  DO  K=1,MORDR
                      DS=DS*RR
                      XPCTR(K)=XPCTR(K)+DS
                      ENDDO
                  ENDIF
              ENDDO
          DO  K= 0,MORDR
              XPCTR(K)= XPCTR(K)/(2.d0*RH)
              ENDDO
        ENDIF
      XPTKE= XPTKE*RH/BFCT
      IF(MORDR.LT.0) GO TO 99
      DMR= 0.d0
      DO  K=0,MORDR
          XPCTR(K)=XPCTR(K)*RH
          DMR= DMR+ DM(K)*XPCTR(K)
          ENDDO
      IF((LXPCT.GT.0).AND.(IABS(LXPCT).LE.2)) THEN
          IF(EPR.LE.VLIM) WRITE(6,600) KV,JR,EPR,DMR,XPTKE
          IF(EPR.GT.VLIM) WRITE(6,602) KV,JR,EPR,DMR,XPTKE,GAMA
          IF(IABS(IRFN).LE.9) WRITE(6,604) (K,XPCTR(K),K=1,MORDR)
          IF(IPNCH.GE.1) WRITE(7,701) KV,JR,EPR,GAMA,XPTKE,DMR,
     1                                        (XPCTR(K),K=1,MORDR)
          ENDIF
      IF(ITRY.GT.19) GO TO 99
c** If appropriate, iteratively correct RREF value till distance
c  coordinate expectation value is identically zero.
      IF(IRFN.EQ.-1) THEN
c** For Dunham expansion parameter, define revised function here
          RREF=XPCTR(1)
          DRT=RREF
          WRITE(6,603) ITRY,DRT,RREF
          DO  I= 1,NPP
              RFN(I)= RFN(I)/RREF - 1.D0
              ENDDO
          ITRY=99
          GO TO 2
          ENDIF  
c** For Surkus-type expansion parameter, define revised function 
      ITRY=ITRY+1
      IF(ITRY.EQ.1) THEN
          RXPCT= XPCTR(1)
          RREF= 0.D0
          DRT= RXPCT
        ELSE
          DER= -IRFN/(2.d0*RREF)
          DRT= -XPCTR(1)/DER
        ENDIF
      RREF=RREF+DRT
      WRITE(6,603) ITRY,DRT,RREF
c** Redefine Surkus-type distance variable RFN using new RREF 
      RMINN= RMIN- RH
      DO  I= 1,NPP
          RR= RMINN+ I*RH
          RFN(I)= (RR**IRFN - RREF**IRFN)/(RR**IRFN + RREF**IRFN)
          ENDDO
      IF(DABS(DRT/RREF).GE.1.D-12) GO TO 2
   99 RETURN
  600 FORMAT(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,
     1  '   <KE>=',F11.3)
  602 FORMAT(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,
     1 '   <KE>=',F11.3/'   Tunneling predissociation  Width(FWHM)=',
     2 G13.6,'    <X**',I2,'>=',F13.8)
  604 FORMAT((8x,3('   <X**',I2,'>=',F13.8:)))
  603 FORMAT(' On iteration #',I2,'  change RREF by',1PD10.2,
     1  '  to   RREF=',0PF13.10,' [Angstroms]')
  701 FORMAT(2I4,F11.3,G11.4,F11.3,3(F12.7)/(5X,6F12.7))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE MATXEL(KV1,JROT1,IOMEG1,EO1,KV2,JROT2,IOMEG2,IRFN,EO2,
     1  NBEG,NEND,LXPCT,MORDR,DM,RH,WF1,WF2,RFN)
c** Subroutine to calculate matrix elements of powers of the distance
c  coordinate between vib. eigenfunction WF1(i) for v=KV1, J=JROT1 of
c  potential-1 & WF2(I), corresponding to KV2 & JROT2 of potentl.-2
      INTEGER I,J,IOMEG1,IOMEG2,IOMUP,IOMLW,IRFN,JROT1,JROT2,JUP,JLW,
     1 KV1,KV2,KVUP,KVLW,LXPCT,NBEG,NEND,MORDR
      REAL*8  ZMAT(0:20),WF1(NEND),WF2(NEND),RFN(NEND),DM(0:MORDR)
      REAL*8  AEINST,DEG,DME,DSM,EO1,EO2,ELW,FCF,FREQ,OMUP,RH,RI,SJ,
     1    ZJUP
      CHARACTER*1  DJ(-3:3)
      DATA DJ/'N','O','P','Q','R','S','T'/
      ZMAT(0)= 0.D0
      IF(MORDR.GE.1) THEN
          DO  J=1,MORDR
              ZMAT(J)=0.D0
              ENDDO
          ENDIF
      IF(IRFN.NE.-4) THEN
c** For regular power series or function matrix elements ...
          DO  I=NBEG,NEND
              DSM=WF2(I)*WF1(I)
              ZMAT(0)=ZMAT(0)+DSM
              RI= RFN(I)
              IF(MORDR.GE.1) THEN
                  DO  J=1,MORDR
                      DSM=DSM*RI
                      ZMAT(J)=ZMAT(J)+DSM
                      ENDDO
                  ENDIF
              ENDDO
        ELSE
c** For partial derivative matrix elements ...
          DO  I=NBEG+1,NEND-1
              DSM=WF1(I)*(WF2(I+1)- WF2(I-1))
              ZMAT(0)=ZMAT(0)+DSM
              RI= RFN(I)
              IF(MORDR.GE.1) THEN
                  DO  J=1,MORDR
                      DSM=DSM*RI
                      ZMAT(J)=ZMAT(J)+DSM
                      ENDDO
                  ENDIF
              ENDDO
          DO  J= 0,MORDR
              ZMAT(J)= ZMAT(J)/(2.d0*RH)
              ENDDO
        ENDIF
      DME=0.D0
      FCF= (ZMAT(0)*RH)**2
      IF(MORDR.GE.0) THEN
          DO  J=0,MORDR
              ZMAT(J)=ZMAT(J)*RH
              DME=DME+DM(J)*ZMAT(J)
              ENDDO
          ENDIF
      FREQ= EO2-EO1
      ELW= DMIN1(EO1,EO2)
c** Now calculate the Honl-London Factor for the particular transition
c   Factors updated as per Hansson & Watson JMS (2005).
      SJ= 0.D0
      KVUP= KV1
      KVLW= KV2
      JUP= JROT1
      JLW= JROT2
      IOMUP= IOMEG1
      IOMLW= IOMEG2
      IF(EO2.GT.EO1) THEN
          KVUP= KV2
          KVLW= KV1
          JUP= JROT2
          JLW= JROT1 
          IOMUP= IOMEG2
          IOMLW= IOMEG1
          ENDIF
      ZJUP= JUP
      OMUP= IOMUP
      DEG= 2*JUP+ 1
      IF((JLW.LT.IOMLW).OR.(JUP.LT.IOMUP)) GO TO 50
      IF(IOMUP.EQ.IOMLW) THEN
c** Factors for  DELTA(LAMBDA) = 0  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (ZJUP+ OMUP)*(JUP- IOMUP)/ZJUP
          IF((JUP.EQ.JLW).AND.(JUP.GT.0)) 
     1                       SJ= DEG*OMUP**2/(ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1)) SJ= (ZJUP+1.D0+OMUP)*(JUP+1-IOMUP)/
     1                                                     (ZJUP+1.D0)
          ENDIF
      IF(IOMUP.EQ.(IOMLW+1)) THEN
c** Factors for  DELTA(LAMBDA) = +1  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (ZJUP+OMUP)*(JUP-1+IOMUP)/(2.D0*ZJUP)
          IF((JUP.EQ.JLW).AND.(JUP.GT.0)) 
     1       SJ= (ZJUP+OMUP)*(JUP+1-IOMUP)*DEG/(2.D0*ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1)) 
     1       SJ= (JUP+1-IOMUP)*(ZJUP+2.D0-OMUP)/(2.D0*ZJUP+2.D0)
          ENDIF 
      IF(IOMUP.LT.IOMLW) THEN
c** Factors for  DELTA(LAMBDA) = -1  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (JUP-IOMUP)*(JUP-1-IOMUP)/(2.D0*ZJUP)
          IF((JUP.EQ.JLW).AND.(JUP.GT.0)) 
     1      SJ= (JUP-IOMUP)*(ZJUP+1.D0+OMUP)*DEG/(2.D0*ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1)) 
     1           SJ= (ZJUP+1.D0+OMUP)*(ZJUP+2.D0+OMUP)/(2.D0*ZJUP+2.D0)
          ENDIF 
c... finally, include Hansson-Watson  w0/w1  term in Honl-London factor
      IF((MIN(IOMUP,IOMLW).EQ.0).and.(IOMUP.NE.IOMLW)) SJ= SJ+SJ
c
c** For FREQ in  cm-1  and dipole moment in  debye , AEINST is the
c  absolute Einstein radiative emission rate (s-1) , using the
c  rotational intensity factors for sigma-sigma transitions.
   50 CONTINUE
      AEINST = DABS(3.1361891D-7 *DABS(FREQ)**3*DME**2 * SJ/DEG)
      IF(LXPCT.GT.0) THEN
          WRITE(6,600) KV1,JROT1,EO1,KV2,JROT2,EO2
          IF(IABS(IRFN).LE.9) WRITE(6,602) (J,ZMAT(J),J= 0,MORDR)
          WRITE(6,604) FCF,DME,FREQ,AEINST
          WRITE(6,606)
          ENDIF
      IF((IABS(LXPCT).EQ.4).OR.(IABS(LXPCT).EQ.5).AND.(SJ.GT.0.D0)) THEN
          IF(IABS(JUP-JLW).LE.3) WRITE(8,801) DJ(JUP-JLW),JLW,KVUP,
     1                                    KVLW,ELW,FREQ,AEINST,FCF,DME
c... Special printout for Hui/LeRoy N2 Quadrupole paper [JCP 1XX (2007)]
ccc       E00= 1175.7693d0
cc        WRITE(11,811) -FREQ,KVUP,JUP,KVLW,JLW,-FREQ,ELW-FREQ-E00,
cc   1                                      ELW-E00,DME**2
cc811 FORMAT(F12.4,2I4,I6,I4,3f12.4,1PD15.6)
          IF(IABS(JUP-JLW).GT.3) WRITE(8,802) JUP-JLW,JLW,KVUP,
     1                                    KVLW,ELW,FREQ,AEINST,FCF,DME
          ENDIF
      IF(IABS(LXPCT).GE.5) 
c    1         WRITE(7,701) KVUP,JUP,KVLW,JLW,(ZMAT(J),J=0,MORDR)
     1         WRITE(7,701) KVUP,JUP,KVLW,JLW,FREQ,(ZMAT(J),J=0,MORDR)
      RETURN
  600 FORMAT(' Coupling   E(v=',I3,', J=',I3,')=',F12.4,'   to   E(v=',
     1 I3,', J=',I3,')=',F12.4)
  602 FORMAT(5x,'Moment matrix elements:',2('   <X**',I2,'>=',F14.10:),
     1  1x/(3x,3('   <X**',I2,'>=',F14.10:),1x))
  604 FORMAT(' FCF=',1PD11.4,'   <M>=',D12.5,'   dE=',0PF12.4,
     1  '   A(Einst)=',1PD11.4,' s-1')
  606 FORMAT(1X,79('+'))
  701 FORMAT(4I4,F12.4,4F12.8:/(4X,6F12.8))
c 701 FORMAT(4I4,6F12.8:/(16X,6F12.8))
  801 FORMAT(1x,A1,'(',I3,')  ',I3,' -',I3,F10.2,F12.4,3(1PD13.5))
  802 FORMAT(i2,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE CDJOEL(EO,NBEG,NEND,BvWN,RH,WARN,V,WF0,RM2,RCNST)
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
c               RH    is the integration stepsize (in units [Ang])
c               WARN  is an integer flag: > 0 print internal warnings,
c               V(i)  is the effective potential (including centrifugal
c                     term if calculation performed at  J > 0) in 
c                     'internal' units, including the factor  RH**2/BvWN
c               RM2(i) is the array  1/(distance**2) in units [1/Ang**2]
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
      PARAMETER (NDMINT= 80001)
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
      REAL*8 V(NEND),WF0(NEND),RM2(NEND),P(NDMINT),WF1(NDMINT),
     1                                            WF2(NDMINT),RCNST(7)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,RH,RHSQ,
     1  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,
     2  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,
     3  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,
     4  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33
c
      IF(NEND.GT.NDMINT) THEN
          WRITE(6,602) NEND,NDMINT
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      RHSQ = RH*RH
      DV = RHSQ/12.D0
      E= EO*RHSQ/BvWN
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
      R2IN = R2IN*RH
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
      V1 = V(NEND) - E
      V2 = V(NEND-1) - E
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*(RM2(NEND) - R2IN)*WF0(NEND)
          G2 = (RM2(NEND-1) - R2IN)*WF0(NEND-1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NEND) - R2IN)*WF1(NEND)
     1                                                - DVV*WF0(NEND))
          G2 = (RM2(NEND-1) - R2IN)*WF1(NEND-1) - DVV*WF0(NEND-1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NEND) - R2IN)*WF2(NEND)
     1                                - DVV*WF1(NEND) - HVV*WF0(NEND))
          G2 = (RM2(NEND-1) - R2IN)*WF2(NEND-1) - DVV*WF1(NEND-1)
     1                                               - HVV*WF0(NEND-1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
c** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          IF(IPASS.EQ.1) G3 = (RM2(M) - R2IN)*WF0(M)
          IF(IPASS.EQ.2) G3 = (RM2(M) - R2IN)*WF1(M) - DVV*WF0(M)
          IF(IPASS.EQ.3) G3 = (RM2(M) - R2IN)*WF2(M) - DVV*WF1(M) 
     1                                                    - HVV*WF0(M)
          V3 = V(M) - E
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
      V1 = V(NBEG) - E
      V2 = V(NBEG+1) - E
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*(RM2(NBEG) - R2IN)*WF0(NBEG)
          G2 = (RM2(NBEG+1) - R2IN)*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NBEG) - R2IN)*WF1(NBEG)
     1                                                - DVV*WF0(NEND))
          G2 = (RM2(NBEG+1) - R2IN)*WF1(NBEG+1) - DVV*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) - DV*((RM2(NBEG) - R2IN)*WF2(NBEG)
     1                                - DVV*WF1(NEND) - HVV*WF0(NEND))
          G2 = (RM2(NBEG+1) - R2IN)*WF2(NBEG+1) - DVV*WF1(NBEG+1)
     2                                               - HVV*WF0(NBEG+1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
c** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + RHSQ*G2 + V2*P2
          P0 = WF0(I)
          IF(IPASS.EQ.1) G3 = (RM2(I) - R2IN)*P0
          IF(IPASS.EQ.2) G3 = (RM2(I)-R2IN)*WF1(I) - DVV*P0
          IF(IPASS.EQ.3) G3 = (RM2(I)-R2IN)*WF2(I) - DVV*WF1(I) - HVV*P0
          V3 = V(I) - E
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          AR = AR + P0*P3
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
          AR = AR + PI*P0
          ENDDO
      OV = AR*RH
      DO  I = NBEG,NEND
          P0 = WF0(I)
c ... and project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*RM2(I)
          IF(IPASS.EQ.1) THEN
c** Now - on first pass accumulate integrals for Dv and Hv
              WF1(I) = PI
              OV01 = OV01 + PI*P0
              OV11 = OV11 + PI*PI
              PER01 = PER01 + PIF*P0
              PER11 = PER11 + PI*PIF
            ELSEIF(IPASS.EQ.2) THEN
c ... and on next pass, accumulate integrals for Lv and Mv
              WF2(I) = PI
              P1 = WF1(I)
              OV02 = OV02 + PI*P0
              OV12 = OV12 + PI*P1
              OV22 = OV22 + PI*PI
              PER02 = PER02 + PIF*P0
              PER12 = PER12 + PIF*P1
              PER22 = PER22 + PI*PIF
            ELSEIF(IPASS.EQ.3) THEN
c ... and on next pass, accumulate integrals for Nv and Ov
              P1 = WF1(I)
              P2 = WF2(I)
              OV03 = OV03 + PI*P0
              OV13 = OV13 + PI*P1
              OV23 = OV23 + PI*P2
              OV33 = OV33 + PI*PI
              PER03 = PER03 + PIF*P0
              PER13 = PER13 + PIF*P1
              PER23 = PER23 + PIF*P2
              PER33 = PER33 + PIF*PI
            ENDIF
          ENDDO
      IF(IPASS.EQ.1) THEN
          DVV = RH*PER01
          HVV = RH*(PER11 - R2IN*OV11)
          IPASS = 2
          RCNST(2) = DVV*BvWN
          RCNST(3) = HVV*BvWn
          GO TO 10
        ELSEIF(IPASS.EQ.2) THEN
          HV2 = RH*PER02*BvWN
          LVV = RH*(PER12 - R2IN*OV12 - DVV*OV11)
          MVV = RH*(PER22 - R2IN*OV22 - 2.D0*DVV*OV12 - HVV*OV11)
          IPASS = 3
          RCNST(4) = LVV*BvWN
          RCNST(5) = MVV*BvWN
          GO TO 10
        ELSEIF(IPASS.EQ.3) THEN
          LV2 = RH*PER03*BvWN
          MV2 = RH*(PER13 - R2IN*OV13 - DVV*OV12 - HVV*OV11)*BvWN
          NVV = RH*(PER23 - R2IN*OV23 - DVV*(OV13 + OV22) 
     1                                     - 2.D0*HVV*OV12 - LVV*OV11)
          OVV = RH*(PER33 - R2IN*OV33 - 2.D0*DVV*OV23 
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

c***********************************************************************
      SUBROUTINE MASSES(IAN,IMN,NAME,GELGS,GNS,MASS,ABUND)
c***********************************************************************
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
c                          COPYRIGHT 2005
c** By R.J. Le Roy (with assistance from G.T. Kraemer & J.Y. Seto).
c                      Last modified  1 June 2005
c***********************************************************************
      REAL*8 zm(123,0:10),mass,ab(123,10),abund
      INTEGER i,ian,imn,gel(123),nmn(123),mn(123,10),ns2(123,10),
     1        gelgs, gns
      CHARACTER*2 NAME,AT(123)
c
      DATA  at(1),gel(1),nmn(1),(mn(1,i),i=1,3)/' H',2,3,1,2,3/
      DATA  (zm(1,i),i=0,3)/1.00794d0, 1.00782503207d0, 2.0141017778d0,
     1                      3.0160492777d0/
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
     1    170.93632358d0, 171.9363815d0, 172.9382108d0, 173.9388621d0,
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
      SUBROUTINE ALF(NDP,RMIN,RH,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,NCN,GV,
     1                                                 INNODE,INNR,IWR) 
c***********************************************************************
c   Version 1.0 dated July 6, 2003. Cosmetic formatting changes 11/04/07
c-----------------------------------------------------------------------
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from the first vibrational
c   level (v=0) to a user specified level (v=KVMAX) or the highest
c   allowed vibrational level of a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++   COPYRIGHT 1998 - 2003  by  Jenning Seto and Robert J. Le Roy   +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Please inform me of any bugs, by phone at: (519)888-4567, ext. 4051 +
c++++++++ by e-mail to: leroy@uwaterloo.ca , or write me at: ++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Based on the automatic level finding routine from LEVEL 6.0 written 
c   by Robert J. Le Roy
c** Uses the Schrodinger solver subroutine SCHRQ.
c
c** On entry:
c    NDP    is the number of datapoints used for the potential.
c    RMIN   is the inner radial distance of the potential (ang).
c    RH     is the meshvalue (ang).
c           NDP, RMIN, and RH define the radial range over which the
c           potential is defined.
c    V(i)   is the scaled input potential (cm-1).
c           The scaling factor BFCT is (2*mu/hbar^2)*RH^2.
c    VLIM   is the potential asymptote (cm-1).
c    KVMAX  is v for the highest vibrational level we wish to find.
c    AFLAG  is rot.quantum J for the (centrifugally distorted) potential
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    NCN    is the near dissociation limit radial exponential.
c    INNODE specifies whether wave fx. initiation @ RMIN starts with a
c        note (normal case: INNODE > 0) or zero slope (when INNODE.le.0)
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   is vib.quantum number for the highest allowed vibrational 
c            level found (may be less than the input value of KVMAX).
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : Subroutine found all levels to v=KVMAX(input).
c             = -1 : KVMAX larger than number of allowed levels.
c             = -2 : Initial trial energy is unusable.
c             = -3 : Calculated trial energy is unusable.
c             = -4 : Cannot find first vibrational level.
c             = -5 : Calculated trial energy too low.
c             = -6 : Calculated trial energy too high.
c             = -7 : An impossible situation occured.
c             = -8 : Potential found to have a second minimum.
c    GV(v)   contains the vibrational energy levels found for v=0-KVMAX
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          <> 0 : print error & warning descriptions.
c          >  0 : also print intermediate ALF messages.
c    MCO   specifies the level of printing of molecular constants.
c          >  0 : print out vibrational energies to channel-21.
c    INNER specifies wave function matching (& initiation) conditions.
c        .le.0 : Match inward & outward solutions at outermost well t.p.
c          > 0 : Match at innermost well inner turning point
c        For most normal cases set INNER = 0,  but ......
c            To find "inner-well-dominated" solutions of an asymmetric
c            double minimum potential, set  INNER > 0.
c    LPRWF specifies option of printing out generated wavefunction
c          > 0 : print wave function every LPRWF-th  point.
c          < 0 : compactly write to channel-7 every |LPRWF|-th wave
c                function value.
c          A lead "card" identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The dimensioning parameters must be consistant with the sizes of the
c   arrays used in the calling program.
c
c    NVIBMX  is the maximum number of vibrational levels considered.
c            Note: NVIBMX should be larger than KVMAX.
c
      INTEGER NVIBMX
      PARAMETER (NVIBMX= 400)
c
c** NF counts levels found in automatic search option
c
c** OWL holds the vibrational levels that are contained in the outer
c   well.
c** IWL holds the vibrational levels that are contained in the inner
c   well (if present).
c
      INTEGER NDP,KVMAX,NCN,KV,AFLAG,NF,NBEG,NEND,INNR(0:KVMAX),IWR,
     1  I,IZPE,IVDIF,IVCOR,IQT,IEG,LTRY,TRIALE,AWO,MCO,INNODE,INNER,
     2  NoPRWF,JROT,NPMIN, NPMAX, NIWL,IWL(0:NVIBMX),NOWL,OWL(0:NVIBMX)
c
      REAL*8 RMIN,RMAX,RH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,LHIE,LLOE,
     1  BZ,BFCT,PW,PWI,GAMA,VMIN,VMAX,RE,PMAX,VDMV,VDL,VDU,DRAVOUT,
     2  DRAVIN,AO,VD,GV(0:KVMAX),RAVG(0:NVIBMX),VPMIN(10),RPMIN(10),
     3  VPMAX(10),RPMAX(10)
c
      DATA AWO/-1/,MCO/0/,NoPRWF/0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check that the array dimensions are adequate.
      IF(KVMAX.GT.NVIBMX) THEN
          WRITE(6,610)
          WRITE(6,613) KVMAX, NVIBMX
          STOP
          ENDIF
c
c** Initialize level counters for each well.
      DO  I= 0,KVMAX
          INNR(I)= -1
          IWL(I)= 0
          OWL(I)= 0
          END DO
c
c** Initialize the remaining variables and flags.
      NF= 0
      NIWL= 0
      NOWL= 0
      KV= 0
      INNER= 0
      LTRY= 0
      LHIE= 1.0d99
      CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** Store rotational quantum number.
      JROT= AFLAG
c
c** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
c  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920d0
      BFCT= BZ*RH*RH
c
c** RMAX is the outer radial distance over which potential is defined. 
      RMAX= RMIN + DBLE(NDP-1)*RH
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
      NPMIN= 0
      DO  I= 2,NDP-1
          IF((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
              NPMIN= NPMIN + 1
              RPMIN(NPMIN)= RMIN + DBLE(I-1)*RH
              VPMIN(NPMIN)= V(I) / BFCT
              IF(NPMIN.EQ.10) GOTO 100
              ENDIF
          END DO
c** If NO minimum can be found, then print a warning and stop.
  100 IF(NPMIN.EQ.0) THEN
          IF(V(2).LE.V(1)) THEN
              WRITE(6,614) JROT
              KVMAX= -1
              RETURN
              ENDIF
          NPMIN= 1
          VPMIN(NPMIN)= V(1)/BFCT
          RPMIN(NPMIN)= RMIN
          WRITE(6,618) VPMIN(1),RMIN
          ENDIF
c** If more than two minima are found, then print a warning and hope ...
      IF(NPMIN.GT.2) THEN
          WRITE(6,605)
          WRITE(6,615) NPMIN, 'minima'
c         STOP
          ENDIF
      NPMAX= 0
c** Locate the potential maxima (if it exists).
      DO  I= 2,NDP-1
          IF((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
              NPMAX= NPMAX + 1
              RPMAX(NPMAX)= RMIN + DBLE(I-1)*RH
              VPMAX(NPMAX)= V(I) / BFCT
              IF(NPMAX.EQ.10) GOTO 150
              ENDIF
          END DO
c
  150 IF(NPMAX.EQ.0) THEN
c** If no maxima were found, then set energy maximum to be the value 
c   at the end of the radial range.
          NPMAX= 1
          RPMAX(NPMAX)= RMAX
          VPMAX(NPMAX)= V(NDP) / BFCT
          ENDIF
c
c** If more than three maxima found, then print a warning and hope ...
      IF(NPMAX.GT.3) THEN
          WRITE(6,605)
          WRITE(6,615) NPMAX, 'maxima'
c         STOP
          ENDIF
c
c** If there is no barrier to dissociation, then set the
c   final VPMAX to be the value at the end of the range. [huh ????]
      IF(RPMAX(NPMAX).LT.RPMIN(NPMIN)) THEN
          NPMAX= NPMAX + 1
          RPMAX(NPMAX)= RMAX
          VPMAX(NPMAX)= V(NDP) / BFCT
          ENDIF
c
c** If innermost maximum occurs before innermost minimum, the potential 
c   turns over in short range region and should not be used.  
c   Print a warning and STOP.
      IF(RPMAX(1).LT.RPMIN(1)) THEN
          WRITE(6,610)
          WRITE(6,616) RPMAX(1)
          STOP
          ENDIF
c
c** Now find the absolute potential minimum.
      VMIN= VPMIN(1)
      RE= RPMIN(1)
      DO  I= 2,NPMIN
          IF(VMIN.GT.VPMIN(I)) THEN
              VMIN= VPMIN(I)
              RE= RPMIN(I)
              ENDIF
          END DO 
c
c** Now find the absolute potential maximum.
      VMAX= VPMAX(1)
      DO  I= 2,NPMAX
          IF(VMAX.LT.VPMAX(I)) VMAX= VPMAX(I)
          END DO 
c
c** If the absolute potential maximum is lower than the absolute
c   potential minimum, then print out an error statement and quit.
      IF(VMAX.LE.VMIN) THEN
          WRITE(6,610)
          WRITE(6,617)
          STOP
          ENDIF
c
c** Otherwise, print out potential extrema count, if desired.
      IF(AWO.GT.0) THEN
          WRITE(6,650) NPMIN, VMIN
          WRITE(6,651) NPMAX, VMAX
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculate 2*NCN/(NCN - 2) for use when calculating trial energies. 
      PW= 20.0d0
      IF((NCN.GT.0).AND.(NCN.NE.2)) PW= 2.d0*DBLE(NCN)/(DBLE(NCN)-2.d0)
      IF(VMAX.GT.VLIM) PW= 2.d0
      PWI= 1.d0/PW
c
c** Use Lennard-Jones estimation of zero point energy to determine the
c   initial trial energy.
c                            _____________________________
c            vD + 0.5 = ao \/ZMU * De * Re^2 / 16.85762920
c                 De = A (vD - v)^3 = A (vD + 0.5)^3
c              E(v=0) = VMIN + A [(vD + 0.5)^3 - vD^3]
c** Choose AO to have a value of 0.25.
c
      AO= 0.25d0
      VD= AO*DSQRT(BZ*(VMAX-VMIN))
      VD= VD*RE - 0.5d0
      AO= (VMAX-VMIN)*(1.d0 - (VD/(VD+0.5d0))**3)
      EO= VMIN + AO
c
      IF(MCO.GE.1) THEN
c** If desired, write out SCHRQ control information.
          WRITE(21,2100)
          WRITE(21,2110) RMIN, RMAX, RH, BZ, ZMU
          WRITE(21,2111) EPS
          WRITE(21,2112)
          WRITE(21,2101)
          ENDIF
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c
   10 IF(AWO.GT.0) THEN
          IF(INNER.EQ.0) WRITE(6,601)
          IF(INNER.EQ.1) WRITE(6,603)
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,RMIN,RH,NDP,
     1                              NBEG,NEND,INNODE,INNER,IWR,NoPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a level higher than NF has already been found and the level found
c   not the correct level, ignore all errors and continue with the 
c   binary search.
c
      IF((KV.GE.0).AND.(KV.LT.NF).AND.((IVCOR.GT.0).OR.(IEG.GT.3))) THEN
          KV= NF + 1
          EO= VMAX + 1.0d0
          IF(IVCOR.GT.0) THEN
              LLOE= (LHIE + LLOE) / 2.0d0
            ELSE
              LLOE= GV(NF-1)
              LHIE= VMAX
              IVCOR= 1
            ENDIF
          TRIALE= 0
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The SCHRQ error condition is KV < 0.
c   There are three possible situations to consider:
c     EO > VMAX : Trial energy greater than potential maximum
c     NF = 0 : Looking for the first vibrational level (v = 0)
c     NF > 0 : Looking for the other vibrational levels (v > 0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the case when the next trial energy is higher than the potential
c   maximum, try one last ditch attempt to find the highest bound level
c   (quasi or otherwise) in the potential.
c
      IF((KV.LT.0).AND.(EO.GT.VMAX)) THEN
          IF(IVCOR.LT.20) THEN
              KV= NF + 1
              IF(IVCOR.EQ.0) THEN
c??????
                  if(nf.eq.0) then
                      lloe= GV(0)
                    else
                      LLOE= GV(NF-1)
                    endif
c??????           LLOE= GV(NF-1)
c??????
                  LHIE= VMAX
                  IVCOR= 1
                  ENDIF
              ENDIF
          IF(LTRY.LT.1) THEN
              LTRY= 1
              KV= 999
              EO= VMAX - 1.0d-2
c
c** If unsuccessful, then print out a warning and exit.
            ELSE
              IF(AWO.NE.0) THEN
                  WRITE(6,605)
                  WRITE(6,606) NF, EO, VMAX
                  ENDIF
              AFLAG= -1
              GOTO 200
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find the first vibrational level (v=0), then double the
c   zero point energy estimation (AO):   E(v=0) = VMIN + IQT*AO
c
        ELSEIF((KV.LT.0).AND.(NF.EQ.0)) THEN
          IF(IQT.GT.1) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,611)
                  WRITE(6,620) IQT, EO
                  ENDIF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero-point level.
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -2
                  GOTO 200
                ENDIF
              ENDIF
          IQT= IQT + 1
          EO= VMIN + DBLE(IQT)*AO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find other vibrational levels (v > 0) then switch to
c   use of differences for estimating spacing.
c
        ELSEIF((KV.LT.0).AND.(NF.GT.0)) THEN
          IF(IVDIF.GT.0) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,612)
                  WRITE(6,621) NF,IVDIF
                ENDIF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -3
                  GOTO 200
                ENDIF
              ENDIF
          IVDIF= 1
          IF(INNER.EQ.0) THEN
              CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
            ELSE
              CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If first level found isn't v=0, try up to 3 times to 'harmonically'
c   estimate improved trial ground state energy.
c           E(v=0)= E(v=KV) - (E(v=KV) - VMIN)/(1 + KV/2)
c
        ELSEIF((KV.GT.0).AND.(NF.EQ.0)) THEN
          IF(IZPE.GT.3) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,611)
                  WRITE(6,622) IZPE,GV(0),KV,EO
                  ENDIF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -4
                  GOTO 200
                ENDIF
              ENDIF
          IZPE= IZPE + 1
          EO= EO - (EO-VMIN)/(1.d0+0.5d0/KV)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the next three cases, KV >= 0 and NF > 0.
c** If the calculated vibrational level is less than the next expected
c   level, then the estimated trial energy is too low.
c** Perhaps the difference in energy between vibrational levels v and
c   v-1 is much greater than the energy between levels v-1 and v-2. 
c                 E(v) - E(v-1) >> E(v-1) - E(v-2)
c   In which case (most likely a potential with a shelf), try twice to
c   estimate a higher trial energy.
c
c   E(v)= E(v-1) + (1+IEG/2) * (2*(E(v-1)-E(v-2)) - (E(v-2)-E(v-3)))
c
        ELSEIF(KV.LT.NF) THEN
          IF(IEG.GT.3) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,612)
                  WRITE(6,623) NF, KV
                  ENDIF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
                ELSE
                  AFLAG= -5
                  GOTO 200
                ENDIF
              ENDIF
          IEG= IEG + 1
c
c** If a second minimum is present, then the next vibrational level may
c   be in the inner well. If so, use the inner well vibrational levels
c   to estimate the next trial energy.
          IF(INNER.EQ.0) THEN
              CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
            ELSE
              CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
            ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is the next expected level, ...
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSEIF(KV.EQ.NF) THEN
          GV(NF)= EO
          INNR(NF)= INNER
          LTRY= 0
          TRIALE= 1
          LLOE= EO
          LHIE= VMAX
          CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c-----------------------------------------------------------------------
c
c** To ease confusion when using a potential with two minima, use
c  <R>=RAVG  to keep track of which is which.
c
   20     RAVG(NF)= 0.5d0*(SWF(NBEG)*(RMIN+DBLE(NBEG-1)*RH)
     1                        + SWF(NEND)*(RMIN+DBLE(NEND-1)*RH))
          DO  I= NBEG+1,NEND-1
              RAVG(NF)= RAVG(NF) + SWF(I)*(RMIN+DBLE(I-1)*RH)
              END DO
          RAVG(NF)= RAVG(NF)*RH
c
c** Double check that the calculated level is in fact located in the
c   correct well. This can be done (for v.ne.0) by comparing <R> for
c   the new level and with previous values in each well. If the
c   difference is greater than 1.5 times the difference in the other
c   well, assume the calculated level is probably in the wrong well.
c
          IF(NOWL.GT.0) THEN
              DRAVOUT= DABS(RAVG(NF) - RAVG(OWL(NOWL-1)))
            ELSE
              DRAVOUT= 9999.9d0
            ENDIF
          IF(NIWL.GT.0) THEN
              DRAVIN= DABS(RAVG(NF) - RAVG(IWL(NIWL-1)))
            ELSE
              DRAVIN= 9999.9d0
            ENDIF
          IF(INNER.EQ.0) THEN
              IF((NOWL.GT.0).AND.(DRAVOUT.GT.(1.5d0*DRAVIN))) THEN
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Inner',NIWL,GV(NF)-VMIN,RAVG(NF)
                  IWL(NIWL)= NF
                  NIWL= NIWL + 1
                ELSE
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Outer',NOWL,GV(NF)-VMIN,RAVG(NF)
                  OWL(NOWL)= NF
                  NOWL= NOWL + 1
                ENDIF
            ELSE
              IF((NIWL.GT.0).AND.(DRAVIN.GT.(1.5d0*DRAVOUT))) THEN
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Outer',NOWL,GV(NF)-VMIN,RAVG(NF)
                  OWL(NOWL)= NF
                  NOWL= NOWL + 1
                ELSE
                  IF(MCO.GE.1) 
     1             WRITE(21,2113) NF,'Inner',NIWL,GV(NF)-VMIN,RAVG(NF)
                  IWL(NIWL)= NF
                  NIWL= NIWL + 1
                ENDIF
              INNER= 0
            ENDIF
c
c** Look for the next uncalculated level.
          NF= NF + 1
          IF(NF.LE.KVMAX) THEN
              IF(INNR(NF).GE.0) THEN
c ... if this next level was found earlier, its INNER > initial -1 value
c  so skip up one, but first calculate its Bv on the way.
                  GOTO 20
                  ENDIF
              ENDIF
c-----------------------------------------------------------------------
c** Now estimate trial energy for next higher vibrational energy level
c   by using the Near-Dissociation Theory result that:
c                  (binding energy)**((NCN-2)/(2*NCN))
c   is (at least locally) linear in vibrational quantum number.
c
          IF(NF.EQ.1) THEN
              VDMV= 0.5d0/(((VMAX-VMIN)/(VMAX-GV(0)))**PWI - 1.d0)
            ELSE
              VDMV= 1.d0/(((VMAX-GV(NF-2))/(VMAX-GV(NF-1)))**PWI- 1.d0)
            ENDIF
c
c** If unable to calculate the next trial energy, see if all of the
c   desired levels have been calculated. If not then turn on the warning
c   flag and quit, otherwise print out success message and quit.
c
          IF((VDMV.LT.1.d0).AND.(NCN.GT.2)) THEN
              IF(NF.LE.KVMAX) THEN
                  AFLAG= -1
                  WRITE(6,640) JROT, NF - 1 + VDMV
                ELSEIF(AWO.GT.0) THEN
                  WRITE(6,630) KVMAX
                ENDIF
              GOTO 200
              ENDIF
c
c** Now calculate the next trial energy.
c
          EO= VMAX - (VMAX-GV(NF-1))*(1.d0-1.d0/VDMV)**PW
c
c** However, if the level is above the dissociation limit (for
c   potentials with barriers) then use differences to calculate the
c   next trial energy.
c
          IF(EO.GT.VMAX) CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is higher then the next expected
c   level, then try once to interpolate harmonically for the missed
c   level:      E(v)= E(v-1) + (E(KV) - E(v-1)) / 2
c
        ELSEIF(KV.GT.NF) THEN
c
c** Record vibrational level (if haven't already) for posterity.
c
          IF((KV.LE.KVMAX).AND.(EO.LT.VMAX)) THEN
              IF(INNR(KV).EQ.-1) THEN
                  GV(KV)= EO
                  INNR(KV)= INNER
                  ENDIF
              ENDIF
          IF(IVCOR.GT.19) THEN
              IF(AWO.NE.0) THEN
                  WRITE(6,610)
                  WRITE(6,612)
                  WRITE(6,624) IVCOR,KV,EO,(NF-1),GV(NF-1)
                  ENDIF
c
c** If interpolation fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   missing level. 
              IF(INNER.EQ.0) THEN
                  INNER= 1
                  LLOE= GV(NF-1)
                  LHIE= VMAX
                  CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print a warning message and exit 
                ELSE
                  AFLAG= -6
                  GOTO 200
                ENDIF
              ENDIF
c
c** Use NDE theory to predict eigenvalue for the missing level.
          IF(IVCOR.EQ.0) THEN
              VDU= (VMAX-EO)**PWI
              VDL= (VMAX-GV(NF-1))**PWI
              EO= VMAX - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
c         IF((NPMIN.EQ.1).AND.(IVCOR.EQ.0)) THEN
c             VDU= (VPMAX(1)-EO)**PWI
c             VDL= (VPMAX(1)-GV(OWL(NOWL-1)))**PWI
c             EO= VPMAX(1)- (VDL+ (VDU - VDL)/ DBLE(KV - NF + 1))**PW
c           ELSEIF((((INNER.EQ.1).AND.(NIWL.GT.0)).OR.
c    1          ((INNER.EQ.0).AND.(NOWL.EQ.0))).AND.(IVCOR.EQ.0)) THEN
c             VDU= (VPMAX(1)-EO)**PWI
c             VDL= (VPMAX(1)-GV(IWL(NIWL-1)))**PWI
c             EO= VPMAX(1)- (VDL+ (VDU - VDL)/ DBLE(KV - NF + 1))**PW
c           ELSEIF((((INNER.EQ.0).AND.(NOWL.GT.0)).OR.
c     1         ((INNER.EQ.1).AND.(NIWL.EQ.0))).AND.(IVCOR.EQ.0)) THEN
c             VDU= (VPMAX(2)-EO)**PWI
c             VDL= (VPMAX(2)-GV(OWL(NOWL-1)))**PWI
c             EO= VPMAX(2)- (VDL+ (VDU - VDL)/ DBLE(KV - NF + 1))**PW
c
c** Otherwise use a binary search.
c
            ELSE
              IF((TRIALE.EQ.1).AND.(EO.LT.LHIE)) THEN
                  LHIE= EO
                  IVCOR= 0
                ELSEIF(TRIALE.EQ.1) THEN
                  LHIE= (LHIE + LLOE) / 2.0d0
                ENDIF
              TRIALE= 1
              EO= (LHIE + LLOE) / 2.0d0
            ENDIF
          IVCOR= IVCOR + 1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If an unknown case occurs (quite impossible but don't quote me on
c   it) then write out an error message and exit.
        ELSE
          IF(AWO.NE.0) THEN
              WRITE(6,610)
              WRITE(6,666) KV,NF
              ENDIF
          AFLAG= -7
          GOTO 200
        ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set KV to the next uncalculated vibrational level to be found unless
c   looking for the highest vibrational level.
c
      IF(KV.NE.999) KV= NF
c
c** If still haven't found all of the vibrational levels then
c   look for the next vibrational level.
c
      IF((KV.LE.KVMAX).OR.(KV.EQ.999)) GOTO 10
c
c** Otherwise, print out a message saying that all is well.
c
      IF((KV.GT.KVMAX).AND.(AWO.GT.0)) WRITE(6,630) KVMAX
c
c** If the potential has levels in a second minimum, then print out a
c   list of those levels to channel-21 if desired.
c
      IF((NIWL.GT.0).AND.(NOWL.GT.0)) THEN
          IF(MCO.GE.1) WRITE(21,2114) NIWL, NOWL
          IF(AWO.NE.0) THEN
              WRITE(6,605)
              WRITE(6,607)
              ENDIF
          AFLAG= -8
          ENDIF
c
  200 IF(AFLAG.LT.0) THEN
c** If unable to find all KVMAX+1 levels requested, then return KVMAX as
c  v for the highest vibrational level actually found, and print out the
c  the energy of that level.
          KVMAX= NF - 1
          IF(AWO.NE.0) WRITE(6,626) KVMAX, GV(KVMAX)
          ENDIF
      IF(MCO.GE.1) THEN
c         WRITE(21,2100)
c         DO  I= 0,NF-1
c             IF(INNR(I).EQ.0)WRITE(21,2113)I,'Outer',I,GV(I)-VMIN,RAVG(I)
c             IF(INNR(I).EQ.1)WRITE(21,2113)I,'Inner',I,GV(I)-VMIN,RAVG(I)
c             END DO
          WRITE(21,2100)
          ENDIF
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(/' Solve by matching inward and outward solutions at the'/
     1  5x,'outermost wave function maximum, S(max), where  R = RR(M)')
  603 FORMAT(/' Solve by matching inward and outward solutions at the'/
     1  5x,'innermost turning point   R1 = R(M)')
  605 FORMAT(/'  *** ALF WARNING ***')
  606 FORMAT(4x,'Next estimated trial energy  E(v=',I3,') =',G15.8/4X,
     1  'lies above potential maximum  VMAX =',G15.8)
  607 FORMAT(4X,'Potential found to have a second minimum.')
  610 FORMAT(/'  *** ALF ERROR ***')
  611 FORMAT(4X,'Attempt to find zero point level fails!')
  612 FORMAT(4X,'Attempt to find next higher vibrational level fails!')
  613 FORMAT(4X,'Number of vib levels requested=',i4,' exceeds internal 
     1ALF array dimension  NVIBMX=',i4)
  614 FORMAT(/'  *** ALF ERROR ***   Unable to find a potential minimum
     1 for   J=',i4)
  615 FORMAT(4X,'There are',I3,'  potential ',A6,' in this potential. St
     1op searching after 10.')
  616 FORMAT(4X,'The potential turns over in the short range region at R
     1 = ',G15.8)
  617 FORMAT(4X,'VMAX =',G15.8,' found to be less than VMIN =',G15.8)
  618 FORMAT(/'  ALF  finds onee potential minimum of',1PD15.7,
     1  '  at  R(1)=',0Pf9.6)
  620 FORMAT(4X,'Use of energy ',I1,'0% up the potential well (E =',
     1G15.8,')'/4X,' fails to produce a viable vibrational eigenstate.')
  621 FORMAT(4X,'Use of differences to estimate the energy for the next'
     1/4X,' vibrational level (v=',I3,') failed after',I3,'  attempt.')
  622 FORMAT(4X,'After',I3,' tries to harmonically estimate the zero-poi
     1nt energy,'/4X,' initial trial energy',G15.8,'   had yielded   E(v
     2=',I3,') =',G15.8)
  623 FORMAT(4X,'Expecting to find level (v=',I3,') but found level (v='
     1,I3,')')
  624 FORMAT(4X,'After',I3,' tries, failed to interpolate trial energy b
     1etween'/4X,'E(v=',I3,') =',G15.8,'   and   E(v=',I3,') =',G15.8)
  626 FORMAT(4X,'The highest calculated level is  E(v=',I3,') =',G15.8)
  630 FORMAT(/' ALF successfully finds all vibrational levels up to   v=
     1 KVMAX=',I3)
  640 FORMAT(/' ALF finds all  J=',i3,'  vib. levels below  vD=',F7.3,
     1  '  estimated by N-D theory')
  650 FORMAT(/' There were',I3,'  potential minima found with the absolu
     1te minimum'/4X,'VMIN =',G15.8,'  cm-1.')
  651 FORMAT(/' There were',I3,'  potential maxima found with the absolu
     1te maximum'/4X,'VMAX =',G15.8,'  cm-1.')
  666 FORMAT(4X,'Undefined case for automatic search.'/,4X,'Values of KV
     1 =',I3,'  and NF =',I3)
 2100 FORMAT(/1X,39('=='))
 2101 FORMAT(/1X,39('--'))
 2110 FORMAT(/' Limits and increment of integration (in Angstroms):'
     1 /'    RMIN =',F6.3,'    RMAX =',F7.3,'    RH =',F9.6,
     2 //' Generate    BZ =',G19.12,' ((1/cm-1)(1/Angstroms**2))'
     3 /' from ZMU:',F15.11,' (amu)')
 2111 FORMAT(/' Eigenvalue convergence criterion is   EPS =',G11.4,'(cm-
     11)')
 2112 FORMAT(/' Calculating properties of the potential described above.
     1 '/' Use Airy function at 3-rd turning point as outer boundary'
     2 /' condition for quasibound levels.')
 2113 FORMAT(' v=',I3,4X,'v(',A5,')=',I3,4X,'Gv=',F16.9,4X,'Bv=',F16.12)
 2114 FORMAT(/' Found',I4,' level(s) in the inner well and',I4,' level(s
     1) in the outer well.')
      END
c***********************************************************************
      SUBROUTINE INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c***********************************************************************
c** This subroutine reinitializes the condition flags when considering a
c   new case (found next vibrational level or finding level in inner
c   well - INNER = 1).
c
c** On entry and exit:
c    IQT     Case when KV < 0 and NF = 0
c            determines the value used for the initial trial energy.
c    IVDIF   Case when KV < 0 and NF > 0
c            is the flag denoting the use of differences to calculate
c            trial energies.
c    IZPE    Case when KV > 0 and NF = 0
c            is the number of times the zero point energy (v = 0) has
c            been estimated harmonically.
c    IEG     Case when KV < NF and NF > 0
c            are the number of times that a larger trial energy is used
c            to find the next level.
c    IVCOR   Case when KV > NF and NF > 0
c            are the number of times that a smaller trial energy is used
c            to find the next level.
c
      INTEGER IZPE,IVDIF,IVCOR,IQT,IEG
c
      IQT= 1
      IVDIF= 0
      IZPE= 0
      IEG= 0
      IVCOR= 0
      RETURN
      END
c***********************************************************************
      SUBROUTINE DTENG(IEG,NF,NVEL,VEL,NVIBMX,VMIN,GV,EO)
c***********************************************************************
c** This subroutine calculates the next trial energy using differences.
c
c** On entry:
c    IEG     factor by which a larger trial energy should be calculated:
c              NVEL = 2 : Increase correction by increments of 25%
c              NVEL > 2 : Increase correction by increments of 50%
c    NF      is the highest calculated vibrational level.
c    NVEL    is the number of levels found in the potential well.
c    VEL(v)  keeps track of all levels in the potential well.
c    NVIBMX  is the maximum number of vibrational levels (dimension).
c    VMIN    is the absolute value of the potential minimum (cm-1).
c    GV(v)   contains the vibrational energy level spacings
c            and rotational constants for each level (cm-1).
c
c** On exit:
c    EO      is the calculated trial energy.
c
      INTEGER IEG,NF,NVEL,NVIBMX,VEL(0:NVIBMX)
c
      REAL*8 VMIN,GV(0:NVIBMX),EO
c
c** If determining the first (non-zero point energy) level in the well,
c   then use the last determined level in the other well plus a larger
c   than harmonic correction that becomes smaller with each new
c   iteration:      E(v=0)= E(v=NF-1) + (E(v=NF-1)-VMIN)/(NF-1+IEG/4) 
c
      IF(NVEL.EQ.0) THEN
        EO= GV(NF-1) + (GV(NF-1) - VMIN)/(NF - 1 + 0.25d0*DBLE(IEG))
c
c** Try to get v= 1 using smaller-than-harmonic spacing.
c
c              E(v=1)= E(v=0) + (1.0+IEG/4)*(E(v=0)-VMIN)
c
      ELSEIF(NVEL.EQ.1) THEN
        EO= GV(VEL(0)) + (1.d0+DBLE(IEG)*0.25d0)*(GV(VEL(0))-VMIN)
c
c** Try to get v= 2 using a sequentially increasing correction.
c
c             E(v=2)= E(v=1) + (0.8+IEG/4)*(E(v=1)-E(v=0))
c
      ELSEIF(NVEL.EQ.2) THEN
        EO= GV(VEL(1))+ (0.8d0+DBLE(IEG)*0.25d0)*(GV(VEL(1))-GV(VEL(0)))
c
c** Try to get v > 2 using a sequentially increasing correction.
c
c        E(v)= E(v-1) + (1.0+IEG)*(2.0*E(v-1)-3.0*E(v-2)+E(v-3))
c
      ELSE
        EO= GV(VEL(NVEL-1)) + (1.d0+DBLE(IEG))
     1    *(2.d0*GV(VEL(NVEL-1))-3.0d0*GV(VEL(NVEL-2))+GV(VEL(NVEL-3)))
      ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c***** R.J. Le Roy  subroutine SCHRQ, last updated  8 April 2007 *******
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                COPYRIGHT 2007  by  Robert J. Le Roy                  +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** SCHRQ solves radial Schrodinger equation in dimensionless form
c  d2WF/dR2 = - (E-V(R))*WF(R) ,  where WF(I) is the wave function.
c** Integrate by Numerov method over N mesh points with increment
c  H=RH across range beginning at RMIN .
c** Input trial energy EO, eigenvalue convergence criterion EEPS
c  potential asymptote VLIM, and all returned energies (EO, GAMA & VMAX)
c  have units (cm-1).
c** On entry, the input potential V(I) must include the centrifugal
c  term and the factor:  'BFCT'=2*mu*(2*pi*RH/hPLANCK)**2  (1/cm-1) ,
c  which is also internally incorporated into EO, VLIM & EEPS.
c* Note that these reduced quantities (& the internal eigenvalue E)
c  contain a factor of the squared integration increment  RH**2 .
c  This saves arithmetic work in the innermost loop of the algorithm.
c** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.85762920 (1/cm-1)
c** INNODE > 0  specifies that wavefx. initiates at RMIN with a node 
c     (normal default case);  INNODE.le.0  spoecifies  zero slope  at
c     RMIN (for finding symmetric eigenfunctions of symmetric potential
c     with potential mid-point @ RMIN).
c** INNER specifies wave function matching condition: INNER = 0  makes
c     matching of inward & outward solutions occur at outermost turning
c     point;  INNER > 0 makes matching occur at innermost turning point.
c * Normally use  INNER=0 ,  but to find inner-well levels of double 
c     minimum potential, set  INNER > 0 .
c----------------------------------------------------------------------
      SUBROUTINE SCHRQ(KV,JROT,EO,GAMA,VMAX,VLIM,V,WF,BFCT,EEPS,RMIN,
     1                          RH,N,NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c----------------------------------------------------------------------
c** Output vibrational quantum number KV, eigenvalue EO, normalized
c  wave function WF(I), and range, NBEG .le. I .le. NEND  over
c  which WF(I) is defined. *** Have set  WF(I)=0  outside this range.
c* (NBEG,NEND), defined by requiring  abs(WF(I)) < RATST=1.D-9  outside.
c** If(LPRWF.gt.0) print wavefunction WF(I) every LPRWF-th point.
c* If(LPRWF.lt.0) "punch" (i.e., WRITE(10,XXX)) every |LPRWF|-th point
c  of the wave function on disk starting at R(NBEG) with step size
c  of  IPSIQ=|LPRWF|*RH. 
c** For energies above the potential asymptote VLIM, locate quasibound
c  levels using Airy function boundary condition and return the level
c  width GAMA and barrier height VMAX, as well as EO.
c** ERROR condition on return is  KV < 0 ; usually KV=-1, but return
c  KV=-2 if error appears to arise from too low trial energy.
c** If(IWR.ne.0) print error & warning descriptions
c  If (IWR.gt.0) also print final eigenvalues & node count.
c  If (IWR.ge.2) also show end-of-range wave function amplitudes
c  If (IWR.ge.3) print also intermediate trial eigenvalues, etc.
c** If input KV.ge.998 , tries to find highest bound level, and
c  trial energy should be only slightly less than VLIM.
c** If input KV < -10 , use log-derivative outer boundary condition at
c  mesh point |KV| , based on incoming value of wave function WF(|KV|)
c  and of the wavefunction derivative at that point, SPNEND, which is
c  brought in as WF(|KV|-1).  For a hard wall condition at mesh point
c  |KV|, set WF(|KV|)=0 and WF(|KV|-1)= -1 before entry.
c----------------------------------------------------------------------
c++ "SCHRQ" calls subroutineas "QBOUND" and "WIDTH", and the latter
c++ calls "LEVQAD" .
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IBEGIN,ICOR,IJ,IJK,INNODE,INNER,IPSID,IQTST,IT,
     1         ITER,ITP1,ITP1P,ITP3,IWR,J,JJ,J1,J2,JPSIQ,JQTST,JROT,
     2         KKV,KV,KVIN,LPRWF,M,MS,MSAVE,
     3         N,NBEG,NBEGB,NBEG2,NDN,NEND,NENDCH,NLINES,NPR
      REAL*8  BFCT,DE,DEP,DEPRN,DF,DOLD,DSOC,
     2        E,EEPS,EO,EPS,F,FX,GAMA,GI,GN,H,HT,PROD,PPROD,
     3        RATIN,RATOUT,RATST,RH,RINC,RMIN,RMINN,RR,RSTT,RWR(20),
     4        WF(N),SB,SI,SM,SN,SNEND,SPNEND,SRTGI,SRTGN,SWR(20),
     5        V(N),VLIM,VMAX,VMX,VPR,
     6        WKBTST,XEND,XPR,XPW,DXPW,Y1,Y2,Y3,YIN,YM,YOUT
      DATA RATST/1.D-9/,XPW/20.72d0/
      DATA NDN/15/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DXPW= (XPW+ 2.30d0)/NDN
      ICOR= 0
      KVIN= KV
      KV= -1
      RMINN= RMIN-RH
      GAMA= 0.d0
      VMAX= VLIM
      VMX= VMAX*BFCT
      H= RH
      HT= 1.d0/12.D+0
      E= EO*BFCT
      EPS= EEPS*BFCT
      DSOC= VLIM*BFCT
      DE= 0.d0
      RATIN= 0.d0
      RATOUT= 0.d0
      IF(IWR.GT.2) THEN
          IF(KVIN.GE.998) then
              WRITE(6,610) EO
            ELSE
              WRITE(6,601) KVIN,JROT,EO,INNER
            ENDIF
          WRITE(6,602)
        ENDIF
      NEND= N
      IF(KVIN.LT.-10) THEN
          NEND= -KVIN
          SNEND= WF(NEND)
          SPNEND= WF(NEND-1)
          ENDIF
      JQTST = 0
c** Start iterative loop; try to converge for up to 15 iterations.
      DO 90 IT= 1,15
          ITER= IT
          IF(INNER.GT.0) GO TO 38
   10     IF(KVIN.LT.-10) THEN
c** If desired, (KVIN < -10) outer boundary set at NEND=|KVIN| and 
c  initialize wavefunction with log-derivative condition based on value
c  WF(NEND) & derivative SPNEND at that mesh point (brought in in CALL)
              GN= V(NEND)-E
              GI= V(NEND-1)-E
              SB= SNEND
              SI= SB*(1.d0+ 0.5d0*GN)- RH*SPNEND
              GO TO 24
              END IF
          IF(E.GE.DSOC) THEN
c** For quasibound levels, initialize wave function in "QBOUND"
              CALL QBOUND(KVIN,JROT,E,EO,VMX,DSOC,V,RMIN,H,GN,GI,
     1                                 SB,SI,N,ITP3,IWR,IQTST,BFCT,IT)
              NEND= ITP3
              VMAX= VMX/BFCT
              IF(IQTST.GT.0) GO TO 24
              IF(IQTST.LT.0) THEN
                  JQTST = JQTST+IQTST
                  IF((JQTST.LE.-2).OR.(VMAX.LT.VLIM)) GO TO 999
c** Try up to once to find level using trial value just below maximum
                  EO = VMAX-0.1D0
                  E = EO*BFCT
                  GO TO 90
                  ENDIF
              GO TO 20
              ENDIF
c** For  E < DSOC  begin inward integration by using JWKB to estimate
c  optimum (minimum) inward starting point which will still give 
c  RATOUT < RATST = exp(-XPW) (ca. 1.d-9) [not needed after 1'st 2 ITER]
          IF(ITER.LE.2) THEN
              NEND= N
c ... first do rough inward search for outermost turning point
              DO  M= N,1,-NDN
                  MS= M
                  GI= V(M)- E
                  IF(GI.LE.0.D0) GO TO 12
                  GN= GI
                  ENDDO
              IF(IWR.NE.0) WRITE(6,611) JROT,EO
              GO TO 999
   12         IF(MS.GE.N) GO TO 998
              FX= GN/(GI-GN)
              SM= 0.5d0*(1.d0+ FX)*DSQRT(GN)
              MS= MS+ 2*NDN
              IF(MS.GE.N) GO TO 20
c ... now integrate exponent till JWKB wave fx. would be negligible
              DO  M= MS,N,NDN
                  NEND= M
                  SM= SM+ DSQRT(V(M)- E)
                  IF(SM.GT.DXPW) GO TO 18
                  ENDDO
   18         IF(NEND.LT.N) NEND= NEND+ NDN
              ENDIF
c** For truly bound state initialize wave function as 1-st order WKB
c   solution increasing inward
   20     GN= V(NEND)- E
          GI= V(NEND-1)- E
          MS= NEND-1
          IF(GI.LT.0.d0) GO TO 998
          SRTGN= DSQRT(GN)
          SRTGI= DSQRT(GI)
          SB= 1.d0
          SI= SB*DSQRT(SRTGN/SRTGI)*DEXP((SRTGN+SRTGI)*0.5d0)
          IF(SB.GT.SI) THEN
c WOOPS - JWKB gives inward DEcreasing solution, so initialize with node
              IF(IWR.NE.0) WRITE(6,618) JROT,EO,SB/SI
              SI= 1.d0
              SB= 0.d0
              ENDIF
   24     M= NEND-1
          Y1= (1.d0-HT*GN)*SB
          Y2= (1.d0-HT*GI)*SI
          WF(NEND)= SB
          WF(NEND-1)= SI
          MS= NEND
          NENDCH= NEND
          IBEGIN= 3
          IF(INNER.GT.0) IBEGIN= ITP1+2
c** Actual inward integration loop starts here
          DO  I= IBEGIN,NEND
              M= M-1
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(M)-E
              SB= SI
              SI= Y3/(1.d0-HT*GI)
              WF(M)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically
c  forbidden region where  (V(I) .gt. E)
                  SI= 1.d0/SI
                  DO  J= M,MS
                      WF(J)= WF(J)*SI
                      ENDDO
                  NENDCH= MS
                  MS= M
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SB= SB*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
c** Test for outermost maximum of wave function.
c ... old matching condition - turning point works OK & is simpler.
cc            IF((INNER.EQ.0).AND.(SI.LE.SB)) GO TO 32
c** Test for outer well turning point 
              IF((INNER.EQ.0).AND.(GI.lt.0.d0)) GO TO 32
              ENDDO
          IF(INNER.EQ.0) THEN
c** Error mode ... inward propagation finds no turning point
              KV= -2
              IF(IWR.NE.0) WRITE(6,616) KV,JROT,EO
              GO TO 999
              ENDIF
c** Scale outer part of wave function before proceding
   32     SI= 1.d0/SI
          MSAVE= M
          RR= RMINN+MSAVE*H
          YIN= Y1*SI
          RATOUT= WF(NEND)*SI
          NEND= NENDCH
          DO  J= MSAVE,NEND
              WF(J)= WF(J)*SI
              ENDDO
          IF(INNER.GT.0) GO TO 70
c-------------------------------------------------------------------
c** Set up to prepare for outward integration **********************
   38     NBEG= 1
          IF(INNODE.LE.0) THEN
c** Option to initialize with zero slope at beginning of the range
              SB= 1.d0
              GN= V(1)-E
              Y1= SB*(1.d0-HT*GN)
              Y2= Y1+GN*SB*0.5d0
              GI= V(2)-E
              SI= Y2/(1.d0-HT*GI)
            ELSE
c** Initialize outward integration with a node at beginning of range
   40         GN= V(NBEG)-E
              IF(GN.GT.10.D0) THEN
c** If potential has [V(1)-E] so high that H is (locally) much too
c  large, then shift inner starting point outward.
                  NBEG= NBEG+1
                  IF(NBEG.LT.N) GO TO 40
                  IF(IWR.NE.0) WRITE(6,613)
                  GO TO 999
                  ENDIF
              IF((ITER.LE.1).AND.(IWR.NE.0)) THEN
                  IF(NBEG.GT.1) WRITE(6,609) JROT,EO,NBEG
                  IF(GN.LE.0.d0) WRITE(6,604) JROT,EO,NBEG,V(NBEG)/BFCT
                  ENDIF
c** Initialize outward wave function with a node:  WF(NBEG) = 0.
              SB= 0.d0
              SI= 1.d0
              GI= V(NBEG+1)-E
              Y1= SB*(1.d0- HT*GN)
              Y2= SI*(1.d0- HT*GI)
            ENDIF
c
          WF(NBEG)= SB
          WF(NBEG+1)= SI
          NBEGB= NBEG
          NBEG2= NBEG+2
          IF(INNER.GT.0) MSAVE= N
c** Actual outward integration loops start here
          DO  I= NBEG2,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GE.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically forbidden
c  region where  V(I) .gt. E
                  SI= 1.d0/SI
                  NBEG= NBEGB
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  NBEGB= I
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
              ITP1= I
c** Exit from this loop at onset of classically allowed region
              IF(GI.LE.0.d0) GO TO 52
              ENDDO
          MS= MSAVE
          IF((INNER.EQ.0).AND.(GN.LE.0.d0)) GO TO 60
          IF(IWR.NE.0) WRITE(6,612) KVIN,JROT,EO,MSAVE
          GO TO 999
   52     ITP1P= ITP1+1
          MS= ITP1
          IF(INNER.GT.0) GO TO 60
          DO  I= ITP1P,MSAVE
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(I)-E
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GT.1.D+17) THEN
c** Renormalize to prevent overflow of  WF(I) , as needed.
                  SI= 1.d0/SI
                  NBEG= NBEGB
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  NBEGB= I
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
              ENDDO
          MS= MSAVE
c** Finished outward integration.  Normalize w.r.t. WF(MSAVE)
   60     SI= 1.d0/SI
          YOUT= Y1*SI
          YM= Y2*SI
          RATIN= WF(NBEG+1)*SI
          DO  I= NBEG,MS
              WF(I)= WF(I)*SI
              ENDDO
          IF(INNER.GT.0) GO TO 10
c----- Finished numerical integration ... now correct trial energy
c** DF*H  is the integral of  (WF(I))**2 dR
   70     DF= 0.d0
          DO  J= NBEG,NEND
              DF= DF+WF(J)**2
              ENDDO
c** Add edge correction to DF assuming wave function dies off as simple
c  exponential past R(NEND);  matters only if WF(NEND) unusually large.
          IF((E.LE.DSOC).AND.(WF(NEND).NE.0)) THEN
              IF((KVIN.GE.-10).AND.(WF(NEND-1)/WF(NEND).GT.1.d0))
     1              DF= DF+ WF(NEND)**2/(2.d0*DLOG(WF(NEND-1)/WF(NEND)))
              ENDIF
          F= (-YOUT-YIN+2.d0*YM+GI)
          DOLD= DE
          IF(DABS(F).LE.1.D+30) THEN
              DE= F/DF
            ELSE
              F= 9.9D+30
              DF= F
              DE= DABS(0.01D+0 *(DSOC-E))
            ENDIF
          IF(IWR.GT.2) THEN
              DEPRN = DE/BFCT
              XEND= RMINN+NEND*H
c** RATIN & RATOUT  are wave fx. amplitude at inner/outer ends of range
c  relative to its value at outermost extremum.
              WRITE(6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,
     1                                                  XEND,NBEG,ITP1
              ENDIF
c** Test trial eigenvalue for convergence
          IF(DABS(DE).LE.DABS(EPS)) GO TO 100
          E= E+DE
c** KV.ge.998  Option ... Search for highest bound level.  Adjust new
c  trial energy downward if it would have been above dissociation.
          IF((KVIN.GE.998).AND.(E.GT.VMX)) E= VMX- 2.d0*(VMX-E+DE)
          EO= E/BFCT
          IF((IT.GT.4).AND.(DABS(DE).GE.DABS(DOLD)).AND.
     1                                       ((DOLD*DE).LE.0.d0)) THEN
c** Adjust energy increment if having convergence difficulties.  Not
c  usually needed except for some quasibounds extremely near  VMAX .
              ICOR= ICOR+1
              DEP= DE/BFCT
              IF(IWR.NE.0) WRITE(6,617) IT,DEP
              DE= 0.5d0*DE
              E= E-DE
              EO= E/BFCT
              ENDIF
   90     CONTINUE
c** End of iterative loop which searches for eigenvalue ************
c-------------------------------------------------------------------*
c** Convergence fails, so return in error condition
      E= E-DE
      EO= E/BFCT
      DEPRN= DE/BFCT
      IF(IWR.NE.0) WRITE(6,620) KVIN,JROT,ITER,DEPRN
      GO TO 999
  100 IF(IWR.NE.0) THEN
          IF(IWR.GE.3) WRITE(6,619)
          IF((DABS(RATIN).GT.RATST).AND.(INNODE.GT.0)
     1                  .AND.(RMIN.GT.0.d0)) WRITE(6,614) JROT,EO,RATIN
          IF((E.LT.DSOC).AND.(DABS(RATOUT).GT.RATST)) THEN
              WKBTST=0.5d0*DABS(V(NEND)-V(NEND-1))/DSQRT((V(NEND)-E)**3)
              IF(WKBTST.GT.1.d-3)WRITE(6,615)JROT,EO,RATOUT,RATST,WKBTST
              ENDIF
          ENDIF
      KKV = 0
c** Perform node count on converged solution
      PROD= WF(ITP1)*WF(ITP1-1)
      J1= ITP1+1
      J2= NEND-1
      DO  J= J1, J2
          PPROD= PROD
          PROD= WF(J)*WF(J-1)
          IF((PPROD.LE.0.d0).AND.(PROD.GT.0.d0)) KKV= KKV+1
          ENDDO
      KV = KKV
c** Normalize & find interval (NBEG,NEND) where WF(I) is non-negligible
      SN= 1.d0/DSQRT(H*DF)
      DO  I= NBEG,NEND
          WF(I)= WF(I)*SN
          ENDDO
      IF(ITP1.LE.1) GO TO 122
      J= ITP1P
      DO  I= 1,ITP1
          J= J-1
          IF(DABS(WF(J)).LT.RATST) GO TO 119
          ENDDO
  119 NBEG= J
      IF(NBEG.LE.1) GO TO 122
      J= J-1
      DO  I= 1,J
          WF(I)= 0.d0
          ENDDO
  122 IF(KVIN.GE.-10) THEN
c** For "non-wall" cases, move NEND inward to where wavefunction 
c  "non-negligible"
          J= NEND-1
          DO  I= NBEG,NEND
              IF(DABS(WF(J)).GT.RATST) GO TO 126
              J= J-1
              ENDDO
  126     NEND= J+1
          END IF
      IF(NEND.LT.N) THEN
c** Zero out wavefunction array at distances past NEND
          DO  I= NEND+1,N
              WF(I)= 0.d0
              ENDDO
          ENDIF
      IF(LPRWF.LT.0) THEN
c** If desired, write every |LPRWF|-th point of the wave function 
c  to a file on channel-10, starting at the NBEG-th mesh point.
          JPSIQ= -LPRWF
          NPR= 1+(NEND-NBEG)/JPSIQ
          RINC= RH*JPSIQ
          RSTT= RMINN+NBEG*RH
c** Write every JPSIQ-th point of the wave function for level  v=KV
c  J=JROT , beginning at mesh point NBEG & distance RSTT where
c  the NPR values written separated by mesh step RINC=JPSIQ*RH
          WRITE(10,701) KV,JROT,EO,NPR,RSTT,RINC,NBEG,JPSIQ
          WRITE(10,702) (RMINN+I*RH,WF(I),I=NBEG,NEND,JPSIQ)
          GO TO 140
          ENDIF
c** Print solutions every  LPRWF-th  point, 6 to a line, in columns.
      IF(LPRWF.GT.0) THEN
          NLINES= ((1+(NEND-NBEG)/LPRWF)+3)/4
          IPSID= LPRWF*NLINES
          WRITE(6,605) KV,JROT,EO
          DO  J= 1,NLINES
              JJ= NBEG+(J-1)*LPRWF
              IJK= 0
              DO  IJ= JJ,NEND,IPSID
                  IJK= IJK+1
                  RWR(IJK)= RMINN+IJ*H
                  SWR(IJK)= WF(IJ)
                  ENDDO
              WRITE(6,606) (RWR(I),SWR(I),I= 1,IJK)
              ENDDO
          ENDIF
  140 IF(IWR.EQ.1) WRITE(6,607) KV,JROT,EO
cc
cc    IF(IWR.NE.0) WRITE(6,699)  rminn+itp1*rh,eO,rminn+msave*rh,eo
cc699 format('    & turning points:',2(f8.5,f11.4))
cc
      IF(IWR.GE.2) WRITE(6,607) KV,JROT,EO,ITER,RR,RATIN,RATOUT
c** For quasibound levels, calculate width in subroutine "WIDTH"
      IF((E.GT.DSOC).AND.(KVIN.GT.-10)) CALL WIDTH(KV,JROT,E,EO,DSOC,
     1  V,WF,VMX,RMIN,H,BFCT,IWR,ITP1,ITP3,INNER,N,GAMA)
      RETURN
c** ERROR condition if  E.gt.V(R)  at outer end of integration range.
  998 XPR= RMINN+MS*H
      VPR= V(MS)/BFCT
      IF(IWR.NE.0) WRITE(6,608) EO,MS,VPR,XPR,IT
c** Return in error mode
  999 KV= -1
      RETURN
  601 FORMAT(/' Solve for  v=',I3,'   J=',I3,'   ETRIAL=',1PD15.7,
     1   '  INNER=',i2,'   WF(1st) WF(NEND)' )
  602 FORMAT(' ITER    ETRIAL',8X,'F(E)      DF(E)     D(E)',
     1 5X,'M    R(M)  /WF(M)   /WF(M)  R(NEND) NBEG ITP1'/
     2  1X,96('-'))
  603 FORMAT(I4,1PD15.7,3D10.2,I5,0PF7.3,1P2D9.1,0PF8.2,I4,I5)
  604 FORMAT('   NOTE:  for  J=',I3,'   EO=',F12.4,' .ge. V(',i3,')=',
     1  F12.4)
  605 FORMAT(/' Solution of radial Schr. equation for   E(v=',I3,',J=',
     1  I3,') =',F15.7/2x,4('    R(I)   WF(I)   ')/2X,38('--') )
  606 FORMAT(2X,4(F8.3,F11.7))
  607 FORMAT('E(v=',I3,',J=',I3,')=',F11.4,1x,I3,' Iterations',
     1  '   R(M)=',F6.3,'  WF(NBEG)/WF(M)=',1PD8.1/
     2  57x,'WF(NEND)/WF(M)=',D8.1)
  608 FORMAT(' *** SCHRQ Error:  E=',F9.2,' > V(',I5,')=',F9.2,
     1  '  at  Rmax=',F6.2,'  for  IT=',I2)
  609 FORMAT(' *** For  J=',I3,'   E=',1PD15.7,"  integration can't",
     1 ' start till past mesh'/37x,'point',I5,',  so RMIN smaller than n
     2eeded')
  610 FORMAT(/' Attempt to find the highest bound level starting from',
     1 '   ETRIAL =',1PD9.2)
  611 FORMAT(' *** SCHRQ inward search at   J=',i3,'   E=',f11.2,
     1  ' finds no classical region')
  612 FORMAT(/' *** ERROR *** for   v =',I3,'   J =',I3,'   E =',
     1  F12.4,'  Innermost turning point not found by   M = MSAVE =',I5)
  613 FORMAT(/' *** ERROR in potential array ... V(I) everywhere',
     1 ' too big to integrate with given  increment')
  614 FORMAT(' *** CAUTION *** For  J=',I3,'  E=',G15.8/16x,
     1 'WF(first)/WF(Max)=',D9.2,'  suggests  RMIN  may be too large')
  615 FORMAT(' ** CAUTION ** For  J=',I3,'  E=',1PD13.6,
     1 '  WF(NEND)/WF(Max)=',D8.1,' >',D8.1/4X,'& initialization ',
     2 'quality test ',1PD8.1,' > 1.D-3   so RMAX may be too small')
  616 FORMAT(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',G14.7,
     1  ':  inward propagation finds no turning point ... Energy too low
     2 or potential too weak' )
  617 FORMAT(' *** SCHRQ has a convergence problem, so for  IT=',I2,
     1 '  cut  DE=',1PD10.2,'  in HALF' )
  618 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  JWKB start gives  SB/SI=',
     1  1PD10.3,'  so use a node.')
  619 FORMAT(1X,96('-'))
  620 FORMAT(' *** CAUTION for  v=',I3,'  J=',I3,"  SCHRQ doesn't conver
     1ge by  ITER=",I2,'  DE=',1PD9.2)
  701 FORMAT(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave funct
     1ion at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,
     2  '   NBEG=',I4,'   |LPRWF|=',I3)
  702 FORMAT((1X,4(f9.4,f10.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c*******************************************************************
      SUBROUTINE QBOUND(KV,JROT,E,EO,VMX,DSOC,V,RMIN,H,GB,GI,SB,SI,N,
     1  ITP3,IWR,IQTST,BFCT,IT)
c*******************************************************************
c** Subroutine to initialize quasibound level wave function as Airy
c  function at third turning point (if possible). For the theory see 
c  J.Chem.Phys. 54, 5114 (1971),  J.Chem.Phys. 69, 3622-31 (1978) 
c----------------------------------------------------------------------
c** IQTST  is error flag. *** If (IQTST.lt.0) initialization fails
c  so eigenvalue calculation aborts *** (IQTST.gt.0) for successful
c  Airy function initialization. *** (IQTST=0) if Airy function
c  initialization prevented because 3-rd turning point beyond
c  range, so that WKB initialization is used.
c----------------------------------------------------------------------
      INTEGER I,II,IQTST,IT,ITP3,IWR,J,JROT,K,KV,N
      REAL*8  A1,A2,A13,A23,BFCT,
     1        C1A,C2A,DF,DSOC,E,EO,FBA,FIA,FJ,GB,GBA,GI,GIA,H,
     2        RMIN,RMINN,SB,SI,SL,V(N),VMX,VMXPR,XJ1
      DATA C1A/0.355028053887817D0/,C2A/0.258819403792807D0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IQTST=1
      RMINN=RMIN-H
c** Start by searching for third turning point.
      J=N
      IF(V(N).GT.E) GO TO 22
      DO  I=2,N
          J=J-1
          IF(V(J).GT.E) GO TO 10
          ENDDO
      GO TO 14
c** Check that there is a classically allowed region inside this point
c  and determine height of barrier maximum.
   10 II=J
      VMX=DSOC
      DO  I=2,J
          II=II-1
          IF(V(II).LE.E) GO TO 16
          IF(V(II).GT.VMX) VMX=V(II)
          ENDDO
c** Energy too high ... find no more than one turning point.
   14 XJ1=RMINN+J*H
c ... Search outward for barrier height to facilitate energy correction
      IF(J.EQ.1) J= 2
      K=J-1
      DO  I=J,N
          IF(V(I).GT.V(K)) GO TO 120
          K=I
          ENDDO
      VMX=V(K)
      GO TO 130
  120 K=K+2
      J=K-1
      DO  I=K,N
          IF(V(I).LT.V(J)) GO TO 126
          J=I
          ENDDO
  126 VMX=V(J)
  130 VMXPR=VMX/BFCT
      IF(IWR.NE.0) WRITE(6,608) JROT,EO,VMXPR,XJ1
      ITP3= J
      IQTST=-1
      GO TO 100
   16 ITP3= J+1
c** ITP3 is the first mesh point outside classically forbidden region
      GB=V(ITP3)-E
      GI=V(ITP3-1)-E
      FJ=GI/(GI-GB)
c** Treat quasibound levels as bound using outer boundary condition
c  of Airy function at third turning point ... as discussed by
c  R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
      SL=(GI-GB)**(1.d0/3.d0)/H
      IF((SL*H).LT.1.d0) THEN
          A1=GI/(SL*H)**2
          A2=GB/(SL*H)**2
          A13=A1*A1*A1
          A23=A2*A2*A2
          FIA= 1.d0+ A13*(A13*(A13+72.D0)+2160.D0)/12960.D0
          GIA=A1+A1*A13*(A13*(A13+90.D0)+3780.D0)/45360.D0
          FBA= 1.d0+ A23*(A23*(A23+72.D0)+2160.D0)/12960.D0
          GBA=A2+A2*A23*(A23*(A23+90.D0)+3780.D0)/45360.D0
c** Airy function  Bi(X)  at points straddling 3-rd turning point
          SI=C1A*FIA+C2A*GIA
          SB=C1A*FBA+C2A*GBA
          GO TO 100
          ENDIF
c** If Airy function expansion unreliable, use zero slope at third
c  turning point as quasibound outer boundary condition.
      DF=GI-GB
      SI= 1.d0+ DF*FJ**3/6.d0
      SB= 1.d0 -DF*(1.d0- FJ)**3/6.d0
      IF(IWR.NE.0) WRITE(6,606) KV,JROT,EO,IT
      GO TO 100
c** If 3-rd turning point beyond range start with WKB wave function
c  at end of range.
   22 IF(IWR.NE.0) WRITE(6,607) JROT,EO
      ITP3= N
      IQTST=0
      GB=V(ITP3)-E
      GI=V(ITP3-1)-E
      VMX=V(ITP3)
      II=ITP3
      DO  I=2,ITP3
          II=II-1
          IF(V(II).LT.VMX) GO TO 100
          VMX=V(II)
          ENDDO
      IF(IWR.NE.0) WRITE(6,604)
c** End of quasibound level initialization schemes.
      IQTST=-9
  100 RETURN
  604 FORMAT(" **** QBOUND doesn't work ... no classically allowed regio
     1n accessible at this energy.")
  606 FORMAT(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',1PD13.6,
     1 '   IT=',I2/5x,'Airy initialization unstable so use  zero slope',
     2 'at  R(3-rd)' )
  607 FORMAT(' *** For  J=',I3,'  E=',F9.2,
     1  '  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
  608 FORMAT(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,
     1  '  find onee turn point:  R=',F6.2)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
c** Subroutine to calculates quasibound level tunneling lifetime/width
c** For relevant theory see Le Roy & Liu [J.Chem.Phys.69,3622-31(1978)]
c  and Connor & Smith [Mol.Phys. 43, 397 (1981)] and Huang & Le Roy 
c  [J.Chem.Phys. 119, 7398 (2003); Erratum, ibid, 127, xxxx (2007)]
c** Final level width calculation from Eq.(4.5) of Connor & Smith.
c  Rearranged slightly for consistency with PotFit derivatives 9/05/02
c-----------------------------------------------------------------------
      SUBROUTINE WIDTH(KV,JROT,E,EO,DSOC,V,S,VMX,RMIN,H,BFCT,IWR,ITP1,
     1  ITP3,INNER,N,GAMA)
c++ "WIDTH" calls subroutine "LEVQAD" ++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IMM,INNER,IRM,ITP1,ITP1P,ITP1P1,ITP2,ITP2M,ITP2M2,
     1         ITP2P1,ITP2P2,ITP3,IWR,JROT,KV,KVI,KVO,
     2         M,M2,N,NN,NST
      REAL*8  ANS1,ANS2,ARG,BFCT,COR,
     1        D1,D2,D3,DFI,DSGB,DSGN,DSOC,DWEB,OMEGJC,
     2        E,EO,EMSC,EMV,G1,G2,G3,GA,GAMA,GAMALG,
     3        H,H2,HBW,HBWB,PI,PMX,RMIN,RMINN,RMX,RT,RT1,RT2,
     4        S(N),SM,TAU,TAULG,TI,TUN0,U1,U2,V(N),VMAX,VMX,
     7        XJ,XX
      CHARACTER*5 LWELL(2)
      DATA PI/3.141592653589793D0/
      DATA LWELL/'INNER','OUTER'/
      RMINN= RMIN- H
      H2= H*H
c** ITP1 is first mesh point to right of innermost turning point.
   40 ITP1P= ITP1+ 1
      ITP1P1= ITP1P+ 1
      IRM= ITP1- 1
c** Calculate JWKB tunneling probability from quadrature over barrier
c** First must locate 2-nd turning point.
      DO  I= ITP1P1,ITP3
          ITP2= I
          IF(V(I).GT.E) GO TO 202
          ENDDO
      GAMA= 0.d0
      GO TO 250
  202 ITP2P1= ITP2+ 1
      ITP2P2= ITP2+ 2
c** ITP2M is the last mesh point before the 2-nd turning point.
      ITP2M= ITP2- 1
      ITP2M2= ITP2- 2
      G1= V(ITP2M)- E
      G2= V(ITP2)- E
      GA= V(ITP2P1)- E
c** Quadrature over barrier starts here.
      CALL LEVQAD(G1,G2,GA,H,RT,ANS1,ANS2)
      SM= ANS2/H
      IF(GA.LT.0.d0) GO TO 218
      SM= SM+ 0.5d0*DSQRT(GA)
      PMX= VMX
      M2= ITP2P2
  204 DO  I=M2,ITP3
          M= I
          GA= V(I)- E
          IF(V(I).GT.PMX) PMX=V(I)
          IF(GA.LT.0.d0) GO TO 210
          SM= SM+ DSQRT(GA)
          ENDDO
      IF(V(M).GT.V(M-1)) THEN
          IF(IWR.NE.0) WRITE(6,602) KV,JROT
          GO TO 250
          ENDIF
      RMX= RMINN+ M*H
      U1= DSQRT(GA/(V(M)- DSOC))
      U2= DSQRT((E- DSOC)/(V(M)- DSOC))
      SM= SM- 0.5d0*DSQRT(GA)+ (DLOG((1.d0+U1)/U2)-U1)*RMX*
     1                                             DSQRT(V(M)- DSOC)/H
      XJ= (DSQRT(1.d0+ 4.d0*(V(M)-DSOC)*(RMX/H)**2)- 1.d0)*0.5d0
      IF(IWR.NE.0) WRITE(6,603) JROT,EO,XJ,RMX
      GO TO 218
  210 IF(M.LT.ITP3) THEN
c** If encounter a double-humped barrier, take care here.
          IF(IWR.NE.0) WRITE(6,609) KV,JROT,EO,M
          KVO= 0
          DSGN= DSIGN(1.d0,S(M-1))
c** Find the effective quantum number for the outer well
          DO  I= M,ITP3
              DSGB= DSGN
              DSGN= DSIGN(1.d0,S(I))
              IF((DSGN*DSGB).LT.0.d0) KVO=KVO+1
              ENDDO
          KVI= KV- KVO
          IF(INNER.EQ.0) THEN
c** For levels of outer well, get correct width by changing ITP1
              ITP1= M
              IF(IWR.GT.0) WRITE(6,610) KVO,LWELL(2)
              GO TO 40
              ENDIF
          IF(IWR.GT.0) WRITE(6,610) KVI,LWELL(1)
c** For "inner-well" levels, locate outer barrier
          DO  I= M,ITP3
              M2= I
              GA= V(I)- E
              IF(GA.GE.0.d0) GO TO 204
              ENDDO
          GO TO 218
          ENDIF 
      G3= V(M-2)- E
      G2= V(M-1)- E
      CALL LEVQAD(GA,G2,G3,H,RT,ANS1,ANS2)
      SM= SM- 0.5d0*DSQRT(G3)-DSQRT(G2) + ANS2/H
  218 EMSC= -SM/PI
      IF(INNER.GT.0) VMX= PMX
      VMAX= VMX/BFCT
c** Tunneling factors calculated here ** TUN0 is simple WKB result
c  as in Child's eqs.(57c) & (59).
c .....  EPSRJ= -2.* PI* EMSC 
      TUN0= 0.5d0*DEXP(2.d0*PI*EMSC)
c ... for permeability calculate Connor-Smith's Eq.(3.7) \omega=OMEGJC 
      OMEGJC= DSQRT(1.d0+ 2.d0*TUN0) - 1.d0
c ... alternate calculation to give better precision for small TUN0
      IF(TUN0.LT.1.d-5) OMEGJC= TUN0*(1.d0-0.5d0*TUN0*(1.d0-TUN0))
      OMEGJC= 4.d0*OMEGJC/(OMEGJC + 2.d0)
c** Quadrature for JWKB calculation of vibrational spacing in well HBW
      D1= E- V(IRM)
      D2= E- V(ITP1)
      D3= E- V(ITP1P)
      CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      RT1= RT
      SM= ANS1/H
      IF(D3.LT.0.d0) GO TO 228
      SM= SM+ 0.5d0/DSQRT(D3)
      DO  I= ITP1P1,ITP2M2
          IMM= I
          EMV= E- V(I)
          IF(EMV.LT.0.d0) GO TO 222
          SM= SM+ 1.d0/DSQRT(EMV)
          ENDDO
      D3= E- V(ITP2M2)
      D2= E- V(ITP2M)
      D1= E- V(ITP2)
      GO TO 226
c** If encounter a double-minimum well, take care here.
  222 D1= EMV
      D2= E- V(IMM-1)
      D3= E- V(IMM-2)
      IF(IWR.NE.0) WRITE(6,605) KV,JROT,EO
  226 CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      RT2=RT
      SM=SM-0.5d0/DSQRT(D3) + ANS1/H
c** Get HBW in same energy units (1/cm) associated with BFCT
  228 HBW=2.d0*PI/(BFCT*SM)
c** HBW fix up suggested by Child uses his eqs.(48)&(62) for HBW
c** Derivative of complex gamma function argument calculated as
c  per eq.(6.1.27) in Abramowitz and Stegun.
      NST= INT(DABS(EMSC)*1.D2)
      NST= MAX0(NST,4)
      ARG= -1.963510026021423d0
      DO  I= 0,NST
          NN= I
          XX= I + 0.5d0
          TI= 1.d0/(XX*((XX/EMSC)**2 + 1.d0))
          ARG= ARG+TI
          IF(DABS(TI).LT.1.D-10) GO TO 233
          ENDDO
c ... and use continuum approximation for tail of summation (???)
  233 COR= 0.5d0*(EMSC/(NN+1.d0))**2
      ARG= ARG+ COR- COR**2
c** Now use WKL's Weber fx. approx for (?) derivative of barrier integral ..
      DWEB= (EO-VMAX)*BFCT/(H2*EMSC)
      DFI= (DLOG(DABS(EMSC)) - ARG)*BFCT/(H2*DWEB)
      HBWB= 1.d0/(1.d0/HBW + DFI/(2.d0*PI))
c** Width from formula (4.5) of  Connor & Smith, Mol.Phys.43,397(1981)
c [neglect time delay integral past barrier in their Eq.(4.16)].
      IF(EMSC.GT.-25.D0) THEN
          GAMA= (HBWB/(2.d0*PI))* OMEGJC
          TAU= 0.D0
          IF(GAMA.GT.1.D-60) TAU= 5.308837457D-12/GAMA
c** GAM0 = TUN0*HBW/PI  is the simple WKB width GAMMA(0) discussed by
c  Le Roy & Liu in J.C.P.69,3622(1978).
          IF(IWR.GT.0) WRITE(6,601) TAU,GAMA,HBWB,VMAX
        ELSE
          GAMALG= DLOG10(HBWB/(2.d0*PI))+2.d0*PI*EMSC/2.302585093D0
          TAULG= DLOG10(5.308837457D-12)-GAMALG
          IF(IWR.GT.0) WRITE(6,611) TAULG,GAMALG,HBWB,VMAX
        ENDIF
  250 RETURN
  601 FORMAT('    Lifetime=',1PD10.3,'(s)   Width=',D10.3,'   dG/dv=',
     1 0PF7.2,'   V(max)=',F9.2)
  602 FORMAT(' *** WARNING ***  For   v =',I3,'   J =',I3,'   cannot cal
     1culate width since barrier maximum beyond range')
  603 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) beyond range so tu
     1nneling calculation uses'/8X,'pure centrifugal potential with  J(a
     2pp)=',F7.2,'  for  R > R(max)=',F7.2)
  605 FORMAT(' **** CAUTION *** Width estimate only qualitative, as have
     1 a double-minimum well for   E(v=',I3,', J=',I3,')=',F15.7/15X,
     2 'a more stable result may be obtained by searching for the quasib
     3ound levels using option: INNER > 0 .')
  609 FORMAT(' *** CAUTION - Permeability estimate not exact as have a d
     1ouble-humped barrier:  E(v=',I3,', J=',I3,') =',G15.8,I6)
  610 FORMAT(16X,'(NOTE: this has the node count of a   v=',I3,2X,A5,
     1 '-well level')
  611 FORMAT(12X,'Log10(lifetime/sec)=',F10.5,' ;   Log10(width/cm-1)=',
     1 F10.5,'   Spacing=',G12.5,'   V(max)=',G14.7,'(cm-1)')
      END
c**********************************************************************
      SUBROUTINE LEVQAD(Y1,Y2,Y3,H,RT,ANS1,ANS2)
c** Subroutine "LEVQAD" fits quadratic  Y = A + B*X + C*X**2  through
c  function values  Y1, Y2, Y3  at equally spaced points separated by
c  distance H, where  Y1 < 0  and (Y2,Y3 .ge.0), locates the function
c  zero (at RT, relative to  X1 < X2 = 0) between points X1 & X2, and
c  evaluates the integral from RT to R3 of   1/sqrt(Y)  , called
c  ANS1, and the integral (same range) of  sqrt(Y) , which is ANS2
c** Alternately, if Y1 & Y3 both  < 0  and only the middle point
c  Y2.ge.0 ,   fit the points to:  Y = A - B*(X-X0)**2 , locate the
c  turning points between which  Y(X) > 0  and evaluate these integrals
c  on this interval.  *************************************************
c----------------------------------------------------------------------
      REAL*8  A,ANS1,ANS2,B,C,CQ,H,HPI,R1,R2,RCQ,RR,RT,SL3,SLT,
     1        X0,Y1,Y2,Y3,ZT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DATA HPI/1.570796326794896D0/
      IF((Y1.GE.0).OR.(Y2.LT.0)) GO TO 99
      IF(Y3.LT.0.d0) GO TO 50
c** Here treat case where both 'Y2' & 'Y3' are positive
      IF(DABS((Y2-Y1)/(Y3-Y2) -1.D0).LT.1.d-10) THEN
c ... special case of true (to 1/10^10) linearity ...
          RT= -H*Y2/(Y2-Y1)
          ANS1= 2.d0*(H-RT)/DSQRT(Y3)
          ANS2= ANS1*Y3/3.D0
          RETURN
          ENDIF
      C= (Y3-2.d0*Y2+Y1)/(2.d0*H*H)
      B= (Y3-Y2)/H-C*H
      A= Y2
      CQ= B**2- 4.d0*A*C
      RCQ= DSQRT(CQ)
      R1= (-B-RCQ)/(2.d0*C)
      R2= R1+ RCQ/C
      IF((R2.LE.0.d0).AND.(R2.GE.-H)) RT=R2
      IF((R1.LE.0.d0).AND.(R1.GE.-H)) RT=R1
      SL3= 2.d0*C*H+B
      SLT= 2.d0*C*RT+B
      IF(C.LT.0.d0) GO TO 10
      ANS1= DLOG((2.d0*DSQRT(C*Y3)+SL3)/SLT)/DSQRT(C)
      GO TO 20
   10 ANS1= -(DASIN(SL3/RCQ)- DSIGN(HPI,SLT))/DSQRT(-C)
   20 ANS2= (SL3*DSQRT(Y3)- CQ*ANS1/2.d0)/(4.d0*C)
      IF(RT.GE.H) WRITE(6,601) H,R1,R2
  601 FORMAT(' *** CAUTION *** in LEVQAD, turning point not between poin
     1ts 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)
      RETURN
c** Here treat case when only 'Y2' is non-negative
   50 RR= (Y2-Y1)/(Y2-Y3)
      X0= H*(RR-1.d0)/((RR+1.d0)*2.d0)
      B= (Y2-Y1)/(H*(2.d0*X0+H))
      A= Y2+ B*X0**2
      ZT= DSQRT(A/B)
      RT= X0- ZT
      ANS1= 2.d0*HPI/DSQRT(B)
      ANS2= ANS1*A*0.5d0
      RETURN
   99 WRITE(6,602) Y1,Y2
  602 FORMAT(' *** ERROR in LEVQAD *** No turning point between 1-st two
     1 points as   Y1=',D10.3,'   Y2=',D10.3)
      ANS1= 0.d0
      ANS2= 0.d0
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PREPOT(LNPT,IAN1,IAN2,IMN1,IMN2,NPP,OMEGA,RR,RM2,VLIM,
     1                                                        VV,NCN)
c** Driver subroutine of package to read parameters and/or generate
c  values of a potential V(I) at the NPP input distances RR(I).
c====================== Version of  8 April 2007 ========================
c**** Subroutine Input:
c----------------------
c  LNPT  is an integer specifying the operational mode:
c      *  LNPT > 0  : for a new case for which all potential-defining
c                     parameters are read in & a description printed
c      *  LNPT.le.0 : if potential points are to be generated in exactly
c                     the same manner as on preceding call, but at
c                     different distances RR(I) (no reads or writes)
c  IAN1 & IAN2 are the atomic numbers and IMN1 & IMN2 the mass numbers
c        of atoms #1 & 2, used (if needed) to specify isotope masses for
c        calculating adiabatic and/or non-adiabatic BOB correction fx.
c  NPP (integer) is the number of input distances  RR(i) (in Angstroms)
c        at which potential values  VV(i) (in cm-1) are to be generated
c  RR  (real array) is set of NPP distances where potential calculated
c  RM2 (real array) on input is the (centrifugal) array of  1/RR(i)**2
c----------------------
c**** Subroutine Output:
c----------------------
c  OMEGA the (integer) elelectronic angular momentum projection q.no.
c  VLIM (cm-1) is the absolute energy at the potential asymptote
c  VV (real array) is the set of function values generated (in cm-1)
c  RM2 values returned may (if appropriate) be modified to include BOB
c      corrections to the (centrifugal) potential  1/RR(i)**2
c  NCN is an integer power defining the asymptotically-dominant 
c      inverse-power long-range potential tail:  CNN/R**NCN 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Calls GENINT (which calls PLYINTRP, SPLINT & SPLINE) ,  or POTGEN ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set maximum array dimension for the input function values to be
c  interpolated over & extrapolated beyong
      INTEGER NTPMX
      PARAMETER (NTPMX= 1600) 
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LPPOT,LWR,
     1  NCN,NLIN,NPP,NROW,NTP,NUSE, OMEGA
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP),
     1  XI(NTPMX),YI(NTPMX),RWR(20),VWR(20),VWRB(3),D1V(3),D1VB(3),
     2  D2V(3),CNN
c
c** Save variables needed for 'subsequent' LNPT.le.0 calls
      SAVE ILR,IR2,LPPOT,NTP,NUSE
      SAVE CNN,VSHIFT,XI,YI
c
      DATA VWRB/3*0.D0/,D1VB/3*0.D0/
      LPPOT= 0
c
      IF(LNPT.GT.0) THEN
c** If NTP > 0 :  define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
c   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
c** If LPPOT > 0 : at every |LPPOT|-th point, print potential and 
c      derivatives-by-differences. ***  If  LPPOT < 0  write potential
c      at every |LPPOT|-th point to channel-8 in a compact format **
c  OMEGA  the (integer) total elextronic angular momentum projection
c            quantum number (required for proper rotational intensities)
c** VLIM (cm-1) is the energy associated with the potential asymptote.
c-----------------------------------------------------------------------
          READ(5,*) NTP, LPPOT, OMEGA, VLIM
c-----------------------------------------------------------------------
          WRITE(6,600) OMEGA,VLIM
          IF(NTP.GT.0) THEN
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP (read above) is number of turning points (XI,YI) to be read in.
c** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
c    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
c    interpolate with cubic spline instead of local polynomials.
c** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
c   This may help if interpolation has trouble on steep repulsive wall.
c** ILR specifies how to extrapolate beyond largest input distance XI(i)
c  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
c  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
c* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c* If ILR > 3 : successive higher power terms differ by factor  1/R
c-----------------------------------------------------------------------
              READ(5,*) NUSE, IR2, ILR, NCN, CNN
c-----------------------------------------------------------------------
              IF(NTP.GT.NTPMX) THEN
c** If interpolation being requested, but the number of input points 
c   exceeds the array size, print a warning and stop.
                  IF(NUSE.GE.0) THEN
                      WRITE(6,602) NTP,NTPMX
                      STOP
                      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++
c** IF NTP > NTPMX, and NUSE < 0  read in the full final array of  NTP
c  mesh points {RR(i),VV(i)} in Angst., cm-1.
                  NPP= NTP
c-----------------------------------------------------------------------
                  READ(5,*) VSHIFT
                  READ(5,*) (RR(I),VV(I),RM2(I),I= 1, NPP)
c-----------------------------------------------------------------------
                  DO  I= 1, NPP
                      RM2(I)= (1.d0 + RM2(I))/RR(I)**2
                      VV(I)= VV(I) + VSHIFT
                      ENDDO
                  RMIN= RR(1)
                  RH= RR(2)- RR(1)
                  WRITE(6,626) NPP,RMIN,RR(NPP),RH
  626 FORMAT(/'  Potential and centrifugal BOB functions defined by',
     1  I6,'-point input arrays'/5x,' on the range ',f9.6,'  to',f11.6,
     2  '[Angst.]  with mesh   RH=',F9.6/'     and no interpolation')
                  WRITE(6,612) VSHIFT
                  GOTO 20
                  ENDIF
              IF(NUSE.GT.0) WRITE(6,604) NUSE,NTP
              IF(NUSE.LE.0) WRITE(6,606) NTP
              IF(IR2.GT.0) WRITE(6,608)
              IF((ILR.GT.1).AND.(DABS(CNN).GT.0.D0))WRITE(6,610)CNN,NCN
c** Read in turning points to be interpolated over
c** RFACT & EFACT are factors required to convert units of input turning
c       points (XI,YI) to Angstroms & cm-1, respectively (may = 1.d0)
c** Turning points (XI,YI) must be ordered with increasing XI(I)
c** Energy VSHIFT (cm-1) is added to the input potential points to
c   make their absolute energy consistent with VLIM (often VSHIFT=Te).
c-----------------------------------------------------------------------
              READ(5,*) RFACT, EFACT, VSHIFT
              READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-----------------------------------------------------------------------
              WRITE(6,612) VSHIFT, RFACT, EFACT
              NROW= (NTP+2)/3
              DO  J= 1,NROW
                  IF(EFACT.LE.10.D0) THEN
                      WRITE(6,614) (XI(I),YI(I),I= J,NTP,NROW)
                    ELSE
                      WRITE(6,616) (XI(I),YI(I),I= J,NTP,NROW)
                    ENDIF
                  ENDDO
                  WRITE(6,624)
              DO  I= 1,NTP
                  YI(I)= YI(I)*EFACT+ VSHIFT
                  XI(I)= XI(I)*RFACT
                  ENDDO
              IF(IR2.GT.0) THEN
                  DO  I= 1,NTP
                      YI(I)= YI(I)*XI(I)**2
                      ENDDO
                  ENDIF
              IF((DABS(YI(NTP)-YI(NTP-1)).LE.0).AND.
     1                              (XI(NTP).LT.RR(NPP))) WRITE(6,618)
              ENDIF
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NTP.GT.0) THEN
          CALL GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                                        NCN,CNN)
        ELSE
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If (NTP.le.0) PREPOT uses subroutine POTGEN to generate a fully
c  analytic potential defined by the following read-in parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
c  and for all cases except IPOTL=2, the potential asymptote energy is
c  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
c  potential minimum and  DSCM  the leading (quadratic) potential coeft.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR, NSR & NCMM are integers cwcharacterizing the chosen potential
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass 
c      dependent Born-Oppenheimer breakdown corrections will be included
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NCN) potential.
c** If IPOTL=2  use Seto's modification of Surkus' GPEF expansion in
c       z = [R**MPAR - Re**MPAR]/[a*R**MPAR + b*Re**MPAR] where 
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NCN=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NCN=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NCN=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus MPAR < 0 case:  z(MPAR,a,b)= z(|MPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c      with exponent factor "beta" defined as a power series of order
c      (NVARB-1) in  y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR)
c      with NVARB coefficients PARM(i).    [!! MPAR .ge.1 !!]
c    * For conventional "simple" Morse potential,  NVARB=1 & MPAR dummy
c*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter 
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  generate an MLR potential [Mol.Phys. 105, 691 (2007)]
c      If MPAR > 0  exponent parameter defined in terms of a polynomial
c           of order (NVARB-2) with the NVARB coefficients  PARM(j).
c           in y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR), 
c           and long-range defined by NCN inverse-power terms
c      If MPAR = 0  exponent polynomial variable is  2*y_{1}= y{O-T}
c      If MPAR < 0  exponent polynomial variable is  y_{|MPAR|}
c      If MPAR.le.0  exponent polynomial connected to limiting inverse-
c           power potential exponent by exponential switching function 
c           with parameters  Asw= PARM(NVARB-1)  and  RSW= PARM(NVARB).
c** If IPOTL=5  generate a Double-Exponential Long-Range (DELR) 
c       potential [JCP 119, 7398 (2003)] with additive long-range part
c       defined by a sum of NCMM damped inverse-power terms, & exponent
c       polynomial radial variable defined by parameter MPAR (=p)
c** If IPOTL=6  generate generalized HFD({m_i},i=1,NCMM) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c** If IPOTL=7  use Tiemann-type polynomial potential attached to an
c     inverse-power long-range tail and an 1/R^{12} (or exponential)
c     inner wall. 
c----------------------------------------------------------------------
c++     READ(5,*) IPOTL, MPAR, NSR, NCMM, NVARB, IBOB, DSCM, REQ
c++     IF(IPOTL.GE.4) READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
c++     IF((IPOTL.EQ.4).OR.(IPOT.EQ.7)) READ(5,*) (MMLR(I),CMM(I),I= 1,MPAR)
c++     IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
c++     IF(IBOB.GT.0) THEN
c++         READ(5,*) MN1R, MN2R, PAD, MAD, NU1, NU2, PNA, NT1, NT2
c++         IF(PAD.GT.0) THEN
c++              IF(NU1.GE.0) READ(5,*) U1INF, (U1(I), I=0,NU1)
c++              IF(NU2.GE.0) READ(5,*) U2INF, (U2(I), I=0,NU2)
c++              IF(NT1.GE.0) READ(5,*) T1INF, (T1(I), I=0,NT1)
c++              IF(NT2.GE.0) READ(5,*) T2INF, (T2(I), I=0,NT2)
c++           ELSE
c++              IF(NU1.GE.0) READ(5,*) U1INF,(U1(I), I=0,NU1), Aad1, Rad1
c++              IF(NU2.GE.0) READ(5,*) U2INF,(U2(I), I=0,NU2), Aad2, Rad2
c++              IF(NT1.GE.0) READ(5,*) T1INF,(T1(I), I=0,NT1), Ana1, Rna1
c++              IF(NT2.GE.0) READ(5,*) T2INF,(T2(I), I=0,NT2), Ana2, Rna2
c++           ENDIF
c++         ENDIF
c++     ENDIF
c-----------------------------------------------------------------------
          NCN= 99
          CALL POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,
     1                                                        NCN,CNN)
        ENDIF
   20 IF(LPPOT.NE.0) THEN
c** If desired, on the first pass (i.e. if LNPT > 0) print the potential
          RH= RR(2)-RR(1)
          INPTS= IABS(LPPOT)
          IF(LPPOT.LT.0) THEN
c** Option to write resulting function compactly to channel-8. 
              RMIN= RR(1)
              NLIN= NPP/INPTS+ 1
              WRITE(8,800) NLIN,VLIM
              WRITE(8,802) (RR(I),VV(I),I= 1,NPP,INPTS)
            ELSE
c** Option to print potential & its 1-st three derivatives, the latter
c  calculated by differences, assuming equally spaced RR(I) values.
              WRITE(6,620)
              NLIN= NPP/(2*INPTS)+1
              RH= INPTS*RH
              DO  I= 1,NLIN
                  LWR= 1+ INPTS*(I-1)
                  DO  J= 1,2
                      JWR= LWR+(J-1)*NLIN*INPTS
                      IF(JWR.LE.NPP) THEN
                          RWR(J)= RR(JWR)
                          VWR(J)= VV(JWR)
                          D1V(J)= (VWR(J)-VWRB(J))/RH
                          VWRB(J)= VWR(J)
                          D2V(J)= (D1V(J)-D1VB(J))/RH
                          D1VB(J)= D1V(J)
                        ELSE
                          RWR(J)= 0.d0
                          VWR(J)= 0.d0
                        ENDIF
                      IF(I.LE.2) THEN
                          D2V(J)= 0.d0
                          IF(I.EQ.1) D1V(J)= 0.d0
                          ENDIF
                      ENDDO
                  WRITE(6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J= 1,2)
                  ENDDO
            ENDIF
          ENDIF
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' State has  OMEGA=',i2,'   and energy asymptote:   Y(lim)=
     1',F12.4,'(cm-1)')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'
     1 ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  604 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  606 FORMAT(' Perform cubic spline interpolation over the',I5,
     1  ' input points' )
  608 FORMAT(' Interpolation actually performed over modified input arra
     1y:   Y(I) * r(I)**2')
  610 FORMAT( ' Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(r)  =  Y(lim) - (',D16.7,')/r**',I2)
  612 FORMAT(' To make input points Y(i) consistent with  Y(lim),  add'
     1 ,'  Y(shift)=',F12.4:/' Scale input points:  (distance)*',
     2 1PD16.9,'  &  (energy)*',D16.9/13x,'to get required internal unit
     3s  [Angstroms & cm-1 for potentials]'/
     4  3('      r(i)         Y(i)  ')/3(3X,11('--')))
  614 FORMAT((3(F13.8,F12.4)))
  616 FORMAT((3(F12.6,F13.8)))
  618 FORMAT(/' !!! CAUTION !!! Last two mesh point  YI  values are equa
     1l'/17x,'so extrapolation to large  r  will be unreliable !!!'/)
  620 FORMAT(/'  Function and first 2 derivatives by differences'/
     1  2('     r       Y(r)     d1Y/dr1    d2Y/dr2')/2(2X,19('--')))
  622 FORMAT(2(0PF8.3,F11.3,1PD11.3,D10.2))
c 622 FORMAT(2(0PF7.2,F12.5,1PD11.3,D10.2))
  624 FORMAT(1x,38('--'))
  800 FORMAT(I7,' function values with asymptotic value:',F14.6)
  802 FORMAT((1X,3(F12.8,F14.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GENINT(LNPT,NPP,XX,YY,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                                        NCN,CNN)
c** GENINT produces a smooth function YY(i) at the NPP input distances
c  XX(i) by performing numerical interpolation over the range of the 
c  NTP input function values YI(j) at the distances XI(j), and using
c  analytic functions to extrapolate beyond their range to with an
c  exponential at short range and a form specified by ILR, NCN & CNN
c** ILR specifies how to extrapolate beyond largest given turning pts
c   If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c   If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c   If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c* If(ILR.ge.2) fit last turning points to:  VLIM - sum(of ILR
c  inverse-power terms beginning with  1/R**NCN). *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  ((cm-1)(Angstroms)**'NCN').
c  If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c  If ILR > 3 : this factor is  1/R .
c=== Calls subroutines PLYINTRP, SPLINT & SPLINE ==================
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,J,IFXCN,IDER,IR2,ILR,ISR,LNPT,MBEG,MFIN,MINNER,
     1  NN,NPP,NUSE,NUST,NORD,NCN,NCN2,NCN4,NTP,
     2  IMX1,NMX,JR2,JMAX,MI(10),MF(10)
      REAL*8  ASR,BSR,CSR,ALR,BLR,CLR,DCSR,ADCSR,PDCSR,VRAT,
     1  DX1,DX2,DX3,EX1,EX2,EX3,CNN,VLIM,X1,X2,X3,Y1,Y2,Y3,
     1  XX(NPP),YY(NPP),XI(NTP),YI(NTP),XJ(20),YJ(20),DUMM(20)
c
      SAVE ASR,BSR,CSR,ISR,ALR,BLR,CLR,IMX1,NMX,JR2,JMAX
c
      NUST= NUSE/2
      IF(NUSE.LE.0) NUST= 2
      IDER= 0
c** Determine if/where need to begin extrapolation beyond input data
c  XX(MI(J))  is the 1-st mesh point past turning point  XI(J) .
c  XX(MF(J))  is the last mesh point before turning point  XI(NTP+1-J)
      DO 6 J = 1,NUST
          MI(J)= 1
          MF(J)= 0
          DO  I= 1,NPP
              IF(XX(I).LE.XI(J)) MI(J)= I+ 1
              IF(XX(I).GE.XI(NTP+1-J)) GO TO 6
              MF(J)= I
              ENDDO
    6     CONTINUE
      IF(NUST.LT.2) THEN
          MFIN= MI(1)-1
        ELSE
          MFIN= MI(2)-1
        ENDIF
      IF(LNPT.GT.0) THEN
c-----------------------------------------------------------------------
c** For a new case determine analytic functions for extrapolating beyond
c  the range of the input points (if necessary) on this or later calls.
c** Try to fit three innermost turning points to  V(R)=A+B*DEXP(-C*R).
c** If unsatisfactory, extrapolate inward with inverse power function
          IF(IR2.LE.0) THEN
              DO  I= 1,4
                  YJ(I)= YI(I)
                  ENDDO
            ELSE
              DO  I= 1,4
                  YJ(I)= YI(I)/XI(I)**2
                  ENDDO
            ENDIF
          X1= XI(1)
          X2= XI(2)
          X3= XI(3)
          Y1= YJ(1)
          Y2= YJ(2)
          Y3= YJ(3)
          IF((Y1-Y2)*(Y2-Y3).LE.0.d0) THEN
c** If 3 innermost points not monotonic, use A+B/X inward extrapoln.
              ISR= 0
              WRITE(6,600)
            ELSE
c** Use cubic through innermost points to get initial trial exponent
c  from ratio of derivatives,  Y''/Y'
              IDER= 2
              ISR= 4
              CALL PLYINTRP(XI,YJ,ISR,X2,XJ,ISR,IDER)
              CSR= XJ(3)/XJ(2)
              DCSR= DABS(CSR*X2)
              IF(DCSR.GT.1.5D+2) THEN
c** If exponential causes overflows, use inverse power inward extrapoln.
                  ISR= 0
                  WRITE(6,602) CSR
                  GO TO 20
                  ENDIF
c** Prepare parameters for inward exponential extrapolation
              VRAT= (Y3- Y2)/(Y1- Y2)
              DX1= X1- X2
              DX3= X3- X2
              EX2= 1.D0
              ADCSR= 1.d99
c** Now iterate (with actual point) to get exact exponent coefficient 
              DO  J= 1,15
                  PDCSR= ADCSR
                  EX1= DEXP( CSR*DX1)
                  EX3= DEXP( CSR*DX3)
                  DCSR= (VRAT- (EX3- EX2)/(EX1- EX2)) /
     1   ((X3*EX3- X2 - (X1*EX1- X2)*(EX3-EX2)/(EX1- EX2))/(EX1- EX2))
                  ADCSR= ABS(DCSR)
                  IF((ADCSR.GT.PDCSR).AND.(ADCSR.LT.1.d-8)) GO TO 12
                  IF(ADCSR.LT.1.d-12) GO TO 12
                  CSR= CSR+ DCSR 
                  ENDDO
              WRITE(6,604) DCSR
   12         BSR= (Y1-Y2)/(EX1-EX2)
              ASR= Y2-BSR*EX2
              BSR= BSR*DEXP(-CSR*X2)
              WRITE(6,606) X2,ASR,BSR,CSR
            ENDIF
   20     IF(ISR.LE.0) THEN
              IF((X1*X2).LE.0.d0) THEN
c** If 1'st two mesh points of opposite sign, extrapolate linearly
                  ISR= -1
                  ASR= Y2
                  BSR= (Y2- Y1)/(X2- X1)
                  CSR= X2
                  WRITE(6,608) X2,ASR,BSR,CSR
                ELSE
c** For inward extrapolation as inverse power through 1'st two points ..
                  BSR= (Y1-Y2)* X1*X2/(X2- X1)
                  ASR= Y1-BSR/X1
                  CSR= X2
                  WRITE(6,610) X2,ASR,BSR
                ENDIF
              ENDIF
          ENDIF
  600 FORMAT('  ** CAUTION ** Exponential inward extrapolation fails'/
     1 16x,'since first 3 points not monotonic, ... so ...')
  602 FORMAT(' *** CAUTION ** inward extrapolation exponent coefficient
     1   C=',D12.4/10x,'could cause overflows, ... so ...')
  604 FORMAT(' *** CAUTION ** after 15 tries inward extrap. exponent coe
     1fft change is',1PD9.1)
  606 FORMAT(' Extrapolate to   X .le.',F7.4,'  with'/'   Y=',F13.3,
     1  SP,1PD15.6,' * exp(',SS,D13.6,'*X)')
  608 FORMAT(' Extrapolate to   X .le.',F8.4,'   with'/'   Y=',F13.3,
     1  SP,1PD16.7,' * [X - (',SS,F8.4,')]')
  610 FORMAT(' Extrapolate to  X .le.',F8.4,'   with   Y=',F12.3,
     1  SP,1PD15.6,')/X**1')
c
      IF(MFIN.GT.0) THEN
c** If needed, calculate function in inner extrapolation region
          IF(ISR.GT.0) THEN
c ... either as an exponential
              DO  I= 1,MFIN
                  EX1= CSR*XX(I)
                  IF(DABS(EX1).GT.1.D+2) EX1= 1.D+2*DSIGN(1.d0,EX1)
                  YY(I)= ASR+BSR*DEXP(EX1)
                  ENDDO
            ELSEIF(ISR.EQ.0) THEN
c ... or if that fails, as an inverse power
              DO  I= 1,MFIN
                  YY(I)= ASR+BSR/XX(I)
                  ENDDO
            ELSEIF(ISR.LT.0) THEN
c ... or if X changes sign, extrapolate inward linearly
              DO  I= 1,MFIN
                  YY(I)= ASR+ BSR*(XX(I)- CSR)
                  ENDDO
            ENDIF
          ENDIF
c** End of inward extrapolation procedure
c-----------------------------------------------------------------------
      MINNER= MFIN
      IF(NUST.GT.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near inner end of range
          DO  J= 3,NUST
              NORD= 2*(J-1)
              MBEG= MI(J-1)
              MFIN= MI(J)-1
              IF(MFIN.GE.MBEG) THEN
                  DO  I=  MBEG,MFIN
                      CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                      YY(I)= DUMM(1)
                      ENDDO
                  ENDIF
              ENDDO
          ENDIF
c** Main interpolation step begins here
c=======================================================================
      MBEG= MI(NUST)
      MFIN= MF(NUST)
      IF(MFIN.GE.MBEG) THEN
          IF(NUSE.LE.0) THEN
c** Either ... use cubic spline for main interpolation step
              CALL SPLINT(LNPT,NTP,XI,YI,MBEG,MFIN,XX,YY)
            ELSE
c ... or use piecewise polynomials for main interpolation step
              DO  I= MBEG,MFIN
                  CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NUSE,IDER)
                  YY(I)= DUMM(1)
                  ENDDO
            ENDIF
          ENDIF
      IF(MFIN.LT.NPP) THEN
          IF(NUST.LE.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near outer end of range
              MBEG= MF(NUST)+1
            ELSE
              NN= NUST-2
              DO  J= 1,NN
                  NORD= 2*(NUST-J)
                  MBEG= MF(NUST-J+1)+1
                  MFIN= MF(NUST-J)
                  IF(MFIN.GE.MBEG) THEN
                      DO  I= MBEG,MFIN
                          CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                          YY(I)= DUMM(1)
                          ENDDO
                      END IF
                  ENDDO
            ENDIF
          ENDIF
      MBEG= MFIN+1
      IF((MFIN.GT.MINNER).AND.(IR2.GT.0)) THEN
c** In (IR2.gt.0) option, now remove X**2 from the interpolated function
          DO  I= MINNER+1,MFIN
              YY(I)= YY(I)/XX(I)**2
              ENDDO
          ENDIF
c** Print test of smoothness at join with analytic inward extrapolation
c     IF(LNPT.GT.0) THEN
c         MST= MAX0(MINNER-4,1)
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         IF(MFN.GT.MFIN) MFN= MFIN
c         IF(MINNER.GT.0) WRITE(6,611) X2,((XX(I),YY(I),I= J,MFN,3),
c    1        J= MST,MST+2)
c 611 FORMAT('     Verify smoothness of inner join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
c         ENDIF
c-----------------------------------------------------------------------
c** To extrapolate potential beyond range of given turning points ...
      IF(LNPT.GT.0) THEN
c** On first entry, calculate things needed for extrapolation constants
          Y1= YI(NTP)
          Y2= YI(NTP-1)
          Y3= YI(NTP-2)
          X1= XI(NTP)
          X2= XI(NTP-1)
          X3= XI(NTP-2)
          IF(IR2.GT.0) THEN
              Y1= Y1/X1**2
              Y2= Y2/X2**2
              Y3= Y3/X3**2
              ENDIF
          ENDIF
c** Check inverse-power tail power ...
      IF(NCN.LE.0) NCN= 6
      IF(ILR.LT.0) THEN
          IF(LNPT.GT.0) THEN
C** For  ILR.lt.0  use  Y = VLIM - ALR * exp[-CLR*(X - BLR)**2]
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              BLR= (X1+X2 - (X2+X3)*EX1/EX2)/(2.d0- 2.d0*EX1/EX2)
              CLR= -EX1/(X1+X2-2.d0*BLR)
              ALR= (VLIM-Y1)*DEXP(CLR*(X1-BLR)**2)
              WRITE(6,614) X2,VLIM,ALR,CLR,BLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*DEXP(-CLR*(XX(I) - BLR)**2)
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.0) THEN
c** For ILR.le.0  use  Y = VLIM - ALR * X**p * exp(-CLR*X)
          IF(LNPT.GT.0) THEN
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              DX1= DLOG(X1/X2)/(X1-X2)
              DX2= DLOG(X2/X3)/(X2-X3)
              BLR= (EX1-EX2)/(DX1-DX2)
              CLR= BLR*DX1- EX1
              ALR= (VLIM-Y1)* DEXP(CLR*X1)/X1**BLR 
              WRITE(6,616) X2,VLIM,ALR,BLR,CLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*XX(I)**BLR *DEXP(-CLR*XX(I))
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
   50 IF(ILR.EQ.1) THEN
c** For  ILR=1 ,  use     Y = VLIM + ALR/X**BLR
          IF(LNPT.GT.0) THEN
              BLR= DLOG((VLIM-Y2)/(VLIM-Y1))/DLOG(X1/X2)
              ALR= (Y1- VLIM)*X1**BLR
              NCN= NINT(BLR)
              WRITE(6,618) X2,VLIM,ALR,BLR,NCN
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM+ ALR/XX(I)**BLR
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
c** Set constants for long-range extrapolation
      IFXCN= 0
      IF((CNN.GT.0.d0).OR.(CNN.LT.0.d0)) IFXCN= 1
      NCN2= NCN+2
      IF(ILR.EQ.2) THEN
c** For ILR=2 ,  use   Y = VLIM - CNN/X**NCN - BSR/X**(NCN+2)
c*  If CNN held fixed need ILR > 2  to prevent discontinuity
          IF(LNPT.GT.0) THEN
              IF(IFXCN.LE.0) THEN
                  CNN= ((VLIM-Y1)*X1**NCN2 -
     1                 (VLIM-Y2)*X2**NCN2)/(X1**2-X2**2)
                  ENDIF
              ALR= CNN
              BLR= (VLIM-Y1)*X1**NCN2 - CNN*X1**2
              WRITE(6,620) X2,VLIM,CNN,NCN,BLR,NCN2
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM-(ALR+BLR/XX(I)**2)/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.3) THEN
c** For ILR=3 , use   Y = VLIM - (CN + CN2/X**2 + CN4/X**4)/X**NCN
          IF(LNPT.GT.0) THEN
              NCN4= NCN+4
              IF(IFXCN.GT.0) THEN
                  ALR= CNN
                  BLR= (((VLIM-Y1)*X1**NCN-ALR)*X1**4-((VLIM-Y2)
     1                     *X2**NCN-ALR)*X2**4)/(X1**2-X2**2)
                  CLR= ((VLIM-Y1)*X1**NCN-ALR-BLR/X1**2)*X1**4
                ELSE
                  EX1= X1**2
                  EX2= X2**2
                  EX3= X3**2
                  DX1= (VLIM-Y1)*X1**NCN4
                  DX2= (VLIM-Y2)*X2**NCN4
                  DX3= (VLIM-Y3)*X3**NCN4
                  BLR= (DX1-DX2)/(EX1-EX2)
                  ALR= (BLR-(DX2-DX3)/(EX2-EX3))/(EX1-EX3)
                  BLR= BLR-ALR*(EX1+EX2)
                  CLR= DX1-(ALR*EX1+BLR)*EX1
                ENDIF
              WRITE(6,622) X2,VLIM,ALR,NCN,BLR,NCN2,CLR,NCN4
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  EX2= 1.d0/XX(I)**2
                  YY(I)= VLIM-(ALR+EX2*(BLR+EX2*CLR))/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.GE.4) THEN
c** For ILR.ge.4,   Y = VLIM-SUM(BB(K)/X**K) , (K=NCN,NMX=NCN+ILR-1)
          IF(LNPT.GT.0) THEN
              IF(NCN.LE.0) NCN= 1
              IMX1= ILR-1
              NMX= NCN+IMX1
              JR2= 0
              IF(IR2.GT.0) JR2= 2
              IDER= 0
              JMAX= ILR
              IF(IFXCN.GT.0) JMAX= IMX1
              WRITE(6,624) X2,ILR,NCN,VLIM
              IF(IFXCN.GT.0) WRITE(6,626) NCN,CNN
              ENDIF
c** Actually extrapolate with polynomial fitted to the last JMAX
c  values of  (VLIM - YI(I))*XI(I)**NMX  , & then convert back to  YY(I).
          IF(MBEG.LE.NPP) THEN
              J= NTP- JMAX
              DO  I= 1,JMAX
                  J= J+1
                  XJ(I)= XI(J)
                  YJ(I)= (VLIM-YI(J)/XI(J)**JR2)*XI(J)**NMX
                  IF(IFXCN.GT.0) YJ(I)= YJ(I)- CNN*XI(J)**IMX1
                  ENDDO
              DO  I= MBEG,NPP
                  CALL PLYINTRP(XJ,YJ,JMAX,XX(I),DUMM,JMAX,IDER)
                  YY(I)= DUMM(1)
                  IF(IFXCN.GT.0) YY(I)= YY(I)+ CNN*XX(I)**IMX1
                  YY(I)= VLIM-YY(I)/XX(I)**NMX
                  ENDDO
              ENDIF
          ENDIF
c** Finished extrapolation section.
   90 CONTINUE
c** Test smoothness at outer join to analytic extrapolation function
c     IF((LNPT.GT.0).AND.(MBEG.LE.NPP)) THEN
c         MST= MBEG-5
c         IF(MST.LT.1) MST= 1
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         WRITE(6,627) X2,((XX(I),YY(I),I= J,MFN,3),J= MST,MST+2)
c         ENDIF
c 627 FORMAT('     Verify smoothness of outer join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
      RETURN
  612 FORMAT('  *** BUT *** since exponent has positive coefficient, swi
     1tch form ...')
  614 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1  F12.4,' - (',1PD13.6,') * exp{-',0PF10.6,' * (r -',F9.6,')**2}')
  616 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1 F12.4,' - (',1PD13.6,') * r**',0PF10.6,'  * exp{-(',F11.6,'*r)}')
  618 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,SP,1PD15.6,'/X**(',SS,D13.6,')] ,  yielding   NCN=',I3)
  620 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',SS,I1,']')
  622 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/
     1  '   Y=',F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',
     2  SS,I1,SP,D14.6,'/X**',SS,I2,']')
  624 FORMAT(' Function for  X .GE.',F7.3,'  generated by',I3,
     1 '-point inverse-power interpolation'/'   with leading term  1/r**
     2',I1,'  relative to dissociation limit   YLIM=',F11.3)
  626 FORMAT('   and (dimensionless) leading coefficient fixed as   C',
     1  I1,'=',G15.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PLYINTRP(XI,YI,NPT,RR,C,NCFT,IDER)
c* From the NPT known mesh points (XI,YI) ,given in order of increasing
c  or decreasing XI(I), select the NCFT points (XJ,YJ) surrounding the 
c  given point RR, and by fitting an (NCFT-1)-th degree polynomial through
c  them, interpolate to find the function CC(1) and its first IDER 
c  derivatives (CC(I+1),I=1,IDER) evaluated at RR.
c* Adapted by  R.J. Le Roy  from algorithm #416,Comm.A.C.M.;  27/02/1988
c=======================================================================
      INTEGER  I,J,K,I1,I2,IFC,IM,IDER,J1,NH,NPT,NCFT
      REAL*8  RR,XX,XI(NPT),YI(NPT),C(NCFT),XJ(20),YJ(20)
c
      IF((NCFT.GT.20).OR.(NCFT.GT.NPT)) GO TO 101
      NH= NCFT/2
c** First locate the known mesh points (XJ,YJ) bracketing RR
      I1= 1
      I2= NCFT
      IF(NCFT.NE.NPT) THEN
          IF(XI(NPT).LE.XI(1)) THEN
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).LT.RR) GO TO 20
                  ENDDO
            ELSE
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).GT.RR) GO TO 20
                  ENDDO
            ENDIF
   20     I1= IM-NH
          IF(I1.LE.0) I1= 1
          I2= I1+NCFT-1
          IF(I2.GT.NPT) THEN
              I2= NPT
              I1= I2-NCFT+1
              ENDIF
          ENDIF
      J= 0
      DO  I= I1,I2
          J= J+1
          XJ(J)= XI(I)-RR
          YJ(J)= YI(I)
          ENDDO
c** Now determine polynomial coefficients C(I).
      DO  I= 2,NCFT
          I1= I-1
          K= I1+1
          DO  J= 1,I1
              K= K-1
              YJ(K)= (YJ(K+1)-YJ(K))/(XJ(I)-XJ(K))
              ENDDO
          ENDDO
      C(1)= YJ(1)
      DO  I= 2,NCFT
          XX= XJ(I)
          C(I)= C(I-1)
          IF(I.NE.2) THEN
              I1= I-1
              K= I1+1
              DO  J= 2,I1
                  K= K-1
                  C(K)= -XX*C(K)+C(K-1)
                  ENDDO
              ENDIF
          C(1)= YJ(I)-XX*C(1)
          ENDDO
c** Finally, convert polynomial coefficients to derivatives at RR.
      IFC= 1
      IF(IDER.GE.NCFT) IDER= NCFT-1
      IF(IDER.LE.1) GO TO 99
      DO  I= 2,IDER
          J= I+1
          IFC= IFC*I
          C(J)= C(J)*IFC
          ENDDO
      IF(J.LT.NCFT) THEN
          J1= J+1
          DO  I= J1,NCFT
              C(I)= 0.D+0
              ENDDO
          ENDIF
   99 RETURN
  101 WRITE(6,601) NCFT,NCFT,NPT
      STOP
  601 FORMAT(/' *** Dimensioning ERROR in PLYINTRP :  either   (NCFT=',
     1  I2,' .GT. 20)   or   (NCFT=',I2,' .GT. NPT=',I3,')')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE SPLINT(LNPT,NTP,R1,V1,MBEG,MEND,XX,YY)
c** Subroutine to generate (if LNPT.ge.0) 4*NTP coefficients CSP(J)
c  of a cubic spline passing through the NTP points (R1(J),V1(J))
c  and to then calculate values of the resulting function YY(I) at the
c  entering abscissae values XX(I) for  I=MBEG to MEND.
c** If LNPT < 0 , generate function values at the given XX(I) using
c  the coefficients CSP(J) obtained and SAVEd on a preceding call.
c** Assumes both R1(J) & XX(I) are monotonic increasing.
c+++++ Calls only subroutines SPLINE and PLYINTRP ++++++++++++++++++++++
c=======================================================================
      INTEGER MAXSP
      PARAMETER (MAXSP=6400)
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
          CALL PLYINTRP(R1(I1ST),V1(I1ST),NIPT,R1(NTP),CSP,NIPT,IDER)
          TTMP= CSP(2)
          CALL PLYINTRP(R1,V1,NIPT,R1(1),CSP,NIPT,IDER)
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
  250 IER= NER
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

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
     1 '/r**',i1,'  with  C_',i1,'=',1PD14.7/7x,'defined by   A_sw=',
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

c***********************************************************************
      SUBROUTINE dampF(R,RHOd,m,Dm,Dmp,Dmpp,IDF)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c  'Dmpp' derivatives w.r.t. R of the chosen potential damping function
c   IDF = 1 :  Tang-Toennies damping function [JCP 80,3726(1984)]
c   IDF = 2 :  damped-dispersion damping function [MolP 52,763(1984)]
c   IDF = 3 : extend IDF=2 by multiplying by Scoles' global damping fx.
c  all assuming  R  in units Angst. 
c***********************************************************************
      INTEGER i,m,IDF
      REAL*8 R,RHOd,RHOdR,Dm,Dmp,Dmpp,fd,fdp,fdpp,Coef,XPOW
      REAL*8 sum(20)
c
      RHOdR= RHOd*R
      IF(IDF.EQ.1) THEN
c** Tang-Toennies damping function: numerical factor 3.16 makes these
c  functions equivalent to Scoles functions, and allows the  RHOd  
c  scaling factor to have the same significance in both cases.
          RHOdR= 3.16d0*RHOdR
          XPOW = DEXP(-RHOdR)
          DO  i= 1,m
              sum(i) = 1.0d0
              END DO
          fd = 1.0d0
          DO  i= 1,m
              fd = fd*RHOdR/i
              sum(i) = sum(i) + fd
              ENDDO
          Dm = 1.0d0-XPOW*sum(m)
          Dmp = 3.16d0*RHOd*XPOW*(sum(m) - sum(m-1))
          Dmpp = 3.16d0*RHOd**2*XPOW*(-sum(m) + 2*sum(m-1) - sum(m-2))
          ENDIF
c
      IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
c** Douketis et al. [MolP 52, 763 (1984)]  damping function.
          XPOW= DEXP(-3.97d0*RHOdR/m-0.39d0*(RHOdR)**2/DSQRT(m*1.0d0))
ccc       XPOW= DEXP(-3.968424881d0*RHOdR/m
ccc  1                      - 0.389246070d0*(RHOdR)**2/DSQRT(m*1.0d0))
          Dm = (1.0d0 - XPOW)**m
          Coef = 3.97d0*RHOd/m + 0.78d0*RHOd*RHOdR/DSQRT(m*1.0d0)
ccc       Coef = 3.968424881d0*RHOd/m 
ccc  1                      + 0.7784921396d0*RHOd*RHOdR/DSQRT(m*1.0d0)
          Dmp = m*Coef*XPOW*(1.0d0 - XPOW)**(m-1)
          Dmpp = XPOW*(1.0d0 - XPOW)**(m-1)
     1                                 *(0.78d0*SQRT(m*1.0d0)*RHOd**2
ccc  1                          *(0.7784921396d0*SQRT(m*1.0d0)*RHOd**2 
     2              - m*Coef**2*(1.0d0-m*XPOW)/(1.0d0 - XPOW))
          ENDIF
c
      IF(IDF.EQ.3) THEN
c** If desired, also apply Scoles' global scaling/damping factor
          XPOW= DEXP(-1.474d0*RHOdR)
          Coef= (RHOdR)**1.68d0
          fd= Coef*XPOW
          fdp= RHOd*fd*(1.474d0 - 1.68/(RHOdR))
          fdpp= RHOd**2 *fd*(-2.172676d0 + 
     1                       (4.95264d0 - 1.1424d0/(RHOdR))/(RHOdR))
          fd= 1.d0 - fd
          Dmpp= Dmpp*fd + 2.d0*Dmp*fdp + Dm*fdpp
          Dmp= Dmp*fd + Dm*fdp
          Dm= Dm*fd
          END IF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

