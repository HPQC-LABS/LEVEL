c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c** This 'Block' Data Utility routine  that governs array dimensioning
c  in program  LEVEL16 must reside with the name 'arrsizes.h' in the 
c  same directory containing the FORTRAN file(s) for this Program when 
c  it is being compiled, **OR** be incorporated into the program 
c  wherever the statement 'INCLUDE arrsizes.h' appears !!
c-----------------------------------------------------------------------
      INTEGER NDIMR, NVIBMX, NTPMX, MAXSP, MORDRMX, RORDR, NbetaMX,
     1                                            LMAX, NBOBmx, NCMMAX
c** NDIMR  is maximum size of PEC, wavefx, and various radial arrary
      PARAMETER (NDIMR= 250001)
c** NVIBMX  is the maximum no. vibrational levels, or rotational sublevel
c       for a given 'v' whose energies may be generated and stored
      PARAMETER (NVIBMX= 400)
c** NTPMX  is maximum no. of PEC or TMF points that may be read-in and 
c   interplated over; MAXSP = no. cubic spline cfts for these NTPMX pts.
      PARAMETER (NTPMX= 2000, MAXSP=4*NTPMX)
c** RORDR is maximum order of rot. constants generated for each vib level
      PARAMETER (RORDR = 7)
c** MORDRMX is maximum polynomial order for TMF or martix element argument
      PARAMETER (MORDRMX = 20)
c** NbetaMX  is the largest no. PEC exponent polynomial parameter
      PARAMETER (NbetaMX  = 50, LMAX= NbetaMX)
c** NBOBmx  is the largest no. of BOB expansion parameters
      PARAMETER (NBOBmx  = 20)
c** NCMMax  is max. no. long-range inverse-power PEC coeffts. allowed
      PARAMETER (NCMMax= 20)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
