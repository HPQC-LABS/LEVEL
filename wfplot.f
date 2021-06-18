c** Program to read LEVEL-"punched" potential and wave function arrays
c  and combine them onto Potential/Wavefunction plot
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION RR(1000),YY(1000)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ This is a SUN-Unix specific subroutine call, required to customize
c  input/output so that "standard" (e.g., VAX) read(n, ) write(m, )
c  will work OK.
      call ioinit(.true.,.false.,.false.,'fort',.false.)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C*** Read potential array, additive energy shift EAD and multiplicative
c  scaling factor ESC
c-----------------------------------------------------------------------
    2	READ(7,*,end=99) NPPT,EAD,ESC,ESH
        if(nppt.le.0) GO TO 15
	READ(7,*) (RR(I),YY(I),I=1,NPPT)
c-----------------------------------------------------------------------
	DO 10 I=1,NPPT
   10	YY(I)= (YY(I)+ EAD)*ESC+ ESH
	WRITE(8,801) NPPT,EAD,ESC,(RR(I),YY(I),I=1,NPPT)
  801 FORMAT('0',I6,'-point Potential shifted by'G15.7,' and scaled by'
     1  ,D11.3/(f6.3,f9.4,4(f7.3,f9.4)))
c-----------------------------------------------------------------------
   15	READ(7,*,END=99) NPW,EW
	IF(NPW.LE.0) GO TO 2
	READ(7,*) (RR(I),YY(I),I=1,NPW)
c-----------------------------------------------------------------------
	EFX= (EW+EAD)*ESC+ ESH
	DO 20 I=1,NPW
   20 	YY(I)= YY(I)+EFX
	WRITE(8,802) NPW,EFX,(RR(I),YY(I),I=1,NPW)
  802 FORMAT('ZONE'/(f6.3,f9.5,4(f7.3,f9.5)))
	GO TO 15
   99	STOP
	END
