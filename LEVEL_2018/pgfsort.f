      Program SORT
      implicit real*8 (a-h,o-z)
      dimension iv(500),ij(500),evj(500),iwv(50),iwj(50),ewvj(50),
     1  iilev(200),ivlev(200),ijlev(200),elev(200)
c** Program to prepare data set for LEVEL by selecting all possible 
c  transitions whose frequencies lie in the range  EMIN to EMAX (cm-1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++ First call a SUN-Unix specific subroutine call, required to 
c  customize input/output so that "standard" (e.g., VAX) read(n, ) 
c  write(m, ) will work OK.
      call ioinit(.true.,.false.,.false.,'FORT',.false.)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	read(5,*) EMIN,EMAX
  	write(6,605) EMIN,EMAX
  605   format('0  Data set for transitions with frequencies in the ',
     1 'range    EMIN=',f8.2,'   to   EMAX =',f9.2,'(cm-1)'/)
	do 10 i=1,1000
	read(5,*,end=12) iv(i),ij(i),evj(i)
   10 nvj= i
   12   nlev= 0
	ntot= 0
	nvjused= 0
	do 90 ilev=2,nvj
c** Preliminary sort to find which initial state levels can emit in the
c  frequency range EMIN to EMAX subject to the chosen selection rules.
	k= 0
	do 20 j=1,ilev-1
c** Apply the P and R rotational selection rule.
	if(iabs(ij(ilev)-ij(j)).ne.1) go to 20
	diff= evj(ilev)-evj(j)
c** Here, count final state levels yielding transition frequencies in
c  the range  EMIN to EMAX (cm-1)
	if((diff.gt.EMAX).or.(diff.lt.EMIN)) go to 20
	k= k+1
   20 continue
	if(k.le.0) go to 26
	nlev= nlev+1
        nvjused= nvjused+1
	iilev(nlev)= ilev
	ivlev(nlev)= iv(ilev)
	ijlev(nlev)= ij(ilev)
	elev(nlev)= evj(ilev)
   26	if((ilev.eq.nvj).or.(nlev.eq.200)) then
	write(6,602) nlev,(ivlev(j),ijlev(j),elev(j),j=1,nlev)
  602   format('0',i5,' 4  0 0  0 -1 0  0.d0'/(4(2I4,f10.2)))
	do 50 i=1,nlev
	ii= iilev(i)
	k= 0
	do 40 j=1,ii-1
c** Apply the P and R rotational selection rule.
	if(iabs(ij(ii)-ij(j)).ne.1) go to 40
	diff= evj(ii)-evj(j)
c** Now, select final state levels yielding transition frequencies in the
c  range  EMIN to EMAX (cm-1)
	if((diff.gt.EMAX).or.(diff.lt.EMIN)) go to 40
	k= k+1
	iwv(k)= iv(j)
	iwj(k)= ij(j)
	ewvj(k)= evj(j)
   40 continue
	if(k.le.0) go to 50
	ntot= ntot+k
	write(6,601) k,(iwv(n),iwj(n),ewvj(n),n=1,k)
  601 format( i6/(1x,4(2i4,f10.2)))
  50	continue
	nlev=0
	write(6,604)
  604  	format(////)
	end if
   90 continue
	write(6,603) nvj,nvjused,ntot
  603	format('0 Total number of (v,J) levels considered is    nvj =',
     1 i4/'  Number of levels contributing emission within the specified
     2 frequency range is   nvj(used) =',i4/'  Total number of matrix el
     3ements in the specified frequency range is   ntot =',i6)
	stop
	end
