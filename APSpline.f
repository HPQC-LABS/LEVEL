c***********************************************************************
c** Asen Pashov's subroutines for constructing spline functions and
c   their derivatives.
      double precision function Scalc(x,m,n,XGRID,rKL,LMAX)
c** At the position 'x', Scalc is returned as the value of the m'th 
c  of the 'n' Sm(x) function defining a natural cubic spline through the
c  mesh points located at  x= XGRID(x_i), for i=1,n.  LMAX specifies the 
c  maximum number of mesh points x= XGRID(x_i) allowed by the calling program
c---------------------------------------------------------------------
      INTEGER  LMAX,I,K,KK,M,N
      REAL*8  x,y1,y2,XGRID(LMAX),rKL(LMAX,LMAX)
      k= 0
      kk= 0
      do i=2,n
c... select interval
          if ((x.gt.XGRID(i-1)).and.(x.le.XGRID(i)))  k=i
          end do
      if (x.lt.XGRID(1)) then
          k=2
          kk=1
          end if
      if (x.gt.XGRID(n)) then
          k=n
          kk=1
          end if
      if(x.eq.XGRID(1)) k=2
      y1=XGRID(k-1)
      y2=XGRID(k)
      Scalc= 0.d0
      IF(kk.eq.0) 
     1    Scalc= rKL(m,k)*((y1-x)*(((y1-x)/(y1-y2))**2-1)/6)*(y1-y2)
     2         + rKL(m,k-1)*((x-y2)*(((x-y2)/(y1-y2))**2-1)/6)*(y1-y2)
      IF(k.EQ.m) Scalc= Scalc + (y1-x)/(y1-y2)
      IF(k-1.EQ.m) Scalc= Scalc + (x-y2)/(y1-y2)
c... Asen's original coding ...
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)+
cc   +   C(x,y1,y2)*rKL(m,k)+D(x,y1,y2)*rKL(m,k-1)
cc       else
cc       Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)
cc     A=(x1-z)/(x1-x2)
cc     B=(z-x2)/(x1-x2)
cc     C=((x1-z)*(((x1-z)/(x1-x2))**2-1)/6)*(x1-x2)
cc     D=((z-x2)*(((z-x2)/(x1-x2))**2-1)/6)*(x1-x2)
c... Asen's original coding ...
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Sprime(x,m,n,XGRID,rKL,LMAX)
c** At the position 'x', evaluate the derivative w.r.t. x of the m'th 
c  Sm(x) function contributing the definition of the the natural cubic
c  spline defined by function values at the  n  points  XGRID(i) [i=1,n]
      INTEGER i,k,kk,m,n,LMAX
      REAL*8 x,del,y1,y2,XGRID(LMAX),rKL(LMAX,LMAX)
      k=0
      kk=0
      do i=2,n
          if((x.gt.XGRID(i-1)).and.(x.le.XGRID(i)))  k=i
          enddo
      if(x.lt.XGRID(1)) then
          k=2
          kk=1
          end if
      if (x.gt.XGRID(n)) then
          k=n
          kk=1
          end if
      if (x.eq.XGRID(1)) k=2
      y1=XGRID(k-1)
      y2=XGRID(k)
      del=y1-y2
      Sprime= 0.d0
      if(kk.eq.0) Sprime= (del-3.d0*(y1-x)**2/del)*rKL(m,k)/6.d0 +
     1                        (3.d0*(x-y2)**2/del-del)*rKL(m,k-1)/6.d0
      IF(k-1.eq.m) Sprime= Sprime + 1.d0/del 
      IF(k.eq.m) Sprime= Sprime - 1.d0/del 
ccc     if(kk.eq.0) then
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del+
ccc  +                    (del-3*(y1-x)**2/del)*rKL(m,k)/6+
ccc  +                    (3*(x-y2)**2/del-del)*rKL(m,k-1)/6
ccc       else
ccc         Sprim=ndirac(k-1,m)/del-ndirac(k,m)/del
ccc       end if
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef(NGRID,XGRID,rKL)   
c** Call this subroutine with list of the 'NGRID' spline x_i values in 
c   array 'XGRID' with maximum dimension 'LMAX', and it will return the 
c   LMAX x LMAX  array of 'rKL' coefficients used for generating the 
c   'NGRID' S_{NGRID}(x) spline coefficient functions
c----------------- Based on nespl subroutine ---------------------------
c** CAUTION .. must dimension internal arrays B, INDX & vv @ compilation
      INCLUDE 'arrsizes.h'                !! needed only to define  LMAX
c***--------------------------------------------------------------------
      INTEGER I,J,NGRID,INDX(1:LMAX)
      REAL*8 XGRID(LMAX),rKL(LMAX,LMAX),B(LMAX,LMAX),vv(LMAX), d
c ...  note vv dimensioned here, but only used in   ludcmp !!
      DO  i= 1,LMAX
          DO  j= 1,LMAX
              rKL(i,j)= 0.d0
              B(i,j)= 0.d0
              ENDDO
          ENDDO
      rKL(1,1)= (XGRID(3)-XGRID(1))/3.d0
      rKL(1,2)= (XGRID(3)-XGRID(2))/6.d0
      do i= 2,NGRID-3
          rKL(i,i-1)= (XGRID(i+1)-XGRID(i))/6.d0
          rKL(i,i)= (XGRID(i+2)-XGRID(i))/3.d0
          rKL(i,i+1)= (XGRID(i+2)-XGRID(i+1))/6.d0
          end do
      rKL(NGRID-2,NGRID-3)= (XGRID(NGRID-1)-XGRID(NGRID-2))/6.d0
      rKL(NGRID-2,NGRID-2)= (XGRID(NGRID)-XGRID(NGRID-2))/3.d0  
      do i= 1,NGRID-2
          B(i,i)= 1.d0/(XGRID(i+1)-XGRID(i))
          B(i,i+1)= -1.d0/(XGRID(i+2)-XGRID(i+1))-1.d0/
     1                                           (XGRID(i+1)-XGRID(i))
          B(i,i+2)= 1.d0/(XGRID(i+2)-XGRID(i+1))
          end do  
      call ludcmp(rKL,NGRID-2,LMAX,indx,vv,d)
      do i= 1,NGRID 
          call lubksb(rKL,NGRID-2,LMAX,indx,B(1,i))
          end do 
      do i= 1,NGRID-2
          do j= 1,NGRID
              rKL(j,i+1)= B(i,j)
              end do
          end do 
      do i= 1,NGRID
          rKL(i,1)= 0.0d0
          rKL(i,NGRID)= 0.0d0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(rKL,NGRID,LMAX,indx,vv,d)
      INTEGER NGRID,LMAX,indx(LMAX),NMAX,i,imax,j,k
      double precision d,rKL(LMAX,LMAX),vv(LMAX),TINY,aamax,dum,sum
      PARAMETER (TINY= 1.0e-20)
      d= 1.d0
      do  i= 1,NGRID
          aamax= 0.d0
          do  j= 1,NGRID
              if (abs(rKL(i,j)).gt.aamax) aamax= abs(rKL(i,j))
              enddo
          if (aamax.eq.0.) WRITE(6,*) 'singular matrix in ludcmp'
          vv(i)= 1.d0/aamax
          enddo
      do  j= 1,NGRID
          do  i= 1,j-1
              sum= rKL(i,j)
              do  k= 1,i-1
                  sum= sum-rKL(i,k)*rKL(k,j)
                  enddo
              rKL(i,j)= sum
              enddo
          aamax= 0.d0
          do  i= j,NGRID
              sum= rKL(i,j)
              do  k= 1,j-1
                  sum= sum-rKL(i,k)*rKL(k,j)
                  enddo
              rKL(i,j)= sum
              dum= vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax= i
                  aamax= dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k= 1,NGRID
                  dum= rKL(imax,k)
                  rKL(imax,k)= rKL(j,k)
                  rKL(j,k)= dum
                  enddo
              d= -d
              vv(imax)= vv(j)
              endif
          indx(j)= imax
          if(rKL(j,j).eq.0.)rKL(j,j)= TINY
              if(j.ne.NGRID)then
                  dum= 1.d0/rKL(j,j)
                  do  i= j+1,NGRID
                      rKL(i,j)= rKL(i,j)*dum
                      enddo
                  endif
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE lubksb(rKL,NGRID,LMAX,indx,b)
      INTEGER i,ii,j,ll, NGRID,LMAX,indx(LMAX)
      double precision rKL(LMAX,LMAX),b(LMAX), sum
      ii= 0
      do  i= 1,NGRID
          ll= indx(i)
          sum= b(ll)
          b(ll)= b(i)
          if (ii.ne.0)then
              do  j= ii,i-1
                  sum= sum-rKL(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii= i
            endif
          b(i)= sum
          enddo
      do  i= NGRID,1,-1
          sum= b(i)
          do  j= i+1,NGRID
              sum= sum-rKL(i,j)*b(j)
              enddo
          b(i)= sum/rKL(i,i)
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

