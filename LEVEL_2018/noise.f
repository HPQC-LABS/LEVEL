c** Read input smooth vib-rot energies and use random number generator
c  to add noise
c
      implicit real*8 (a-h,o-z)
      character*4 title(20)
      read (5,601) (title(i),i=1,20)
      write(6,601) (title(i),i=1,20)
      read (5,601) (title(i),i=1,20)
      write(6,601) (title(i),i=1,20)
  601 format(20A4)
      call g05cbf(22)
   20 read(5,*) iv,ij,e,ue
      ue= 0.0005 + ue
      e1 = e+ 0.003*(g05caf(1.d0)-0.5d0)
      de= e1-e
      write(6,602) iv,ij,e1,ue
  602 format(2i8,f13.4,e15.4,f13.4)
      go to 20
   99 stop
      end
