      subroutine rayfit(ix,ixz,np,fx,ft,res,ctrl,ierr)
      implicit none

      integer np
      real ix(np)
      real ixz(np)
      real fx, ft
      real res
      integer ctrl,ierr
      
      real a(2,2),b(2),w(2)
      integer i,j
      real det,dx

      do i = 1,2
         b(i) = 0.
         do j = 1,2
            a(i,j) = 0.
         enddo
      enddo
      if (np .lt. 2) then
         ierr = 1
         return
      endif
      do i = 1,np
         b(1) = b(1) + ix(i) * ixz(i)
         b(2) = b(2) + ix(i)
         a(1,1) = a(1,1) + ixz(i)**2
         a(1,2) = a(1,2) + ixz(i)
         a(2,2) = a(2,2) + 1.0
      enddo
      a(2,1) = a(1,2)
      det = a(1,1)*a(2,2)-a(2,1)*a(1,2)
      if (det .eq. 0) then 
         ierr = 1
         return
      endif
      w(1) = a(2,2)/det * b(1) - a(1,2)*b(2)/det
      w(2) = -a(2,1)/det *b(1) + a(1,1)/det*b(2)
      fx = w(2)
      ft = w(1)
      do i = 1,np
         dx = abs(ix(i) - (ft*ixz(i)+fx))
         res = res + dx**2
      enddo
      res = res/np
      ierr = 0
      return
      end
