      real function pol1(val,prm)
      implicit none
      real val
      real prm(2)
      call poln(val,pol1,prm,1)
      return
      end

      real function pol2(val,prm)
      implicit none
      real val
      real prm(3)
      call poln(val,pol2,prm,2)
      return
      end

      real function pol3(val,prm)
      implicit none
      real val
      real prm(4)
      call poln(val,pol3,prm,3)
      return
      end


      subroutine poln(raw,cal,prm,n) 
      implicit none
      integer i,n
      real raw,cal,prm
      dimension prm(n+1)
      cal = 0.
      do i = 1,n+1
         cal = cal + prm(n)*raw**(n-1)
      enddo
      return
      end

      logical function gate1d(x,xmin,xmax,ctrl)
      implicit none
c  ctrl=0  xmin.gt.x.lt.xmax
c  ctrl=1  xmin.ge.x.le.xmax
      real x,xmin,xmax
      integer ctrl
      gate1d = .FALSE.
      if (ctrl.eq.1) then
         if ((x.ge.xmin).and.(x.le.xmax)) then
            gate1d = .TRUE.
         endif
      else
         if ((x.gt.xmin).and.(x.lt.xmax)) then
            gate1d = .TRUE.
         endif
      endif
      return
      end
