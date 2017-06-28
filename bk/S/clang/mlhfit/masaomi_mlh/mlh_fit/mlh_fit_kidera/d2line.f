************************************************************************
      subroutine d2func(x,aa,f2a,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30),f2a(30,30)

      kerr = 0

      f2a(1,1) = 0.d0
      f2a(1,2) = 0.d0
      f2a(2,1) = 0.d0
      f2a(2,2) = 0.d0

      return
      end
