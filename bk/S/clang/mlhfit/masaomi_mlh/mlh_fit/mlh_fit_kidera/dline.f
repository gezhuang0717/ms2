************************************************************************
      subroutine dfunc(x,aa,fa,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30),fa(30)

      kerr = 0

      fa(1) = 1.d0
      fa(2) = x

      return
      end
