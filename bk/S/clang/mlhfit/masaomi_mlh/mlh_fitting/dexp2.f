c************************************************************************             
      subroutine dfunc(x,aa,fa,ncef,kerr)
c************************************************************************                  
      implicit double precision(a-h,o-z)
      dimension aa(30),fa(30)

      kerr = 0

      arg = - x/aa(2)

      fa(1) = (1/aa(2)) * dexp(arg)
      fa(2) = ( -aa(1)/(aa(2))**2 + aa(1)*x/(aa(2))**3 ) * dexp(arg)
      fa(3) = 1.d0
      fa(4) = x

      return
      end