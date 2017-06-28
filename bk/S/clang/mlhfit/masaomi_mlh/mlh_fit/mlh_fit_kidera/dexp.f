************************************************************************
      subroutine dfunc(x,aa,fa,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30),fa(30)

      kerr = 0

      s=x-aa(3)
      s2=s*s
      arg = - 0.5d0*((x-aa(3))**2.d0)*((aa(2))**-2.d0)
      fnc=dexp(arg)
      fa(1) = fnc
      fa(2) = aa(1)*s2*((aa(2))**-3.d0)*fnc
      fa(3) = aa(1)*((aa(2))**-2.d0)*s*fnc
      fa(4)=1.d0
      fa(5)=x

      return
      end
