************************************************************************
      double precision function func(x,aa,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30)

      kerr = 0

      arg = - 0.5d0*((x-aa(3))**2.d0)*((aa(2))**-2.d0)
      func = aa(1)*dexp(arg)
      func=func+aa(4)
      func=func+aa(5)*x

      return
      end
