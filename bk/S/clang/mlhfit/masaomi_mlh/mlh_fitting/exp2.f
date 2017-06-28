c************************************************************************                 
      double precision function func(x,aa,ncef,kerr)
c************************************************************************                 
      implicit double precision(a-h,o-z)
      dimension aa(30)

      kerr = 0

      arg  = - x/aa(2)
      func = (aa(1)/aa(2)) * dexp(arg) + aa(3) + aa(4) * x

      return
      end
