************************************************************************
      double precision function func(x,aa,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30)

      kerr = 0

      func = aa(1)+aa(2)*x

      return
      end
