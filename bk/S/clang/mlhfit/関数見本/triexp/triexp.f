************************************************************************      double precision function func(x,aa,ncef,kerr)************************************************************************      implicit double precision(a-h,o-z)      dimension aa(30)      kerr = 0      func = aa(1) * dexp( aa(2) * x ) +      >       aa(3) * dexp( aa(4) * x ) +      >       aa(5) * dexp( aa(6) * x )       return      end