      double precision function func(x,aa,ncef,kerr)
      implicit double precision (A-H,O-Z)
      dimension aa(30)
      kerr=0
      func=aa(1)*dexp(aa(2)*x) ! 2 parameter exponential
c      func=aa(1)*dexp(aa(2)*(x+aa(3))) ! 3 parameter exponential
      return
      end
