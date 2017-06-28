      subroutine dfunc(x,aa,fa,ncef,kerr)
      implicit double precision (A-H,O-Z)
      dimension aa(30),fa(30)
      kerr=0     
      func=dexp(aa(2)*x)
      fa(1)=func
      fa(2)=aa(1)*x*func
c 3 parameter
c      func=dexp(aa(2)*(x+aa(3)))
c      fa(1)=func
c      fa(2)=aa(1)*(x+aa(3))*func
c      fa(3)=aa(1)*aa(2)*func

      return
      end
