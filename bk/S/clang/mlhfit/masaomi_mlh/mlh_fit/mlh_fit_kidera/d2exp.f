************************************************************************
      subroutine d2func(x,aa,f2a,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30),f2a(30,30)

      kerr = 0

      s=x-aa(3)
      s2=s**2.d0
      arg=-0.5d0*((x-aa(3))**2.d0)*((aa(2))**-2.d0)
      fnc=dexp(arg)

      f2a(1,1) = 0.d0
      f2a(1,2) = s2*((aa(2))**-3.d0)*fnc
      f2a(1,3) = ((aa(2))**-2.d0)*s*fnc

      f2a(2,1) = f2a(1,2)
      f2a(2,2) = (-3+s2*((aa(2))*-2.d0))
     >           *aa(1)*s2*((aa(2))**-4.d0)*fnc
      f2a(2,3) = (-2+s2*((aa(2))*-2.d0))
     >           *aa(1)*s*((aa(2))**-3.d0)*fnc
      f2a(3,1) = f2a(1,3)
      f2a(3,2) = f2a(2,3)
      f2a(3,3) = (-1+s2*((aa(2))**-2.d0))
     >           *aa(1)*((aa(2))**-2.d0)*fnc

 1000 do 1001 j=1,5
         f2a(j,4) = 0.d0
         f2a(4,j) = 0.d0
         f2a(j,5) = 0.d0
         f2a(5,j) = 0.d0
 1001 continue

      return
      end
