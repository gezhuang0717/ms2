      subroutine d2func(x,aa,f2a,ncef,kerr)

      implicit double precision(a-h,o-z)
      dimension aa(30),f2a(30,30)

      kerr = 0

      arg = - x/aa(2)

      f2a(1,1) = 0.d0
      f2a(1,2) = ( -1/(aa(2))**2 + x/(aa(2))**3 ) * dexp(arg)
      f2a(1,3) = 0.d0
      f2a(1,4) = 0.d0

      f2a(2,1) = f2a(1,2)
      f2a(2,2) = ( ( 2*aa(1)/(aa(2))**3 - 3*aa(1)*x/(aa(2))**4 )
     >                  +( -aa(1)/(aa(2))**2 + aa(1)*x/(aa(2))**3) )
     >           * dexp(arg)
      f2a(2,3) = 0.d0
      f2a(2,4) = 0.d0

      f2a(3,1) = 0.d0
      f2a(3,2) = 0.d0
      f2a(3,3) = 0.d0
      f2a(3,4) = 0.d0

      f2a(4,1) = 0.d0
      f2a(4,2) = 0.d0
      f2a(4,3) = 0.d0
      f2a(4,4) = 0.d0

      return
      end
