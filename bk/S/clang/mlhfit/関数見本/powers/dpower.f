************************************************************************      subroutine dfunc(x,aa,fa,ncef,kerr)************************************************************************      implicit double precision(a-h,o-z)      dimension aa(30),fa(30)      kerr = 0      fa(1) = 1.d0      if (ncef .eq. 1) return 1000 do 1001 j=2,ncef      fa(j) = fa(j-1) * x 1001 continue      return      end