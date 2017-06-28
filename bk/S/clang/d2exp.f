      subroutine d2func(x,aa,f2a,ncef,kerr)
      implicit double precision (A-H,O-Z)
      dimension aa(30),f2a(30,30)
      kerr=0
c 2 parameter
      arg=aa(2)*x
      func=dexp(arg)
      f2a(1,1)=0.d0
      f2a(1,2)=x*func
      f2a(2,1)=f2a(1,2)
      f2a(2,2)=aa(1)*x*x*func
c 3 parameters
c      s1=x+aa(3)
c      s2=s1**2.0d0
c      arg=aa(2)*s1
c      func=dexp(arg)
c      f2a(1,1)=0.d0
c      f2a(1,2)=s1*func
c      f2a(1,3)=aa(2)*func
c      f2a(2,1)=s1*func
c      f2a(2,2)=aa(1)*s2*func
c      f2a(2,3)=aa(1)*func*(1.d0+aa(2)*s1)
c      f2a(3,1)=aa(2)*func
c      f2a(3,2)=aa(1)*func*(aa(2)*x+aa(2)*aa(3)+1.d0)
c      f2a(3,3)=aa(1)*(aa(2)**2.d0)*func
c
c      do j=1,5
c      f2a(j,4)=0.d0
c      f2a(4,j)=0.d0
c      f2a(j,5)=0.d0
c      f2a(5,j)=0.d0
c      enddo
 
      return
      end 






