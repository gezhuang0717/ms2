      subroutine mlhfit(XX,YY,NDATA,AA,SA,NCEF) 

      implicit double precision (A-H,O-Z)
      integer NDAT
      dimension XX(1000),YY(1000),ET(1000),TT(1000),FITC(1000)
      dimension ZT(1000),SZ(1000),AA(NCEF),SA(NCEF),SE(1000)
      external func,dfunc,d2func

      DATA TT/1000*1.D0/
      DATA ET/1000*1.D0/
      DATA SE/1000*0.D0/
      DATA ZT/1000*1.D0/
      DATA SZ/1000*0.D0/

c      do i=1,NDATA
c      write(*,*) XX(i), YY(i)
c      enddo

c      write(*,*) XX(NDATA), YY(NDATA)
      NETA=0
      KERR=0
      MODE=2 ! second order differential 

c      write(*,*) AA(1), AA(2)

c      write(*,*) func(XX(1),AA,NCEF,KERR) ! don't remove
      a=func(XX(1),AA,NCEF,KERR) ! don't remove

        call CFDQR2(XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
     &              func,dfunc,d2func,ZT,SZ,CC,SC,NCEFB,FNCTB,MODE,KERR)

      RETURN
      END
