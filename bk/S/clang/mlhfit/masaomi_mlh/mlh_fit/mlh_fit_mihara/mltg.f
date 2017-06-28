c=================================
c     Multi-component Gaussian
c=================================
      double precision function fmltg(x,AA,NCEF)
      implicit double precision(a-h,o-z)
      dimension AA(30),c(10,3),g(10)
      COMMON/user/pi,DX,fr
      COMMON /PARAM/AAI(30),II(30),NCOEF,ncmp,ic

      fmltg = 0.0d0

      do 10 i = 1 , ncmp
         do 100 j = 1 , 3
            c(i,j) = aa(3*(i-1)+j)
 100     continue
         if(c(i,2).eq.0.0) then
            g(i) = 0.0d0
         else
            g(i) = c(i,1)/(2.d0*pi)**0.5/c(i,2)
     &            *dexp(-.5d0*((x-c(i,3))/c(i,2))**2)
         endif
         fmltg = fmltg + g(i)
 10   continue

      return
      end

c=================================
c     response function (gaussian + gaussian)
c=================================
      double precision function resp(x,y,AA,NCEF)
      implicit double precision(a-h,o-z)
      dimension AA(30)
      COMMON/user/pi,DX,fr

      pi = dacos(-1.d0)

      f = 1.d0
cc      f = 1.1 d0
c     a1 = 168039.5137
c     a2 = 0.1976351676 * f
c     a3 = 23.7236603
c     b1 = 35234.84288
c     b2 = 0.2311057893 * f
c     b3 = 23.87983774

      a1 = 165743.
      a2 = 0.200232 * f
      a3 = 23.7964
      b1 = 20171.5
      b2 = 0.241728 * f
      b3 = 24.1728

      c = 1./(2.*pi)**0.5/(a1+b1)
      gr = (a1*a3+b1*b3)/(a1+b1)

      ta = a3 - gr
      tb = b3 - gr

c     dt_dp = -0.0102431 + 2.*5.65684e-6*y - 3.*3.0277e-9*y**2
c    &       + 4.*8.2773e-12*y**3
      dt_dp = -0.0108133 + 2.*6.70958e-6*y - 3.*4.0599e-9*y**2

      pa = ta / dt_dp
      pb = tb / dt_dp

      sa = a2 / dabs(dt_dp)
      sb = b2 / dabs(dt_dp)

      ca = c * a1 / sa
      cb = c * b1 / sb

      g1 = ca * exp(-0.5*(x-y-pa)**2/sa**2)
      g2 = cb * exp(-0.5*(x-y-pb)**2/sb**2)

c     write(*,*) pa,pb,sa,sb
c     write(*,*) ca,cb

      resp = (g1 + g2) * fmltg(y,aa,ncef)

      return
      end

c=================================
c     Convoluting resolution
c=================================
      DOUBLE PRECISION FUNCTION fmltig(X,AA0,NCEF,KERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AA(30),AA0(30)
      COMMON/user/pi,DX,fr
      COMMON /PARAM/AAI(30),II(30),NCOEF,ncmp,ic

      EXTERNAL resp

      KERR = 0
      pi = acos(-1.0)

      DO 100 I = 1,NCEF
         AA(II(I)) = AA0(I)
 100  CONTINUE
      IF(NCOEF.GT.NCEF) THEN
         DO 101 I = NCEF+1,NCOEF
            AA(II(I)) = AAI(II(I))
 101     CONTINUE
      ENDIF
      if(ic.ne.0) then                    !ic = 1 : center are common
         do 10 i = 1,ncmp
            aa(3*(i-1)+3) = aa(3)
 10      continue
      endif
      DO 1 I = 1,NCOEF
         AA0(I) = AA(II(I))
 1    CONTINUE

      call integ(x,AA,-400.,400.,resp,ss,NCEF,200,n1)

      fmltig = SS * dx

c      write(*,*) n1,n2,n3

      RETURN
      END

c=================================
c     dmltig （一回微分）
c=================================
      SUBROUTINE dmltig(X,AA,DAA,FA,NCEF,KERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AA(30),FA(30),DAA(30)

      KERR = 0

      DO 10 I = 1, NCEF
         aj = AA(I)
         AA(I) = aj + DAA(I)
         fu = fmltig(X,AA,NCEF,KERR)
         AA(I) = aj - DAA(I)
         fl = fmltig(X,AA,NCEF,KERR)
         FA(I) = 0.5 * (fu - fl)/DAA(I)
         AA(I) = aj
 10   CONTINUE

c     do 11 i = 1,NCEF
c        write(*,*) I,II(I),FA(I)
c11   continue

      RETURN
      END

c=================================
c     d2mltig （二回微分）
c=================================
      SUBROUTINE d2mltig(X,AA,F2A,NCEF,KERR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AA(30),F2A(30,30)

      KERR = 0

      DO 1 I = 1,NCEF
         DO 10 J = 1,NCEF
            F2A(I,J) = 0.0
 10      CONTINUE
 1    CONTINUE

      RETURN
      END

c==========================================
c     SUBROUTINE FOR SIMPSON INTEGRATION
c==========================================
      SUBROUTINE INTEG(X,AA,RI,RF,FUNC,AINTG,NCEF,NI,N1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AA(30)
      COMMON/user/pi,DX,fr
      COMMON /PARAM/AAI(30),II(30),NCOEF,ncmp,ic

c      SMPERR = 1.0D-4
      SMPERR = 1.0D-6
      BINTG = 0.0
      RANGE = RF - RI

      SUMIF = FUNC(X,RI,AA,NCEF)+FUNC(X,RF,AA,NCEF)
 
c      N1 = 60
      N1 = NI
C         .......     N1*2 IS NUMBER OF INTEGRATION RANGE CUT

      DH1 = RANGE/DFLOAT(N1)/2.0
      DH2 = DH1*2.0

      Y = RI + DH2
      SUME = 0.0
      DO 44 I = 1, N1-1
        SUME = SUME + FUNC(X,Y,AA,NCEF)
        Y = Y + DH2
 44   CONTINUE

 1    CONTINUE

      Y = RI + DH1
      SUMO = 0.0
      DO 10 JJ = 1, N1
        SUMO = SUMO + FUNC(X,Y,AA,NCEF)
        Y = Y + DH2
 10   CONTINUE

      AINTG = DH1/3.0 * (SUMIF + 2.0*SUME + 4.0*SUMO)
      IF(AINTG.EQ.0.0) GO TO 2
      D = (AINTG - BINTG)/AINTG
      IF(DABS(D)-SMPERR) 3,3,2

      goto 3

 2    BINTG = AINTG
      SUME = SUME + SUMO
      DH2 = DH1
      DH1 = DH1/2.0
      N1 = N1*2

c      IF(N1.GT.200) GOTO 3
      IF(N1.GT.20000) GOTO 3

      GO TO 1

 3    RETURN
      END
