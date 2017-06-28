      FUNCTION FUNC1(X)
      COMMON/PAWPAR/PAR(5)

c      PAR(1) = constant
c      PAR(2) = mean
c      PAR(3) = sigma
c      PAR(4) = slope prm1
c      PAR(5) = offset

      real sqrt2pi
      data sqrt2pi/2.5067588/

      FUNC1 = (PAR(1)/sqrt2pi/PAR(3))*EXP(-(X-PAR(2))**2/2./PAR(3)**2)
     &       + PAR(4)*X +PAR(5)
      RETURN
      END

