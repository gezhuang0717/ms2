      FUNCTION FUNC8(X)
      COMMON/PAWPAR/PAR(5)

c      PAR(1) = constant
c      PAR(2) = mean
c      PAR(3) = constant
c      PAR(4) = mean
c      PAR(5) = sigma

      real sqrt2pi
      data sqrt2pi/2.5067588/

      FUNC8 =  (PAR(1)/sqrt2pi/PAR(5))*EXP(-(X-PAR(2))**2/2./PAR(5)**2)
     &       + (PAR(3)/sqrt2pi/PAR(5))*EXP(-(X-PAR(4))**2/2./PAR(5)**2)
      RETURN
      END

