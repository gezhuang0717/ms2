      FUNCTION FUNC3(X)
      COMMON/PAWPAR/PAR(11)

c      PAR(1) = constant
c      PAR(2) = mean
c      PAR(3) = sigma
c      PAR(4) = constant
c      PAR(5) = mean
c      PAR(6) = sigma
c      PAR(7) = constant
c      PAR(8) = mean
c      PAR(9) = sigma
c      PAR(10) = slope prm1
c      PAR(11) = offset

      real sqrt2pi
      data sqrt2pi/2.5067588/

      FUNC3 =  (PAR(1)/sqrt2pi/PAR(3))*EXP(-(X-PAR(2))**2/2./PAR(3)**2)
     &       + (PAR(4)/sqrt2pi/PAR(6))*EXP(-(X-PAR(5))**2/2./PAR(6)**2)
     &       + (PAR(7)/sqrt2pi/PAR(9))*EXP(-(X-PAR(8))**2/2./PAR(9)**2)
     &       + PAR(10)*X + PAR(11)
      RETURN
      END

