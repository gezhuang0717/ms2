      FUNCTION FUNCPT(X)
      COMMON/PAWPAR/PAR(6)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCPT = (PAR(1)*bin/sqrt2pi/60.3)*EXP(-(X-4181.)**2/2./60.3**2)
     &        + (PAR(2)*bin/sqrt2pi/45.1)*EXP(-(X-4423.)**2/2./45.1**2)
     &        + (PAR(3)*bin/sqrt2pi/49.3)*EXP(-(X-4673.)**2/2./49.3**2)
     &        + (PAR(4)*bin/sqrt2pi/62.9)*EXP(-(X-4953.)**2/2./62.9**2)
     &        + EXP(PAR(5)-PAR(6)*X)
      RETURN
      END

