      FUNCTION FUNCPT2(X)
      COMMON/PAWPAR/PAR(6)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCPT2 = (PAR(1)*bin/sqrt2pi/37.0)*EXP(-(X-2972.)**2/2./37.0**2)
     &        + (PAR(2)*bin/sqrt2pi/70.4)*EXP(-(X-3389.)**2/2./70.4**2)
     &        + (PAR(3)*bin/sqrt2pi/88.5)*EXP(-(X-3754.)**2/2./88.5**2)
     &        + (PAR(4)*bin/sqrt2pi/39.4)*EXP(-(X-3911.)**2/2./35.2**2)
     &        + EXP(PAR(5)-PAR(6)*X)
      RETURN
      END

