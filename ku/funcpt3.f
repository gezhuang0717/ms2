      FUNCTION FUNCPT3(X)
      COMMON/PAWPAR/PAR(7)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCPT3 = (PAR(1)*bin/sqrt2pi/58.0)*EXP(-(X-1945.)**2/2./58.0**2)
     &        + (PAR(2)*bin/sqrt2pi/38.8)*EXP(-(X-2208.)**2/2./38.8**2)
     &        + (PAR(3)*bin/sqrt2pi/58.0)*EXP(-(X-2406.)**2/2./58.0**2)
     &        + (PAR(4)*bin/sqrt2pi/38.2)*EXP(-(X-2701.)**2/2./38.2**2)
     &        + (PAR(5)*bin/sqrt2pi/37.0)*EXP(-(X-2972.)**2/2./37.0**2)
     &        + EXP(PAR(6)-PAR(7)*X)
      RETURN
      END

