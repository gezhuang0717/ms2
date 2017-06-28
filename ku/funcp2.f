      FUNCTION FUNCP2(X)
      COMMON/PAWPAR/PAR(18)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCP2 =  (PAR(1)*bin/sqrt2pi/49.4)*EXP(-(X-617.)**2/2./49.4**2)
     &        + (PAR(2)*bin/sqrt2pi/70.6)*EXP(-(X-856.)**2/2./70.6**2)
     &        + (PAR(3)*bin/sqrt2pi/62.5)*EXP(-(X-1091.)**2/2./62.5**2)
     &        + (PAR(4)*bin/sqrt2pi/53.1)*EXP(-(X-1437.)**2/2./53.1**2)
     &        + (PAR(5)*bin/sqrt2pi/47.9)*EXP(-(X-1652.)**2/2./47.9**2)
     &        + (PAR(6)*bin/sqrt2pi/58.0)*EXP(-(X-1940.)**2/2./58.0**2)
     &        + (PAR(7)*bin/sqrt2pi/38.8)*EXP(-(X-2408.)**2/2./38.8**2)
     &        + (PAR(8)*bin/sqrt2pi/38.3)*EXP(-(X-2701.)**2/2./38.3**2)
     &        + (PAR(9)*bin/sqrt2pi/37.0)*EXP(-(X-2975.)**2/2./37.0**2)
     &        + (PAR(10)*bin/sqrt2pi/60.1)*EXP(-(X-3364.)**2/2./60.1**2)
     &        + (PAR(11)*bin/sqrt2pi/88.5)*EXP(-(X-3754.)**2/2./88.5**2)
     &        + (PAR(12)*bin/sqrt2pi/35.2)*EXP(-(X-3911.)**2/2./35.2**2)
     &        + (PAR(13)*bin/sqrt2pi/55.1)*EXP(-(X-4186.)**2/2./55.1**2)
     &        + (PAR(14)*bin/sqrt2pi/35.1)*EXP(-(X-4423.)**2/2./35.1**2)
     &        + (PAR(15)*bin/sqrt2pi/49.3)*EXP(-(X-4673.)**2/2./49.3**2)
     &        + (PAR(16)*bin/sqrt2pi/54.0)*EXP(-(X-4956.)**2/2./54.0**2)
     &        + EXP(PAR(17)-PAR(18)*X)
      RETURN
      END

