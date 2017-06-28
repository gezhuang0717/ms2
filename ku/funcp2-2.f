      FUNCTION FUNCP2(X)
      COMMON/PAWPAR/PAR(16)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCP2 =  (PAR(1) *bin/sqrt2pi/80.0)*EXP(-(X-601.)**2/2./88.0**2)
     &        + (PAR(2) *bin/sqrt2pi/106.)*EXP(-(X-858.)**2/2./106.**2)
     &        + (PAR(3) *bin/sqrt2pi/76.3)*EXP(-(X-1090)**2/2./76.3**2)
     &        + (PAR(4) *bin/sqrt2pi/55.8)*EXP(-(X-1439)**2/2./55.8**2)
     &        + (PAR(5) *bin/sqrt2pi/52.7)*EXP(-(X-1659)**2/2./52.7**2)
     &        + (PAR(6) *bin/sqrt2pi/88.8)*EXP(-(X-1940)**2/2./88.8**2)
     &        + (PAR(7) *bin/sqrt2pi/68.3)*EXP(-(X-2416)**2/2./68.3**2)
     &        + (PAR(8) *bin/sqrt2pi/48.4)*EXP(-(X-2699)**2/2./48.4**2)
     &        + (PAR(9) *bin/sqrt2pi/42.8)*EXP(-(X-2976)**2/2./42.8**2)
     &        + (PAR(10)*bin/sqrt2pi/77.4)*EXP(-(X-3389)**2/2./77.4**2)
     &        + (PAR(11)*bin/sqrt2pi/58.5)*EXP(-(X-3754)**2/2./58.5**2)
     &        + (PAR(12)*bin/sqrt2pi/39.4)*EXP(-(X-3915)**2/2./39.4**2)
     &        + (PAR(13)*bin/sqrt2pi/71.1)*EXP(-(X-4189)**2/2./71.1**2)
     &        + (PAR(14)*bin/sqrt2pi/46.9)*EXP(-(X-4414)**2/2./46.9**2)
     &        + EXP(PAR(15)-PAR(16)*X)
      RETURN
      END

