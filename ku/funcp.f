      FUNCTION FUNCP2(X)
      COMMON/PAWPAR/PAR(16)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCP2 =  (PAR(1) *bin/sqrt2pi/88.8)*EXP(-(X-599.)**2/2./88.8**2)
     &        + (PAR(2) *bin/sqrt2pi/88.4)*EXP(-(X-867.)**2/2./88.4**2)
     &        + (PAR(3) *bin/sqrt2pi/65.1)*EXP(-(X-1093)**2/2./65.1**2)
     &        + (PAR(4) *bin/sqrt2pi/49.7)*EXP(-(X-1444)**2/2./49.7**2)
     &        + (PAR(5) *bin/sqrt2pi/52.0)*EXP(-(X-1650)**2/2./52.0**2)
     &        + (PAR(6) *bin/sqrt2pi/58.0)*EXP(-(X-1945)**2/2./58.0**2)
     &        + (PAR(7) *bin/sqrt2pi/52.7)*EXP(-(X-2409)**2/2./52.7**2)
     &        + (PAR(8) *bin/sqrt2pi/41.3)*EXP(-(X-2706)**2/2./41.3**2)
     &        + (PAR(9) *bin/sqrt2pi/45.6)*EXP(-(X-2978)**2/2./45.6**2)
     &        + (PAR(10)*bin/sqrt2pi/122.5)*EXP(-(X-3385)**2/2./122.5**2)
     &        + (PAR(11)*bin/sqrt2pi/69.8)*EXP(-(X-3800)**2/2./69.8**2)
     &        + (PAR(12)*bin/sqrt2pi/31.2)*EXP(-(X-3912)**2/2./31.2**2)
     &        + (PAR(13)*bin/sqrt2pi/49.9)*EXP(-(X-4183)**2/2./49.9**2)
     &        + (PAR(14)*bin/sqrt2pi/96.0)*EXP(-(X-4454)**2/2./96.0**2)
     &        + EXP(PAR(15)-PAR(16)*X)
      RETURN
      END

