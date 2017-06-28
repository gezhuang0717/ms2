      FUNCTION FUNCX(X)
      COMMON/PAWPAR/PAR(16)

      real sqrt2pi
      real bin
      data sqrt2pi/2.5067588/
      data bin/20./

      FUNCX =  (PAR(1)*bin/sqrt2pi/49.4)*EXP(-(X-617.0)**2/2./49.4**2)
     &       + (PAR(2)*bin/sqrt2pi/63.6)*EXP(-(X-856.4)**2/2./63.6**2)
     &       + (PAR(3)*bin/sqrt2pi/57.5)*EXP(-(X-1093)**2/2./57.5**2)
     &       + (PAR(4)*bin/sqrt2pi/50.1)*EXP(-(X-1437)**2/2./50.1**2)
     &       + (PAR(5)*bin/sqrt2pi/44.9)*EXP(-(X-1654)**2/2./44.9**2)
     &       + (PAR(6)*bin/sqrt2pi/41.6)*EXP(-(X-1945)**2/2./41.6**2)
     &       + (PAR(7)*bin/sqrt2pi/38.8)*EXP(-(X-2406)**2/2./38.8**2)
     &       + (PAR(8)*bin/sqrt2pi/38.3)*EXP(-(X-2701)**2/2./38.3**2)
     &       + (PAR(9)*bin/sqrt2pi/37.0)*EXP(-(X-2972)**2/2./37.0**2)
     &       + (PAR(10)*bin/sqrt2pi/98.7)*EXP(-(X-3364)**2/2./98.7**2)
     &       + (PAR(11)*bin/sqrt2pi/99.9)*EXP(-(X-3765)**2/2./99.9**2)
     &       + (PAR(12)*bin/sqrt2pi/35.2)*EXP(-(X-3911)**2/2./35.2**2)
     &       + (PAR(13)*bin/sqrt2pi/60.3)*EXP(-(X-4181)**2/2./60.3**2)
     &       + (PAR(14)*bin/sqrt2pi/35.1)*EXP(-(X-4423)**2/2./35.1**2)
     &       + PAR(15)*(EXP(1/(X*PAR(16))-1))
      RETURN
      END

