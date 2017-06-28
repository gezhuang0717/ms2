      FUNCTION FUNCR(X)
      COMMON/PAWPAR/PAR(4)

c      PAR(1) = amplitude
c      PAR(2) = g-factor
c      PAR(3) = ini-phase
c      PAR(4) = const
c
      real pi
      data pi/3.1415926/

      real b0
c      data b0/0.259225/
      data b0/0.5/

      real sqrt2pi
      data sqrt2pi/2.5067588/

      real larmor
      data larmor/0.047895/

      FUNCR= PAR(1)*SIN(PAR(2)*b0*larmor*X*2. + PAR(3)) + PAR(4)
      RETURN
      END
