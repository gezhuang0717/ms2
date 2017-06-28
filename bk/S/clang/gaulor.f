      Real*8 FUNCTION GAULOR(SIG,SIG2,IX)
C      RANDUM NUMBER GENERATOR OF GAUSSIAN DISTRIBUTIO
C      EXP(-(X/SIG)**2/2) : CENTER OF DISTRIBUTION IS 0!
C      IA : SEED OF RANDUM NUMBER
C      VL AND VU RESTRICT THE RANGE OF RANDUM NUMBERS
C      VL : LOWER LIMIT, VU : UPPER LIMIT
      Implicit Double Precision (A-H,O-Z)
*      Real*8 X,VL,VU

*      write(6,*)'sig in gausran : ',sig
      VL = -5.*sig2
      VU =  5.*sig2
 10   X=URAND2(IX)*(VL-VU)+VU
      Y=URAND2(IX)
c      GUS=(X/SIG)**2/2.
      GUS=0.095*EXP(-(X/SIG2)**2/2)+0.905*(1.d0/(X**2+SIG**2))*(SIG**2)

      IF(GUS.LT.Y) GO TO 10
c      write(6,*)X,GUS
      GAULOR=X
      RETURN
      END

