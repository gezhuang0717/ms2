      Real*8 FUNCTION GAUSRAN(SIG,IX)
C      RANDUM NUMBER GENERATOR OF GAUSSIAN DISTRIBUTIO
C      EXP(-(X/SIG)**2/2) : CENTER OF DISTRIBUTION IS 0!
C      IA : SEED OF RANDUM NUMBER
C      VL AND VU RESTRICT THE RANGE OF RANDUM NUMBERS
C      VL : LOWER LIMIT, VU : UPPER LIMIT
      Implicit Double Precision (A-H,O-Z)
*      Real*8 X,VL,VU

*      write(6,*)'sig in gausran : ',sig

      VL = -5.*sig
      VU =  5.*sig
 10   X=URAND2(IX)*(VL-VU)+VU
      Y=URAND2(IX)
      GUS=(X/SIG)**2/2.
      GUS=EXP(-GUS)
c      write(6,*)X,Y,GUS,SIG

      IF(GUS.LT.Y) GO TO 10
      GAUSRAN=X
      RETURN
      END

