      Real*8 FUNCTION SOGRAN(sig,IX)
C      RANDUM NUMBER GENERATOR OF GAUSSIAN DISTRIBUTIO
C      EXP(-(X/SIG)**2/2) : CENTER OF DISTRIBUTION IS 0!
C      IA : SEED OF RANDUM NUMBER
C      VL AND VU RESTRICT THE RANGE OF RANDUM NUMBERS
C      VL : LOWER LIMIT, VU : UPPER LIMIT
      Implicit Double Precision (A-H,O-Z)
*      Real*8 X,VL,VU

*      write(6,*)'sig in gausran : ',sig
      VL = 50.
      VU = 150.
 10   X=URAND2(IX)*(VL-VU)+VU
      Y=URAND2(IX)
c      GUS=(X/SIG)**2/2.
      GUS=(-1.9718+7565.87*EXP(-((X-97.5938)/16.9141)**2)
     &          +10928*EXP(-((X-90.6404)/9.58397)**2))/17353.
c&          +4408*EXP(-((X-0.0839)/0.1)**2))

      IF(GUS.LT.Y) GO TO 10
c      write(6,*)X,Y,GUS

      SOGRAN=(X-95.)/24.
      RETURN
      END

