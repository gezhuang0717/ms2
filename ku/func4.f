      FUNCTION FUNC4(X)
      COMMON/PAWPAR/PAR(3)

c      PAR(1) = life1
c      PAR(2) = const1
c      PAR(3) = bg

   real  LN 
   data  LN/0.693147/

      FUNC4 = PAR(2)*EXP(-X/PAR(1)*LN) + PAR(3)
      RETURN
      END
