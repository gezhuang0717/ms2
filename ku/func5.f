      FUNCTION FUNC5(X)
      COMMON/PAWPAR/PAR(5)

c      PAR(1) = life1
c      PAR(2) = const1
c      PAR(3) = life2
c      PAR(4) = const4
c      PAR(5) = bg

   real  LN 
   data  LN/0.693147/

      FUNC5 = PAR(2)*EXP(-X/PAR(1)*LN) + PAR(4)*EXP(-X/PAR(3)*LN)+PAR(5)
      RETURN
      END
