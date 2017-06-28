      FUNCTION FUNCE(X)
      COMMON/PAWPAR/PAR(3)

c      PAR(1) = constant
c      PAR(2) = life
c      PAR(3) = offset
      FUNCE= PAR(1)*EXP(-X/(PAR(2)) + PAR(3)
      RETURN
      END
