c for fitting by x** function @ANAPAW
      REAL FUNCTION FUNC(X,Y)   !for a 2-Dim histogram
      COMMON/PAWPAR/PAR(1)
      FUNC=PAR(1)*X+1
      END
