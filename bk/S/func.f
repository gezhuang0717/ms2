c for fitting by x** function @ ANAPAW                                         
c  
      REAL FUNCTION FUNC(X,Y) !for a 2-Dim histogram                          
c   
      COMMON/PAWPAR/PAR(1)

      FUNC = 1. - PAR(1)*X 

      END
