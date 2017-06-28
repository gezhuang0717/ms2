C===========================================                                              
C     \203\206\203E\223x\212\326\220\224\202\314\221\316\220\224\202\314\214v\216Z        
C===========================================                                              
        DOUBLE PRECISION FUNCTION FLVAL2(YY,FIT,NDATA,INFL)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DIMENSION YY(2000),FIT(2000)

        INFL=0
        FLVAL2=0.D0

 4000      DO 4001 N=1,NDATA
c           WRITE(*,*) N,FIT(N)                                                           
           IF(FIT(N) .LE. 0.D0) GOTO 4020

           FLVAL2=FLVAL2+YY(N)*DLOG(FIT(N))-FIT(N)
 4001         CONTINUE
        write(*,*) '   FLVAL2 = ',FLVAL2
c        write(5,*) '   FLVAL2 = ',FLVAL2                                                 
        RETURN

 4020      WRITE(*,*)
        WRITE(5,*)

        INFL=1

        FLVAL2=0.
        RETURN
        END
