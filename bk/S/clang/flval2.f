C===========================================
C     �楦�ٴؿ����п��η׻�
C===========================================
        DOUBLE PRECISION FUNCTION FLVAL2(YY,FIT,NDATA,INFL)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DIMENSION YY(2000),FIT(2000)

        INFL=0
        FLVAL2=0.D0

 4000   DO 4001 N=1,NDATA
c           WRITE(*,*) N,FIT(N)
ccc           IF(FIT(N) .LE. 0.D0) GOTO 4020

           FLVAL2=FLVAL2+YY(N)*DLOG(FIT(N))-FIT(N)
 4001   CONTINUE
cc        write(*,*) '   FLVAL2 = ',FLVAL2
c        write(5,*) '   FLVAL2 = ',FLVAL2
        RETURN

ccc 4020   WRITE(*,9000)
ccc        WRITE(5,9000)
ccc* 9000   FORMAT (1H ,"�ؿ��ͤ�����ޤ�������ͤǤ�")
ccc 9000   FORMAT (1H ,"Parameter is less than 0")

        INFL=1

        FLVAL2=0.
        RETURN
        END
