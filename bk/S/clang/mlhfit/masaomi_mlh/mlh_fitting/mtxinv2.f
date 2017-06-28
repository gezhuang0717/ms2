C================================================                                         
C    Gauss Jordan\226@\202\311\202\346\202\351\213t\215s\227\361\202\314\214v\216Z        
C================================================                                         
        SUBROUTINE MTXINV(ARRAY,NC,DET)

        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION ARRAY(30,30),JA(30),KA(30)

        DET=1.D0

 2000      DO 2001 L=1,NC

        AMAX=0.D0

 2020      DO 2021 J=L,NC
 2040               DO 2041 K=L,NC
              IF(DABS(AMAX)-DABS(ARRAY(J,K)) .GT. 0.D0) GOTO 2041
              AMAX=ARRAY(J,K)
              JA(L)=J
              KA(L)=K
 2041               CONTINUE
 2021                  CONTINUE

        IF (AMAX .NE. 0.D0) GOTO 2100
        DET=0.D0
        RETURN

 2100      J=JA(L)
        IF(J .EQ. L) GOTO 2200
 2120      DO 2121 K=1,NC
           SAVE=ARRAY(L,K)
           ARRAY(L,K)=ARRAY(J,K)
           ARRAY(J,K)=-SAVE
 2121         CONTINUE

 2200            K=KA(L)
        IF(K .EQ. L) GOTO 2300
 2220      DO 2221 J=1,NC
           SAVE=ARRAY(J,L)
           ARRAY(J,L)=ARRAY(J,K)
           ARRAY(J,K)=-SAVE
 2221         CONTINUE

 2300            DO 2301 J=1,NC
           IF (J .EQ. L) GOTO 2301
           ARRAY(J,L)=-ARRAY(J,L)/AMAX
 2301         CONTINUE

 2320            DO 2321 J=1,NC
           IF(J .EQ. L) GOTO 2321
           DO 2341 K=1,NC
              IF (K .EQ. L) GOTO 2341
              ARRAY(J,K)=ARRAY(J,K)+ARRAY(J,L)*ARRAY(L,K)
 2341               CONTINUE
 2321                  CONTINUE

 2360                     DO 2361 K=1,NC
           IF(K .EQ. L) GOTO 2361
           ARRAY(L,K)=ARRAY(L,K)/AMAX
 2361         CONTINUE

        ARRAY(L,L)=1.D0/AMAX
        DET=DET*AMAX

 2001      CONTINUE

 2500         DO 2501 L=1,NC
           M=NC-L+1
           K=JA(M)
           IF (K .LE. M) GOTO 2540
 2520            DO 2521 J=1,NC
              SAVE=ARRAY(J,M)
              ARRAY(J,M)=-ARRAY(J,K)
              ARRAY(J,K)=SAVE
 2521               CONTINUE


 2540                     J=KA(M)
           IF (J .LE. M) GOTO 2501
 2560            DO 2561 K=1,NC
              SAVE=ARRAY(M,K)
              ARRAY(M,K)=-ARRAY(J,K)
              ARRAY(J,K)=SAVE
 2561               CONTINUE
 2501                  CONTINUE

        RETURN
        END
