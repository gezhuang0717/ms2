C============================================
C     修正グラムシュミット法による
C               行列のQR変換
C============================================

        SUBROUTINE QRTRN(AQ,RR,JNX,NDATA,NCEF,NRANK)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION AQ(2030,30),RR(30,30),JNX(30),RN(30)

        EPSM=1.D-10

 8000   DO 8001 J=1,NCEF
           JNX(J)=J
           RN(J)=1.D0
 8020      DO 8021 K=1,NCEF
              RR(J,K)=0.D0
 8021      CONTINUE
 8001   CONTINUE

        NRANK=0

 8100   DO 8101 J=1,NCEF
           AMAX=0.D0
 8120      DO 8121 K1=J,NCEF
              K=JNX(K1)
              RK=0.D0

 8140         DO 8141 N=1,NDATA
                 RK=RK+AQ(N,K)*AQ(N,K)
 8141         CONTINUE

              RK=DSQRT(RK)
              RA=RK/RN(K)

              IF(J .NE. 1) GOTO 8160
              RN(K)=RK

 8160         IF (AMAX .GE. RA .AND. J .LT. NCEF) GOTO 8121
              AMAX=RA
              RM=RK
              M=K1
 8121      CONTINUE

           MSAVE=JNX(J)
           JNX(J)=JNX(M)
           JNX(M)=MSAVE

           M=JNX(J)
           RR(J,M)=RM

           IF(AMAX .LE. EPSM) GOTO 8500

           NRANK = NRANK+1

 8200      DO 8201 N=1,NDATA
              AQ(N,M)=AQ(N,M)/RR(J,M)
 8201      CONTINUE

           IF (J .EQ. NCEF) RETURN

 8220      DO 8221 K1=J+1,NCEF
              K=JNX(K1)

              DO 8241 N=1,NDATA
                 RR(J,K)=RR(J,K)+AQ(N,M)*AQ(N,K)
 8241         CONTINUE

 8260         DO 8261 N=1,NDATA
                 AQ(N,K)=AQ(N,K)-AQ(N,M)*RR(J,K)
 8261         CONTINUE

 8221      CONTINUE

 8101   CONTINUE

       RETURN

 8500  DO 8501 K=NRANK+1,NCEF
          KK=JNX(K)
 8520     DO 8521 L=1,NCEF
             RR(L,KK)=0.D0
 8521     CONTINUE
 8501  CONTINUE

       RETURN
       END
