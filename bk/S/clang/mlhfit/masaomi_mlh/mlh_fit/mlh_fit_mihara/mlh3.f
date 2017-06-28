C================================================
C     Maximum Likelihood法によるFitting Programm
C================================================

        SUBROUTINE CFDQR2(XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,FNCTD,
     >                    FNDRV,F2DRV,ZT,SZ,CC,SC,NCEFB,FNCTB,MODE,KERR)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XX(500),YY(500),ET(500),TT(500),FITC(500),SE(500)
        DIMENSION ZT(500),SZ(500),YD(510),AA(30),DA(30),SA(30),FA(30)
        DIMENSION CC(30),SC(30),FC(30),QQ(500,30),QX(510,30),DD(30)
        DIMENSION RR(30,30),JNX(30),F2A(30,30),DFDY(500,30),DFDA(30,30)
        DIMENSION BB(30)
        DIMENSION QV(30)
        DIMENSION QZ(500,30)
        DIMENSION DAA(30),DBB(30),DCC(30)
        COMMON /PARAM/AAI(30),II(30),NCOEF,ncmp,ic

C     このProgramは”データ解析”の本に載っているものです。

	EPSM=1.D-10
	EPSF=1.D-6
	PI=3.141592653589793D0

	KERR=0
	NITER1=0
	FLMDA=0.D0

 1000	DO 1001 J=1,NCEF
           BB(J)=AA(J)
 1001	CONTINUE

	FCRD=0.D0
        RTPI=.5D0*DLOG(2.D0*PI)
 1010   DO 1011 N= 1,NDATA
           IF(YY(N) .LE. 1.D0) GOTO 1040
           Y=YY(N)
           IF (Y .LE. 50.D0) GOTO 1020
C===========================
C     Stirlingの公式
C===========================        
           FCRD=FCRD+RTPI+(Y+.5D0)*DLOG(Y)-Y+1.D0/(12.D0*Y)
           GOTO 1040

 1020      FCRD=FCRD+DLOG(Y)
           Y=Y-1.D0
           IF (Y .GT. 1.D0) GOTO 1020

 1040      FITC(N)=FNCTD(XX(N),AA,NCEF,KERR)*TT(N)*ET(N)
           IF (KERR .GE. 5) RETURN
c          WRITE(*,*) N,XX(N),YY(N),FITC(N)
 1011   CONTINUE

        IF (NCEFB .EQ. 0) GOTO 1080

 1060   DO 1061 N=1,NDATA
           FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN
           FITC(N)=FITC(N)+FUNCA*TT(N)*ZT(N)
 1061   CONTINUE

 1080   FLV1=FLVAL2(YY,FITC,NDATA,INFL)
        IF(INFL .EQ. 0) GOTO 1100
        KERR=3
        GOTO 3820

C===========================
C     外側のループ
C===========================

 1100   NITER2=0
        NITER1=NITER1+1

 1200   DO 1201 J=1,NCEF
           DD(J)=0.D0
 1220      DO 1221 K=1,NCEF
              RR(J,K)=0.D0
 1221      CONTINUE
 1201   CONTINUE

        do 1299 i=1,NCEF
           DAA(i) = 1.0d-6 * AA(i)
           if(DAA(i).eq.0.0) DAA(i)=1.0d-10
 1299   continue

 1300   DO 1301 N=1,NDATA

           CALL FNDRV(XX(N),AA,DAA,FA,NCEF,KERR)
c           write(*,*) FA
           IF (KERR .GE. 5)RETURN

           RPSI=1.D0/DSQRT(FITC(N))
           YD(N)=(YY(N)-FITC(N))*RPSI

 1320      DO 1321 J=1,NCEF
              QQ(N,J)=TT(N)*ET(N)*FA(J)*RPSI
              DD(J)=DD(J)+QQ(N,J)**2
 1321      CONTINUE
 1301   CONTINUE

 1340   DO 1341 J=1,NCEF
           DD(J)=1.D0/DSQRT(DD(J))
 1341   CONTINUE

 1360   DO 1361 N=1,NDATA
 1380      DO 1381 J=1,NCEF
              QQ(N,J)=QQ(N,J)*DD(J)
 1381      CONTINUE
 1361   CONTINUE

 1400   DO 1401 J=1,NCEF
 1420      DO 1421 K=J,NCEF

 1440         DO 1441 N=1,NDATA
                 RR(J,K)=RR(J,K)+QQ(N,J)*QQ(N,K)
 1441         CONTINUE
 1421      CONTINUE
           IF (J .EQ. NCEF) GOTO 1470

 1460      DO 1461 K=J+1,NCEF
              RR(K,J)=RR(J,K)
 1461      CONTINUE
 1470      CONTINUE
 1401   CONTINUE

        CALL MTXINV(RR,NCEF,DET)
        IF (DET .EQ. 0.D0) GOTO 3900
C=================================
C     λCの算出
C=================================
        TRACE=0.D0
 1480   DO 1481 J=1,NCEF
           TRACE=TRACE+RR(J,J)
 1481   CONTINUE

        FLMDC=1.D0/TRACE
C=================================
C     内側のループ
C=================================
 2000   DO 2001 N=1,NDATA
 2020      DO 2021 J=1,NCEF
              QX(N,J)=QQ(N,J)
 2021      CONTINUE
 2001   CONTINUE

 2040   DO 2041 J=1,NCEF
           DA(J)=0.D0
 2041   CONTINUE

        NITER2=NITER2+1
C=================================
C     修正マルカールの方法
C=================================

 2100   DO 2101 N=NDATA+1,NDATA+NCEF
           K=N-NDATA
           YD(N)=0.D0
 2120      DO 2121 J=1,NCEF
              IF(J .NE. K) GOTO 2140
              QX(N,J)=DSQRT(FLMDA)
              GOTO 2121
 2140         QX(N,J)=0.D0
 2121      CONTINUE
 2101   CONTINUE

C=================================
C     QR変換を呼び出す
C=================================

 2200   CALL QRTRN(QX,RR,JNX,NDATA+NCEF,NCEF,NRANK)
C=================================
C     後退代入法による計算
C=================================
 2300   DO 2301 J=NRANK,1,-1
           JJ=JNX(J)
           QV(JJ)=0.D0
 2320      DO 2321 N=1,NDATA
              QV(JJ)=QV(JJ)+QX(N,JJ)*YD(N)
 2321      CONTINUE
           ZV=QV(JJ)
           IF (J .EQ. NRANK) GOTO 2360
 2340      DO 2341 K=J+1,NRANK
              KK=JNX(K)
              ZV=ZV-RR(J,KK)*DA(KK)
 2341      CONTINUE

 2360      DA(JJ)=ZV/RR(J,JJ)

 2301   CONTINUE

        DLLN=0.D0
 2400   DO 2401 J=1,NRANK
           JJ=JNX(J)
           DLLN=DLLN+QV(JJ)**2
           DA(JJ)=DA(JJ)*DD(JJ)
           BB(JJ)=AA(JJ)+DA(JJ)
           write(*,*) JJ,BB(JJ)
           write(5,*) JJ,BB(JJ)
 2401   CONTINUE

 2440   DO 2441 N=1,NDATA
           FITC(N)=FNCTD(XX(N),BB,NCEF,KERR)*TT(N)*ET(N)
c           write(*,*) XX(N),FITC(N)
           IF (KERR .GE. 5) RETURN
 2441   CONTINUE
        IF (NCEFB .EQ. 0) GOTO 2500

 2460   DO 2461 N=1,NDATA
           FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN
           FITC(N)=FITC(N)+FUNCA*TT(N)*ZT(N)
 2461   CONTINUE

 2500   FLV2=FLVAL2(YY,FITC,NDATA,INFL)

        IF (INFL .EQ. 0) GOTO 2600
 2520   FLMDA=DMAX1(FLMDA*2.D0,FLMDC)
        IF (NITER2  .LE. 100) GOTO 2000
 2540   WRITE(7,7020)
        WRITE(5,7020)
 7020   FORMAT(1H0,"内側のループが100回を越えました")
        KERR=2
        RETURN

C===============================
C     係数の誤差の計算
C===============================

 2600   DO 2601 J=1,NCEF
           SA(J)=0.D0
 2620      DO 2621 K=1,NCEF
              DFDA(J,K)=0.D0
 2621      CONTINUE
           IF (NCEFB .EQ. 0) GOTO 2601
 2640      DO 2641 L=1,NCEFB
              RR(J,L)=0.D0
 2641      CONTINUE
 2601   CONTINUE

        do 2999 i=1,NCEF
           if(SA(i).eq.0.0) then
              DBB(i) = 1.0d-6 * BB(i)
c              write(*,*) ' aaaa '
              if(DBB(i).eq.0.0) DBB(i)=1.0d-10
              goto 2999
           endif
           DBB(i) = 1.0d-6 * SA(i)
 2999   continue

 3000   DO 3001 N=1,NDATA
           FUNC=FNCTD(XX(N),BB,NCEF,KERR)
           IF (KERR .GE. 5) RETURN
           FACT =TT(N)*(YY(N)-FITC(N))/FITC(N)
           ETF=TT(N)*ET(N)/FITC(N)
           ETX=(TT(N)/FITC(N))**2*ET(N)*YY(N)

           CALL FNDRV(XX(N),BB,DBB,FA,NCEF,KERR)
           IF (KERR .GE. 5) RETURN
           IF (MODE .EQ. 1) GOTO 3020
           CALL F2DRV(XX(N),BB,F2A,NCEF,KERR)
           IF (KERR .GE. 5) RETURN

 3020      DO 3021 J=1,NCEF
              DFDY(N,J)=FA(J)*ETF

 3040         DO 3041 K=J,NCEF
                 DFDA(J,K)=DFDA(J,K)-ETF**2*YY(N)*FA(J)*FA(K)
                 IF (MODE .EQ. 1) GOTO 3041
                 DFDA(J,K)=DFDA(J,K)+FACT*ET(N)*F2A(J,K)
 3041         CONTINUE

              IF (J .EQ. NCEF) GOTO 3021
 3060         DO 3061 K=J+1,NCEF
                 DFDA(K,J)=DFDA(J,K)
 3061         CONTINUE
 3021      CONTINUE

           IF (NETA .EQ. 0) GOTO 3100
 3080      DO 3081 J=1,NCEF
              QX(N,J)=(FACT-ETF**2*YY(N)*FUNC/ET(N))*FA(J)
 3081      CONTINUE

 3100      IF(NCEFB .EQ. 0) GOTO 3001

        do 3101 i=1,NCEF
           if(SA(i).eq.0.0) then
              DCC(i) = 1.0d-6 * CC(i)
c              write(*,*) ' bbbb '
              if(DCC(i).eq.0.0) DCC(i)=1.0d-10
              goto 3101
           endif
           DCC(i) = 1.0d-6 * SA(i)
 3101   continue

           CALL FNDRV(XX(N),CC,DCC,FC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN

 3120      DO 3121 J=1,NCEF
 3140         DO 3141 L=1,NCEFB
                 RR(J,L)=RR(J,L)-ETX*FA(J)*FC(L)*ZT(N)
 3141         CONTINUE

              IF (NETA .EQ. 1) GOTO 3121
              FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
              IF (KERR .GE. 5) RETURN
              QZ(N,J)=-ETX*FUNCA*FA(J)

 3121      CONTINUE

 3001   CONTINUE

C===========================================
C     行列Dの逆行列の計算
C===========================================

        CALL MTXINV(DFDA,NCEF,DET)
        IF (DET .EQ. 0.D0) GOTO 3900

 3200   DO 3201 J=1,NCEF

 3220      DO 3221 N=1,NDATA
              DADY=0.D0
 3240         DO 3241 K=1,NCEF
                 DADY=DADY-DFDA(J,K)*DFDY(N,K)
 3241         CONTINUE
              SA(J)=SA(J)+DADY**2*FITC(N)

              IF (NETA .EQ. 0) GOTO 3221
              DADE = 0.D0
 3260         DO 3261 K=1,NCEF
                 DADE=DADE-DFDA(J,K)*QX(N,K)
 3261         CONTINUE
              SA(J)=SA(J)+(DADE*SE(N))**2
 3221      CONTINUE

           IF (NCEFB .EQ. 0) GOTO 3400
 3300      DO 3301 L=1,NCEFB
              DADC=0.D0
 3320         DO 3321 K=1,NCEF
                 DADC=DADC-DFDA(J,K)*RR(K,L)
 3321         CONTINUE
              SA(J)=SA(J)+(DADC*SC(L))**2
 3301      CONTINUE

           IF(NETA .EQ. 1) GOTO 3400
 3340      DO 3341 N=1,NDATA
              DADZ=0.D0
 3360         DO 3361 K=1,NCEF
                 DADZ=DADZ-DFDA(J,K)*QZ(N,K)
 3361         CONTINUE
              SA(J)=SA(J)+(DADZ*SZ(N))**2
 3341      CONTINUE
 
 3400      SA(J)=DSQRT(SA(J))

 3201   CONTINUE

C========================================
C     収束の判定
C========================================

 3500   DLTR=FLV2-FLV1
        IF (DABS(DLTR) .GT. DABS(FLV2-FCRD)*EPSF) GOTO 3600

 3520   DO 3521 J=1,NCEF
           EPSA=DMAX1(SA(J)*1.D-6,DABS(BB(JJ))*EPSM)
           IF(DABS(DA(J)) .GT. EPSA) GOTO 3600
 3521   CONTINUE
        GOTO 3800

 3600   RAT=DLTR/DLLN
        IF (RAT .GE. 0.D0) GOTO 3700
        FU=DMAX1(2.D0,DMIN1(2.D0-RAT,1.D1))
        FLMDA=FLMDA*FU
        IF(FLMDA .EQ. 0.D0) FLMDA=FLMDC*FU*.5D0
        IF(NITER2 .LE. 100)THEN 
           GOTO 2000 
        ELSE 
           GOTO 2540 
        ENDIF

 3700   IF(RAT .LT. .75D0) GOTO 3720
        FLMDA=.5D0*FLMDA
        IF(FLMDA .LT. FLMDC) FLMDA=0.D0
        GOTO 3740
 3720   IF(RAT .GE. .25D0) GOTO 3740
        FLMDA=DMAX1(FLMDA*2.D0,FLMDC)

 3740   DO 3741 J=1,NCEF
           AA(J)=BB(J)
 3741   CONTINUE
        FLV1=FLV2
        IF(NITER1 .LE. 100) GOTO 1100
        WRITE (*,7040)
        WRITE (5,7040)
 7040   FORMAT (1H0,"外側のループが100回を越えました")
        KERR=2
        RETURN

 3800   DO 3801 J=1,NCEF
           AA(J)=BB(J)
 3801   CONTINUE
        IF(FLMDA .LE. 10.D0*FLMDC) RETURN

        WRITE(*,9100)
        WRITE(5,9100)
 9100   FORMAT(1H0,"λが規定値の10倍を越えたまま収束しました")
        KERR=1
        RETURN

 3820   WRITE(*,9120)
        WRITE(5,9120)
 9120   FORMAT(1H0,"初期の設定が不適当です")
        RETURN

 3900   WRITE(*,9140)
        WRITE(5,9140)
 9140   FORMAT(1H0,"行列式の値がゼロになりました")
        KERR=4
        RETURN
        END

C===========================================
C     ユウ度関数の対数の計算
C===========================================
        DOUBLE PRECISION FUNCTION FLVAL2(YY,FIT,NDATA,INFL)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        
        DIMENSION YY(500),FIT(500)

        INFL=0
        FLVAL2=0.D0

 4000   DO 4001 N=1,NDATA
c           WRITE(*,*) N,FIT(N)
           IF(FIT(N) .LE. 0.D0) GOTO 4020

           FLVAL2=FLVAL2+YY(N)*DLOG(FIT(N))-FIT(N)
 4001   CONTINUE
        write(*,*) FLVAL2
        write(5,*) FLVAL2
        RETURN

 4020   WRITE(*,9000)
        WRITE(5,9000)
 9000   FORMAT (1H0,"関数値がゼロまたは負の値です")

        INFL=1

        FLVAL2=0.
        RETURN
        END
C================================================
C    Gauss Jordan法による逆行列の計算
C================================================
        SUBROUTINE MTXINV(ARRAY,NC,DET)

        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION ARRAY(30,30),JA(30),KA(30)

        DET=1.D0

 2000   DO 2001 L=1,NC

        AMAX=0.D0

 2020   DO 2021 J=L,NC
 2040      DO 2041 K=L,NC
              IF(DABS(AMAX)-DABS(ARRAY(J,K)) .GT. 0.D0) GOTO 2041
              AMAX=ARRAY(J,K)
              JA(L)=J
              KA(L)=K
 2041      CONTINUE
 2021   CONTINUE

        IF (AMAX .NE. 0.D0) GOTO 2100
        DET=0.D0
        RETURN

 2100   J=JA(L)
        IF(J .EQ. L) GOTO 2200
 2120   DO 2121 K=1,NC
           SAVE=ARRAY(L,K)
           ARRAY(L,K)=ARRAY(J,K)
           ARRAY(J,K)=-SAVE
 2121   CONTINUE

 2200   K=KA(L)
        IF(K .EQ. L) GOTO 2300
 2220   DO 2221 J=1,NC
           SAVE=ARRAY(J,L)
           ARRAY(J,L)=ARRAY(J,K)
           ARRAY(J,K)=-SAVE
 2221   CONTINUE

 2300   DO 2301 J=1,NC
           IF (J .EQ. L) GOTO 2301
           ARRAY(J,L)=-ARRAY(J,L)/AMAX
 2301   CONTINUE

 2320   DO 2321 J=1,NC
           IF(J .EQ. L) GOTO 2321
           DO 2341 K=1,NC
              IF (K .EQ. L) GOTO 2341
              ARRAY(J,K)=ARRAY(J,K)+ARRAY(J,L)*ARRAY(L,K)
 2341      CONTINUE
 2321   CONTINUE

 2360   DO 2361 K=1,NC
           IF(K .EQ. L) GOTO 2361
           ARRAY(L,K)=ARRAY(L,K)/AMAX
 2361   CONTINUE
        
        ARRAY(L,L)=1.D0/AMAX
        DET=DET*AMAX

 2001   CONTINUE

 2500   DO 2501 L=1,NC
           M=NC-L+1
           K=JA(M)
           IF (K .LE. M) GOTO 2540
 2520      DO 2521 J=1,NC
              SAVE=ARRAY(J,M)
              ARRAY(J,M)=-ARRAY(J,K)
              ARRAY(J,K)=SAVE
 2521      CONTINUE
           

 2540      J=KA(M)
           IF (J .LE. M) GOTO 2501
 2560      DO 2561 K=1,NC
              SAVE=ARRAY(M,K)
              ARRAY(M,K)=-ARRAY(J,K)
              ARRAY(J,K)=SAVE
 2561      CONTINUE
 2501   CONTINUE
        
        RETURN
        END

C============================================
C     修正グラムシュミット法による
C               行列のQR変換
C============================================

        SUBROUTINE QRTRN(AQ,RR,JNX,NDATA,NCEF,NRANK)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION AQ(510,30),RR(30,30),JNX(30),RN(30)

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
