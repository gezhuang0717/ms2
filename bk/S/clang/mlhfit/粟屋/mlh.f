C========================================================
C     Maximum Likelihood法による
C             fitting program
C     (粟屋 隆氏著  学会出版センター刊
C      ’データ解析’の中に、このプログラムはあります。）
C========================================================

C********************************************************
C     ＋＋このプログラムの使い方＋＋
C      サブルーチン'CFDQR'を呼び出すだけで、Fitting の結果
C     と、その誤差が返されます。
C     あらかじめ、XX,YY,ET,TT,SE,NDATA,NCEF,NETA,MODE
C     と、合わせる関数と、その一回微分、二回微分を
C     を入力しておき、'CALL CFDQR'で、
C     AA(N)にFitting 結果(Nは整数）
C     SA(N)にその誤差   が返されます。
C     
C     今のところ、Fitできる関数は
C     べき級数、EXP関数、GAUSSIANの三種ですが、簡単にFitting関数は
C     定義できます。
C     これは、実際のプログラムを見て下さい。
C********************************************************     

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      DIMENSION XX(46),YY(46),ET(46),TT(46),SE(46),ZT(46)
      DIMENSION IDAT(46),SZ(46),AA(3),SA(3),CC(3),SC(3),NY(46),A1(3)

      EXTERNAL FEXPN,DREXP,D2EXP
     
C      DATA NY/17,15,11,14,7,13,9,7,6,5,
C     >5,4,9,3,3,1,2,2,4,2,
C     >2,5,3,2,2,2,0,0,0,2,
C     >2,3,2,4,0,1,1,1,2,3,
C     >2,0,0,1,0,2,0,1,0,1/
 
C      DATA NY /40,29,21,22,21,24,12,15,11,17,18,
C     >18,8,17,10,15,14,14,13,4,10,14,13,15,17,10,9,
C     >9,7,12,9,14,12,14,15,12,5,9,7,10,13,8,13,10,10,
C     >8,10,13,11,17,16,9,9,7,10,13,17,11,7,11,15,7,9,
C     >10,9,8,7,15,11,9,10,18,16,7,8/

C      DATA NY /22,31,25,18,19,13,14,19,16,14,13,12,8,12,9,
C     >11,6,7,6,5,8,8,13,11,4,13,10,9,6,17,10,8,9,6,7,8,6,
C     >9,5,7,5,10,9,4,5,7,3,5,7,4,6,6,6,3,8,3,10,7,4,6,6,4,6,
C     >10,6,5,11,7,8,6,6,7,9,4,9,5/ 

C========================================================
C     ファイルの読み込み（ここではDATAは整数型）
C========================================================

      OPEN (1,FILE='TestData',STATUS='UNKNOWN')
      DO 10 I=1,46
         READ(1,*) IDAT(I)
 10   CONTINUE

C********************************************************
C     初期設定
C********************************************************

C========================================================
C     フィッティング係数の初期値
C              （絶対値が大きいと収束しない！）
C========================================================    
 
      DATA A1/1,-0.3,15/

C========================================================
C     各チャンネルの計測時間
C               （ここでは固定）
C========================================================

      DATA TT/46*1.D0/

C========================================================
C     各チャンネルに置ける計測効率
C               （ここでは固定）
C========================================================

      DATA ET/46*1.D0/

C========================================================
C     係数効率の誤差
C               （ここではすべてゼロとした）
C========================================================

      DATA SE/46*0.D0/

C========================================================
C     バックグラウンドのチャンネルごとの係数効率
C           およびその誤差
C           （バックグラウンドの解析をしない場合は
C           ダミーパラメーターであるが、他の変数（TT,ET等）
C           と形を合わせなければならない。）
C========================================================

      DATA ZT/46*1.D0/
      DATA SZ/46*0.D0/

C      OPEN (*,FILE='PRN2')

C========================================================
C     整数型で読み込んだDATAを
C             倍精度型に変換
C               （文字変数の最初の文字に注意）
C========================================================

 1000 DO 1001 N = 1,46
         XX(N)=N
C        YY(N)=NY(N)
         YY(N)=IDAT(N)
 1001 CONTINUE

C========================================================
C     DATAの個数を NDATAに代入
C========================================================

      NDATA = 46

C========================================================
C     未定変数の個数をNCEFに代入
C       ここでは
C        f(x)=A1*exp(A2*t)+A3  で三個
C    （あわせる関数の形はCALL CFDQR文の中で定義（後述））
C========================================================

      NCEF = 3

C========================================================
C     係数効率を考えるかどうか
C     NETA = 0  <- 考えない。（SE等、係数誤差のパラメータは
C                                  ダミーとなる）
C          = 1  <- RIの効率を考える。
C          = 2  <- バックグラウンドの効率まで考える。
C     （詳しいことは’データ解析’を参照）
C========================================================

      NETA = 0

C========================================================
C     エラーフラグ
C========================================================

      KERR = 0

 1100 DO 1101 MODE = 1,2

 1120    DO 1121 J = 1,3
            AA(J)=A1(J)
 1121    CONTINUE

         IF (MODE .EQ. 1) THEN
            WRITE(*,7010)
 7010       FORMAT (1H0,'＊＊ 一回微分係数のみ ＊＊' /)
         ELSE
            WRITE(*,7020)
 7020       FORMAT (1H0,'＊＊ 二回微分係数も入れる ＊＊'/)
         ENDIF

         CALL CFDQR2 (XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
     >        FEXPN,DREXP,D2EXP,ZT,SZ,CC,SC,0,FEXPN,MODE,KERR)

 1140    DO 1141 J=1,3
            WRITE (*,7000) J,AA(J),SA(J)
 7000       FORMAT (1H,'AA(',I1,')=',D12.5,'±',D12.5)
 1141    CONTINUE
   
 1101 CONTINUE
  
      STOP
      END
      
C================================================
C     Maximum Likelihood法によるFitting Programm
C================================================

        SUBROUTINE CFDQR2(XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,FNCTD,
     >                    FNDRV,F2DRV,ZT,SZ,CC,SC,NCEFB,FNCTB,MODE,KERR)
C        SUBROUTINE CFDQR2(XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,FLINE,DLINE,D2LIN,ZT,SZ,CC,SC,NCEFB,FNCTB,MODE,KERR)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XX(500),YY(500),ET(500),TT(500),FITC(500),SE(500)
        DIMENSION ZT(500),SZ(500),YD(510),AA(10),DA(10),SA(10),FA(10)
        DIMENSION CC(10),SC(10),FC(10),QQ(500,10),QX(510,10),DD(10)
        DIMENSION RR(10,10),JNX(10),F2A(10,10),DFDY(500,10),DFDA(10,10)
        DIMENSION BB(10)
        DIMENSION QV(10)
        DIMENSION QZ(500,10)
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
           WRITE(*,*) FITC(N)
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

 1300   DO 1301 N=1,NDATA

           CALL FNDRV(XX(N),AA,FA,NCEF,KERR)
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
 2401   CONTINUE

 2440   DO 2441 N=1,NDATA
           FITC(N)=FNCTD(XX(N),BB,NCEF,KERR)*TT(N)*ET(N)
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

 3000   DO 3001 N=1,NDATA
           FUNC=FNCTD(XX(N),BB,NCEF,KERR)
           IF (KERR .GE. 5) RETURN
           FACT =TT(N)*(YY(N)-FITC(N))/FITC(N)
           ETF=TT(N)*ET(N)/FITC(N)
           ETX=(TT(N)/FITC(N))**2*ET(N)*YY(N)

           CALL FNDRV(XX(N),BB,FA,NCEF,KERR)
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
           CALL FNDRV(XX(N),CC,FC,NCEFB,KERR)
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
 7040   FORMAT (1H0,"外側のループが100回を越えました")
        KERR=2
        RETURN

 3800   DO 3801 J=1,NCEF
           AA(J)=BB(J)
 3801   CONTINUE
        IF(FLMDA .LE. 10.D0*FLMDC) RETURN

        WRITE(*,9100)
 9100   FORMAT(1H0,"λが規定値の10倍を越えたまま収束しました")
        KERR=1
        RETURN

 3820   WRITE(*,9120)
 9120   FORMAT(1H0,"初期の設定が不適当です")
        RETURN

 3900   WRITE(*,9140)
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

           IF(FIT(N) .LE. 0.D0) GOTO 4020

           FLVAL2=FLVAL2+YY(N)*DLOG(FIT(N))-FIT(N)
 4001   CONTINUE 
        RETURN

 4020   WRITE(7,9000)
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
        DIMENSION ARRAY(10,10),JA(10),KA(10)

        DET=1.D0

 2000   DO 2001 L=1,NC

        AMAX=0.D0

 2020   DO 2021 J=L,NC
 2040      DO 2041 K=L,NC
              IF(DABS(AMAX)-DABS(ARRAY(J,K)) .GT. 0.D0) GOTO 2041
              AMAX=ARRAY(J,K)
              JA(L)=J
             