C========================================================It's OK ?
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
C********************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XX(2000),YY(2000),ET(2000),TT(2000),SE(2000),ZT(2000)
      DIMENSION SZ(2000),AA(30),SA(30),CC(30),SC(30)
      DIMENSION DAT(2000)
      dimension AAA(30),SAA(30)
c      dimension aai(30),II(30)
      CHARACTER DATFN*50,DATFNC*30,ANS,fcond*30,fres*30,
     &          ftblt*30,ftbldp*30

      common/PARAM/aai(30),II(30),ncoef

      EXTERNAL func,dfunc,d2func

C========================================================
C     ファイルの読み込み（ここではDATAは整数型）
C========================================================

      fcond = 'fcond.dat'
      OPEN (1,FILE=fcond,STATUS='UNKNOWN')! open initial condition file

      READ(1,'(A30)') fres                    ! read result file name
c      OPEN (5,FILE=fres,STATUS='UNKNOWN')     ! open result file
      OPEN (5,FILE=fres,STATUS='NEW')     ! open result file
      read(1,'(a30)') datfnc
      write(*,*) 'output file : ',fres,datfnc
      write(5,*) 'output file : ',fres,datfnc

c-- read data ---------------------------------------------------------
      READ(1,'(A30)') DATFN                       ! read data file name
      write(*,'(a,a)') 'Data file name : ',DATFN
      write(5,'(a,a)') 'Data file name : ',DATFN
      OPEN (2,FILE=DATFN,STATUS='UNKNOWN')        ! open data file
      i = 1                                       ! read data file
 1    read(2,*,end=10) dat(i)
      write(*,*) i,dat(i)
      i = i + 1
      if(dat(i).le.0.0) dat(i)=0.0
      goto 1
 10   CONTINUE
      ndat = i - 1                          ! number of points of data
      write(*,*) 'Number of points of data = ',ndat
      write(5,*) 'Number of points of data = ',ndat
      close(2)
c---------------------------------------------------------------------
c********************************************************
C     初期設定
c********************************************************
      read(1,*) from,to                          ! read range of x axis
      CMP = (TO - FROM)/(DFLOAT(ndat)-1.d0)
      read(1,*) ncoef,ncef
      write(*,*) 'Number of parameters = ',ncoef
      write(5,*) 'Number of parameters = ',ncoef
      write(*,*) 'Number of free parameters = ',ncef
      write(5,*) 'Number of free parameters = ',ncef
      read(1,*) il,ih                    ! read fitting region in ch
      write(*,*) 'Fitting region : ',il,' -> ',ih
      write(5,*) 'Fitting region : ',il,' -> ',ih
      read(1,*) (II(i),i=1,ncoef)    ! read order of parameters
      read(1,*) (AAI(i),i=1,ncoef)   ! read initial value of parameters
      write(*,*) 'Initial value of free parameters'
      write(5,*) 'Initial value of free parameters'
      do 31 i = 1,ncef
         write(*,1000) II(i),AAI(II(i))
         write(5,1000) II(i),AAI(II(i))
 31   continue
      if(ncef.eq.ncoef) goto 34
      write(*,*) 'Initial value of fixed parameters'
      write(5,*) 'Initial value of fixed parameters'
      do 32 i = ncef+1,ncoef
         write(*,1000) II(i),AAI(II(i))
         write(5,1000) II(i),AAI(II(i))
 32   continue
 1000 FORMAT (2H  ,'AAI(',I2,') = ',D17.10)
 34   DO 33 I = 1,NCOEF
         AA(I) = AAI(II(I))
 33   CONTINUE

      DATA TT/2000*1.D0/
      DATA ET/2000*1.D0/
      DATA SE/2000*0.D0/
      DATA ZT/2000*1.D0/
      DATA SZ/2000*0.D0/

c========================================================
C     DATAを XX,YYに代入
c========================================================
      NDATA = IH - IL + 1
      DO 22 N = 1,NDATA
         XX(N) = DFLOAT(IL+N-1)
         XX(N) = FROM + CMP * (XX(N)-1.)
         YY(N) = DAT(IL+N-1)
 22   CONTINUE

C========================================================
C     計測効率を考えるかどうか
C     NETA = 0  <- 考えない。（SE等、係数誤差のパラメータは
C                                  ダミーとなる）
C          = 1  <- RIの効率を考える。
C          = 2  <- バックグラウンドの効率まで考える。
C     （詳しいことは’データ解析’を参照）
C========================================================
      NETA = 0
c========================================================
C     エラーフラグ
c========================================================
      KERR = 0

c 1100 DO 1101 MODE = 1,2
         MODE = 2

         IF (MODE .EQ. 1) THEN
            WRITE(*,7010)
            WRITE(5,7010)
 7010       FORMAT (1H ,'＊＊ 一回微分係数のみ ＊＊' /)
         ELSE
            WRITE(*,7020)
            WRITE(5,7020)
 7020       FORMAT (1H ,'＊＊ 二回微分係数も入れる ＊＊'/)
         ENDIF

         CALL CFDQR2 (XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
     >         func,dfunc,d2func,ZT,SZ,CC,SC,0,func,MODE,KERR)

         WRITE(*,*) 'Fitting region(ch) : ',IL,' --> ',IH
         WRITE(5,*) 'Fitting region(ch) : ',IL,' --> ',IH
         WRITE(*,*) 'Fitting region (arb.) : ',XX(1),' --> ',XX(NDATA)
         WRITE(5,*) 'Fitting region (arb.) : ',XX(1),' --> ',XX(NDATA)
         WRITE(*,*) 'Free parameters '
         WRITE(5,*) 'Free parameters '
 1140    DO 1141 J = 1,NCEF
            WRITE (*,7000) II(J),AA(J),SA(J)
            WRITE (5,7000) II(J),AA(J),SA(J)
 7000       FORMAT (2H  ,'AA(',I2,') = ',D17.10,' +- ',D17.10)
 1141    CONTINUE
         IF(NCEF.GE.NCOEF) GOTO 445
         WRITE(*,*) 'Fixed parameters '
         WRITE(5,*) 'Fixed parameters '
         DO 1142 J = NCEF+1,NCOEF
            WRITE (*,7001) II(J),AA(J)
            WRITE (5,7001) II(J),AA(J)
 7001       FORMAT (2H  ,'AA(',I2,') = ',D17.10)
 1142    CONTINUE

 445     CONTINUE

         do 100 i=1,NCOEF
               AAA(II(i)) = AA(i)
               SAA(II(i)) = SA(i)
 100     continue

         chisq = 0.0
         do 443 i=1,ndata
            x = xx(i)
            fff = func(x,aa,ncef,kerr)
            if(yy(i).le.0.) then
               w = 1.d-10
            else
c               w = 1./yy(i)
               w = 1./fff
            endif
            chisq = chisq + (yy(i)-fff)**2*w
 443     continue
         write(*,*) 'chisq = ',chisq
         write(5,*) 'chisq = ',chisq
         chisq = chisq/(dfloat(ndata-ncef))
         write(*,*) 'reduced chisq = ',chisq
         write(5,*) 'reduced chisq = ',chisq

 444     CONTINUE
cc       WRITE(*,*) ' Save fitting curve ? (y or n) '
cc       READ(*,'(A1)') ANS
cc       IF(ANS.EQ.'y') THEN
cc          WRITE(*,*) ' Input file name '
cc          READ(*,'(A30)') DATFNC
c         OPEN (3,FILE=DATFNC,STATUS='UNKNOWN')
         OPEN (3,FILE=DATFNC,STATUS='NEW')

            DO 1143 I = 1,ndat
               X = DFLOAT(I)
               X = FROM + CMP * (X - 1.)
               fff = func(x,aa,ncef,kerr)
cc             write(*,*) x,fff
               WRITE(3,*) fff
 1143       CONTINUE

cc       ELSEIF(ANS.EQ.'n') THEN
cc          GOTO 1101
cc       ELSE
cc          GOTO 444
cc       ENDIF

c 2222    DO 3333 J = 1,NCEF
c            WRITE (*,*) AA(J)
c 3333    CONTINUE

 1101 CONTINUE

 999  stop

      END
C================================================
C     Maximum Likelihood法によるFitting Programm
C================================================

        SUBROUTINE CFDQR2(XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
     >              FNCTD,FNDRV,F2DRV,ZT,SZ,CC,SC,NCEFB,FNCTB,MODE,KERR)

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XX(2000),YY(2000),ET(2000),TT(2000),FITC(2000),
     >            SE(2000)
        DIMENSION ZT(2000),SZ(2000),YD(2030),AA(30),DA(30),SA(30),FA(30)
        DIMENSION BB(30),CC(30),SC(30),FC(30),QQ(2000,30),QX(2030,30)
        DIMENSION DD(30),QV(30),RR(30,30),JNX(30),F2A(30,30)
        DIMENSION DFDY(2000,30),DFDA(30,30),QZ(2000,30)

C     このProgramは”データ解析”の本に載っているものです。

C     XX        横軸のデータ
C     YY        縦軸のデータ
C     ET        計数効率
C     SE        計数効率の誤差
C     TT        測定時間
C     NDATA     データの個数
C
C     FITC      合わせる関数
C
C     FNCTD     合わせる関数の計算
C     FNDRV     係数による関数の１階微分
C     F2DRV     係数による関数の２階微分
C     FNCTB     バックグラウンドの関数の計算
C
C     AA        合わせる関数の計算
C     SA        係数の誤差
C     NCEF      係数の個数
C
C     FA        合わせる関数の計数による１階微分
C     FA2       計数による２階微分
C
C     ZT        バックグラウンドの計数効率
C     SZ        同上の誤差
C     CC        バックグラウンド関数の計数
C     SC        同上の誤差
C     NCEFB     計数の誤差
C
C     NETA      0   計数効率は考えない
C               1   RIの計数効率を考える
C               2   バックグラウンドの計数効率まで考える
C
C     MODE      1   係数の誤差の計算の際１階微分まで考える
C               2   係数の誤差の計算の際２階微分まで考える
C
C     KERR      エラーメッセージ
C                1   λが大きいまま収束（警告のみ）
C                2   繰り返し数オーバー（警告のみ）
C                3   初期値不良
C                4   ベクトルが独立でない
C
C     必要なサブプログラム  FLVAL2,QRTRN,MTXINV
C     必要なサブプログラム（仮の名前）  FNCTD,FNDRV,F2DRV,FNCTB
C
c      write(*,*)'NDATA,NEAT,AA,SA,NCEF = ',NDATA,NETA,AA,SA,NCEF
c      write(*,*)'MODE,KERR = ',MODE,KERR
        EPSM=1.D-10
        EPSF=1.D-6
        PI=3.141592653589793D0

        KERR=0
        NITER1=0
        FLMDA=0.D0

 1000   DO 1001 J=1,NCEF
           BB(J)=AA(J)
 1001   CONTINUE

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

 1040      call coef(aa,ncef,aa)
           FITC(N)=FNCTD(XX(N),AA,NCEF,KERR)*TT(N)*ET(N)
           IF (KERR .GE. 5) RETURN
c          WRITE(*,*) N,XX(N),YY(N),FITC(N)
 1011   CONTINUE

        IF (NCEFB .EQ. 0) GOTO 1080

 1060   DO 1061 N=1,NDATA
           call coef(cc,ncefb,cc)
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

           call coef(aa,ncef,aa)
           CALL FNDRV(XX(N),AA,FA,NCEF,KERR)
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
c      write(*,*)' NRANK = ',nrank
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
c           write(5,*) JJ,BB(JJ)
C 2402       FORMAT (10H          ,'AA(',I2,') = ',D17.10)
 2401   CONTINUE

 2440   DO 2441 N=1,NDATA
           call coef(bb,ncef,bb)
           FITC(N)=FNCTD(XX(N),BB,NCEF,KERR)*TT(N)*ET(N)
c           write(*,*) XX(N),FITC(N)
           IF (KERR .GE. 5) RETURN
 2441   CONTINUE

        IF (NCEFB .EQ. 0) GOTO 2500

 2460   DO 2461 N=1,NDATA
           call coef(cc,ncefb,cc)
           FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN
           FITC(N)=FITC(N)+FUNCA*TT(N)*ZT(N)
 2461   CONTINUE

 2500   FLV2=FLVAL2(YY,FITC,NDATA,INFL)

        IF (INFL .EQ. 0) GOTO 2600
 2520   FLMDA=DMAX1(FLMDA*2.D0,FLMDC)
        IF (NITER2  .LE. 100) GOTO 2000
 2540   WRITE(*,7020)
        WRITE(5,7020)
 7020   FORMAT(1H ,"内側のループが100回を越えました")
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
           call coef(bb,ncef,bb)
           FUNC=FNCTD(XX(N),BB,NCEF,KERR)
           IF (KERR .GE. 5) RETURN
           FACT =TT(N)*(YY(N)-FITC(N))/FITC(N)
           ETF=TT(N)*ET(N)/FITC(N)
           ETX=(TT(N)/FITC(N))**2*ET(N)*YY(N)

           call coef(bb,ncef,bb)
           CALL FNDRV(XX(N),BB,FA,NCEF,KERR)
           IF (KERR .GE. 5) RETURN
           IF (MODE .EQ. 1) GOTO 3020
           call coef(bb,ncef,bb)
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
           call coef(cc,ncefb,cc)
           CALL FNDRV(XX(N),CC,FC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN

 3120      DO 3121 J=1,NCEF
 3140         DO 3141 L=1,NCEFB
                 RR(J,L)=RR(J,L)-ETX*FA(J)*FC(L)*ZT(N)
 3141         CONTINUE

              IF (NETA .EQ. 1) GOTO 3121
              call coef(cc,ncefb,cc)
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
 7040   FORMAT (1H ,"外側のループが100回を越えました")
        KERR=2
        RETURN

 3800   DO 3801 J=1,NCEF
           AA(J)=BB(J)
 3801   CONTINUE
        IF(FLMDA .LE. 10.D0*FLMDC) RETURN

        WRITE(*,9100)
        WRITE(5,9100)
 9100   FORMAT(1H ,"λが規定値の10倍を越えたまま収束しました")
        KERR=1
        RETURN
C
C   KERR = 3
C
 3820   WRITE(*,9120)
        WRITE(5,9120)
 9120   FORMAT(1H ,"初期の設定が不適当です")
        RETURN
C
C   KERR = 4
C
 3900   WRITE(*,9140)
        WRITE(5,9140)
 9140   FORMAT(1H ,"行列式の値がゼロになりました")
        KERR=4
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
C===========================================
C     ユウ度関数の対数の計算
C===========================================
        DOUBLE PRECISION FUNCTION FLVAL2(YY,FIT,NDATA,INFL)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)

        DIMENSION YY(2000),FIT(2000)

        INFL=0
        FLVAL2=0.D0

 4000   DO 4001 N=1,NDATA
c           WRITE(*,*) N,FIT(N)
           IF(FIT(N) .LE. 0.D0) GOTO 4020

           FLVAL2=FLVAL2+YY(N)*DLOG(FIT(N))-FIT(N)
 4001   CONTINUE
        write(*,*) '   FLVAL2 = ',FLVAL2
c        write(5,*) '   FLVAL2 = ',FLVAL2
        RETURN

 4020   WRITE(*,9000)
        WRITE(5,9000)
 9000   FORMAT (1H ,"関数値がゼロまたは負の値です")

        INFL=1

        FLVAL2=0.
        RETURN
        END
c======================================
      subroutine coef(cfi,ncef,cf)
c======================================
      
      implicit double precision(a-h,o-z)
      dimension cf(30),cfi(30)
      common/PARAM/aai(30),II(30),ncoef

      do 1 I = 1,ncef
         cf(II(I)) = cfi(I)
 1    continue
      if(ncoef.gt.ncef) then
         do 101 I = ncef+1,ncoef
            cf(II(I)) = aai(II(I))
 101     continue
      endif

      return
      end

************************************************************************
      double precision function func(x,aa,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30)

      kerr = 0

      arg = - 0.5d0*((x-aa(3))**2.d0)*((aa(2))**-2.d0)
      func = aa(1)*dexp(arg)
      func=func+aa(4)
      func=func+aa(5)*x

      return
      end
************************************************************************
      subroutine dfunc(x,aa,fa,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30),fa(30)

      kerr = 0

      s=x-aa(3)
      s2=s*s
      arg = - 0.5d0*((x-aa(3))**2.d0)*((aa(2))**-2.d0)
      fnc=dexp(arg)
      fa(1) = fnc
      fa(2) = aa(1)*s2*((aa(2))**-3.d0)*fnc
      fa(3) = aa(1)*((aa(2))**-2.d0)*s*fnc
      fa(4)=1.d0
      fa(5)=x

      return
      end
************************************************************************
      subroutine d2func(x,aa,f2a,ncef,kerr)
************************************************************************
      implicit double precision(a-h,o-z)
      dimension aa(30),f2a(30,30)

      kerr = 0

      s=x-aa(3)
      s2=s**2.d0
      arg=-0.5d0*((x-aa(3))**2.d0)*((aa(2))**-2.d0)
      fnc=dexp(arg)

      f2a(1,1) = 0.d0
      f2a(1,2) = s2*((aa(2))**-3.d0)*fnc
      f2a(1,3) = ((aa(2))**-2.d0)*s*fnc

      f2a(2,1) = f2a(1,2)
      f2a(2,2) = (-3+s2*((aa(2))*-2.d0))
     >           *aa(1)*s2*((aa(2))**-4.d0)*fnc
      f2a(2,3) = (-2+s2*((aa(2))*-2.d0))
     >           *aa(1)*s*((aa(2))**-3.d0)*fnc
      f2a(3,1) = f2a(1,3)
      f2a(3,2) = f2a(2,3)
      f2a(3,3) = (-1+s2*((aa(2))**-2.d0))
     >           *aa(1)*((aa(2))**-2.d0)*fnc

 1000 do 1001 j=1,5
         f2a(j,4) = 0.d0
         f2a(4,j) = 0.d0
         f2a(j,5) = 0.d0
         f2a(5,j) = 0.d0
 1001 continue

      return
      end
