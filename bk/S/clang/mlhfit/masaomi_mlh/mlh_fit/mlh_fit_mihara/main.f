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
C********************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      DIMENSION XX(500),YY(500),ET(500),TT(500),SE(500),ZT(500)
      DIMENSION IDAT(500),SZ(500),AA(30),SA(30),CC(30),SC(30),NY(500)
      DIMENSION DAT(500)
      dimension AAA(30),SAA(30)
      CHARACTER DATFN*50,DATFNC*30,ANS,fcond*30,fres*30
      COMMON/user/pi,DX,fr
      COMMON/PARAM/AAI(30),II(30),NCOEF,ncmp,ic

c     EXTERNAL DBLG,DDBLG
c     EXTERNAL FEXPN,DREXP,D2EXP
      EXTERNAL fmltig,dmltig,d2mltig

      pi = dacos(-1.d0)

C========================================================
C     ファイルの読み込み（ここではDATAは整数型）
C========================================================

      fcond = 'fcond.dat'
      OPEN (1,FILE=fcond,STATUS='UNKNOWN')             ! open initial condition file

      READ(1,'(A30)') fres                             ! read result file name
      OPEN (5,FILE=fres,STATUS='new')                  ! open result file

      READ(1,'(A30)') DATFN                            ! read data file name
      write(*,'(a,a)') 'Data file name : ',DATFN
      write(5,'(a,a)') 'Data file name : ',DATFN
      read(1,*) ndat                                   ! number of points of data
      write(*,*) 'Number of points of data = ',ndat
      write(5,*) 'Number of points of data = ',ndat
      OPEN (2,FILE=DATFN,STATUS='UNKNOWN')             ! open data file
      DO 10 I = 1,ndat                                 ! read data file
         READ(2,*) DAT(I)
         if(dat(i).le.0.0) dat(i)=0.0
c        write(*,*) dat(i)
 10   CONTINUE
      close(2)
C********************************************************
C     初期設定
C********************************************************
      read(1,*) from,to                                ! read range of x axis
      CMP = (TO - FROM)/DFLOAT(ndat)
      DX = CMP
      read(1,*) ncmp,ncef,ic                           ! read number of coefficients
      read(1,*) fr                                     ! read resolution correction factor
      ncoef = 3*ncmp
      write(*,*) 'Number of parameters = ',ncoef
      write(5,*) 'Number of parameters = ',ncoef
      write(*,*) 'Number of free parameters = ',ncef
      write(5,*) 'Number of free parameters = ',ncef
      write(*,111) ic
      write(5,111) ic
 111  format (i2,' [0 : independent center, 1 : common center]')
      write(*,*) 'Transverse resolution correction factor = ',fr
      write(5,*) 'Transverse resolution correction factor = ',fr
      read(1,*) il,ih                                 ! read fitting region in ch
      write(*,*) 'Number of Gaussian components = ',ncmp
      write(5,*) 'Number of Gaussian components = ',ncmp
      write(*,*) 'Fitting region : ',il,' -> ',ih
      read(1,*) (II(i),i=1,ncoef)                      ! read order of parameters
      read(1,*) (AAI(i),i=1,ncoef)                     ! read initial value of parameters
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

      DATA TT/500*1.D0/
      DATA ET/500*1.D0/
      DATA SE/500*0.D0/
      DATA ZT/500*1.D0/
      DATA SZ/500*0.D0/

C========================================================
C     DATAを XX,YYに代入
C========================================================
      NDATA = IH - IL + 1
      DO 22 N = 1,NDATA
         XX(N) = DFLOAT(IL+N-1)
         XX(N) = FROM + CMP/2. + CMP * (XX(N)-1.)
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
C========================================================
C     エラーフラグ
C========================================================
      KERR = 0

 1100 DO 1101 MODE = 1,1

         IF (MODE .EQ. 1) THEN
            WRITE(*,7010)
            WRITE(5,7010)
 7010       FORMAT (1H0,'＊＊ 一回微分係数のみ ＊＊' /)
         ELSE
            WRITE(*,7020)
            WRITE(5,7020)
 7020       FORMAT (1H0,'＊＊ 二回微分係数も入れる ＊＊'/)
         ENDIF

c        CALL CFDQR2 (XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
c    >           DBLG,DDBLG,D2DBLG,ZT,SZ,CC,SC,0,DBLG,MODE,KERR)
c        CALL CFDQR2 (XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
c    >           FEXPN,DREXP,D2EXP,ZT,SZ,CC,SC,0,FEXPN,MODE,KERR)
         CALL CFDQR2 (XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
     >           fmltig,dmltig,d2mltig,ZT,SZ,CC,SC,0,fmltig,MODE,KERR)

         WRITE(*,*) 'Fitting region(ch) : ',IL,' --> ',IH
         WRITE(5,*) 'Fitting region(ch) : ',IL,' --> ',IH
         WRITE(*,*) 'Fitting region (MeV/c) : ',XX(1),' --> ',XX(NDATA)
         WRITE(5,*) 'Fitting region (MeV/c) : ',XX(1),' --> ',XX(NDATA)
         WRITE(*,*) 'Free parameters '
         WRITE(5,*) 'Free parameters '
 1140    DO 1141 J = 1,NCEF
            WRITE (*,7000) II(J),AA(J),SA(J)
            WRITE (5,7000) II(J),AA(J),SA(J)
 7000       FORMAT (2H  ,'AA(',I2,') = ',D17.10,' ± ',D17.10)
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
         do 101 i = 1,ncmp
            write(*,*) 'Component : ',i
            write(5,*) 'Component : ',i
            tn = AAA(3*(i-1)+1)
            dtn= SAA(3*(i-1)+1)
            sg = AAA(3*(i-1)+2)
            ds = SAA(3*(i-1)+2)
            ww = sg*2.d0*(2.d0*dlog(2.d0))**.5
            dw = ds*2.d0*(2.d0*dlog(2.d0))**.5
            ce = AAA(3*(i-1)+3)
            dc = SAA(3*(i-1)+3)
            write(*,*) ' Total number = ',tn,' +- ',dtn
            write(5,*) ' Total number = ',tn,' +- ',dtn
            write(*,*) ' Sigma(MeV/c) = ',sg,' +- ',ds
            write(5,*) ' Sigma(MeV/c) = ',sg,' +- ',ds
            write(*,*) ' FWHM(MeV/c) = ',ww,' +- ',dw
            write(5,*) ' FWHM(MeV/c) = ',ww,' +- ',dw
            write(*,*) ' Center(ch) = ',ce,' +- ',dc
            write(5,*) ' Center(ch) = ',ce,' +- ',dc
 101     continue
         chisq = 0.0
         do 443 i=1,ndata
            x = xx(i)
            w = 1./yy(i)
            chisq = chisq + (yy(i)-fmltig(x,aa,ncef,kerr))**2*w
 443     continue
         write(*,*) 'chisq = ',chisq
         write(5,*) 'chisq = ',chisq
         chisq = chisq/(dfloat(ndata-ncef))
         write(*,*) 'reduced chisq = ',chisq
         write(5,*) 'reduced chisq = ',chisq

         close(5)

 444     CONTINUE
         WRITE(*,*) ' Save fitting curve ? (y or n) '
         READ(*,'(A1)') ANS
         IF(ANS.EQ.'y') THEN
            WRITE(*,*) ' Input file name '
            READ(*,'(A30)') DATFNC
            OPEN (3,FILE=DATFNC,STATUS='UNKNOWN')

            DO 1143 I = 1,ndat
               X = DFLOAT(I)
               X = FROM + CMP/2. + CMP * (X - 1.)
               WRITE(3,*) fmltig(X,AA,NCEF,KERR)
c              WRITE(3,*) DBLG(X,AA,NCEF,KERR)
c              WRITE(3,*) FEXPN(X,AA,NCEF,KERR)
 1143       CONTINUE

         ELSEIF(ANS.EQ.'n') THEN
            GOTO 1101
         ELSE
            GOTO 444
         ENDIF

 1101 CONTINUE

      STOP
      END
