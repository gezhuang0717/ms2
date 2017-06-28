C================================================                                         
C     Maximum Likelihood\226@\202\311\202\346\202\351Fitting Programm                     
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

c      write(*,*)'NDATA,NEAT,AA,SA,NCEF = ',NDATA,NETA,AA,SA,NCEF                         
c      write(*,*)'MODE,KERR = ',MODE,KERR                                                 
        EPSM=1.D-10
        EPSF=1.D-6
        PI=3.141592653589793D0

        KERR=0
        NITER1=0
        FLMDA=0.D0

 1000      DO 1001 J=1,NCEF
           BB(J)=AA(J)
 1001         CONTINUE

        FCRD=0.D0
        RTPI=.5D0*DLOG(2.D0*PI)
 1010      DO 1011 N= 1,NDATA
           IF(YY(N) .LE. 1.D0) GOTO 1040
           Y=YY(N)
           IF (Y .LE. 50.D0) GOTO 1020
C===========================                                                              
C     Stirling\202\314\214\366\216\256                                                    
C===========================                                                              
           FCRD=FCRD+RTPI+(Y+.5D0)*DLOG(Y)-Y+1.D0/(12.D0*Y)
           GOTO 1040

 1020            FCRD=FCRD+DLOG(Y)
           Y=Y-1.D0
           IF (Y .GT. 1.D0) GOTO 1020

 1040            call coef(aa,ncef,aa)
           FITC(N)=FNCTD(XX(N),AA,NCEF,KERR)*TT(N)*ET(N)
           IF (KERR .GE. 5) RETURN
c          WRITE(*,*) N,XX(N),YY(N),FITC(N)                                               
 1011         CONTINUE

        IF (NCEFB .EQ. 0) GOTO 1080

 1060      DO 1061 N=1,NDATA
           call coef(cc,ncefb,cc)
           FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN
           FITC(N)=FITC(N)+FUNCA*TT(N)*ZT(N)
 1061         CONTINUE

 1080            FLV1=FLVAL2(YY,FITC,NDATA,INFL)
        IF(INFL .EQ. 0) GOTO 1100
        KERR=3
        GOTO 3820

C===========================                                                              
C     \212O\221\244\202\314\203\213\201[\203v                                             
C===========================                                                              

 1100      NITER2=0
        NITER1=NITER1+1

 1200      DO 1201 J=1,NCEF
           DD(J)=0.D0
 1220            DO 1221 K=1,NCEF
              RR(J,K)=0.D0
 1221               CONTINUE
 1201                  CONTINUE

 1300                     DO 1301 N=1,NDATA

           call coef(aa,ncef,aa)
           CALL FNDRV(XX(N),AA,FA,NCEF,KERR)
c           write(*,*) FA                                                                 
           IF (KERR .GE. 5)RETURN

           RPSI=1.D0/DSQRT(FITC(N))
           RPSI=1.D0/DSQRT(FITC(N))
           YD(N)=(YY(N)-FITC(N))*RPSI

 1320            DO 1321 J=1,NCEF
              QQ(N,J)=TT(N)*ET(N)*FA(J)*RPSI
              DD(J)=DD(J)+QQ(N,J)**2
 1321               CONTINUE
 1301                  CONTINUE

 1340                     DO 1341 J=1,NCEF
           DD(J)=1.D0/DSQRT(DD(J))
 1341         CONTINUE

 1360            DO 1361 N=1,NDATA
 1380                     DO 1381 J=1,NCEF
              QQ(N,J)=QQ(N,J)*DD(J)
 1381               CONTINUE
 1361                  CONTINUE

 1400                     DO 1401 J=1,NCEF
 1420                              DO 1421 K=J,NCEF

 1440                                          DO 1441 N=1,NDATA
                 RR(J,K)=RR(J,K)+QQ(N,J)*QQ(N,K)
 1441                     CONTINUE
 1421                           CONTINUE
           IF (J .EQ. NCEF) GOTO 1470

 1460            DO 1461 K=J+1,NCEF
              RR(K,J)=RR(J,K)
 1461               CONTINUE
 1470                     CONTINUE
 1401                        CONTINUE

        CALL MTXINV(RR,NCEF,DET)
        IF (DET .EQ. 0.D0) GOTO 3900
C=================================                       
C     \203\311C\202\314\216Z\217o                                                         
C=================================                                                        
        TRACE=0.D0
 1480      DO 1481 J=1,NCEF
           TRACE=TRACE+RR(J,J)
 1481         CONTINUE

        FLMDC=1.D0/TRACE
C=================================                                                        
C     \223\340\221\244\202\314\203\213\201[\203v                                          
C=================================                                                        
 2000      DO 2001 N=1,NDATA
 2020               DO 2021 J=1,NCEF
              QX(N,J)=QQ(N,J)
 2021               CONTINUE
 2001                  CONTINUE

 2040                     DO 2041 J=1,NCEF
           DA(J)=0.D0
 2041         CONTINUE

        NITER2=NITER2+1
C=================================                                                        
C     \217C\220\263\203}\203\213\203J\201[\203\213\202\314\225\373\226@                   
C=================================                                                        

 2100      DO 2101 N=NDATA+1,NDATA+NCEF
           K=N-NDATA
           YD(N)=0.D0
 2120            DO 2121 J=1,NCEF
              IF(J .NE. K) GOTO 2140
              QX(N,J)=DSQRT(FLMDA)
              GOTO 2121
 2140                  QX(N,J)=0.D0
 2121                        CONTINUE
 2101                           CONTINUE
C=================================                                                        
C     QR\225\317\212\267\202\360\214\304\202\321\217o\202\267                             
C=================================                                                        

 2200                     CALL QRTRN(QX,RR,JNX,NDATA+NCEF,NCEF,NRANK)
c      write(*,*)' NRANK = ',nrank                                                        
C=================================                                                        
C     \214\343\221\336\221\343\223\374\226@\202\311\202\346\202\351\214v\216Z             
C=================================                                                        
 2300                                 DO 2301 J=NRANK,1,-1
           JJ=JNX(J)
           QV(JJ)=0.D0
 2320            DO 2321 N=1,NDATA
              QV(JJ)=QV(JJ)+QX(N,JJ)*YD(N)
 2321               CONTINUE
           ZV=QV(JJ)
           IF (J .EQ. NRANK) GOTO 2360
 2340            DO 2341 K=J+1,NRANK
              KK=JNX(K)
              ZV=ZV-RR(J,KK)*DA(KK)
 2341               CONTINUE

 2360                     DA(JJ)=ZV/RR(J,JJ)

 2301                        CONTINUE

        DLLN=0.D0
 2400      DO 2401 J=1,NRANK
           JJ=JNX(J)
           DLLN=DLLN+QV(JJ)**2
           DA(JJ)=DA(JJ)*DD(JJ)
           BB(JJ)=AA(JJ)+DA(JJ)
           write(*,*) JJ,BB(JJ)
c           write(5,*) JJ,BB(JJ)                                                          
C 2402       FORMAT (10H          ,'AA(',I2,') = ',D17.10)   
 2401         CONTINUE

 2440            DO 2441 N=1,NDATA
           call coef(bb,ncef,bb)
           FITC(N)=FNCTD(XX(N),BB,NCEF,KERR)*TT(N)*ET(N)
c           write(*,*) XX(N),FITC(N)                                                      
           IF (KERR .GE. 5) RETURN
 2441         CONTINUE

        IF (NCEFB .EQ. 0) GOTO 2500

 2460      DO 2461 N=1,NDATA
           call coef(cc,ncefb,cc)
           FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN
           FITC(N)=FITC(N)+FUNCA*TT(N)*ZT(N)
 2461         CONTINUE

 2500            FLV2=FLVAL2(YY,FITC,NDATA,INFL)

        IF (INFL .EQ. 0) GOTO 2600
 2520      FLMDA=DMAX1(FLMDA*2.D0,FLMDC)
        IF (NITER2  .LE. 100) GOTO 2000
 2540      WRITE(*,*)
        WRITE(5,*)
        KERR=2
        RETURN

C===============================                                                          
C     \214W\220\224\202\314\214\353\215\267\202\314\214v\216Z                             
C===============================                                                          

 2600      DO 2601 J=1,NCEF
           SA(J)=0.D0
 2620            DO 2621 K=1,NCEF
              DFDA(J,K)=0.D0
 2621               CONTINUE
           IF (NCEFB .EQ. 0) GOTO 2601
 2640            DO 2641 L=1,NCEFB
              RR(J,L)=0.D0
 2641               CONTINUE
 2601                  CONTINUE

 3000                     DO 3001 N=1,NDATA
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

 3020            DO 3021 J=1,NCEF
              DFDY(N,J)=FA(J)*ETF

 3040                  DO 3041 K=J,NCEF
                 DFDA(J,K)=DFDA(J,K)-ETF**2*YY(N)*FA(J)*FA(K)
                 IF (MODE .EQ. 1) GOTO 3041
                 DFDA(J,K)=DFDA(J,K)+FACT*ET(N)*F2A(J,K)
 3041                     CONTINUE

              IF (J .EQ. NCEF) GOTO 3021
 3060                  DO 3061 K=J+1,NCEF
                 DFDA(K,J)=DFDA(J,K)
 3061                     CONTINUE
 3021                           CONTINUE

           IF (NETA .EQ. 0) GOTO 3100
 3080            DO 3081 J=1,NCEF
              QX(N,J)=(FACT-ETF**2*YY(N)*FUNC/ET(N))*FA(J)
 3081               CONTINUE

 3100                     IF(NCEFB .EQ. 0) GOTO 3001
           call coef(cc,ncefb,cc)
           CALL FNDRV(XX(N),CC,FC,NCEFB,KERR)
           IF (KERR .GE. 5) RETURN

 3120            DO 3121 J=1,NCEF
 3140                        DO 3141 L=1,NCEFB
                 RR(J,L)=RR(J,L)-ETX*FA(J)*FC(L)*ZT(N)
 3141                     CONTINUE

              IF (NETA .EQ. 1) GOTO 3121
              call coef(cc,ncefb,cc)
              FUNCA=FNCTB(XX(N),CC,NCEFB,KERR)
              IF (KERR .GE. 5) RETURN
              QZ(N,J)=-ETX*FUNCA*FA(J)

 3121               CONTINUE

 3001                  CONTINUE

C===========================================                                              
C     \215s\227\361D\202\314\213t\215s\227\361\202\314\214v\216Z                          
C===========================================                                              

        CALL MTXINV(DFDA,NCEF,DET)
        IF (DET .EQ. 0.D0) GOTO 3900


 3200      DO 3201 J=1,NCEF

 3220               DO 3221 N=1,NDATA
              DADY=0.D0
 3240                  DO 3241 K=1,NCEF
                 DADY=DADY-DFDA(J,K)*DFDY(N,K)
 3241                     CONTINUE
              SA(J)=SA(J)+DADY**2*FITC(N)

              IF (NETA .EQ. 0) GOTO 3221
              DADE = 0.D0
 3260                  DO 3261 K=1,NCEF
                 DADE=DADE-DFDA(J,K)*QX(N,K)
 3261                     CONTINUE
              SA(J)=SA(J)+(DADE*SE(N))**2
 3221               CONTINUE

           IF (NCEFB .EQ. 0) GOTO 3400
 3300            DO 3301 L=1,NCEFB
              DADC=0.D0
 3320                  DO 3321 K=1,NCEF
                 DADC=DADC-DFDA(J,K)*RR(K,L)
 3321                     CONTINUE
              SA(J)=SA(J)+(DADC*SC(L))**2
 3301               CONTINUE

           IF(NETA .EQ. 1) GOTO 3400
 3340            DO 3341 N=1,NDATA
              DADZ=0.D0
 3360                  DO 3361 K=1,NCEF
                 DADZ=DADZ-DFDA(J,K)*QZ(N,K)
 3361                     CONTINUE
              SA(J)=SA(J)+(DADZ*SZ(N))**2
 3341               CONTINUE

 3400                     SA(J)=DSQRT(SA(J))

 3201                        CONTINUE

C========================================                                                 
C     \216\373\221\251\202\314\224\273\222\350                                            
C========================================                                                 

 3500                           DLTR=FLV2-FLV1
        IF (DABS(DLTR) .GT. DABS(FLV2-FCRD)*EPSF) GOTO 3600

 3520      DO 3521 J=1,NCEF
           EPSA=DMAX1(SA(J)*1.D-6,DABS(BB(JJ))*EPSM)
           IF(DABS(DA(J)) .GT. EPSA) GOTO 3600
 3521         CONTINUE
        GOTO 3800

 3600      RAT=DLTR/DLLN
        IF (RAT .GE. 0.D0) GOTO 3700
        FU=DMAX1(2.D0,DMIN1(2.D0-RAT,1.D1))
        FLMDA=FLMDA*FU
        IF(FLMDA .EQ. 0.D0) FLMDA=FLMDC*FU*.5D0
        IF(NITER2 .LE. 100)THEN
           GOTO 2000
        ELSE
           GOTO 2540
        ENDIF

 3700      IF(RAT .LT. .75D0) GOTO 3720
        FLMDA=.5D0*FLMDA
        IF(FLMDA .LT. FLMDC) FLMDA=0.D0
        GOTO 3740
 3720      IF(RAT .GE. .25D0) GOTO 3740
        FLMDA=DMAX1(FLMDA*2.D0,FLMDC)

 3740      DO 3741 J=1,NCEF
           AA(J)=BB(J)
 3741         CONTINUE
        FLV1=FLV2
        IF(NITER1 .LE. 100) GOTO 1100
        WRITE (*,*)
        WRITE (5,*)

        KERR=2
        RETURN

 3800      DO 3801 J=1,NCEF
           AA(J)=BB(J)
 3801         CONTINUE
        IF(FLMDA .LE. 10.D0*FLMDC) RETURN

        WRITE(*,*)
        WRITE(5,*)
        KERR=1
        RETURN
C                                                                                         
c   KERR = 3                                                                              
C                                                                                         
 3820      WRITE(*,*)
       WRITE(5,*)
        RETURN
C                                                                                         
C   KERR = 4                                                                              
C                                                                                         
 3900      WRITE(*,*)
        WRITE(5,*)
        KERR=4
        RETURN
        END




