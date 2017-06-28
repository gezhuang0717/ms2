      program main
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XX(2000),YY(2000),ET(2000),TT(2000),SE(2000),ZT(2000)
      DIMENSION SZ(2000),AA(30),SA(30),CC(30),SC(30)
      DIMENSION DAT(2000)
c      dimension aai(30),II(30)                                                   
      CHARACTER DATFN*50,DATFNC*30,ANS,fcond*30,fres*30,
     &          ftblt*30,ftbldp*30

      common/PARAM/aai(30),II(30),ncoef

      EXTERNAL func,dfunc,d2func

      fcond = 'fcond2.dat'
      OPEN (1,FILE=fcond,STATUS='UNKNOWN')! open initial condition file             

      READ(1,'(A30)') fres                    ! read result file name             
      OPEN (5,FILE=fres,STATUS='UNKNOWN')     ! open result file                 
c      OPEN (5,FILE=fres,STATUS='NEW')     ! open result file                      
      read(1,'(a30)') datfnc
      write(*,*) 'output file : ',fres,datfnc
      write(5,*) 'output file : ',fres,datfnc

c-- read data -----------------------------------------------------------          
      READ(1,'(A30)') DATFN                       ! read data file name           
      write(*,'(a,a)') 'Data file name : ',DATFN
      write(5,'(a,a)') 'Data file name : ',DATFN
      OPEN (2,FILE=DATFN,STATUS='UNKNOWN')        ! open data file                 
      i = 1                                       ! read data file                

 1        read(2,*,end=10) dat(i)
c      write(*,*) i,dat(i)                                                        
      i = i + 1
      if(dat(i).le.0.0) dat(i)=0.0
      goto 1
 10      CONTINUE
      ndat = i - 1                          ! number of points of data            
      write(*,*) 'Number of points of data = ',ndat
      write(*,*) 'Number of points of data = ',ndat
      write(5,*) 'Number of points of data = ',ndat
      close(2)
c---------------------------------------------------------------------             


      read(1,*) from,to                      ! read range of x axis                

      CMP = (TO - FROM)/(DFLOAT(ndat)-1.d0)
      read(1,*) ncoef,ncef
      write(*,*) 'Number of parameters = ',ncoef
      write(5,*) 'Number of parameters = ',ncoef
      write(*,*) 'Number of free parameters = ',ncef
      write(5,*) 'Number of free parameters = ',ncef
      read(1,*) il,ih                 ! read fitting region in ch                  

      write(*,*) 'Fitting region : ',il,' -> ',ih
      write(5,*) 'Fitting region : ',il,' -> ',ih
      read(1,*) (II(i),i=1,ncoef)  ! read order of parameters                      
      read(1,*) (AAI(i),i=1,ncoef) ! read initial value of parameters              
      write(*,*) 'Initial value of free parameters'
      write(5,*) 'Initial value of free parameters'
      do 31 i = 1,ncef
         write(*,1000) II(i),AAI(II(i))
         write(5,1000) II(i),AAI(II(i))
 31         continue
      if(ncef.eq.ncoef) goto 34
      write(*,*) 'Initial value of fixed parameters'
      write(5,*) 'Initial value of fixed parameters'
      do 32 i = ncef+1,ncoef
         write(*,1000) II(i),AAI(II(i))
         write(5,1000) II(i),AAI(II(i))
 32         continue
 1000        FORMAT (2H  ,'AAI(',I2,') = ',D17.10)
 34             DO 33 I = 1,NCOEF
         AA(I) = AAI(II(I))
 33         CONTINUE

      DATA TT/2000*1.D0/
      DATA ET/2000*1.D0/
      DATA SE/2000*0.D0/
      DATA ZT/2000*1.D0/
      DATA SZ/2000*0.D0/

      NDATA = IH - IL + 1
      DO 22 N = 1,NDATA
         XX(N) = DFLOAT(IL+N-1)
         XX(N) = FROM + CMP * (XX(N)-1.)
         YY(N) = DAT(IL+N-1)
 22         CONTINUE

      NETA = 0
      KERR = 0

 1100  DO 1101 MODE = 1,2

         IF (MODE .EQ. 1) THEN
            WRITE(*,*)
            WRITE(5,*)
         ELSE
            WRITE(*,*)
            WRITE(5,*)

         ENDIF

         CALL CFDQR2 (XX,YY,TT,ET,SE,NDATA,NETA,AA,SA,NCEF,
     >         func,dfunc,d2func,ZT,SZ,CC,SC,0,func,MODE,KERR)

         WRITE(*,*) 'Fitting region(ch) : ',IL,' --> ',IH
         WRITE(5,*) 'Fitting region(ch) : ',IL,' --> ',IH
         WRITE(*,*) 'Fitting region (arb.) : ',XX(1),' --> ',XX(NDATA)
         WRITE(5,*) 'Fitting region (arb.) : ',XX(1),' --> ',XX(NDATA)
         WRITE(*,*) 'Free parameters '
         WRITE(5,*) 'Free parameters '
 1140        DO 1141 J = 1,NCEF
            WRITE (*,7000) II(J),AA(J),SA(J)
            WRITE (5,7000) II(J),AA(J),SA(J)
 7000              FORMAT (2H  ,'AA(',I2,') = ',D17.10,' +- ',D17.10)
 1141                  CONTINUE
         IF(NCEF.GE.NCOEF) GOTO 445
         WRITE(*,*) 'Fixed parameters '
         WRITE(5,*) 'Fixed parameters '
         DO 1142 J = NCEF+1,NCOEF
            WRITE (*,7001) II(J),AA(J)
            WRITE (5,7001) II(J),AA(J)
 7001              FORMAT (2H  ,'AA(',I2,') = ',D17.10)
 1142                  CONTINUE

 445                        CONTINUE

         chisq = 0.0
         do 443 i=1,ndata
            x = xx(i)
            if(yy(i).le.0.) then
               w = 1.d-10
            else
               w = 1./yy(i)
            endif
            fff = func(x,aa,ncef,kerr)
            chisq = chisq + (yy(i)-fff)**2*w
 443             continue
         write(*,*) 'chisq = ',chisq
         write(5,*) 'chisq = ',chisq
         chisq = chisq/(dfloat(ndata-ncef))
         write(*,*) 'reduced chisq = ',chisq
         write(5,*) 'reduced chisq = ',chisq

 444          CONTINUE
cc       WRITE(*,*) ' Save fitting curve ? (y or n) '                               
cc       READ(*,'(A1)') ANS                                                        
cc       IF(ANS.EQ.'y') THEN                                                       
cc          WRITE(*,*) ' Input file name '             
cc          READ(*,'(A30)') DATFNC                                                
         OPEN (3,FILE=DATFNC,STATUS='UNKNOWN')                                   
c         OPEN (3,FILE=DATFNC,STATUS='NEW')

            DO 1143 I = 1,ndat
               X = DFLOAT(I)
               X = FROM + CMP * (X - 1.)
               fff = func(x,aa,ncef,kerr)
cc             write(*,*) x,fff                                                   
               WRITE(3,*) fff
 1143                 CONTINUE

cc       ELSEIF(ANS.EQ.'n') THEN                                                  
cc          GOTO 1101                                                              
cc       ELSE                                                                      
cc          GOTO 444                                                              
cc       ENDIF                                                                     

c 2222    DO 3333 J = 1,NCEF                                                     

c            WRITE (*,*) AA(J)                                                    
c 3333    CONTINUE                                                                

 1101                  CONTINUE
ccc                                                                               

      write(*,*) ' '

      do 100 i=1,NCOEF
        write(*,*) AA(i)
 100      continue

      write(*,*) ' '

      do 101 i=1,NCOEF
        write(*,*) SA(i)
 101      continue

 999        stop

      END


