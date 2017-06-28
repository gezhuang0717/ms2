      Subroutine SIM(Ai,Zi,isotp,amp,thtgt,fl,ex,count,beame,eout)
c modified in 2015/03/04 by s.suzuki
c Ai:Mass number, Zi:Atomic number, isotp:non-reacted events(1),inelastic events(2)
c amp:Slope of inelastic-event probabaility function 
c thtgt:Thickness of the reaction target, fl:fluctuation of the beam
c ex:excitation energy for the inelastic event, count:Number of particles used in program
c beame:Incident beam energy, eout:Output energy per particle
c
c
c      Subroutine SIM(Apara,icall)
c
c     Simulation program for fitting of tail part
c                                                   
c
c     2002.10.25 M. Fukuda
c
*     Followings are considered.
*       + energy spread according to the position where reaction occurred
*         in the reaction target
*       + energy-loss fluctuation of the primary beam in the reaction target
*       + energy-loss fluctuation of fragments in the reaction target
*       + energy spread by reaction (Goldhaber)
*       + energy-loss fluctuation of all in the CsI counter
*       + thickness ununiformity of CsI counter
*       + energy resolution of CsI counter
*       + tof resolution 
*       + energy resolution of plastic counter
*       + energy resolution of Si counter
*
C      CALCULATE THE ENERGY OF A PARTICLE(P,D...) AFTER PASSING THROUGH
C      MATERIALS
C        ***  Precise Si stopping power included
C        ***  Material Input Method improved
C        ***  Precise C density included
C        ***  This program should be compiled and linked with
C                 eloss.f and rrange.f
C        ***  Eliminated output to for55.dat
C        ***  Excutable is "echange"
      Implicit Real*8 (A-H,O-Z)
      PARAMETER(SIG0=90.,NMASS=931.5,NMASS2=867692.25)

      COMMON /gat2D/XP(20,50),YP(20,50),NPNT(20)
*      COMMON /RDATA/EN,IEN,ENUM,SNUM,LNUM
      COMMON /RAW/RAW,IDATA
*      COMMON /FILE/NFILE,RFILE,DABLK
      common /eventcnt/evenm
      common /regscl/Elow,Ehigh,xbin,nbin,ILdata,IHdata

      integer count,ninel,Ai,Zi
      Real*8 beame,eout(count),thtgt,fl,ex
      

      DIMENSION Apara(10)
      DIMENSION E(20),P(20),MATHC(19),MATNAM(19),THICK(19),MATDEN(23)
      DIMENSION ERANGE(20),numat(19),DDE(20),matnam2(20),mathc2(19)
      dimension Neve(22),AIADJ(23)
      REAL*8 MATHC,matden,mathc2
      CHARACTER*8 MATNAM,MTNM(23),matnam2,matno(25),namat(23)
      CHARACTER AFILE*30,SFILE*30,DFILE*30,MFILE*30,PFILE*30
      Real RAW(1024)
      INTEGER IDATA*2(1024)
      integer spec1,spec2
      parameter(clv=0.2998)
      COMMON /ORDR/II(10)
      COMMON /RANGEO/ ERNG,RRNG,REGY,THMAT
      common /addone/spec1(256),spec2(256,256)
      DATA MTNM/'AIR','PLASTIC','BE','C','NAF','AL','SI','KCL','CU',
     &          'AG','PB','N2','XE','NI','MO','XEN2','AR','HE','AU',
     &          'BAF2','TI','CSI','NAI'/
      DATA matno/'AIR','BE','C','PLASTIC','WATER','MYLAR','NAF','AL',
     >           'SI','KCL','NI','CU','AG','MO','PB','N2','XE','XEN2',
     >           'AR','HE','AU','BAF2','TI','CSI','NAI'/
c      DATA MATDEN/0.001205,1.032,1.848,2.267,2.56,2.70,2.33,0.,8.96,
      DATA MATDEN/0.001205,1.032,1.848,1.79,2.56,2.70,2.33,0.,8.96,
     &		  10.5,11.35,0.00125,0.005896,8.902,10.22,0.00473,
     &            1.6617e-3,0.0001663,19.32,4.89,4.54,4.51,3.67/

      DATA NAMAT /'AIR','PLASTIC','C','NAF','AL','SI','KCL','CU',
     + 'AG','PB','BE','N2','XE','NI','MO','XEN2','AR','HE','AU',
     + 'BAF2','TI','CSI','NAI'/
      DATA AIADJ/85.8,63.4,79.,135.4,164.0,172.,184.0,317.,487.,
     +          793.,60.,91.0,758.0,304.50,438.82,591.3,210.,24.5874,
     +          797.,388.,247.4,554.,431./
      SAVE

1	CONTINUE
*        write(6,*)' start sim',icall
	do i=1,20
	  e(i)=0.
	  p(i)=0.
	enddo
        afile='test.dat'
        OPEN (3,FILE=afile,STATUS='unknown')

ccc        write(6,*) Ai,Zi,amp,thtgt,count,beame,isotp
ccc        return

cc      PFILE = 'Mg24.dat'
cc      AFILE = 'upper.ana'
cc      SFILE = '../pawdata/sim.rzf'
cc      MFILE = 'reso.dat'
cc      DFILE = 'upper.dat'

c        print '(I6)', count
      
ccc      if(ICALL.eq.1) then
ccc
ccc      OPEN(3,FILE=PFILE,FORM='FORMATTED',
ccc     >     err=902,STATUS='OLD')
ccc      read(3,*) Z0
ccc      read(3,*) A0
cccc      read(3,*) beame
ccc      read(3,*) beame2
ccc      read(3,*) LAYER
ccc      read(3,*) ntgt
ccc      read(3,*) thtgt
ccc      read(3,*) numsi1
ccc      read(3,*) sigde
ccc      read(3,*) sigsi1
ccc      read(3,*) sigsi2
ccc      read(3,*) sigsi3
ccc      read(3,*) sigsi4
ccc      read(3,*) signai
ccc*      close(3)
cccc      close(2)
ccc
ccc*      endif
c      icall = icall + 1

cc      Z0=11
cc      A0=30
      A0=dble(Ai)
      Z0=dble(Zi)

      LAYER=2
      ntgt=3
cc      thtgt=10.
c for wedge shaped target
c      thtgt=9.5+0.00983871*xpos
c for plate target 3.6 g/cm2
c      thtgt=19.48
      numsi1=2

      sigde=0.00000001
      sigsi1=0.00000001
      sigsi2=0.009327  
      sigsi3=0.009165  
      sigsi4=0.00759   
      signai=0.004 

c inelastic events
      ninel=1000

      IX = 999
      Nis=2
      Neve(1)=100000
      Neve(1) = count
      Neve(2)=count
c      Neve(2)=1
c      Neve(2)=ninel
c      Neve(2)=1000000

!4000
      Neve(3)=0
!120
      Neve(4)=10000
!740     
      Neve(5)=120
      Neve(6)=10
      Neve(7)=10
      Neve(8)=840 
      Neve(9)=1335
      Neve(10)=105
      Neve(11)=185
      Neve(12)=700 
      Neve(13)=535 
      Neve(14)=20 
      Neve(15)=150
      Neve(16)=590 
      Neve(17)=1040
      Neve(18)=920 
      Neve(19)=437 
      Neve(20)=408 
      Neve(21)=810 
      Neve(22)=382 


*      call INIT
      do k=1,256
         spec1(k)=0
      enddo

* Energy Loss in Reaction Target
      N = ntgt
      matnam(1)=matno(N)
      MATHC(1)=thtgt
      mathc(1)=mathc(1)/10.d0
      matnam2(1)=matnam(1)
      mathc2(1)=mathc(1)
      numat(1)=0
      do j=1,23
         if(matnam(1).eq.mtnm(j)) numat(1)=j
         if(matnam(1).eq.namat(j)) namtgt=j
      enddo
*      write(6,*)'before first ELOSS'
*      write(6,*)thtgt,mathc(1)
*      write(6,*)'beame,Z0,A0,matnam(1),mathc(1)',beame,z0,a0,
*     &           matnam(1),mathc(1)
      detgt = ELOSS(beame,Z0,A0,MATNAM(1),MATHC(1))
*      detgt = ELOSS(einc,Z0,A0,MATNAM(1),MATHC(1))

*      write(6,*)' energy loss in target : ',detgt

*      type *,' '
*      TYPE *,' Incident Particle, Z, A, E in MeV/nucleon'
*      ACCEPT *,ZI,AI,E(1)

*------------------------------------
*      do 777 isotp = 1,Nis
c      isotp = 1
c      isotp = 2
 2    if(isotp.eq.1)then
         ZI=Z0
         AI=A0
      elseif(isotp.eq.2)then
         ZI=Z0
         AI=A0
      elseif(isotp.eq.3)then
         ZI=Z0
         AI=A0+1
      elseif(isotp.eq.4)then
         ZI=Z0
         AI=A0-1
      elseif(isotp.eq.5)then
         ZI=Z0
         AI=A0-2
      elseif(isotp.eq.6)then
         ZI=Z0
         AI=A0-3
      elseif(isotp.eq.7)then    
         ZI=Z0-1
         AI=A0
      elseif(isotp.eq.8)then
         ZI=Z0-1
         AI=A0-1
      elseif(isotp.eq.9)then
         ZI=Z0-1
         AI=A0-2
      elseif(isotp.eq.10)then
         ZI=Z0-1
         AI=A0-4
      elseif(isotp.eq.11)then
         ZI=Z0-2
         AI=A0-2
      elseif(isotp.eq.12)then
         ZI=Z0-2
         AI=A0-3
      elseif(isotp.eq.13)then
         ZI=Z0-2
         AI=A0-5
      elseif(isotp.eq.14)then
         ZI=Z0-3
         AI=A0-3
      elseif(isotp.eq.15)then
         ZI=Z0-3
         AI=A0-4
      elseif(isotp.eq.16)then
         ZI=Z0-3
         AI=A0-5
      elseif(isotp.eq.17)then
         ZI=Z0-3
         AI=A0-6
      elseif(isotp.eq.18)then
         ZI=Z0-4
         AI=A0-8
      elseif(isotp.eq.19)then
         ZI=Z0-4
         AI=A0-9
      elseif(isotp.eq.20)then
         ZI=Z0-5
         AI=A0-9
      elseif(isotp.eq.21)then
         ZI=Z0-5
         AI=A0-10
      elseif(isotp.eq.22)then
         ZI=Z0-5
         AI=A0-11
      endif
      
!
ccc      if(zi.le.0.0) go to 999
!

*      write(6,*)' isotp : ',isotp

      if(AI.ne.A0 .and. AI.lt.A0) then
         sigmom=SIG0*dsqrt(AI*(A0-AI)/(A0-1))
      elseif(AI.ne.A0 .and. AI.gt.A0) then
         sigmom=SIG0*dsqrt(A0*(AI-A0)/(AI-1))   
      endif

*=========== start event loop ==========
      do 333 iii=1,Neve(isotp)

c Energy spread of the beam
c       sigeinc = beame*0.0021
       sigeinc = beame*fl
       beamefl = gausran(sigeinc,IX)
       einc = beame+beamefl

* Reaction Position in the Target
      if(isotp.eq.1) then
         xx = 1.
      else
         xx = urand2(IX)
      endif
      thtgtr = thtgt*(1.-xx)
      deprj = ELOSS(einc,Z0,A0,MATNAM(1),thtgt*xx/10.)

      if(isotp.eq.2 .and. AI.eq.A0 .and. ZI.eq.Z0 ) then
         a = amp
         efrgi  = (einc-deprj+(dlog(1.d0-urand2(IX)))/a-ex/A0)
c Energy loss in the rest of target thickness, for inelastic events
         efrg = efrgi - Eloss(efrgi,Z0,A0,MATNAM(1),thtgtr/10.)-0.2
c for day2 exp. -0.2[MeV/u] (Energy loss in F5PL2 after the target) 
c  first exited state(12c 4.439,9Be 1.665,27Al 0.884)    

      elseif(AI.ne.A0 .and. AI.gt.A0) then
         deprj = ELOSS(einc,Z0,A0,MATNAM(1),thtgt*xx/10.)
         efrg = (A0/AI)**2*(einc - deprj) ! energy per nucleon

      else
         efrg = einc - deprj ! energy per nucleon for non-reacted events
      endif
      
      if(efrg.lt.0.0000001) then
         efrg = 0.
      else
         efrg = efrg
      endif


* Caluculation of energy-loss fluctuation
      pinz=AIADJ(namtgt)*1.e-6  ! Ionization Potential in MeV
      x=dlog(4.*0.511*einc/931.5/pinz) ! logarizumic term
      ss=0.511/(931.5*ai)/x*deprj*ai*(2.*einc-deprj)*ai
      ss=dsqrt(ss)/ai
      if(ss.eq.0)then
         sfl = 0.
      else
         sfl=gausran(ss,IX)
      endif
      efrg = efrg + sfl


* Energy spread by reaction (Goldhaber)
      if(AI.ne.A0) then
         pfl=gausran(sigmom,IX)
         pfrg=AI*dsqrt((NMASS+efrg)**2-NMASS2) + pfl
         efrg=dsqrt((pfrg/AI)**2+NMASS2) - NMASS
      endif


cccccccccccccccccccccccccccccccccccccccccccccc
cc original
cc      if(isotp.eq.1) then
c      if(isotp.eq.1.or.isotp.eq.2) then
c         xx = 1.
c      elseif(isotp.eq.3) then
c 77      xxp = urand2(IX)
c         reratep=0.94296+0.052103*xxp+0.0048666*xxp**2
c
c         xx = urand2(IX)
c         yy=xx
c
c         if(yy.le.reratep) then
c            xx = xxp
c         else
c            goto 77
c         endif
c
c      else
c         xx = urand2(IX)
cc      write(6,*)'********',deprj,matnam(1),thtgtr,'***********'
c      endif
c
c      thtgtr = thtgt * (1.-xx)
c
c
c
cccc modified by s.s
cccc       sigeinc = beame*0.00283
c       sigeinc = beame*0.0021
c       beamefl = gausran(sigeinc,IX)
cc       beamefl=0.0
c       einc = beame+beamefl
c      if(isotp.eq.2 .and. AI.eq.A0 .and. ZI.eq.Z0 ) then
c      
c         deprj = ELOSS(einc,Z0,A0,MATNAM(1),thtgt*xx/10.)
cc      a=1.3
cc         a = Apara(II(1))
c         a = amp
cc s.s
cc         efrg  = (einc-deprj+(dlog(1.d0-urand2(IX)))/a-4.439/A0)
c         efrg  = (einc-deprj+(dlog(1.d0-urand2(IX)))/a-ex/A0)
c
cc  first exited state(12c 4.439,9Be 1.665,27Al 0.884)    
c
c*         write(6,*)' efrag',efrg
c*         if(iii.eq.1)write(6,*)' iii,a,efrg ',iii,a,efrg
c*         efrg  = (einc-deprj+(dlog(1.-urand2(IX)))/a-0.07033)
cc      efrg = einc-deprj
cc      efrg  = (einc-(einc-deprj)*urand2(IX)*urand2(IX)-deprj)
c
c      elseif(AI.ne.A0 .and. AI.gt.A0) then
c         deprj = ELOSS(einc,Z0,A0,MATNAM(1),thtgt*xx/10.)
c         efrg = (A0/AI)**2*(einc - deprj) ! energy per nucleon
c
c      else
c         deprj = ELOSS(einc,Z0,A0,MATNAM(1),thtgt*xx/10.)
c         efrg = einc - deprj ! energy per nucleon
c      endif
c      
c      if(efrg.lt.0.0000001) then
c         efrg = 0.
c      else
c         efrg = efrg
c      endif
c
c
c
c
c
c* Caluculation of energy-loss fluctuation
c      pinz=AIADJ(namtgt)*1.e-6  ! Ionization Potential in MeV
cc      write(6,*)' pinz : ',pinz
c      x=dlog(4.*0.511*einc/931.5/pinz) ! logarizumic term
c      ss=0.511/(931.5*ai)/x*deprj*ai*(2.*einc-deprj)*ai
c      ss=dsqrt(ss)/ai
c      if(ss.eq.0)then
c         sfl = 0.
c      else
c         sfl=gausran(ss,IX)
c      endif
ccc      write(*,*) 'sfl=',sfl
c      efrg = efrg + sfl
c
c
c* Energy spread by reaction (Goldhaber)
c      if(AI.ne.A0) then
c
c*         write(6,*)'efrg : ',efrg
c         pfl=gausran(sigmom,IX)
c*         pfl = 0.
c         pfrg=AI*dsqrt((NMASS+efrg)**2-NMASS2) + pfl
c         efrg=dsqrt((pfrg/AI)**2+NMASS2) - NMASS
c*         write(6,*)'sigmom,efrg : ',sigmom,efrg
c      endif
cccccccccccccccccccccccccccccccccccccccccccccc


 1010 FORMAT(A8)
 1001 FORMAT(I2)

c      if(efrg.le.0.) then
c         deltae = 0.
c         tof=100000.
*         tof1=1000000.
*         tof2=1000000.
c      goto 33
c      endif
         
      e(1)=efrg
c added by s.s
c      write(6,*) efrg
      eout(iii)=efrg
c      write(3,*) efrg-212.35

      go to 111


c      close(3)
      DO 20 I=1,LAYER
*      TYPE 1030,I
*      type *,' 1: Air       2: Be        3: C        4: Pla.Scint.'
*      type *,' 5: Water     6: Mylar     7: NaF      8: Al'
*      type *,' 9: Si       10: KCl      11: Ni      12: Cu'
*      type *,'13: Ag       14: Mo       15: Pb      16: N2'
*      type *,'17: Xe       18: XeN2     19: Ar      20: He'
*      type *,'21: Au       22: BaF2     23: Ti      24: CsI'
*      type *,'25: NaI'
   
         if(I.eq.1) then
            N=ntgt       
            if(isotp.eq.1) then
               mathc(i) = 0.
            else
               mathc(i) = thtgtr
            endif
         elseif(i.eq.2) then
            N=25
            mathc(i) = 1000.
         endif

c         else
c           OPEN (3,FILE=DFILE,FORM='FORMATTED',STATUS='OLD')
c           read(3,*) N         
c           read(3,*) mathc(i)
           
       

         MATNAM(I) = matno(N)
*      TYPE *,' Thickness of Material in mm'
*      ACCEPT *,MATHC(I)
*      MATHC(I)=thde
         mathc(i)=mathc(i)/10.
         matnam2(i)=matnam(i)
         mathc2(i)=mathc(i)
         numat(i)=0
         do j=1,25
            if(matnam(i).eq.mtnm(j)) numat(i)=j
         enddo
   
* Calculation of Energy Loss
*         if(iii.eq.1)write(6,*)' Layer : ',i,'   Material : ',matnam(i)
*         if(iii.eq.3)write(6,*)' Layer : ',i,'   Material : ',matnam(i)
*         if(iii.eq.4)write(6,*)' Layer : ',i,'   Material : ',matnam(i)
         if(E(I)-ELOSS(E(I),ZI,AI,MATNAM(I),MATHC(I)).lt.0.0001) then
*            write(6,*)' inside if'
            do il=i,LAYER
               e(il+1)=0.
               de=de+e(il)
               dde(il)=e(il)
            enddo
            goto 18
         endif
       
!         if(iii.eq.2)write(6,*)10*mathc(i),e(i),matnam(i)
       
         E(I+1)=E(I)-ELOSS(E(I),ZI,AI,MATNAM(I),MATHC(I))
         DE=DE+ELOSS(E(I),ZI,AI,MATNAM(I),MATHC(I))
         DDE(I)=ELOSS(E(I),ZI,AI,MATNAM(I),MATHC(I))
 18      ERANGE(I)=ERNG
         THICK(I)=THMAT
         IF(E(I+1).LE.0.0) GO TO 30
*         if(iii.eq.1)write(6,*)E(I),dde(i),zi,ai,matnam(i),mathc(i)        
 30   CONTINUE
 20   CONTINUE
   

*      write(6,*)'de1,2 = ',dde(1),dde(2)
*      write(6,*)' e15= ',e(15)
*         write(6,*)'iii : ',iii
*         write(6,*)'DDE(1),sigdeI,X',DDE(1),sigde,IX
*         deltae=FLUCT(DDE(1),sigde,IX) ! sigde => fraction
         if(DDE(numsi1).lt.1e-5) then
            deltae1 = 0.
         else
            deltae1 = FLUCT(DDE(numsi1),sigde,IX)*AI
            deltae1 = FLUCT(deltae1,sigsi1,IX) ! counter resolution
            deltae1 = deltae1
         endif

         if(DDE(numsi2).lt.1e-5) then
            deltae2 = 0.
         else
            deltae2 = FLUCT(DDE(numsi2),sigde,IX)*AI
            deltae2 = FLUCT(deltae2,sigsi2,IX) ! counter resolution
            deltae2 = deltae2
         endif

         if(DDE(numsi3).lt.1e-5) then
            deltae3 = 0.
         else
            deltae3 = FLUCT(DDE(numsi3),sigde,IX)*AI
            deltae3 = FLUCT(deltae3,sigsi3,IX) ! counter resolution
            deltae3 = deltae3
         endif

         if(DDE(numsi4).lt.1e-5) then
            deltae4 = 0.
         else
            deltae4 = FLUCT(DDE(numsi4),sigde,IX)*AI
            deltae4 = FLUCT(deltae4,sigsi4,IX) ! counter resolution
            deltae4 = deltae4
         endif

         if(DDE(numnai).lt.1e-5) then
            deltae5 = 0.
         else
            deltae5 = FLUCT(DDE(numnai),sigde,IX)*AI
            deltae5 = FLUCT(deltae5,signai,IX) ! counter resolution
            deltae5 = deltae5
         endif
         if(deltae5.gt.0)then
            deltae = deltae1
cdeltae2+deltae3+deltae4
         else
            deltae = deltae1
         endif
c        write(6,*)deltae
****************************************************************
*     Add One 
****************************************************************
*         write(6,*)' addone',deltae
ccc         call Add1(1,deltae)

cc         print '(f6.2)',  deltae

 
 111       continue

 333  continue
*=========== end event loop =================
       return

 777  continue

cc      write(6,*)' Event loop end'
*Output
*      open(1,FILE='spec2.dat',status='UNKNOWN')
*      do jj=1,256
*         write(1,'(256I2)')(spec2(jj,ii),jj=1,256)
*         write(1,'(256I6)')(spec2(jj,ii),ii=1,256)
*      enddo
*      close(1)

      DO 40 I=1,LAYER
         PMOM=(E(I)+931.5)**2-867692.25
         P(I)=DSQRT(PMOM)
   40 CONTINUE
 1020 FORMAT(/10X,20(1H*),3X,F10.5,1X,'MeV/nucleon',3X,F10.5,1X,
     +      'MeV/c/nucleon'/
     +     / 20X,' Remaining Range = ', F10.5,'g/cm**2'/
     +      /11X,A8/15X,F10.5,'g/cm**2',5X,
     +      F10.5,' cm'/)
      I=LAYER+1
      IF(E(I).GT.0.0) THEN
         PMOM=(E(I)+931.5)**2-867892.25
         P(I)=DSQRT(PMOM)
      ELSE
         P(I)=0.0
      END IF

*      write(6,*)' end of sim'

*	TYPE 1505,INT(ZI),INT(AI),E(1)
 1505	FORMAT(/'  INCIDENT PARTICLE (Z,A) : (',I3,',',I3,')',
     &         '    E =',F8.4,' MeV/u')
*	type 1500
 1500	format(/'  LAYER   MATERIAL   THICKNESS     RRANGE     RRANGE    
     & OENERGY','     NENERGY')
*	TYPE 1510
 1510	FORMAT( '                        (um)      (mg/cm2)     (um)     
     & (MeV/u)','     (MeV/u)') 
*	DO 500 I=1,LAYER
* 500	type 1520,I,MATNAM2(I),MATHC2(I)*1E4,ERANGE(I)*1E3,
*     &            ERANGE(I)/MATDEN(NUMAT(I))*1E4,E(I),E(I+1)
 1520	FORMAT(2X,I3,6X,A8,2X,5(F9.2,2X))
*	type 1530
 1530	format(//'  LAYER     ELOSS      ELOSS')
*        type 1540
 1540	format(  '           (MeV/u)     (MeV)')
*	DO 600 I=1,LAYER
* 600	TYPE 1550,I,DDE(I),DDE(I)*AI 
 1550	FORMAT(2X,I3,3X,2(F9.2,2X))
*	type *,' '
*      TYPE *,'                  TOTAL ENERGY LOSS = ',DE,'   MeV/u'
        return
!
*        DE=0.
*	go to 1
!
ccc 999      STOP
 902     write(6,*)'ana file open error'
      END subroutine SIM



******************************
      SUBROUTINE INIT
C     1��anafile���顢bin����DATA�κ����͡��Ǿ��ͤ��ɤ߹��ߡ�
C        ����(PAWC)��˥ҥ��ȥ������ΰ����ݤ��롣
C     2��anafile���顢gate�˴ؤ���������Ф���
C        GID,HID,LOGID�ˤ��ξ����񤭹��ࡣ
******************************
      Implicit Real*8 (A-H,O-Z)
      COMMON /gat2D/XP(20,50),YP(20,50),NPNT(20)
      COMMON /RAW/RAW,IDATA
      COMMON /PARAM/GNUM,GNUM2,LOGNUM,HNUM,NTU
      REAL*8 GID(999,4),GID2(20,4)
      INTEGER HID(999,4),LOGID(999,4)
      LOGICAL CUTID(999)
      CHARACTER COMG*4,TITLE*30
      REAL RAW(1024)
      INTEGER IDATA*2(1024)
      REAL*8 NUM(9)
      INTEGER HNUM,GNUM,GNUM2,LOGNUM,NTU
C      REAL XMAX,XMIN,YMAX,YMIN,WGT
*      INTEGER XBIN,YBIN
      REAL*8 WGT
      HNUM = 1
      GNUM = 1
      GNUM2 = 1
      LOGNUM = 1
      NTU = 0
      
*      write(6,*)' before read in INIT'
* 1    READ (4,10,END=999) COMG,(NUM(I),I=1,9,1),TITLE
*      write(6,*)'COMG,NUM,TITLE',COMG,NUM,TITLE
*      IF (COMG(1:4) .EQ. 'end') THEN
*         go to 999
*      ELSE IF (COMG(1:4) .EQ. 'gate') THEN
*         GO TO 100
*      ELSE IF (COMG(1:3) .EQ. 'and') THEN
*         GO TO 170
*      ELSE IF (COMG(1:2) .EQ. 'or') THEN
*         GO TO 150
*      ELSE IF (COMG(1:4) .EQ. 'hst1') THEN
*         GO TO 200
*      ELSE IF (COMG(1:4) .EQ. 'hst2') THEN
*         GO TO 300
*      ELSE IF (COMG(1:3) .EQ. 'not') THEN
*         GO TO 190
*      else if (comg(1:4) .eq. 'gat2') then
*         go to 400
*      ELSE
*         GO TO 1   
*      ENDIF      

      
 10   FORMAT (A4,F6.0,F6.0,F6.0,F6.0,F6.0,F6.0
c    >     ,F6.0,F6.0,F6.0,A10)
     >     ,F6.0,F6.0,F6.0,A30)

*******************************************************
C     Read-In Gate Conditions
*******************************************************
 
****	Gate	****
* 100  CONTINUE
*      DO 101 K=1,4
*         GID(GNUM,K) = NUM(K)
* 101  CONTINUE
*      GNUM = GNUM + 1
*      GO TO 1

****  2D-Gate   ****
* 400  continue
*      do 401 K=1,4
* 401     GID2(GNUM2,K) = NUM(K)
*      GNUM2 = GNUM2 + 1
*      go to 1
      
****	Or	****
* 150  CONTINUE
*      LOGID(LOGNUM,1) = 1
*      DO 151 K=1,3
*         LOGID(LOGNUM,K+1) = INT(NUM(K))
* 151  CONTINUE
*      LOGNUM = LOGNUM + 1
*      GO TO 1

****	And	**** 
* 170  CONTINUE
*      LOGID(LOGNUM,1) = 2
*      DO 171 K=1,3
*         LOGID(LOGNUM,K+1) = INT(NUM(K))
* 171  CONTINUE
*      LOGNUM = LOGNUM + 1
*      GO TO 1

****     Not    ****
* 190  CONTINUE
*      LOGID(LOGNUM,1) = 3
*      DO K=1,2
*      LOGID(LOGNUM,K+1) = INT(NUM(K))
*      ENDDO
*      LOGNUM = LOGNUM + 1
*      GO TO 1

******************************************************
C     Making 1-D Histogram
******************************************************
 
* 200  XBIN = INT(NUM(3))
C      XMIN = NUM(4)
C      XMAX = NUM(5)
*      HID(HNUM,1) = 1
*      HID(HNUM,2) = INT(NUM(1))    ! gate #
*      HID(HNUM,3) = INT(NUM(2))    ! data #
*      HID(HNUM,4) = 0
c      write(6,*)' init ',num
* 210  HNUM = HNUM + 1
*      GO TO 1

******************************************************
C     Making 2-D Histogram
******************************************************

* 300  XBIN = INT(NUM(3))
*      YBIN = INT(NUM(7))
C      XMIN = NUM(4)
C      XMAX = NUM(5)
C      YMIN = NUM(8)
C      YMAX = NUM(9)
*      HID(HNUM,1) = 2
*      HID(HNUM,2) = INT(NUM(1))  ! gate #
*      HID(HNUM,3) = INT(NUM(2))  ! x-spec #
*      HID(HNUM,4) = INT(NUM(6))  ! y-spec #
* 310  HNUM = HNUM + 1
*      GO TO 1

 999  CONTINUE
*      WRITE(*,*) 'Histogram initialization Completed!'
*      write(*,*) ' '
*      write(6,*) 'HNUM ',HNUM
c      WRITE(3,*) 'Histogram initialization Completed!'
c      write(3,*) ' '

      RETURN

      END         



      Real*8 Function FLUCT(value,sig,IX)
*
*
*     Fluctuate given value with given width sigma
*
      Implicit Real*8 (A-H,O-Z)
*      Real*8 sigma

      sigma=value*sig
      delta=GAUSRAN(sigma,IX)
*      write(6,*)'delta,sigma',delta,sigma
      
      FLUCT = value + delta

      if(FLUCT.lt.0.) FLUCT = 1.e-4
*      write(6,*)'FLUCT',FLUCT
      return
      end


      Subroutine ADD1(i,value1)

      Implicit Real*8 (A-H,O-Z)
      common /addone/spec1(256),spec2(256,256)
      common /regscl/Elow,Ehigh,xbin,nbin,ILdata,IHdata
      integer xch,ych,spec1,spec2

      xch = (value1-Elow)/xbin + 1
*      write(6,*)'Ehigh,Elow,nbin,xbin,xch',Ehigh,Elow,nbin,xbin,xch

      if(xch.gt.256) return
      if(xch.lt.1) return

      spec1(xch) = spec1(xch) + 1

      return
      end


      Subroutine ADD2(i,value1,value2)

      Implicit Real*8 (A-H,O-Z)
      common /addone/spec1(256),spec2(256,256)
      common /regscl/Elow,Ehigh,xbin,nbin,ILdata,IHdata
      integer xch,ych,spec1,spec2

*      xbinn = (xmax-xmin) / 256.
*      ybin = (ymax-ymin) / 256.

*      xch = (value1-xmin)/xbinn
*      ych = (value2-ymin)/ybin

*      write(6,*)xch,ych
*      if(xch.gt.256)xch=256
*      if(ych.gt.256)ych=256
*      if(xch.lt.1)xch=1
*      if(ych.lt.1)ych=1

*      spec2(xch,ych)=spec2(xch,ych)+1

      return
      end
