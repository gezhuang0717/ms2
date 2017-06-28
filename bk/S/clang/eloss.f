      Real*8 FUNCTION ELOSS(E,ZI,AI,MAT,THICK)
C      CALCULATE THE ENERGY LOSS OF NUCLEAR PARTICLE INSIDE THE 
C      MATERIAL (MAT)

C      E: ENERGY OF INCIDENT PARTICLE (IN MEV/NUCLEON) NO MASS ENERGY.
C      ZI: ATOMIC NUMBER OF INCIDENT PARTICLE
C      AI: MASS NUMBER OF INCIDENT PARTICLE
C      THMAT THICKNESS OF MATERIAL (G/CM**2)
C     THICK: THICKNESS OF MATERIAL IN (CM)
C     For the analysis of CH2 target, CH2 information is added.
C     To analise data older than 2006 Oct, 
C     one must use elossold.f

      Implicit Real*8 (A-H,O-Z)
      real*8 tenerg(138),trange(138,20)
      COMMON /RCMMN/ TENERG,TRANGE
      common /rcmmn2/ MIADJ
c      COMMON /RCMMN/ MIADJ,TENERG,TRANGE
      COMMON /RANGEO/ ERNG,RRNG,REGY,THMAT
      DIMENSION NAMAT(24)
      DIMENSION ZM(24),AM(24),AIADJ(24),RHO(24)
      CHARACTER*8 MAT,NAMAT

      DATA NAMAT /'AIR','PLASTIC','C','NAF','AL','SI','KCL','CU',
     + 'AG','PB','BE','N2','XE','NI','MO','XEN2','AR','HE','AU',
     + 'BAF2','TI','CSI','NAI','CH2'/ 
      DATA ZM / 7.2,3.38,6.0,10.0,13.0,14.0,18.0,29.0,47.0,82.0,
     +          4.0,7.0, 54.0,28.0,42.0,44.0,18.0,2.0,79.,24.667,
     +          22.0,54.,32.0,2.66667/
      DATA AM /14.41,6.24,12.01,21.0,26.98,28.086,37.28,63.54,
     +        107.87,207.9,9.01,14.0,131.30,58.69,95.94,101.2,39.948,
     +        4.002602,196.9665,58.441,47.867,129.91,74.95,4.66667/
*      DATA AIADJ/85.8,63.4,79.,135.4,164.0,172.,184.0,317.,487.,
*     +          793.,60.,91.0,758.0,304.50,438.82,591.3,210.,24.5874,
*     +          797.,269.3,247.4,562./
*      DATA AIADJ/85.8,63.4,79.,135.4,164.0,172.,184.0,317.,487.,
      DATA AIADJ/85.8,63.4,79.,135.4,163.0,172.,184.0,317.,487.,
     +          793.,60.,91.0,758.0,304.50,438.82,591.3,210.,24.5874,
     +          797.,388.,247.4,554.,431.,59.6/
c      DATA RHO /0.001205,1.032,2.267,2.56,2.7,2.33,1.98,8.96,10.5,
      DATA RHO /0.001205,1.032,1.79,2.56,2.7,2.33,1.98,8.96,10.5,
     +         11.35,1.848,0.00125,0.005896,8.902,10.22,0.00473,
     +         1.6617E-3,0.0001663,19.32,4.89,4.54,4.51,3.67,0.95/
c     carbon ionization potential was 1.55, which is changed to 2.25
c     how to calculate the ionization potential on NI, MO
c         experimental formula 
c           Iadj = (9.76 + 58.8 * Z**-1.19) * Z
c                                       (in case of Z>= 13)
c         by Bichsel-Uehling & Barkas-von FRIESEN    
      DATA GFAC /0.9317/

C      GFAC IS THE FACTER GIVE THE ACTUAL GAS DENSITY
C      ERNG: RANGE OF PRATICLE
C      RRNG: REMAINING RANGE OF PARTICLE
C      ZM: Z OF MATERIAL
C      AM: MASS OF MATERIAL
C      AIADJ: INIZATION POTENTIAL
C      RHO: DENSITY OF MATERIAL (G/CM**2)
C      REGY: THE PARTICLE ENERGY AFTER MATERIAL

C       *************************************************************

*      write(6,*)'E,ZI,AI,MAT,THICK',E,ZI,AI,MAT,THICK
      DO 10 I=1,24
      IF(MAT.EQ.NAMAT(I)) GO TO 20
   10 CONTINUE
C      TYPE *,'MAT= ',MAT
      write(6,*)'MAT= ',MAT
      K=19
      write(6,*)'GIVE ME ZM,AM,AIADJ,RHO,IGAS,ETA (4F,I2,F)'
      read(5,*)ZM(19),AM(19),AIADJ(19),RHO(19),IGAS,ETA
*      TYPE *,'GIVE ME ZM,AM,AIADJ,RHO,IGAS,ETA (4F,I2,F)'
*      ACCEPT 1000,ZM(19),AM(19),AIADJ(19),RHO(19),IGAS,ETA
 1000 FORMAT(4F10.5,I2,F10.5)
      GO TO 30
   20 CONTINUE
      K=I
      IGAS=0
      IF(K.EQ.1) IGAS=1
      ETA=1.0
      IF(K.EQ.1) ETA=GFAC
   30 CONTINUE

C      TYPE *,E,ZI,AI,ZM(K),AM(K),AIADJ(K),RHO(K),IGAS,ETA
      ERNG=RRANGE(E,ZI,AI,ZM(K),AM(K),AIADJ(K),1,1,1,RHO(K),IGAS,ETA)
      THMAT=THICK*RHO(K)
      RRNG=ERNG-THMAT
      IF(RRNG.GT.0.0) GO TO 40
      REGY=0.0
      GO TO 50
   40 REGY=RNERGY(RRNG,ZI,AI,ZM(K),AM(K),AIADJ(K),1,1,1,RHO(K),IGAS,ETA)
   50 ELOSS=E-REGY
c      write(*,*) E, regy, eloss
      RETURN
      END


