      SUBROUTINE EncBeam(val,nx,ny,naok,
     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata)
c Modified by abey @ 2015.12.01
c-----------------------------------------------------------------------
c ANALYZER 2 : Encbeam
c-----------------------------------------------------------------------
c  for raw data histgrams
c  ID=1 ILC1 CD-PPAC X1, X2 raw
c  ID=2 ILC1 CD-PPAC Y1, Y2 raw
c  ID=3 ILC2 CD-PPAC X1, X2 raw
c  ID=4 ILC2 CD-PPAC Y1, Y2 raw
c  ID=5 ILC2 T-shape L-side raw  ("L" means view from downstream)
c  ID=6 ILC2 T-shape R-side Up raw
c  ID=7 ILC2 T-shape R-side Down raw
c  ID=8 ELC Plastic Q (U,D) raw
c  ID=10 R-MD1, 2 raw 
c  ID=11 R-MD3, 4 raw 
c  ID=12 NaI raw
c  ID=24 R-MD5 raw, dummy  
c  ID=26 S0 Plastic T (L,R) ("L" means view from downstream)
c  ID=27 S0 Plastic Q (L,R)
c  ID=31 F2 Pla L Q,T
c  ID=32 F2 Pla R Q,T
c  ID=33 F3 Pla L Q,T
c  ID=34 F3 PLa R Q,T
c  ID=35 F3 PLa CFD L,R
c  for view true histgrams
c-----------------------------------------------------------------------
c  ID=90 PPAC X, Y (F3, F4, F5,F6, FH9, FH10, S0)
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F31aX F31aY F31bX F31bY F32aX F32aY F32bX F32bY
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             F51aX F51aY F51bX F51bY F52aX F52aY F52bX F52bY
c
c  W# :  21   22    23    24    25    26    27    28    29    30
c             F4X   F4Y   F6X   F6Y   
c
c  W# :  31   32    33    34    35    36    37    38    39    40
c             FH9X  FH9Y  FHXX  FHXY  S0X   S0Y
c
c-----------------------------------------------------------------------
c  ID=201 F3 reconstruction
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F31X  F31Y  F32X  F32Y 
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             F3X         F3Y
c
c  W# :  21   22    23    24    25    26    27    28    29    30
c             F3A         F3B
c
c-----------------------------------------------------------------------
c  ID=202 F5 reconstruction
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F51X  F51Y  F52X  F52Y 
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             F5X         F5Y
c
c  W# :  21   22    23    24    25    26    27    28    29    30
c             F5A         F5B
c
c-----------------------------------------------------------------------
c  ID=203 S0 reconstruction
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             S0X         S0Y
c
c  W# :  21   22    23    24    25    26    27    28    29    30
c             S0A         S0B
c
c-----------------------------------------------------------------------
c  ID=204 Dispersion correction calc. for F6, ILC1, ILC2
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F46         F6ILC1      F6ILC2
c
c-----------------------------------------------------------------------
c  ID=101 ILC1 PPAC position
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             X1-X2 X1+X2 Y1-Y2 Y1-Y2       
c
c  W# :  21   22    23    24    25    26    27    28    29    30
c             ILC1X       ILC1Y
c
c-----------------------------------------------------------------------
c  ID=102 ILC2 PPAC position
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             X1-X2 X1+X2 Y1-Y2 Y1-Y2       
c
c  W# :  21   22    23    24    25    26    27    28    29    30
c             ILC2X       ILC2Y
c
c-----------------------------------------------------------------------
c  ID=103 ILC2 T-shape plastic
c  W# :  1    2      3      4      5      6      7      8      9      10
c        ID   LTcal  RUTcal RDTcal       RTavr.         Tavr
c
c-----------------------------------------------------------------------
c  ID=104 R-MD plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   MD1   MD2   MD3   MD4   MD5
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c                         MD3   MD4
c                         (TC842 ver.)
c-----------------------------------------------------------------------
c  ID=105 ELC plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   UTcal DTcal Uq    DQ
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             Tavr        Qavr        Qcor
c             
c  W# :  21   22    23    24    25    26    27    28    29    30
c             RTOF(FH10)  RTOF(R3)    RTOFcal(10) RTOFcal(R3)
c             
c-----------------------------------------------------------------------
c  ID=106 S0 plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   LTcal RTcal       LQcal RQcal       LT(R3)RT(R3)
c
c  W# :  11   12    13    14    15    16    17    18    19    20
c             Tavr              Qavr              Lavr(R3)
c                
c  W# :  21   22    23    24    25    26    27    28    29    30
c             Tcal                                Lcal(R3)
c                
c-----------------------------------------------------------------------
c  ID=107 NaI Q
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   Qraw        Qcal
c
c-----------------------------------------------------------------------
c  ID=108 Schottky
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   Trig.
c
c-----------------------------------------------------------------------
c  ID=301 F2 plastic
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F2LT  F2RT        F2LQ  F2RQ  
c             (ch)  (ch)        (ch)  (ch)  
c  W# :  11   12    13    14    15    16    17    18    19    20
c             F2LT  F2RT        F2LQ  F2RQ
c             (ns)  (ns)        (ch)  (ch)
c  W# :  21   22    23    24    25    26    27    28    29    30
c             F2Tavr            F2Qavr
c
c-----------------------------------------------------------------------
c  ID=302 F3 plastic
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F3LT  F3RT        F3LQ  F3RQ        F3LTC F3RTC
c             (ch)  (ch)        (ch)  (ch)        (ch)  (ch)
c  W# :  11   12    13    14    15    16    17    18    19    20
c             F3LT  F3RT        F3LQ  F3RQ        F3LTC F3RTC
c             (ns)  (ns)        (ch)  (ch)        (ns)  (ns)
c  W# :  21   22    23    24    25    26    27    28    29    30
c             F3Tavr            F3Qavr            F3TCavr
c             (ns)              (ch)              (ns)
c  W# :  31   32    33    34    35    36    37    38    39    40
c             TOF23       TOFcorF3Qcor            TOF23C      TOF23Ccor 
c
c-----------------------------------------------------------------------
c  encbeam.f
c
c
c
c
c     ihit_min                  : the number of required data (ihit_min0)
c     hitdet(1:nhitdet)         : hit detector id
c     nhitdata(1:ndet)          : the number of data for each id 
c     nhitdet                   : the number of hit detector
c     ndet                      : the number of total detectors
c     ndata                     : the length of rawdata array
      IMPLICIT NONE
      INCLUDE 'analyslogic.fh'
      INCLUDE 'commonprm.fh'

 0    INTEGER nx, ny, ndet, ndata, naok
      INTEGER rawdata(ndata,ndet)
      INTEGER nhitdata(ndet)
      INTEGER hitdet(ndet)
      INTEGER nhitdet
      INTEGER i,j

      REAL    val(nx,ny)

ccccc F3 ccccccc
      REAL f3eX, f3eY
      REAL f3dX, f3dY
      REAL f3dZ, f3Z

ccccc F6 ccccccc
      REAL f5eX, f5eY
      REAL f5dX, f5dY
      REAL f5dZ, f5Z

ccccc S0 ccccccc
c      REAL s0eX, s0eY
      REAL s0dX, s0dY
      REAL s0dZx, s0dZy
      REAL s0Zx, s0Zy
c
      REAL S0TL, S0TR
      REAL S0QL, S0QR

cccc dispersion correction factor cccc
      REAL cf46dp
      REAL cf6il1dp
      REAL cf6il2dp

ccc For R3 part

      REAL RMD1
      REAL RMD2
      REAL RMD3
      REAL RMD4
      REAL RMD5

      REAL ILC2
      REAL NaI_Q1
      REAL NaI_Q2

      REAL F2LT
      REAL F2RT
      REAL F2LQ
      REAL F2RQ

      REAL F2Tavr

      REAL F3LT
      REAL F3RT
      REAL F3LQ
      REAL F3RQ

      REAL F3LTC
      REAL F3RTC

      REAL F3Tavr
      REAL F3TCavr

      REAL S0TOF_FH10
      REAL S0TOF_R3

ccc local variables ccc                                                                  
      real mu,clight,qe,amu,lpass,lring
      real brho_36,brho_R3,TOF_3S,TOF_R3,aoq_36,aoq_R3,Z_F3,Z_ELC
      real Tsum_F3,Tsum_S0,Tsum_ELC,tof_offset_3S,tof_offset_R3
      real Qsum_F3,Qsum_S0,Qsum_ELC,Z_F3_raw,Z_ELC_raw,beta,gamma,dp
      real moq_0,moq_1,T_0,T_1,tof_0,Turn_raw,Turn_0,dpop,Tsum_F2
      real Turn_1,TOF_3SR,beta_R,brho_R,gamma_R,TOF_beta,beta_F6r
      real moq_2,TIF_R3,c,brho_f6,moq_3,TOF_R3_0,TOF_R3_1,moq_00
      real TOF_23,Qsum_F2,tof_offset_23,beta_calib,Tsum_FH10
      real TOF_23_raw,TOF_3S_raw,TOF_R3_raw,moq_b,moq_4,xx
      real beta_ratio,beta_calib2,beta_ratio2
      INTEGER evtnumber
ccccccccccccccccccccccc                  

c      CHARACTER ihitchara*4
      INTEGER ihit_min0

      INTEGER ihit, id

      SAVE ihit_min0
      SAVE evtnumber

      IF (InitENCflag(2)) THEN
c         CALL LOADDALIPRM
c         CALL GETENV('IHIT_MIN0',ihitchara)
c         READ(ihitchara,*) ihit_min0
c         IF (ihit_min0.EQ.1) THEN
c            WRITE(*,*) 'ENCDALI-M : IHIT_MIN0 = ', ihit_min0
c         ELSE
c            ihit_min0 = 3
c            WRITE(*,*) 'ENCDALI-M : IHIT_MIN0 -> ', ihit_min0
c         ENDIF
         InitENCFlag(2) = .FALSE.
         evtnumber = 0
c         DO i=1,186,1
c            write(*,*)i,ped(i)+gain(i)*3840.
c            write(*,*)i,(100.-ped(i))/gain(i),
c     &           int((100.-ped(i))/gain(i)/2.)
c            write(*,*)i, int((100.-ped(i))/gain(i)/2.)
c         ENDDO
      ENDIF

      S0TL = -1000.
      S0TR = -1000.

      naok = 0

      DO ihit = 1,nhitdet
         id = hitdet(ihit)

c         IF (id.GT.nch_dali) CYCLE
c         IF (nhitdata(id).LT.ihit_min0) CYCLE ! ihit_min
c         IF (rawdata(3,id).GT.1) CYCLE

c         write(*,*) fh10tref

         naok = naok + 1

         val(1,naok) = id

c         write(*,*) id, naok

         val(2,naok) = rawdata(1,id)
         val(3,naok) = rawdata(2,id)

cccccccc mistake in map alignment cccccccc

         RMD1 = rawdata(1,10)
         RMD2 = rawdata(2,10)
         RMD3 = rawdata(1,11)
         RMD4 = rawdata(2,11)
         RMD4 = rawdata(1,24)

         NaI_Q1 = rawdata(1,12)
         NaI_Q2 = rawdata(2,12)

         F2LQ = rawdata(1,31)
         F2LT = rawdata(2,31)
         F2RQ = rawdata(1,32)
         F2RT = rawdata(2,32)

         F3LQ = rawdata(1,33)
         F3LT = rawdata(2,33)
         F3RQ = rawdata(1,34)
         F3RT = rawdata(2,34)

         F3LTC = rawdata(1,35)
         F3RTC = rawdata(2,35)


c         DCCT = rawdata(2,28)
ccccccccccccccccccccccccccccccccccccccc

ccccccc for S0 Pla Tsuika ccccccc
         
         do j=26,27
            if(id.eq.j) then
            do i=1,2
c      write(*,*) rawdata(1,26), rawdata(2,26), fh10tref
c      write(*,*) rawdata(1,27), rawdata(2,27), fh10tref
               rawdata(i,j) = rawdata(i,j) - fh10tref + 20000.
c      write(*,*) rawdata(1,26), rawdata(2,26), fh10tref
c      write(*,*) rawdata(1,27), rawdata(2,27), fh10tref
            enddo
         endif
         enddo
         S0TL = rawdata(1,26)
         S0TR = rawdata(2,26)

         S0QL = rawdata(1,27)
         S0QR = rawdata(2,27)


c        write(*,*) S0QL, S0QR
ccccccccccccccccccccccccccccccccc

         
c         IF (val(2,naok).GT.0. .AND. val(3,naok).GT.0.) THEN
c            mult = mult + 1
c         ENDIF

      ENDDO

c      ID=26
c      naok = naok +1
c      val(1,naok) = id
c      val(2,naok) = fh10tref

cccccccccccccccccccccccccccccccccc
ccc all ppac from encppac Part ccc
cccccccccccccccccccccccccccccccccc
      ID = 90
      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = f31ax * (-1.)
      val(3,naok) = f31ay
      val(4,naok) = f31bx * (-1.)
      val(5,naok) = f31by

      val(6,naok) = f32ax * (-1.)
      val(7,naok) = f32ay
      val(8,naok) = f32bx * (-1.)
      val(9,naok) = f32by

      val(12,naok) = f51ax * (-1.)
      val(13,naok) = f51ay
      val(14,naok) = f51bx * (-1.)
      val(15,naok) = f51by

      val(16,naok) = f52ax * (-1.)
      val(17,naok) = f52ay
      val(18,naok) = f52bx * (-1.)
      val(19,naok) = f52by

      val(22,naok) = f4x
c      val(22,naok) = f4x * (-1.)
      val(23,naok) = f4y

      val(24,naok) = f6x
c      val(24,naok) = f6x * (-1.)
      val(25,naok) = f6y

      val(32,naok) = fh9x - 0.04602
      val(33,naok) = fh9y - 0.5314

      val(34,naok) = fh10x * (-1.) - 0.3947
      val(35,naok) = fh10y + 0.08362

      val(36,naok) = s0x * (-1.) - 0.6050
      val(37,naok) = s0y - 0.4585

c      val(42,naok) = f21x * (-1.)
c      val(43,naok) = f21y
c      val(44,naok) = f22x * (-1.)
c      val(45,naok) = f22y

cccccccccccccc
ccc for F2 ccc
cccccccccccccc

      ID = 301

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = F2LT  ! F2 Pla LT raw
      val(3,naok) = F2RT  ! F2 Pla RT raw
      val(5,naok) = F2LQ  ! F2 Pla LQ raw
      val(6,naok) = F2RQ  ! F2 Pla RQ raw

      val(12,naok) = val(2,naok) * 0.0268680 ! F2 PLa LT cal
      val(13,naok) = val(3,naok) * 0.0271937 ! F2 PLa RT cal
      val(15,naok) = val(5,naok) - 0.0000 ! F2 PLa LQ ped
      val(16,naok) = val(6,naok) - 0.0000 ! F2 Pla RQ ped

      val(22,naok) = (val(12,naok) + val(13,naok)) * 0.5 ! F2 Tavr
      val(25,naok) = (val(15,naok) * val(16,naok)) **0.5 ! F2 Qavr

      F2Tavr = val(22,naok)

      val(50,naok) = 0.5*(val(5,naok)+val(6,naok))
      Tsum_F2 = val(22,naok)
      Qsum_F2 = val(50,naok)

cccccccccccccc
ccc for F3 ccc
cccccccccccccc

      ID = 302

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = F3LT  ! F3 Pla LT raw
      val(3,naok) = F3RT  ! F3 Pla RT raw
      val(5,naok) = F3LQ  ! F3 Pla LQ raw
      val(6,naok) = F3RQ  ! F3 Pla RQ raw
      val(8,naok) = F3LTC ! F3 Pla LT(CFD) raw 
      val(9,naok) = F3RTC ! F3 Pla RT(CFD) raw 

      val(12,naok) = val(2,naok) * 0.0271382 ! F3 PLa LT cal
      val(13,naok) = val(3,naok) * 0.0276006 ! F3 PLa RT cal
      val(15,naok) = val(5,naok) - 0.0000 ! F3 PLa LQ ped
      val(16,naok) = val(6,naok) - 0.0000 ! F3 Pla RQ ped
      val(18,naok) = val(8,naok) * 0.0257660 ! F3 PLa LTC cal
      val(19,naok) = val(9,naok) * 0.0275440 ! F3 Pla RTC cal

      val(22,naok) = (val(12,naok) + val(13,naok)) * 0.5 ! F3 Tavr
      val(25,naok) = (val(15,naok) * val(16,naok)) **0.5 ! F3 Qavr
      val(28,naok) = (val(18,naok) + val(19,naok)) * 0.5 ! F3 TCavr

      F3Tavr = val(22,naok)

      val(32,naok) = val(22,naok) - F2Tavr ! TOF23
      val(38,naok) = val(28,naok) - F2Tavr ! TOF23(CFD)

      val(34,naok) = val(32,naok) * 1.0000 ! TOF23 cor. by F6X <-- input
      val(35,naok) = val(25,naok) * 1.0000 ! F3Q cor. by F6X <-- input
      val(40,naok) = val(38,naok) * 1.0000 ! TOF23(CFD) cor. by F6X <-- input

      val(50,naok) = 0.5*(val(15,naok)+val(16,naok))
      Tsum_F3 = val(22,naok)
      Qsum_F3 = val(50,naok)


ccccccccccccccccccc
ccc for F3 part ccc
ccccccccccccccccccc

      ID = 201

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = (f31ax + f31bx) * 0.5  ! F3PPAC1 X
      val(3,naok) = (f31ay + f31by) * 0.5  ! F3PPAC1 Y
      val(4,naok) = (f32ax + f32bx) * 0.5  ! F3PPAC2 X
      val(5,naok) = (f32ay + f32by) * 0.5  ! F3PPAC2 Y
 
      f3eX = val(4,naok)
      f3eY = val(5,naok)
      f3dX = val(4,naok) - val(2,naok)
      f3dY = val(5,naok) - val(3,naok)
      f3dZ = 890.
      f3Z = -890. ! Distance from PPAC2 to F3-focal plane

      val(12,naok) = f3eX + (f3dX/f3dZ)*f3Z ! X @ F3
      val(14,naok) = f3eY + (f3dY/f3dZ)*f3Z ! Y @ F3

      val(22,naok) = atan(f3dX/f3dZ)*1000. ! angle for X @ F3
      val(24,naok) = atan(f3dY/f3dZ)*1000. ! angle for Y @ F3


ccccccccccccccccccc
ccc for F5 part ccc
ccccccccccccccccccc

      ID = 202

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = (f51ax + f51bx) * 0.5  ! F5PPAC1 X
      val(3,naok) = (f51ay + f51by) * 0.5  ! F5PPAC1 Y
      val(4,naok) = (f52ax + f52bx) * 0.5  ! F5PPAC2 X
      val(5,naok) = (f52ay + f52by) * 0.5  ! F5PPAC2 Y
 
      f5eX = val(4,naok)
      f5eY = val(5,naok)
      f5dX = val(4,naok) - val(2,naok)
      f5dY = val(5,naok) - val(3,naok)
      f5dZ = 650.
      f5Z = -250. ! Distance from PPAC2 to F5-focal plane for X

      val(12,naok) = f5eX + (f5dX/f5dZ)*f5Z ! X @ F5
      val(14,naok) = f5eY + (f5dY/f5dZ)*f5Z ! Y @ F5

      val(22,naok) = atan(f5dX/f5dZ)*1000. ! angle for X @ F5
      val(24,naok) = atan(f5dY/f5dZ)*1000. ! angle for Y @ F5

ccccccccccccccccccc
ccc for S0 part ccc
ccccccccccccccccccc

      ID = 203

      naok = naok + 1
      val(1,naok) = id

      s0dX = s0x - fh10x
      s0dY = s0y - fh10y
      s0dZx = 500.
      s0dZy = 500.
      s0Zx =  1321.7 ! Distance from PPAC2 to s0-focal plane for X
      s0Zy =  1313.1 ! Distance from PPAC2 to s0-focal plane for Y

      val(12,naok) = s0X + (s0dX/s0dZx)*s0Zx ! X @ S0
      val(14,naok) = s0Y + (s0dY/s0dZy)*s0Zy ! Y @ S0

      val(22,naok) = atan(s0dX/s0dZx)*1000. ! angle for X @ S0
      val(24,naok) = atan(s0dY/s0dZy)*1000. ! angle for Y @ S0


cccccccccccccccccccccccccccccccccccccc
ccc for dispersion correction part ccc
cccccccccccccccccccccccccccccccccccccc

      ID = 204

      naok = naok + 1
      val(1,naok) = id

      cf46dp = -3.9000 ! (-75.3 / 18.3 mm) <-- input fitting result F4 X vs F6 X
c      cf46dp = -4.0000 ! (-75.3 / 18.3 mm) <-- input fitting result F4 X vs F6 X
      val(2,naok) = f6x - cf46dp * f4x

      cf6il1dp = -0.6000 ! (-47 / 75.3 mm) <-- input fitting result F6 X vs ILC1 X
      val(4,naok) = ilc1x - cf6il1dp * f6x

      cf6il2dp = 0.9500 ! (-72.5 / 75.3 mm) <-- input fitting result F4 X vs ILC2 X
      val(6,naok) = ilc2x - cf6il2dp * f6x

cccccccccccccc
ccc S0 Pla ccc
cccccccccccccc

      ID = 106

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = S0TL * 0.09766 ! ch2ns
      val(3,naok) = S0TR * 0.09766 ! ch2ns

      val(5,naok) = S0QL * 0.09766 ! ch2ns
      val(6,naok) = S0QR * 0.09766 ! ch2ns

      val(12,naok) = (val(2,naok) + val(3,naok)) * 0.5 ! S0 T Avr. @ FH10

      val(8,naok) = FH10_LT ! in nsec DAQ at B3F
      val(9,naok) = FH10_RT ! in nsec DAQ at B3F

      val(18,naok) = (val(8,naok) + val(9,naok)) * 0.5 ! S0 Tavr @ B3F 

      val(22,naok) = val(12,naok) + 0.0000 ! TOF(F3-S0) correction @ FH10

      val(28,naok) = val(18,naok) + 0.0000 ! TOF(F3-S0) correction @ B3F

      S0TOF_FH10 = val(22,naok)
      S0TOF_R3 = val(28,naok)


c      val(50,naok) = 0.5*(S0QL+S0QR)
      val(50,naok) = 0.5*(val(5,naok)+val(6,naok))
      Tsum_S0 = 0.5*(val(2,naok)+val(3,naok))
      Tsum_FH10 = 0.5*(val(8,naok)+val(9,naok))
      Qsum_S0 = val(50,naok)

ccccccccccccccccccccc   R3 Part  cccccccccccccccccccc

cccccccccccccccccccccccc
ccc ILC1 PPAC Calib. ccc
cccccccccccccccccccccccc

      ID = 101

      naok = naok + 1
      val(1,naok) = id
c      val(2,naok) = val(2,1) - 0.0 ! X1 - ped.
c      val(3,naok) = val(3,1) - 0.0 ! X2 - ped.
c      val(4,naok) = val(2,2) - 0.0 ! Y1 - ped.
c      val(5,naok) = val(3,2) - 0.0 ! Y2 - ped.
      val(2,naok) = val(2,1) - 50.63 ! X1 - ped.
      val(3,naok) = val(3,1) - 47.56 ! X2 - ped.
      val(4,naok) = val(2,2) - 48.50 ! Y1 - ped.
      val(5,naok) = val(3,2) - 54.71 ! Y2 - ped.


      val(12,naok) = val(2,naok) - val(3,naok)
      val(13,naok) = val(2,naok) + val(3,naok)
      val(14,naok) = val(4,naok) - val(5,naok)
      val(15,naok) = val(4,naok) + val(5,naok)

      val(17,naok) = val(12,naok) / val(13,naok)
      val(18,naok) = val(14,naok) / val(15,naok)

      val(20,naok)=val(2,naok)/val(13,naok)
      val(21,naok)=val(4,naok)/val(15,naok)

      val(22,naok) = val(17,naok) * (-26.866) + 0.4579 ! ILC1 PPAC X calib
      val(24,naok) = val(18,naok) * (-27.019) + 0.4142 ! ILC1 PPAC Y calib

      ilc1x = val(24,naok)
      ilc1y = val(22,naok)

      val(26,naok) = ilc1x + 3.247 ! ILC1 PPAC X + measure offset
      val(28,naok) = ilc1y + 5.320 ! ILC1 PPAC Y + measure offset 

cccccccccccccccccccccccc
ccc ILC2 PPAC Calib. ccc
cccccccccccccccccccccccc

      ID = 102

      naok = naok + 1
      val(1,naok) = id
c      val(2,naok) = val(2,3) - 0.0 ! X1 - ped.
c      val(3,naok) = val(3,3) - 0.0 ! X2 - ped.
c      val(4,naok) = val(2,4) - 0.0 ! Y1 - ped.
c      val(5,naok) = val(3,4) - 0.0 ! Y2 - ped.
      val(2,naok) = val(2,3) - 53.49 ! X1 - ped.
      val(3,naok) = val(3,3) - 52.79 ! X2 - ped.
      val(4,naok) = val(2,4) - 49.28 ! Y1 - ped.
      val(5,naok) = val(3,4) - 52.51 ! Y2 - ped.

      val(12,naok) = val(2,naok) - val(3,naok)
      val(13,naok) = val(2,naok) + val(3,naok)
      val(14,naok) = val(4,naok) - val(5,naok)
      val(15,naok) = val(4,naok) + val(5,naok)

      val(17,naok) = val(12,naok) / val(13,naok)
      val(18,naok) = val(14,naok) / val(15,naok)
      
      val(22,naok) = (val(17,naok) * (-52.445) + 0.3874) * (-1.) ! ILC2 PPAC X calib
      val(24,naok) = val(18,naok) * (-51.954) + 0.0237 ! ILC2 PPAC Y calib

      ilc2x = val(22,naok) ! ILC2 PPAC X
      ilc2y = val(24,naok) ! ILC2 PPAC Y


c test
      val(50,naok) = val(3,naok)/val(15,naok)
      val(51,naok) = val(3,naok)/val(15,naok)*(-52.445+0.3874)*(-1)

c      write(*,*) ilc2x

cccccccccccccccccc
ccc ILC2 T-Pla ccc
cccccccccccccccccc

      ID = 103

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = val(3,5) * 0.02541 ! L ch2ns  TDC1 0ch
      val(3,naok) = val(3,6) * 0.02560 ! RU ch2ns  TDC1 1ch
      val(4,naok) = val(3,7) * 0.02524 ! RD ch2ns  TDC1 4ch

      val(6,naok) = (val(3,naok) + val(4,naok)) * 0.5 ! R side avr. 
      val(8,naok) = (val(2,naok) + val(6,naok)) * 0.5 ! T avr.

      val(10,naok) = val(2,naok)-val(6,naok) !Tdiff X
      val(11,naok) = val(3,naok)-val(4,naok) !Tdiff Y
      val(12,naok) = 222.9-21.55*val(10,naok) !X from dT
      val(13,naok) = -17.93+12.89*val(11,naok) !Y from dT


      ILC2 = val(8,naok)

cccccccccccccccc
ccc R-MD Pla ccc
cccccccccccccccc

      ID = 104

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = RMD1 * 0.02541 ! R-MD1 ch2ns  TDC2 0ch
      val(3,naok) = RMD2 * 0.02529 ! R-MD2 ch2ns  TDC2 1ch
      val(4,naok) = RMD3 * 0.02545 ! R-MD3 ch2ns  TDC2 2ch
      val(5,naok) = RMD4 * 0.02551 ! R-MD4 ch2ns  TDC2 3ch
      val(6,naok) = RMD5 * 0.02567 ! R-MD5 ch2ns TDC2 4ch

      val(14,naok) = RMD3_T ! in nsec from TC842
      val(15,naok) = RMD4_T ! in nsec from TC842

ccccccccccccccc
ccc ELC Pla ccc
ccccccccccccccc

      ID = 105

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = ELC_UT ! TU in nsec
      val(3,naok) = ELC_DT ! TD in nsec

      val(4,naok) = val(2,8) ! QUraw 
      val(5,naok) = val(3,8) ! QDraw

      val(12,naok) = (val(2,naok) + val(3,naok))*0.5 ! ELC_Tavr
      val(14,naok) = (val(4,naok) * val(5,naok))**0.5 ! ELC_Qavr

      val(16,naok) = 1.0000 * val(14,naok) ! ELC_Q correction by F6 X

      val(22,naok) = val(12,naok) - S0TOF_FH10 ! Ring TOF use FH10 data
      val(24,naok) = val(12,naok) - S0TOF_R3 ! Ring TOF use R3 data

      val(26,naok) = val(22,naok) + 0.0000 ! Ring TOF cal use FH10 data
      val(28,naok) = val(24,naok) + 0.0000 ! Ring TOF cal use R3 data


      val(50,naok) = 0.5*(val(4,naok)+val(5,naok))
      Tsum_ELC = val(12,naok)
      Qsum_ELC = val(50,naok)

cccccccccccc
ccc NaI  ccc
cccccccccccc
c
      ID = 107

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = NaI_Q1 ! QDC3(HOSHIN C009) ch1
c      val(3,naok) = NaI_Q2 ! QDC3(HOSHIN C009) ch2
c
      val(4,naok) = val(2,naok) * 1.000 -0.00 


cccccccccccccccc
ccc Schottky ccc
cccccccccccccccc

      ID = 108

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = Schottky ! RSA Trig. in ns


c ======================================
c
C PID
c
c ======================================
       mu = 1.6605402e-27 ! kg
       clight = 299.792458   ! Mm/s
       qe = 1.60217733E-19   ! C
c        amu = mu * pow(c,2) / e * 1.e-6;  
       amu = 931.4940023 ! MeV/c2
       lpass = 84.244 ! F3 to S0 [m]
       lring = 60.3507  ! nominal curcumference of ring [m]
       brho_36 = 3.8874 ! F3 to S0 [Tm] (center)
c       brho_36 = 3.8891 ! F3 to S0 [Tm] (center)
       brho_R3 = 3.9734 ! R3 main [Tm]

c       moq_0 = 1.998372691 ! 38K [u/C]
c       moq_0 = 1.99819695028 ! 36Ar [u/C]
       moq_0 = 1.9976484 ! 36Ar [u/C] (for bare 36Ar18+)
       moq_00 = 2.0564427 ! 35Cl [u/C] (for bare 35Cl17+)
c      tof_0 = 700.e+3 ! 38K TOF in RING [ns]
c      T_0 = 372. ! 38K [ns]

       tof_offset_23 = 24.45 + 54.15 ! TOF offset F2 to F3
       tof_offset_3S = -654.1 + 537.1 ! TOF offset F3 to S0
c       tof_offset_R3 =  6100.0 ! TOF offset S0 to ELC
       tof_offset_R3 =  511.0 ! TOF offset S0 to ELC

       TOF_23_raw = Tsum_F3 - Tsum_F2 ! TOF from F2 to F3
       TOF_3S_raw = Tsum_S0 - Tsum_F3  ! TOF from F3 to S0
       TOF_R3_raw = Tsum_ELC - Tsum_S0  ! TOF form S0 to ELC
       TOF_23 = Tsum_F3 - Tsum_F2 + tof_offset_23 ! TOF from F2 to F3
       TOF_3S = Tsum_S0 - Tsum_F3 + tof_offset_3S ! TOF from F3 to S0
       TOF_R3 = Tsum_ELC - Tsum_FH10 + tof_offset_R3 ! TOF form S0 to ELC

c       dp = 1.+(f6x/98./100.) ! F6 dis = + 98mm/%
c       dpop = f6x/75. !dp/p [%]
c      dispersin = 97.65 at F6 gfor MS02
       dpop = f6x/97.65 !dp/p [%] F6 disp=98.0 mm/%
       brho_f6 = brho_36*(1.+dpop/100.)
       xx = brho_f6/(moq_0*amu/clight)
       beta_F6r = xx/sqrt(1.+xx**2.)
c     correction for detector thickness in S0 by LISE(35Cl)
c      0.99680 for 3.8874 Tm for 35Cl and 36Ar
       brho_R = 0.99680*brho_f6

c       beta_calib = 1.0162 -0.00092082*TOF_3S
c       beta_calib = 1.0453 -0.00095883*TOF_3S
c       beta_calib = 281.87/(TOF_3S-5.5951)
c       beta_calib = 279.0/(TOF_3S-9.0487)
c      below case for 3.8874 Tm with 4 nuclei
c      a = 281.16, b = 123.90
       beta_calib = 281.16/(TOF_3S_raw - 123.90)
c     below parameter by Yoshi
c       beta_calib = 279.41/(TOF_3S_raw - 127.14)
c      below case for 3.8891 Tm
c       beta_calib = 282.31/(TOF_3S_raw - 121.85)
       beta_ratio = beta_F6r/beta_calib
c      0.000016244 for 3.8874 Tm and in 35Cl 
c       beta_calib2 = beta_calib*(1. - f6x*0.000016244)
c      below case for 3.8891 Tm and in 35Cl 
c       beta_calib2 = beta_calib*(1. - f6x*0.000016376)
c      0.000018771 for 3.8874 Tm and in 36Ar
       beta_calib2 = beta_calib*(1. - f6x*0.000018771)
c      below case for 3.8891 Tm and in 36Ar
c       beta_calib2 = beta_calib*(1. - f6x*0.000015464)
       beta_ratio2 = beta_F6r/beta_calib2

       beta = lpass/TOF_3S/clight*1.e+3
       gamma = 1./sqrt(1.-beta**2.)
       aoq_36 = clight*brho_36*dp/amu/beta/gamma -0.038
cc      aoq_R3 = clight*brho_R3*dp/amu/beta/gamma

       Z_F3_raw = beta*sqrt(Qsum_F3) ! proportional to Z
       Z_ELC_raw = beta*sqrt(Qsum_ELC) ! proportional to Z
       Z_F3 = 0.58684 + 1.1585*Z_F3_raw ! to Z
       Z_ELC = Z_ELC_raw*1.0 ! to Z

c     for #0 nuclei     
       T_0 = TOF_R3/1892.7
c     for #1 nuclei                                                     
       T_1 = TOF_R3/1853.7
c     correction for detector thickness in S0 by LISE(35Cl)
c       TOF_3SR = 0.5192515/0.5177481*TOF_3S
c       beta_R = 0.997265*beta_calib *1.000067
c      below case for 3.8874 Tm 
c     0.99892 for 35Cl, 0.99894 for 36Ar
       beta_R = 0.99767*beta_calib2
c     below para,eter by Yoshi
c       beta_R = 0.99891*beta_calib2
       gamma_R = 1./sqrt(1.-beta_R**2.)
c       beta_R = beta_calib
c       TOF_R3_0 = TOF_R3 - (TOF_beta - 372650.)/beta_R
       TOF_R3_1 = TOF_R3 + 39*lring/beta_R/clight*1000.

c       Turn_raw = TOF_R3*beta_R*clight*gamma_R/lring
       Turn_raw = TOF_R3*beta_R*clight*gamma_R/lring/1000.
c                                ! the number of turn in R3
c       N_turn = TOF_R3/TOF_3S ! proportional to the number of turn
       TOF_beta = TOF_R3*beta_R
       Turn_0 = T_0*beta_R ! proportional to flight length for #0 
       Turn_1 = T_1*beta_R ! proportional to flight length for #1
       TOF_R3_0 = TOF_R3 - (TOF_beta - 372650.)/beta_R
c       TOF_R3_0 = 372650./beta_R

c      moq_1=(moq_0)*(T_1/T_0)*sqrt((1.-beta_R**2)/
c     &  (1.-(T_1*beta_R/T_0)**2.)) 
      moq_1=(moq_0)*(T_1/378.3178)*sqrt((1.-beta_R**2)/
     &  (1.-(T_1*beta_R/378.3176)**2.)) 
c     701288.6 for 35Cl, 731210.5 for 36Ar
      moq_4=(moq_0)*(TOF_R3/701288.6)
     &   *sqrt(((1.-(701288.6/TOF_R3)**2.)
     &   /(moq_0*mu/qe/brho_R*clight*1.e+6)**2)+1.)
c
      moq_b=(moq_0)*(0.53207/beta_R)*sqrt((1.-beta_R**2)/
     &  (1.-(0.53207*beta_R/beta_R)**2.)) 
c
      moq_2=(moq_0)*(TOF_R3/701288.6)*sqrt((1.-beta_R**2)/
     &  (1.-(TOF_R3*beta_R/701288.6)**2.)) 
c     701292.22 by Yoshi
c      moq_2=(moq_0)*(TOF_R3/701292.22)*sqrt((1.-beta_R**2)/
c     &  (1.-(TOF_R3*beta_R/701292.22)**2.)) 
c
c       moq_1=(moq_0)*(TOF_R3/tof_0)
c     &   *sqrt((1.-beta_calib**2)/(1.-(TOF_R3*beta_calib/tof_0)**2.)) 
c      
       moq_3=(moq_0)*(T_1/378.3176)
     &   *sqrt(((1.-(378.3176/T_1)**2.)
     &   /(moq_0*mu/qe/brho_R*clight*1.e+6)**2)+1.)

       evtnumber = evtnumber + 1

c put ito vals
       id = 500
       naok = naok + 1
       val(1,naok) = id
       val(2,naok) = TOF_23
       val(3,naok) = TOF_3S
       val(4,naok) = TOF_R3 !nano second
       val(5,naok) = TOF_R3_0 !nano second
       val(6,naok) = TOF_R3_1
       val(7,naok) = TOF_3S_raw
       val(8,naok) = TOF_R3_raw
      

       id = 501
       naok = naok + 1
       val(1,naok) = id
       val(2,naok) = beta
       val(3,naok) = gamma
       val(4,naok) = aoq_36
       val(5,naok) = dpop
       val(6,naok) = beta_calib
       val(7,naok) = TOF_beta
       val(8,naok) = beta_F6r
       val(9,naok) = beta_ratio
       val(10,naok) = beta_calib2
       val(11,naok) = beta_ratio2


       id = 502
       naok = naok + 1
       val(1,naok) = id
       val(2,naok) = Z_F3_raw
       val(3,naok) = Z_F3

       id = 503
       naok = naok + 1
       val(1,naok) = id
       val(2,naok) = Z_ELC_raw
       val(3,naok) = Z_ELC

       id = 600
       naok = naok + 1
       val(1,naok) = id
       val(2,naok) = Turn_raw
       val(3,naok) = Turn_0
       val(4,naok) = Turn_1
       val(5,naok) = T_0
       val(6,naok) = T_1
       val(7,naok) = moq_1
       val(8,naok) = moq_2
       val(9,naok) = beta_R
       val(10,naok) = brho_f6
       val(11,naok) = brho_R
       val(12,naok) = moq_3
       val(13,naok) = moq_b
       val(14,naok) = moq_4
       val(15,naok) = evtnumber

      RETURN
      END

c ======================================

