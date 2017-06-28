      SUBROUTINE EncBeam(val,nx,ny,naok,
     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata)
c-----------------------------------------------------------------------
c ANALYZER 2 : Encbeam
c-----------------------------------------------------------------------
c  for raw data histgrams
c  ID=1 ILC1 CD-PPAC X1, X2 raw
c  ID=2 ILC1 CD-PPAC Y1, Y2 raw
c  ID=3 ILC2 CD-PPAC X1, X2 raw
c  ID=4 ILC2 CD-PPAC Y1, Y2 raw
c  ID=5 ILC2 T-shape L-side raw  ("L means view from downstream")
c  ID=6 ILC2 T-shape R-side Up raw
c  ID=7 ILC2 T-shape R-side Down raw
c  ID=8 ELC Plastic Q (L,R) raw
c  ID=9 ELC Plastic T (L,R) raw
c  ID=10 R-MD1, 2 raw 
c  ID=11 R-MD3, 4 raw 
c  ID=26 FH10 Plastic T (L,R) (L means view from upstream)
c
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
c             X           Y           ILC1X       ILC1Y
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
c  W# :  11   12     13     14     15     16     17     18     19     20
c
c
c-----------------------------------------------------------------------
c  ID=104 R-MD plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   MD1   MD2   MD3   MD4   MD5
c
c-----------------------------------------------------------------------
c  ID=105 ELC plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   Traw        Tcal
c
c-----------------------------------------------------------------------
c  ID=106 S0 plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   LTcal RTcal L-R         S0Tavr
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
      INTEGER i
c      INTEGER i,j

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
      REAL S0TL, S0TR

cccc dispersion correction factor cccc
      REAL cf46dp
      REAL cf6il1dp
      REAL cf6il2dp

      REAL RMD4

      REAL ILC2
      REAL ELC

      REAL TOF

c      CHARACTER ihitchara*4
      INTEGER ihit_min0

      INTEGER ihit, id

      SAVE ihit_min0

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

c         write(*,*) id

c         IF (id.GT.nch_dali) CYCLE
c         IF (nhitdata(id).LT.ihit_min0) CYCLE ! ihit_min
c         IF (rawdata(3,id).GT.1) CYCLE

c         write(*,*) fh10tref

         naok = naok + 1

         val(1,naok) = id

c         write(*,*) id, naok

         val(2,naok) = rawdata(1,id)
         val(3,naok) = rawdata(2,id)

         RMD4 = rawdata(2,11)

ccccccc for S0 Pla Tsuika ccccccc
         
         if(id.eq.26) then
            do i=1,2
c              write(*,*) rawdata(1,26), rawdata(2,26), fh10tref
               rawdata(i,26) = rawdata(i,26) - fh10tref + 20000.
c              write(*,*) rawdata(1,26), rawdata(2,26), fh10tref
            enddo
         endif
         S0TL = rawdata(1,26)
         S0TR = rawdata(2,26)

c        write(*,*) S0TL, S0TR
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

      val(32,naok) = fh9x
      val(33,naok) = fh9y

      val(34,naok) = fh10x * (-1.)
      val(35,naok) = fh10y

      val(36,naok) = s0x
      val(37,naok) = s0y

      val(42,naok) = f21x * (-1.)
      val(43,naok) = f21y
      val(44,naok) = f22x * (-1.)
      val(45,naok) = f22y

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
      s0dZx = 1590.4
      s0dZy = 1573.2
      s0Zx = -221.7 ! Distance from PPAC2 to s0-focal plane for X
      s0Zy = -213.1 ! Distance from PPAC2 to s0-focal plane for Y

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

      val(4,naok) = val(2,naok) - val(3,naok) ! FH10 L - R for test

      val(6,naok) = (val(2,naok) + val(3,naok)) * 0.5 ! S0 T Avr.

ccccccccccccccccccccc   R3 Part  cccccccccccccccccccc

cccccccccccccccccccccccc
ccc ILC1 PPAC Calib. ccc
cccccccccccccccccccccccc

      ID = 101

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = val(2,1) - 50.58 ! X1 - ped.
      val(3,naok) = val(3,1) - 46.55 ! X2 - ped.
      val(4,naok) = val(2,2) - 48.50 ! Y1 - ped.
      val(5,naok) = val(3,2) - 54.38 ! Y2 - ped.

      val(12,naok) = val(2,naok) - val(3,naok)
      val(13,naok) = val(2,naok) + val(3,naok)
      val(14,naok) = val(4,naok) - val(5,naok)
      val(15,naok) = val(4,naok) + val(5,naok)

      val(17,naok) = val(12,naok) / val(13,naok)
      val(18,naok) = val(14,naok) / val(15,naok)

      val(20,naok)=val(2,naok)/val(13,naok)
      val(21,naok)=val(4,naok)/val(15,naok)

      val(22,naok) = val(17,naok) * (-26.799) + 0.4354 ! ILC1 PPAC X calib
      val(24,naok) = val(18,naok) * (-26.928) - 0.0865 ! ILC1 PPAC Y calib

      ilc1x = val(24,naok)
      ilc1y = val(22,naok)

      val(26,naok) = ilc1x ! ILC1 PPAC X
      val(28,naok) = ilc1y ! ILC1 PPAC Y

cccccccccccccccccccccccc
ccc ILC2 PPAC Calib. ccc
cccccccccccccccccccccccc

      ID = 102

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = val(2,3) - 53.50 ! X1 - ped.
      val(3,naok) = val(3,3) - 52.24 ! X2 - ped.
      val(4,naok) = val(2,4) - 49.32 ! Y1 - ped.
      val(5,naok) = val(3,4) - 51.83 ! Y2 - ped.

      val(12,naok) = val(2,naok) - val(3,naok)
      val(13,naok) = val(2,naok) + val(3,naok)
      val(14,naok) = val(4,naok) - val(5,naok)
      val(15,naok) = val(4,naok) + val(5,naok)

      val(17,naok) = val(12,naok) / val(13,naok)
      val(18,naok) = val(14,naok) / val(15,naok)
      
      val(22,naok) = (val(17,naok) * (-53.045) + 0.3107) * (-1.) ! ILC2 PPAC X calib
      val(24,naok) = val(18,naok) * (-52.269) + 0.0692 ! ILC2 PPAC Y calib

      ilc2x = val(22,naok) ! ILC2 PPAC X
      ilc2y = val(24,naok) ! ILC2 PPAC Y

cccccccccccccccccc
ccc ILC2 T-Pla ccc
cccccccccccccccccc

      ID = 103

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = val(3,5) * 0.02550 ! L ch2ns
      val(3,naok) = val(3,6) * 0.02568 ! RU ch2ns 
      val(4,naok) = val(3,7) * 0.02536 ! RD ch2ns

      val(6,naok) = (val(3,naok) + val(4,naok)) * 0.5 ! R side avr. 
      val(8,naok) = (val(2,naok) + val(6,naok)) * 0.5 ! T avr.

      ILC2 = val(8,naok)

cccccccccccccccc
ccc R-MD Pla ccc
cccccccccccccccc

      ID = 104

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = val(2,10) * 0.02553 ! R-MD1 ch2ns
      val(3,naok) = val(3,10) * 0.02538 ! R-MD2 ch2ns 
      val(4,naok) = val(2,11) * 0.02557 ! R-MD3 ch2ns
      val(5,naok) = RMD4  * 0.02563 ! R-MD4 ch2ns
      val(6,naok) = val(2,24) * 0.02565 ! R-MD5 ch2ns

ccccccccccccccc
ccc ELC Pla ccc
ccccccccccccccc

      ID = 105

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = rawdata(1,9)
      val(4,naok) = val(2,naok) * 0.2308 - 10.815 ! T-coin ch2us(1ms/4000ch)

      ELC = val(4,naok)

      TOF = ELC - (ILC2 / 1000.) ! in us 

      val(6,naok) = TOF

c      write(*,*) TOF
      RETURN
      END

c ======================================

