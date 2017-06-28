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
c  ID=26 S0 Plastic T (L,R) (L means view from upstream)
c  ID=27 NaI raw
c  ID=31 F2 Pla L
c  ID=32 F2 Pla R
c  ID=33 F3 Pla L
c  ID=34 F3 PLa R
c  ID=28 DCCT
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
c        ID   TUraw TDraw QUraw QDraw
c
c-----------------------------------------------------------------------
c  ID=106 S0 plastic TOF
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   LTcal RTcal       LQcal RQcal       S0Tavr
c
c-----------------------------------------------------------------------
c  ID=107 NaI Q
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   Q
c
c-----------------------------------------------------------------------
c  ID=301 F2 plastic
c  W# :  1    2     3     4     5     6     7     8     9     10
c        ID   F2LT  F2RT  F2LQ  F2RQ  F3LT  F3RT  F3LQ  F3RQ
c             (ch)  (ch)  (ch)  (ch)  (ch)  (ch)  (Ch)  (ch)
c  W# :  11   12    13    14    15    16    17    18    19    20
c             F2LT  F2RT  F2LQ  F2RQ  F3LT  F3RT  F3LQ  F3RQ
c             (ns)  (ns)              (ns)  (ns)
c  W# :  21   22    23    24    25    26    27    28    29    30
c             F2Tavr      F2Qavr      F3Tavr      F3Qavr
c
c  W# :  31   32    33    34    35    36    37    38    39    40
c             TOF23
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

      REAL RMD1
      REAL RMD2
      REAL RMD3
      REAL RMD4
      REAL RMD5

      REAL ILC2
c      REAL ELC_UT, ELC_DT
      REAL NaI_Q1
      REAL NaI_Q2

      REAL TOF

      REAL F2LT
      REAL F2RT
      REAL F2LQ
      REAL F2RQ

      REAL F3LT
      REAL F3RT
      REAL F3LQ
      REAL F3RQ

      REAL DCCT

ccc local variables ccc
      real mu,clight,qe,amu,lpass
      real brho_36,brho_R3,TOF_3S,TOF_R3,aoq_36,aoq_R3,Z_F3,Z_ELC
      real Tsum_F3,Tsum_S0,Tsum_ELC,tof_offset_3S,tof_offset_R3
      real Qsum_F3,Qsum_ELC,Z_F3_raw,Z_ELC_raw,beta,gamma,dp
      real moq_0,moq_1,T_0,T_1,tof_0,N_turn
ccccccccccccccccccccccc

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
         F2RQ = rawdata(2,32)

         F3LQ = rawdata(1,33)
         F3LT = rawdata(2,33)
         F3RQ = rawdata(1,34)
         F3RQ = rawdata(2,34)

         DCCT = rawdata(2,28)
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

      val(32,naok) = fh9x
      val(33,naok) = fh9y

      val(34,naok) = fh10x * (-1.)
      val(35,naok) = fh10y

      val(36,naok) = s0x * (-1.)
      val(37,naok) = s0y

c      val(42,naok) = f21x * (-1.)
c      val(43,naok) = f21y
c      val(44,naok) = f22x * (-1.)
c      val(45,naok) = f22y

ccccccccccccccccc
ccc for F2-F3 ccc
ccccccccccccccccc

      ID = 301

      naok = naok + 1
      val(1,naok) = id

      val(2,naok) = F2LT  ! F2 Pla LT raw
      val(3,naok) = F2RT  ! F2 Pla RT raw
      val(4,naok) = F2LQ  ! F2 Pla LQ raw
      val(5,naok) = F2RQ  ! F2 Pla RQ raw

      val(6,naok) = F3LT  ! F3 Pla LT raw
      val(7,naok) = F3RT  ! F3 Pla RT raw
      val(8,naok) = F3LQ  ! F3 Pla LQ raw
      val(9,naok) = F3RQ  ! F3 Pla RQ raw

      val(12,naok) = val(2,naok) * 1.0000 ! F2 PLa LT cal
      val(13,naok) = val(3,naok) * 1.0000 ! F2 PLa RT cal
      val(14,naok) = val(4,naok) - 0.0000 ! F2 PLa LQ ped
      val(15,naok) = val(5,naok) - 0.0000 ! F2 Pla RQ ped

      val(16,naok) = val(6,naok) * 1.0000 ! F3 Pla LT cal
      val(17,naok) = val(7,naok) * 1.0000 ! F3 Pla RT cal
      val(18,naok) = val(8,naok) - 0.0000 ! F3 Pla LQ ped
      val(19,naok) = val(9,naok) - 0.0000 ! F3 Pla RQ ped

      val(22,naok) = (val(12,naok) + val(13,naok)) * 0.5 ! F2 Tavr
      val(24,naok) = (val(14,naok) * val(15,naok)) **0.5 ! F2 Qavr

      val(26,naok) = (val(16,naok) + val(17,naok)) * 0.5 ! F3 Tavr
      val(28,naok) = (val(18,naok) * val(19,naok)) **0.5 ! F3 Qavr
      val(29,naok) = (val(18,naok) + val(19,naok)) *0.5 ! F3 Qavr


      val(32,naok) = val(26,naok) - val(22,naok) + 0.000 ! TOF23 (F3-F2)

      Tsum_F3 = val(26,naok)
      Qsum_F3 = val(29,naok)

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

c      write(*,*) val(2,naok),val(3,naok),val(5,naok),val(6,naok)

      val(8,naok) = (val(2,naok) + val(3,naok)) * 0.5 ! S0 T Avr.

      Tsum_S0 = val(8,naok)

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

      val(26,naok) = ilc1x ! ILC1 PPAC X
      val(28,naok) = ilc1y ! ILC1 PPAC Y

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
      val(5,naok) = RMD4  * 0.02551     ! R-MD4 ch2ns  TDC2 3ch
      val(6,naok) = RMD5 * 0.02567 ! R-MD5 ch2ns TDC2 4ch

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


      val(20,naok) = 0.5*(val(2,naok) + val(3,naok))
      val(21,naok) = 0.5*(val(4,naok) + val(5,naok))
      
      Tsum_ELC = val(20,naok)
      Qsum_ELC = val(21,naok)

c      ELC = val(4,naok)

c      TOF = ELC - (ILC2 / 1000.) ! in us 

c      val(6,naok) = TOF

cccccccccccc
ccc NaI  ccc
cccccccccccc
c
      ID = 107

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = NaI_Q1 ! QDC3(HOSHIN C009) ch1
      val(3,naok) = NaI_Q2 ! QDC3(HOSHIN C009) ch2
c
c      NaI = val(2,naok) * 1.000 -0.00


cccccccccc
cc DCCT cc
cccccccccc

      ID = 108

      naok = naok + 1
      val(1,naok) = id
      val(2,naok) = DCCT

c ======================================
c
C PID
c
c ======================================
       mu = 1.66053892e-27 ! kg
       clight = 299.792458   ! Mm/s
       qe = 1.60217657E-19   ! C
c        amu = mu * pow(c,2) / e * 1.e-6;  
       amu = 931.4940023 ! MeV/c2
       lpass = 84.244 ! F3 to S0 [m]
       brho_36 = 3.9768 ! F3 to F6 [Tm]
       brho_R3 = 3.9734 ! R3 main [Tm]

       moq_0 = 1.998372691 ! 38K [u/C]
c      tof_0 = 700.e+3 ! 38K TOF in RING [ns]
c      T_0 = 372. ! 38K [ns]

c     tof_offset_3S = ! TOF offset F3 to S0
c     tof_offset_R3 = ! TOF offset S0 to ELC

       TOF_3S = Tsum_S0 - Tsum_F3 + tof_offset_3S ! TOF from F3 to S0
       TOF_R3 = Tsum_ELC - Tsum_S0 + tof_offset_R3 ! TOF form S0 to ELC

       dp = 1.+(f6x/75./100.) ! F6 dis = + 75mm/%

       beta = lpass/TOF_3S/clight*1.e+3
       gamma = 1./sqrt(1.-beta**2.)
       aoq_36 = clight*brho_36*dp/amu/beta/gamma
cc      aoq_R3 = clight*brho_R3*dp/amu/beta/gamma

       Z_F3_raw = beta*sqrt(Qsum_F3) ! proportional to Z
       Z_ELC_raw = beta*sqrt(Qsum_ELC) ! proportional to Z
       Z_F3 = Z_F3_raw*1.0 ! to Z
       Z_ELC = Z_ELC_raw*1.0 ! to Z

c      moq_1=(moq_0)*(T_1/T_0)*sqrt((1.-beta**2)/(1.-(T_1*beta/T_0)**2.)) 

       N_turn = TOF_R3/TOF_3S ! proportional to the number of turn

       moq_1=(moq_0)*(TOF_R3/tof_0)
     &   *sqrt((1.-beta**2)/(1.-(TOF_R3*beta/tof_0)**2.)) 

      RETURN
      END

