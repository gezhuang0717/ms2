c     subroutine : encoptics
c     data selection
c     
c   Input
c     data from enctraj
c      
c
c   Output

      subroutine encoptics(
     &     val,nx,ny,naok,
     &     analyzer)
c     
      implicit none
      include '../core/analyslogic.fh'
      include '../traj/trajectory.fh'
      include './optics.fh'
      integer analyzer
      integer nx,ny,naok
      real val(nx,ny)

      real x6c4,x6c7,x7c6,x9c7,xs0c9,xs0c7
      real xs2c7,xs2c9,xs2cs0,xs2c6
      real xs0ppc7,xs0ppc9
      real xs2c, xs2c3
      naok = 0

      if (initencflag(analyzer)) then
         write(6,10) 
     $        ' ANAPAW-M : [encoptics]          analys     =', analyzer
 10      format(1x,a,1x,i4)
         call loadopticsprm
         call srand(time())
         initencflag(analyzer) = .false. 
      endif

c F3 
      naok = naok + 1
      val(1,naok) = 100
c original
c      val(32,naok) = trj(X,F3DC)
c      val(33,naok) = trj(Y,F3DC)
c      val(34,naok) = trj(A,F3DC)
c      val(35,naok) = trj(B,F3DC)
c      val(36,naok) = trj(T,F3DC)
c for SHARAQ06
c      val(32,naok) = trj(X,F3TR)
c      val(33,naok) = trj(Y,F3TR)
c      val(34,naok) = trj(A,F3TR)
c      val(35,naok) = trj(B,F3TR)
c      val(36,naok) = trj(T,F3TR)
c for tetra12apr uHODO (z=-116mm)
      val(22,naok) = trj(X,F3PP)
      val(23,naok) = trj(Y,F3PP)
      val(24,naok) = trj(A,F3PP)
      val(25,naok) = trj(B,F3PP)
      val(26,naok) = trj(T,F3PP)

      val(32,naok) = trj(X,F3DC)
      val(33,naok) = trj(Y,F3DC)
      val(34,naok) = trj(A,F3DC)
      val(35,naok) = trj(B,F3DC)
      val(36,naok) = trj(T,F3DC)

c      write(*,*) "opt",trj(X,F3DC),trj(Y,F3DC),
c     $           trj(A,F3DC),trj(B,F3DC)

c because of DC32 down (2011/6/4)

c F4
      val(42,naok) = trj(X,F4)
      val(43,naok) = trj(Y,F4)
      val(46,naok) = trj(T,F4)

c F5
      val(52,naok) = trj(X,F5)
      val(53,naok) = trj(Y,F5)
      val(54,naok) = trj(A,F5)
      val(55,naok) = trj(B,F5)
      val(56,naok) = trj(T,F5)

c      write(81,*) trj(X,F5),trj(Y,F5),trj(A,F5),trj(B,F5)

c F6
      x6c4 = trj(X,F6)-c46*trj(X,F4) 
      x6c7 = trj(X,F6)-c76*trj(X,FH7) 
      val(62,naok) = trj(X,F6)
      val(63,naok) = trj(Y,F6)
      val(66,naok) = trj(T,F6)
      val(67,naok) = x6c4
      val(68,naok) = x6c7
c FH7
      x7c6 = trj(X,FH7)-trj(X,F6)/c76 
      val(72,naok) = trj(X,FH7)
      val(73,naok) = trj(Y,FH7)
      val(74,naok) = trj(A,FH7)
      val(75,naok) = trj(B,FH7)
      val(76,naok) = trj(T,FH7)
      val(77,naok) = x7c6
c FH91
      x9c7 = trj(X,FH9)-c79*trj(X,FH7) 
      val(92,naok) = trj(X,FH9)
      val(93,naok) = trj(Y,FH9)
      val(96,naok) = trj(T,FH9)
      val(97,naok) = x9c7

c S0
      xs0c7 = trj(X,S0)-c7s0*trj(X,FH7) 
      xs0c9 = trj(X,S0)-c9s0*trj(X,FH9) 
c      xs0c7 = trj(X,S0)-c7s0*trj(X,F6) 
c      write(*,*) trj(X,FH7), trj(X,S0), xs0c7
      val(102,naok) = trj(X,S0)
      val(103,naok) = trj(Y,S0)
      val(104,naok) = trj(A,S0)
      val(105,naok) = trj(B,S0)
      val(106,naok) = trj(T,S0)
      val(107,naok) = xs0c7
      val(108,naok) = xs0c9
c      write(*,*) trj(X,S0),trj(Y,S0),trj(A,S0),trj(B,S0)

c S0 PPAC
      xs0ppc9 = trj(X,S0PP)-c9s0pp*trj(X,FH9) 
      xs0ppc7 = trj(X,S0PP)-c7s0pp*trj(X,FH7) 
      val(142,naok) = trj(X,S0PP)
      val(143,naok) = trj(Y,S0PP)
c      val(144,naok) = trj(A,S0PP)
c      val(145,naok) = trj(B,S0PP)
c      val(146,naok) = trj(T,S0PP)
c      val(147,naok) = xs0ppc7
c      val(148,naok) = xs0ppc9

c S2
      xs2c6  = trj(X,S2)-c6s2*trj(X,F6) 
      xs2c7  = trj(X,S2)-c7s2*trj(X,FH7) 
      xs2c9  = trj(X,S2)-c9s2*trj(X,FH9) 
      xs2cs0 = trj(X,S2)-cs0s2*trj(X,S0) 
c optics run
c      xs2c  = trj(X,S2)-0.5811*trj(A,S2)
c 7Li run
c      xs2c  = trj(X,S2)-0.5392*trj(A,S2)
c 10/9
      xs2c  = trj(X,S2)-cs2f*trj(A,S2)
      xs2c3 = xs2c-0.4557*trj(X,F3DC)
      val(122,naok) = trj(X,S2)
      val(123,naok) = trj(Y,S2)
      val(124,naok) = trj(A,S2)
      val(125,naok) = trj(B,S2)
      val(126,naok) = trj(T,S2)
      val(127,naok) = xs2c7
      val(128,naok) = xs2c9
      val(129,naok) = xs2cs0
      val(130,naok) = xs2c
      val(131,naok) = xs2c3
c      val(132,naok) = trj(X,S2)+0.7585*trj(X,F6)
      val(132,naok) = xs2c6

c for higher-order Aberrations      
      val(201,naok) = trj(A,F3DC)*trj(B,F3DC)*trj(B,F3DC) ! abb
      val(202,naok) = trj(B,F3DC)*trj(B,F3DC)*trj(X,FH9)  ! bbd
      val(203,naok) = trj(A,F3DC)*trj(X,FH9)*trj(X,FH9)   ! add
      val(204,naok) = trj(A,F3DC)*trj(A,F3DC)*trj(X,FH9)  ! aad
c      val(205,naok) =
      val(205,naok) = trj(A,F3DC)*trj(X,FH9)              ! ad
      val(206,naok) = trj(B,F3DC)*trj(X,FH9)              ! bd

      end
      


