       subroutine calcdata
      implicit none
      integer nx/500/,ny/500/,naok(50)
      real val(500,500,50)
      common/ANALYSVALUE/ val,naok
#include "analyslogic.fh"
#include "analysevent.fh"
#include "rawdata.fh"
#include "commonprm.fh"
#include "runstat.fh"

      integer i
      character*132 ihitminchar
      integer catid, analyzer
      integer analy(4),cat(4)
      integer Ttrg
      integer ndataid
      integer dataid(10)

      external encgeneric
      external encrf
      external encppac
      external encspla
c      external encwinds
c      external encbldc31a
c      external encbldc32a
      external  encbldc30a
      external encbldcs20a
      external encbldcs22a
      external encbldcs23a
      external encbldcs20c
      external enc3
c      external enccrdcanode
c      external enccrdccathode
c      external encnyoki
      external enctraj
c      external encpid
      external encoptics
      external encphys
      external enccoin
c      external encwphys
c      external encs0
      external encspla_qdc
      external enctest
      external encdiamond
c      external encdali
      external encuhodoa, encuhodod, encuhodo
c     external encsharaq
c     external encrid, enctracker
c      external encdali, encgamma
      external encsi
      external encstopper

c BigRIPS ppac
      call callanalyzer(1,1,encppac)
c RF
      call callanalyzer(3,3,encrf)
c SHARAQ plastics
      call callanalyzer(21,21,encspla)
c DC31 anode
      call callanalyzer(22,22,encbldc30a)
c DC32 anode
      call callanalyzer(23,23,encbldc30a)
C F3MWDC Tracking
      analyzer=24
      if (analyzerflag(analyzer))
     $  call enc3(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c DC32 cathode
c      call callanalyzer(25,25,encgeneric)
c      call callanalyzer(24,24,encdcc)
c DC61 anode
c      call callanalyzer(25,25,encbldcs22a)
c DC71 anode
      call callanalyzer(28,28,encbldcs20a)
c DC71 cathode
c      call callanalyzer(29,29,encbldcs20c)
c DC72 anode
      call callanalyzer(30,30,encbldcs20a)
c DC72 cathode
c      call callanalyzer(31,31,encbldcs20c)
c DC91 anode
      call callanalyzer(33,33,encbldcs20a)
c DC91 cathode
c      call callanalyzer(34,34,encbldcs20c)
c DCX1 anode
      call callanalyzer(36,36,encbldcs20a)
c DCX1 cathode
c      call callanalyzer(37,37,encbldcs20c)
c DCX2 anode
      call callanalyzer(38,38,encbldcs20a)
c DCS0 anode
      call callanalyzer(39,39,encbldcs23a)
C Tracking by using MWDCs at FHX and S0
C      analyzer=40
C      if (analyzerflag(analyzer))
C     $  call encs0(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c plastic qdc data at fh10
      call callanalyzer(42,42,encspla_qdc)


c ======== added on 2014/5/8 by motoki ========
c DCS2-1 (upstream) anode
      call callanalyzer(45,45,encbldcs22a)
c DCS2-2 (downstream) anode
      call callanalyzer(46,46,encbldcs22a)

c Tracking by using two MWDCs at S2
      analyzer=48
      if (analyzerflag(analyzer))
     $  call encs2(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c S2 Si
      call callanalyzer(41,41,encsi)
c ===========================================

c === Active Stopper
      call callanalyzer(35,35,encstopper)

c neutron wall 
c      call callanalyzer(44,44,encwinds)
c s2 crdc anode
c      call callanalyzer(45,45,enccrdcanode)
c s2 crdc cathode
c      call callanalyzer(46,46,enccrdccathode)
c     s2 nyokinyoki 
c      call callanalyzer(48,48,encnyoki)


      call callanalyzer(49,44,encgeneric)
c     check reference timing for v1190
      call callanalyzer(50,50,encgeneric)
      call callanalyzer(10,10,encdiamond)
c
c F3 uHodoscope
      call callanalyzer(5,5,encuhodod)
      call callanalyzer(6,6,encuhodoa)
c coincidence register
      call callanalyzer(7,7,enccoin)


c     dali
c      call callanalyzer(12,12,encdali)


c RID
c      call callanalyzer(14,14,encrid)


c gamma
c      analyzer=16
c      if (analyzerflag(analyzer))
c     $  call encgamma(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)


c trajectory
      analyzer=15
      if (analyzerflag(analyzer))
     $   call enctraj(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c particle id
c      analyzer=17
c      if (analyzerflag(analyzer))
c     $   call encpid(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c optics
      analyzer=18
      if (analyzerflag(analyzer))
     $   call encoptics(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c physics
      analyzer=47
      if (analyzerflag(analyzer))
     $   call encphys(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c winds-physics
c      analyzer=43
c      if (analyzerflag(analyzer))
c     $   call encwphys(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)
c sharaq?
c      analyzer=19
c      if(analyzerflag(analyzer))
c     $   call encsharaq(val(1,1,analyzer),nx,ny,naok(analyzer),analyzer)

      return
      end

      subroutine callanalyzer(analyzer,catid,enc)
      implicit none
#include "analyslogic.fh"
#include "analysevent.fh"
#include "rawdata.fh" 

      integer nx/500/,ny/500/,naok(50)
      real val(500,500,50)
      common/ANALYSVALUE/ val,naok

      integer analyzer
      integer catid

      if (.not. analyzerflag(analyzer)) goto 1000

      if (rawheadpos(catid).eq.0) then
         write(*,*) ' ANAPAW-E : invalid rawheadpos(catid =',catid,
     $        ') in analyzer =',analyzer
         goto 1000
      endif

      if (hitheadpos(catid).eq.0) then
         write(*,*) ' ANAPAW-E : invalid hitheadpos(catid =',catid,
     $        ') in analyzer =',analyzer
         goto 1000
      endif

      call enc(val(1,1,analyzer),nx,ny,naok(analyzer),
     $     rawdata(rawheadpos(catid)),nhitdata(hitheadpos(catid)),
     $     hitdet(hitheadpos(catid)),
     $     nhitdet(catid),ndet(catid),ndata(catid),analyzer)

 1000 continue

      return
      end
