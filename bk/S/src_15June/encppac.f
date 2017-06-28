c Analyzer 1 : encppac
      subroutine encppac(val,nx,ny,naok,
     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata)
c ----------------------------------------------------------------------------
c ANALYZER 1 : PPAC
c ----------------------------------------------------------------------------
cID = 1: F1PPAC-1
cID = 2: F1PPAC-2
cID = 3: F3PPAC-1A
cID = 4: F3PPAC-1B
cID = 5: F3PPAC-2A
cID = 6: F3PPAC-2B
cID = 7: F2PPAC-1A
cID = 8; F2PPAC-1B
cID = 9: F4PPAC
cID =10: F5PPAC-1A
cID =11: F5PPAC-1B
cID =12: F5PPAC-2A
cID =13: F5PPAC-2B
cID =14: F6PPAC
cID =15: FH9PPAC-A
cID =16: FH9PPAC-B Not use
cID =17: FH10PPAC-A
cID =18: FH10PPAC-B Not use
cID =19: S0PPAC-A
cID =20: S0PPAC-B Not use (FHX Pla)
cID =21: Dummy
cID =22: Dummy
cID =23: Dummy
cID =24: Dummy
c Non standard
cW# : 1      2      3      4      5      6      7      8      9      10
c     ID     TX1raw TX2raw TY1raw TY2raw TAxraw QX1raw QX2raw QY1raw QY2raw
cW# : 11     12     13     14     15     16     17     18     19     20
c     QAraw  TX1    TX2    TY1    TY2    TA     QX1    QX2    QY1    QY2
cW# : 21     22     23     24     25     26     27     28     29     30
c     QA     
cW# : 31     32     33     34     35     36     37     38     39     40
c     X(mm)  Y(mm)  TXdiff TYdiff TXsum  TYsum  QXsum  QYsum  
c
c focus position
c
cID =101: F1
cID =102: F2
cID =103: F3
cID =104: F4
cID =105: F5
cID =106: F6
cID =107: F7
cID =108: F8
cID =109: F9
cID =110: F10
cID =111: F11
cID =112: F8 including F8-3   --> image @ standard F8 focus
cA = x/z, B = y/z, T = avg(A_T) of upstream PPAC
cW# : 1      2      3      4      5      6      7      8      9      10
c     ID     X      Y      A      B      T      ResX   ResY   A(mrad)B(mard)
cW# : 11     12     13     14     15     16     17     18     19     20
c     fflag  nhitX  nhitY  
c
c ZPOS Zposition is defined in ppac.prm
cID =201 - 208 (nLayer)
c  202: DayOne F8 target 
cW# : 1      2      3      4      5      6      7      8      9      10
c     ID     X      Y      A      B                           A(mrad)B(mard)
c     11
c     R
c---------------------------------------------------------------------------
c     nDataMin                  : the number of required data (ihit_min0)
c     hitDet(1:nHitDet)         : hit detector id, = id in mapfile
c     nData(1:lenDet)           : the number of data for each id (ihit)
c     nHitDet                   : the number of hit detector
c     lenData                   : the length of array
c     lenDet                    : the length of array
      implicit none
      include 'analyslogic.fh'
      include 'commonprm.fh'
      include 'valcommon.fh'
      include 'ppac.fh'
      include 'valppac.fh'

      integer nx, ny, ndet, ndata, naok
      integer rawdata(ndata,ndet)
      integer nhitdata(ndet)
      integer hitdet(ndet)
      integer nhitdet
      real val(nx,ny)

      integer iHit

c local
      real ppacx(nPPAC), ppacy(nPPAC), ppact(nPPAC)
      real raw_ppac(10),cal_ppac(20)
      logical iflag_x(nPPAC),iflag_y(nPPAC)
      integer id,i,j
      real temp

c for tracking
      integer ifocus
      integer npx,npy,ctrl,ierr
      integer npmax,npid(nfocus+1)
      parameter(npmax=5)
      real    fitx(npmax),fity(npmax),fitxz(npmax),fityz(npmax)
      integer npxid(npmax),npyid(npmax)
      real    dis
      real    resx(nfocus+1),resy(nfocus+1)
      integer ppac_id(npmax,nfocus+1)
      data npid/2,4,2,1,4,1,1,2/
      data (ppac_id(i,1),i=1,2)/1,2/
      data (ppac_id(i,2),i=1,4)/3,4,5,6/
      data (ppac_id(i,3),i=1,2)/7,8/
      data (ppac_id(i,4),i=1,1)/9/
      data (ppac_id(i,5),i=1,4)/10,11,12,13/
      data (ppac_id(i,6),i=1,1)/14/
      data (ppac_id(i,7),i=1,1)/15/
      data (ppac_id(i,8),i=1,2)/17,19/
      logical nearhit

c for ZPOS
      real zposx,zposy
c      real fh9tref,fh10tref
c      real tref

      If(InitENCflag(1)) then
         call LOADPPACPRM
         InitENCFlag(1) = .FALSE.
      endif

      do i=1,nPPAC
         ppacx(i)=-1000
         ppacy(i)=-1000
         ppact(i)=-1000
      enddo
      do i=1,nfocus+1
         ppac_fx(i) = -1000.
         ppac_fy(i) = -1000.
         ppac_fa(i) = -1000.
         ppac_fb(i) = -1000.
         ppac_ft(i) = -1000.
         ppac_fflag(i) = 0
      enddo
      naok = 0
c      tref = -10000.
      fh9tref = -10000.
      fh10tref = -10000.

      do iHit = 1,nHitDet
         id = hitDet(iHit)
c     if (ihitflag .and. nData(id).lt.nDataMin) cycle 
c     if (nData(id).lt.nDataMin) cycle 
c     write(*,*) 'nHitData(',id,' ) =', nhitdata(id)
c     write(*,*) 'nHitData(id) = ', nHitData(id), '  nData=', nData
c     if (nHitData(id).lt.nData) cycle 
cccccccfor R3 PPAC, minimum nData is 5 (TX1/TX2/TY1/TY2/TA) cccccccc
         do i=1,10
            raw_ppac(i) = rawdata(i,id) + rand() -0.5
         enddo
      
         if(id.eq.15) then
            fh9tref = raw_ppac(10)
         else if(id.eq.17) then
            fh10tref = raw_ppac(10)
         endif
         
         
         if (nHitData(id).lt.5) cycle 
         naok = naok + 1
         val(1,naok)  = id
         do i = 1,10
            raw_ppac(i) = rawdata(i,id)+rand()-0.5
         enddo
c     
c     subtruct Tref here
c     if(id.eq.15 .or. id.eq.17) then
c     tref = raw_ppac(10)
c     endif
         
         if(id.ge.15 .and. id.le.16) then
            do i=1,5
               raw_ppac(i) = raw_ppac(i) - fh9tref + 20000.
            end do
         end if
         if(id.ge.17 .and. id.le.20) then
            do i=1,5
               raw_ppac(i) = raw_ppac(i) - fh10tref + 20000.
            end do
         end if
c         if(id.ge.15 .and. id.le.20) then
c            do i=1,5
c               raw_ppac(i) = raw_ppac(i) - tref + 20000.
c            end do
c         end if
c      write(*,*) fh9tref, fh10tref
         
c     
      do i = 1,10
         val(1+i,naok)  = raw_ppac(i)
c       write(*,*) id,  raw_ppac(1), raw_ppac(2)
      enddo
c     
c     
         if(id.le.14) then
            call calc_ppac(raw_ppac,cal_ppac,
     &           id,iflag_x(id),iflag_y(id), .TRUE.)
         else
c     For FH9, FHX, don't use overflow
            call calc_ppac(raw_ppac,cal_ppac,
     &           id,iflag_x(id),iflag_y(id), .FALSE.)
         end if
c     write(*,*) cal_ppac(11),cal_ppac(12)
         do i = 1,10
            val(11+i,naok) = cal_ppac(i)
         enddo
         do i = 1,8
            val(30+i,naok) = cal_ppac(10+i)

cccccc add for enc beam by abey ccccccc
c for F2PPAC1
            if(id.eq.7) then
               f21x = cal_ppac(11)
               f21y = cal_ppac(12)
c for F2PPAC2
            else if(id.eq.8) then
               f22x = cal_ppac(11)
               f22y = cal_ppac(12)
c for F3PPAC1a
            else if(id.eq.3) then
               f31ax = cal_ppac(11)
               f31ay = cal_ppac(12)
c for F3PPAC1b
            else if(id.eq.4) then
               f31bx = cal_ppac(11)
               f31by = cal_ppac(12)
c for F3PPAC2b
            else if(id.eq.5) then
               f32ax = cal_ppac(11)
               f32ay = cal_ppac(12)
c for F3PPAC2b
            else if(id.eq.6) then
               f32bx = cal_ppac(11)
               f32by = cal_ppac(12)
c for F4PPAC
            else if(id.eq.9) then
               f4x = cal_ppac(11)
               f4y = cal_ppac(12)
c for F5PPAC1a
            else if(id.eq.10) then
               f51ax = cal_ppac(11)
               f51ay = cal_ppac(12)
c for F5PPAC1b
            else if(id.eq.11) then
               f51bx = cal_ppac(11)
               f51by = cal_ppac(12)
c for F5PPAC2a
            else if(id.eq.12) then
               f52ax = cal_ppac(11)
               f52ay = cal_ppac(12)
c for F5PPAC2b
            else if(id.eq.13) then
               f52bx = cal_ppac(11)
               f52by = cal_ppac(12)
c for F6PPAC
            else if(id.eq.14) then
               f6x = cal_ppac(11)
               f6y = cal_ppac(12)
c for FH9
            else if(id.eq.15) then
               fh9x = cal_ppac(11)
               fh9y = cal_ppac(12)
c for FH10
            else if(id.eq.17) then
               fh10x = cal_ppac(11)
               fh10y = cal_ppac(12)
c for S0
            else if(id.eq.19) then
               s0x = cal_ppac(11)
               s0y = cal_ppac(12)
            endif
ccccccccccccccccccccccccccccccccccccccc

         enddo
c 
c            write(*,*) fh10x,fh10y
c
         if (iflag_x(id)) ppacx(id) = cal_ppac(11)
         if (iflag_y(id)) ppacy(id) = cal_ppac(12)
         if (iflag_x(id).or.iflag_y(id)) then
            ppact(id) = cal_ppac(5)
         else
            ppact(id) = -1000.
         endif
      enddo
      
c     PPAC Tracking
      do ifocus =1,nfocus+1
         npx = 0
         npy = 0
         do i = 1,npid(ifocus)
            id = ppac_id(i,ifocus)
            if (iflag_x(id)) then
               npx = npx + 1
               npxid(npx) = id
               fitx(npx) = ppacx(id)
               fitxz(npx) = ppac_zpos(1,id)
            endif
            if (iflag_y(id)) then
               npy = npy + 1
               npyid(npy) = id
               fity(npy) = ppacy(id)
               fityz(npy) = ppac_zpos(2,id)
            endif
         enddo
         select case (npid(ifocus))
         case (1)
            if ((npx.eq.1).and.(npy.eq.1)) then
               ppac_fx(ifocus) = fitx(1)
               ppac_fy(ifocus) = fity(1)
               ppac_ft(ifocus) = ppact(npxid(1))
            endif
         case (2)
            if (npx.eq.2) then
               ppac_fx(ifocus) = 0.5*(fitx(1)+fitx(2))
               dis = fitxz(2)-fitxz(1)
               ppac_fa(ifocus) = (fitx(2)-fitx(1))/dis
               ppac_ft(ifocus) = 0.5*(ppact(npxid(1))+ppact(npxid(2)))
            endif
            if (npy.eq.2) then
               ppac_fy(ifocus) = 0.5*(fity(1)+fity(2))
               dis = fityz(2)-fityz(1)
               ppac_fb(ifocus) = (fity(2)-fity(1))/dis
               ppac_ft(ifocus) = 0.5*(ppact(npyid(1))+ppact(npyid(2)))
            endif
         case default
            nearhit = .false.
            if (npx.eq.2) then
               if ((npxid(1).eq.ppac_id(1,ifocus)).and.
     &            (npxid(2).eq.ppac_id(2,ifocus))) then
                  nearhit = .true.
               endif
               if ((npxid(1).eq.ppac_id(3,ifocus)).and.
     &            (npxid(2).eq.ppac_id(4,ifocus))) then
                  nearhit = .true.
               endif
            endif
            if (.not.nearhit) then
               call rayfit(fitx,fitxz,npx,
     &         ppac_fx(ifocus),ppac_fa(ifocus),resx(ifocus),ctrl,ierr)
               if (ierr.eq.1) then
                  ppac_fx(ifocus) = -1000.
                  ppac_fa(ifocus) = -1000.
               endif
            endif
            nearhit = .false.
            if (npy.eq.2) then
               if ((npyid(1).eq.ppac_id(1,ifocus)).and.
     &            (npyid(2).eq.ppac_id(2,ifocus))) then
                  nearhit = .true.
               endif
               if ((npyid(1).eq.ppac_id(3,ifocus)).and.
     &            (npyid(2).eq.ppac_id(4,ifocus))) then
                  nearhit = .true.
               endif
            endif
            if (.not.nearhit) then
               call rayfit(fity,fityz,npy,
     &         ppac_fy(ifocus),ppac_fb(ifocus),resy(ifocus),ctrl,ierr)
               if (ierr.eq.1) then
                  ppac_fy(ifocus) = -1000.
                  ppac_fb(ifocus) = -1000.
               endif
            endif
            j = 0
            do i = 1,2  ! Upstream PPAC only
               ppac_ft(ifocus) = 0.
               id = ppac_id(i,ifocus)
               if (iflag_x(id).or.iflag_y(id)) then
                  ppac_ft(ifocus) = ppac_ft(ifocus) + ppact(id)
               endif
            enddo
            if (j.eq.0) then
               ppac_ft(ifocus) = -1000.
            else
               ppac_ft(ifocus) = ppac_ft(ifocus)/j
            endif
         end select
         if ((ppac_fx(ifocus).gt.-1000.).and.
     &       (ppac_fy(ifocus).gt.-1000.)) then
            ppac_fflag(ifocus) = 1
         endif
         naok=naok+1
         val(1,naok) = 100 + ifocus
         val(2,naok) = ppac_fx(ifocus)
         val(3,naok) = ppac_fy(ifocus)
         val(4,naok) = ppac_fa(ifocus)
         val(5,naok) = ppac_fb(ifocus)
         val(6,naok) = ppac_ft(ifocus)
         val(7,naok) = resx(ifocus)
         val(8,naok) = resy(ifocus)
         val(9,naok) = atan(ppac_fa(ifocus))*1000.
         val(10,naok) = atan(ppac_fb(ifocus))*1000.
         val(11,naok) = ppac_fflag(ifocus)
         val(12,naok) = npx
         val(13,naok) = npy
      enddo
c Position calc at ZPOS
c     ZPOS(focus number,z_position)
c      write(*,*) up to here
ccccc 2015.06.16 comment out ccccc
c      do i = 1,nLayer
c         if (ZPOS(1,i).gt.0) then
c            ifocus = int(ZPOS(1,i))
c            if (ppac_fflag(ifocus).eq.1) then
c               zposx = ppac_fx(ifocus) + ppac_fa(ifocus)*ZPOS(2,i)
c               zposy = ppac_fy(ifocus) + ppac_fb(ifocus)*ZPOS(2,i)
c            else
c               zposx = -1000.
c               zposy = -1000.
c            endif
c         endif
c         naok = naok + 1
c         val(1,naok) = 200 + i
c         val(2,naok) = zposx
c         val(3,naok) = zposy
c         val(4,naok) = ppac_fa(ifocus)
c         val(5,naok) = ppac_fb(ifocus)
c         val(9,naok) = atan(ppac_fa(ifocus))*1000.
c         val(10,naok) = atan(ppac_fb(ifocus))*1000.
c
c         val(11,naok) = sqrt( zposx**2 + zposy**2 )
c      enddo

      return
      end
