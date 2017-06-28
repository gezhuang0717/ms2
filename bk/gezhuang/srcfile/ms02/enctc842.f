      subroutine enctc842(val,nx,ny,naok,
     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata,analyzer)
c      subroutine enctc842(val,nx,ny,naok,
c     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata,analyzer)

      implicit none
      integer nx, ny, ndet, ndata, naok
      integer rawdata(ndata,ndet)
      integer nhitdata(ndet)
      integer hitdet(ndet)
      integer nhitdet
      integer analyzer
      real val(nx,ny)
      integer i,j,id

      include 'analyslogic.fh'
      include 'commonprm.fh'

      if (initencflag(analyzer)) then
         write(*,*) ' ANAPAW-M : encadc analys =',analyzer
         write(*,*) ' ANAPAW-M : ndet  =',ndet
         write(*,*) ' ANAPAW-M : ndata =',ndata
         initencflag(analyzer) = .false. 
      endif

c      write(*,*) ' ANAPAW-M : analyzer =',analyzer, ' nhitdet =',nhitdet

      naok = 0

      do i=1,nhitdet
         id = hitdet(i)
         naok = naok + 1
         val(1,naok) = id
c        rawdata is converted double to 10ps step integer by TArtDecoderTC842
c         write(*,*) id, rawdata(1,id)
         val(2,naok) = rawdata(1, id) / 100.                 ! in nsec 
         val(3,naok) = rawdata(1, id) / 100. / 1000.         ! in usec
         val(4,naok) = rawdata(1, id) / 100. / 1000. / 1000. ! in msec

      ELC_UT = val(2,1) ! in nsec for ELC-U
      ELC_DT = val(2,2) ! in nsec for ELC-D

      FH10_LT = val(2,3) ! in nsec for FH10-L
      FH10_RT = val(2,4) ! in nsec for FH10-R
      Schottky = val(4,5) ! in msec for Schottky Trig.

      RMD3_T = val(2,6) ! in nsec for R-MD3
      RMD4_T = val(2,7) ! in nsec for R-MD4

      enddo
c      write(*,*) ELC_UT, ELC_DT

 

      return
      end

