      subroutine enctc890(val,nx,ny,naok,
     &     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata,analyzer)
c      subroutine enctc890(val,nx,ny,naok,
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
         val(2,naok) = rawdata(1, id)
         val(3,naok) = rawdata(1, id) * 50. / 1000.  ! in ns
         val(4,naok) = rawdata(1, id) * 50. / 1000. / 1000. ! in us
      enddo

      return
      end

