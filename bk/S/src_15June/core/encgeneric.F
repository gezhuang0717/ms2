c -*- compile-command: "make all" -*-
c last update: 2010/03/19 10:30:59 JST noji@maps:/home/noji/exp/nov09ivsm/users/noji/src/encgeneric.F
c
      subroutine encgeneric(val,nx,ny,naok,
     $     rawdata,nhitdata,hitdet,nhitdet,ndet,ndata,analyzer)

      implicit none
      integer nx, ny, ndet, ndata, naok
      integer rawdata(ndata,ndet)
      integer nhitdata(ndet)
      integer hitdet(ndet)
      integer nhitdet
      integer analyzer
      real val(nx,ny)
      integer i,j,id

#include "analyslogic.fh"

      if (initencflag(analyzer)) then
         write(*,*) ' ANAPAW-M : encgeneric analys =',analyzer
         write(*,*) ' ANAPAW-M : ndet  =',ndet
         write(*,*) ' ANAPAW-M : ndata =',ndata
         initencflag(analyzer) = .false. 
      endif

c      write(*,*) ' ANAPAW-M : analyzer =',analyzer, ' nhitdet =',nhitdet

      naok = 0

      do i=1,nhitdet
         id = hitdet(i)
c         write(*,*) 'id =',hitdet(i),' nhit =',nhitdata(id)
         naok = naok + 1
         val(1,naok) = id
         do j=1,ndata
            val(10+j,naok) = rawdata(j,id)
c            write(*,*) 'raw( ',j,') = ',rawdata(j,id)
         enddo
      enddo

      return
      end

