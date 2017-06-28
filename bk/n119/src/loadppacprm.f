      subroutine loadppacprm
      implicit none
      include 'ppac.fh'
      integer ier
      character*132 ppacprmfile
      integer i, tmp
c      integer dbflag
c      parameter(dbflag=1)
      

      call getenv('PPAC_PRM',ppacprmfile)
      open(unit=80,file=ppacprmfile,status='old',err=1001)

      write(*,*) 'Load PPAC Param nPPAC=', nPPAC

      Call READ_FLT_LIST(80,PPAC_ch2ns,5*nPPAC,ier)
c      Write(*,*)ppac_ch2ns,'ch2ns '
      If(ier.NE.0)Goto 1002
      Call READ_FLT_LIST(80,PPAC_qped,5*nPPAC,ier)
c      Write(*,*)PPAC_qped,'qpedstal '
      If(ier.NE.0)Goto 1002
      Call READ_INT(80,iflag_outside_offset,1,ier)   ! line calib prm
c      Write(*,*) 'iflag_outside_offset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_FLT(80,Xoutoffset,nPPAC,ier)
c      Write(*,*)'Xoutoffset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_FLT(80,Youtoffset,nPPAC,ier)
c      Write(*,*)'Youtoffset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002

      Call READ_INT(80,iflag_inside_offset,1,ier)    ! kumagai-san prm
c      Write(*,*) iflag_inside_offset
c      Write(*,*)'iflag_inside_offset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_FLT_LIST(80,Xinoffset,nPPAC,ier)
c      Write(*,*) Xinoffset
c      Write(*,*)'Xinoffset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_FLT_LIST(80,Yinoffset,nPPAC,ier)
c      Write(*,*) Yinoffset
c      Write(*,*)'Yinoffset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002

      Call READ_FLT_LIST(80,ns2mm,2*nPPAC,ier)         ! kumagai-san prm
c      Write(*,*) ns2mm
c      Write(*,*)'ns2mm'
c      Write(*,*)
      do i = 1, nPPAC
         ns2mm(1,i) = ns2mm(1,i) / 2.
         ns2mm(2,i) = ns2mm(2,i) / 2.
      enddo
      If(ier.NE.0)Goto 1002
      Call READ_FLT_LIST(80,PPAC_zpos,2*nPPAC,ier)
c      Write(*,*) PPAC_zpos
c      Write(*,*)'PPAC_ZPOS'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_INT(80,iflag_geom_offset,1,ier)
c      Write(*,*) iflag_geom_offset
c      Write(*,*)'iflag_geom_offset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_FLT(80,geomoffset,2*nPPAC,ier)
c      Write(*,*) geomoffset
c      Write(*,*)'geomoffset'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
c     flag for txsum gate
      call read_int(80,tmp,1,ier)
      If(ier.NE.0)Goto 1002
      if (tmp.ne.0) then
         iflag_txsumgate = .true.
      else
         iflag_txsumgate = .false.
      endif
      Call READ_FLT_LIST(80,TXsumgate,2*nPPAC,ier)
c     Write(*,*) TXsumgate
c      Write(*,*)'TXsumgate'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
c     flag for tysum gate
      call read_int(80,tmp,1,ier)
      If(ier.NE.0)Goto 1002
      if (tmp.ne.0) then
         iflag_tysumgate = .true.
      else
         iflag_tysumgate = .false.
      endif
      Call READ_FLT_LIST(80,TYsumgate,2*nPPAC,ier)
c     Write(*,*) TYsumgate
c      Write(*,*)'TYsumgate'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      Call READ_FLT(80,ZPOS,2*nLayer,ier)
c      Write(*,*) ZPOS
c      Write(*,*)'ZPOS'
c      Write(*,*)
      If(ier.NE.0)Goto 1002
      write(*,*) ' ANAPAW-M : <loadppac> parameter loaded'
      close(80)
      return
 1001 write(*,*) 'Cannot open file',
     &     ppacprmfile(1:len_trim(ppacprmfile))
      stop
      return
 1002 write(*,*) 'Error while reading parameter'
      stop
      return
      end
