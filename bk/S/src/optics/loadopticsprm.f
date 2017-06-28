      subroutine loadopticsprm

      include './optics.fh'
      integer lun
      character*256 prmfile,prmfile2

      lun  = 80
c     xxcorrection.prm: mom. correction for focus tuning
      prmfile= 'prm/optics/xxcorrelation_mass14may.prm'
      open(unit=lun,file=prmfile,status='old')
      read(lun,*) c46
      read(lun,*) c76
      read(lun,*) c79
      read(lun,*) c7s0
      read(lun,*) c9s0
      read(lun,*) c7s0pp
      read(lun,*) c9s0pp
      read(lun,*) c7s2
      read(lun,*) c9s2
      read(lun,*) cs0s2
      read(lun,*) c6s2
      close(lun)
      lun  = 80
c     corrections.prm: (x|ad) crrection at S2
      prmfile2= 'prm/optics/corrections.prm'
      open(unit=lun,file=prmfile2,status='old')
      read(lun,*) cs2f
      close(lun)
      return
      end
