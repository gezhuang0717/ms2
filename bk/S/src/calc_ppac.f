      subroutine calc_ppac(raw_ppac,cal_ppac,
     -                        iplane,flag_x,flag_y, oflowflag)
      implicit none
      
      include 'ppac.fh'
      include 'util.fh'
      
      real raw_ppac(10)
      real cal_ppac(20)
      integer iplane
      logical flag_x,flag_y
      logical oflowflag
      
      real tx1,tx2,ty1,ty2,ta,txdiff,tydiff,txsum,tysum
      real qx1,qx2,qy1,qy2,qa,qxsum,qysum
      real ppacx,ppacy
      real oflow
c      parameter(oflow = 4095.)
      parameter(oflow = 4094.)  ! modified 12/12 8:11 AM

      integer i,j

      do i = 1,20
         cal_ppac(i) = -10000.
      enddo
      flag_x = .FALSE.
      flag_y = .FALSE.
c      write(*,*) 'PPAC ch2ns'
c      write(*,*) (PPAC_ch2ns(j,iplane),j=1,5)
c      write(*,*) (raw_ppac(j),j=1,5)
      tx1 = raw_ppac(1)* PPAC_ch2ns(2,iplane)
      tx2 = raw_ppac(2)* PPAC_ch2ns(4,iplane)
      ty1 = raw_ppac(3)* PPAC_ch2ns(3,iplane)
      ty2 = raw_ppac(4)* PPAC_ch2ns(5,iplane)
      ta  = raw_ppac(5)* PPAC_ch2ns(1,iplane)

      qx1 = raw_ppac(6)  - PPAC_qped(1,iplane)
      qx2 = raw_ppac(7)  - PPAC_qped(2,iplane)
      qy1 = raw_ppac(8)  - PPAC_qped(3,iplane)
      qy2 = raw_ppac(9)  - PPAC_qped(4,iplane)
      qa  = raw_ppac(10) - PPAC_qped(5,iplane)

      if (gate1d(raw_ppac(1),0.0,oflow,0).and.
     &     gate1d(raw_ppac(2),0.0,oflow,0).and.
     &     gate1d(raw_ppac(5),0.0,oflow,0)) then
         txsum = tx1 + tx2 -2.0*ta
      else
         if(oflowflag .eqv. .FALSE.) then
            txsum = tx1 + tx2 -2.0*ta
         else
            txsum = -10000.
         end if
      endif
      if (gate1d(raw_ppac(3),0.0,oflow,0).and.
     &     gate1d(raw_ppac(4),0.0,oflow,0).and.
     &     gate1d(raw_ppac(5),0.0,oflow,0)) then
         tysum = ty1 + ty2 -2.0*ta
      else
         if(oflowflag .eqv. .FALSE.) then
            tysum = ty1 + ty2 -2.0*ta
         else
            tysum = -10000.
         end if
      endif
      qxsum = sqrt(qx1*qx2)
      qysum = sqrt(qy1*qy2)

      if (gate1d(raw_ppac(1),0.0,oflow,0).and.
     &     gate1d(raw_ppac(2),0.0,oflow,0)) then
         txdiff = tx1 - tx2
      else
         if(oflowflag .eqv. .FALSE.) then
            txdiff = tx1 - tx2
         else
            txdiff = -10000.
         end if
      endif
      
      if (gate1d(raw_ppac(3),0.0,oflow,0).and.
     &     gate1d(raw_ppac(4),0.0,oflow,0)) then
         tydiff = ty1 - ty2
      else
         if(oflowflag .eqv. .FALSE.) then
            tydiff = ty1 - ty2
         else
            tydiff = -10000.
         end if

      endif
      if (iflag_outside_offset.eq.1) then
         txdiff = txdiff - Xoutoffset(iplane)
         tydiff = tydiff - Youtoffset(iplane)
      endif

c Tsum check
      ppacx = -1000.
      ppacy = -1000.
      if ((.not. iflag_txsumgate)
     &     .or. gate1d(txsum,TXsumgate(1,iplane),
     &     TXsumgate(2,iplane),1)) then
         ppacx = txdiff*ns2mm(1,iplane)
         if (iflag_inside_offset.eq.1) then
            ppacx = ppacx - Xinoffset(iplane)
         endif
         if (iflag_geom_offset.eq.1) then
c            if(iplane.eq.35) then
c               ppacx = (-ppacx) - geomoffset(1,iplane)
c            else
               ppacx = ppacx - geomoffset(1,iplane)
c            endif
         endif
c Definition for optics
         ppacx = -1.0*ppacx
         flag_x = .TRUE.
      endif
      if ((.not. iflag_tysumgate)
     &     .or.gate1d(tysum,TYsumgate(1,iplane),
     &     TYsumgate(2,iplane),1)) then
         ppacy = tydiff*ns2mm(2,iplane)
         if (iflag_inside_offset.eq.1) then
            ppacy = ppacy - Yinoffset(iplane)
         endif
         if (iflag_geom_offset.eq.1) then
c            if(iplane.eq.35) then
c               ppacy = (-ppacy) - geomoffset(2,iplane)
c            else
               ppacy = ppacy - geomoffset(2,iplane)
c            endif
         endif
         flag_y = .TRUE.
      endif

      cal_ppac(1) = tx1
      cal_ppac(2) = tx2
      cal_ppac(3) = ty1
      cal_ppac(4) = ty2
      cal_ppac(5) = ta
      cal_ppac(6) = qx1
      cal_ppac(7) = qx2
      cal_ppac(8) = qy1
      cal_ppac(9) = qy2
      cal_ppac(10) = qa

      cal_ppac(11) = ppacx
      cal_ppac(12) = ppacy
      cal_ppac(13) = txdiff
      cal_ppac(14) = tydiff
      cal_ppac(15) = txsum
      cal_ppac(16) = tysum
      cal_ppac(17) = qxsum
      cal_ppac(18) = qysum
      return
      end
