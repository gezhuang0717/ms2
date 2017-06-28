      subroutine sortRawData
      implicit none
      include 'analysevent.fh'
      include 'rawdata.fh'
      integer i,j,seg_modid
c      write(*,*) '<sortRawData> begin'
      	
      do i = 1,maxsegid
        if (nw(i).gt.0) then
         seg_modid = iand(segid(i),x'0000ff')
         if (seg_modid .eq. 0) then
            call decCamac16fixed(segid(i),evtdata(addr(i)),nw(i))
         endif
         if (seg_modid .eq. 9) then
            call decP7166(segid(i),evtdata(addr(i)),nw(i))
         endif
         if (seg_modid .eq. 5) then
            call decAD413A(segid(i),evtdata(addr(i)),nw(i))
         endif
         if (seg_modid .eq. 6) then
c            call decL3377s(segid(i),evtdata(addr(i)),nw(i))
            call decL3377d(segid(i),evtdata(addr(i)),nw(i))
         endif
c        if (seg_modid .eq. 10) then
c            call decCTM405(segid(i),evtdata(addr(i)),nw(i))
c        endif
         if (seg_modid .eq. 21) then ! DALI-A
            call decV792compat(segid(i),evtdata(addr(i)),nw(i))
         endif
         if (seg_modid .eq. 24) then ! DALI-T
            call decV1190(segid(i),evtdata(addr(i)),nw(i))
         endif
        endif
      enddo	
c      write(*,*) '<sortRawData> end'

      return 
      end


      subroutine initRawData
      implicit none
      include 'rawdata.fh'
      
      call fill_int(rawdata,kMaxDataID*kMaxDetID*kMaxCatID,-1000)
      call fill_int(nData,kMaxDetID*kMaxCatID,0)
      call fill_int(nHitDet,kMaxCatID,0)

      return
      end

      subroutine cleanupRawData
      implicit none
      include 'rawdata.fh'
      integer i,j,k,detid
      do i=1,kMaxCatID
         do j=1,nHitDet(i)
            detid = hitDet(j,i)
            nData(detid,i) = 0
            do k=1,nDataMin(i)
               rawdata(k,detid,i) = -1000
            enddo
         enddo
         nHitDet(i) = 0
      enddo
      return
      end
         
            

      
