c======================================                                                   
      subroutine coef(cfi,ncef,cf)
c======================================                                                   

      implicit double precision(a-h,o-z)
      dimension cf(30),cfi(30)
      common/PARAM/aai(30),II(30),ncoef

      do 1 I = 1,ncef
         cf(II(I)) = cfi(I)
 1           continue
      if(ncoef.gt.ncef) then
         do 101 I = ncef+1,ncoef
            cf(II(I)) = aai(II(I))
 101             continue
      endif

      return
      end
