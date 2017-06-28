      subroutine int2char(num,ch)

      integer num
      character*4 ch

      ch = '    '

      write(ch,*) num

      if(num.gt.99)then
        ch='0'//ch(2:4)
      else
        ch='00'//ch(2:3)
      endif

      return
      end

