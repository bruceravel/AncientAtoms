      subroutine upper (str)
c  changes a-z to upper case.  ascii specific
c-   for ascii:  ichar(upper case 'a') =  65
c-               ichar(lower case 'a') =  97
      character*(*)  str
      integer iupa, iloa, iloz, idif
      data    iupa, iloa / 65, 97/
      idif = iloa - iupa
      iloz = iloa + 25
      jlen = max(1, istrln (str) )
      do 10  i = 1, jlen
         ic = ichar (str(i:i))
         if ((ic.ge.iloa).and.(ic.le.iloz))  str(i:i) = char(ic-idif)
   10 continue
      return
c end subroutine upper
      end
