      subroutine untab (string)
c replace tabs with blanks :    tab is ascii dependent
      integer        itab , i, ilen
      parameter      (itab = 9)
      character*(*)  string, tab*1
      tab  = char(itab)
      ilen = max(1, istrln(string))
 10   continue
        i = index(string(:ilen), tab )
        if (i .ne. 0) then
            string(i:i) = ' '
            go to 10
        end if
      return
c end subroutine untab
      end
