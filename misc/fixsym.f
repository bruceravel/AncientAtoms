      subroutine fixsym(sym)

c  returns a word with the first letter capitalized and the remaining
c  letters in lower case
c  this is very useful for writing out atomic symbols

      character*2 toss, sym*(*)

      toss = sym
      call upper(toss(1:1))
      call lower(toss(2:2))
      ii   = istrln(sym)
      if (ii.gt.2) then
          call lower(sym(3:ii))
          sym  = toss(1:1)//toss(2:2)//sym(3:ii)
      else
          sym  = toss(1:1)//toss(2:2)
      endif
      return
c  end subroutine fixsym
      end
