      subroutine triml (string)
c removes leading blanks.
      character*(*)  string
      jlen = istrln(string)
c
c-- all blank and null strings are special cases.
      if (jlen .eq. 0)  return
c-- find first non-blank char
      do 10  i = 1, jlen
         if (string (i:i) .ne. ' ')  goto 20
  10  continue
  20  continue
c-- if i is greater than jlen, no non-blanks were found.
      if (i .gt. jlen)  return
c-- remove the leading blanks.
      string = string (i:)
      return
c end subroutine triml
      end
