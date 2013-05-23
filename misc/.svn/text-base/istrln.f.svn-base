      function istrln (string)
c
c  returns index of last non-blank character.
c  returns zero if string is null or all blank.
      character*(*)  string
c-- if null string or blank string, return length zero.
      istrln = 0
      if (string(1:1).eq.char(0))  return
      if (string.eq.' ')  return
c
c-- find rightmost non-blank, non-null character.
      ilen = len (string)
      do 20  i = ilen, 1, -1
         if ((string (i:i) .ne. ' ') .and.
     $              (string (i:i) .ne. char(0)))  goto 30
   20 continue
   30 istrln = i

      return
c end function istrln
      end
