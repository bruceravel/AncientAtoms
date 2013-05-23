       logical function isnum (string)
c  returns true if string can be a number, else returns false
c  recognizes e and d exponentials, bit is not foolproof
c  to be a number, a string must contain:
c     - only characters in  'de.+-, 1234567890' (case is checked)
c     - no more than one 'd' or 'e'
c     - no more than one '.'
       character*(*)  string, str*70, number*20
       integer        iexp, idec, i, ilen, ier, j, istrln
       real           x
       external       istrln
c  note:  layout of *number* is important: don't change this!!
       data           number   /'de.+-, 1234567890'/
c-
       isnum = .false.
       iexp  = 0
       idec  = 0
       str   = string
       ilen  = max(1, istrln (str) )
       call smcase(str, number )
       do 100  i = 1, ilen
          j = index(number,str(i:i) )
          if (j.le.0)               go to 200
          if((j.eq.1).or.(j.eq.2))  iexp = iexp + 1
          if (j.eq.3)               idec = idec + 1
 100   continue
c  every character in the string was found in  *number*
c  so the string probably is a number
       isnum = .true.
c  but let's do a few more tests:
c    number of exponential and decimal markers
       if (iexp.ge.2) isnum = .false.
       if (idec.ge.2) isnum = .false.
c    read with iostat (this may report an error, but not always)
       read(str,150,iostat=ier)  x
 150   format (bn,f70.0)
       if (ier.ne.0)  isnum = .false.
c  all tests done
 200   continue
       return
c  end logical function isnum
       end
