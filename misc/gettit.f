      subroutine gettit(keywrd,string,title,ntit,stdout)

      character*(*) keywrd, string, title
      character*72  messg, toss
      integer       ntit, i, j
      logical       stdout

      ntit  = ntit + 1
      toss  = string
      call case(keywrd,toss)
      i     = index(toss, keywrd)
      j     = istrln(keywrd)
      title = string(i+j+1:70)
      call triml(title)
      if ( (title(:1).eq.'=') .or. (title(:1).eq.',') ) then
          toss  = title(2:)
          title = toss
          call triml(title)
      endif
      if (.not.stdout) then
          messg = '  title > '//title
          call messag(messg)
      endif

      return
c  end subroutine gettit
      end
