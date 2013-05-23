      subroutine getint(keywrd,string,ivalue,line,ierr)

      character*(*) keywrd, string
      character*72  messg
      integer       ivalue, j, k, line

 400  format(bn,i10)
 430  format(' *** Error at line ', i3, ' of input file while')

      read(string,400,iostat=ierr)ivalue
      if (ierr.ne.0) then
          j = istrln(string)
          k = istrln(keywrd)
          write(messg,430)line
          call messag(messg)
          messg = '     reading '//string(:j)//
     $                ' as an integer for "'// keywrd(:k)//'".-'
          call messag(messg)
          ierr = 3
c           stop
      endif

      return
c  end subroutine getint
      end
