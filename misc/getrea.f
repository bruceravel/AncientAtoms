      subroutine getrea(keywrd,string,value,line,ierr)

      character*(*) keywrd, string
      character*72  messg
      integer       j, k, id, ie, line
      real          value
      logical       isnum

 400  format(bn,f13.0)
 410  format(bn,e19.13)
 420  format(bn,d19.13)
 430  format(' *** Error at line ', i3, ' of input file while')

      if (isnum(string)) then
          call lower(string)
          ie = index(string, 'e')
          id = index(string, 'd')
          if ((id.eq.0).and.(ie.eq.0)) then
              read(string,400)value
          elseif (ie.ne.0) then
              read(string,410)value
          elseif (id.ne.0) then
              read(string,420)value
          endif
      else
          j = istrln(string)
          k = istrln(keywrd)
          write(messg,430)line
          call messag(messg)
          messg = '     reading '//string(:j)//
     $                ' as a real number for "'//keywrd(:k)//'".-'
          call messag(messg)
          ierr = 3
c          stop
      endif

      return
c  end subroutine getrea
      end
