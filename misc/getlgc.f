      subroutine getlgc(keywrd,string,lvalue,line,ierr)

      character*(*) keywrd, string
      character*72  test*2
      logical       lvalue
      integer       line, ierr

      test   = 'ab'
      lvalue = .false.
      call triml(string)
      call case(test,string)
      if ( (string(:1).eq.'t') .or. (string(:1).eq.'y')
     $                         .or. (string(:2).eq.'on') )
     $            lvalue=.true.

      return
c  end subroutine getlog
      end
