      subroutine getatm(nline, ierr, string, elem, x, y, z, tag)
      implicit real(a-h,o-z)
      implicit integer(i-n)
c--------------------------------------------------------------
c  copyright (c) 1998 Bruce Ravel
c--------------------------------------------------------------
c========================================================================
c  This routine parses a line containing information about a unique atom
c  position in a unit cell and returns the information about that site.
c========================================================================
      parameter(nwdx=20)
      character*2   elem, test, check*10
      character*20  tag*10, words(nwdx)
      character*78  messg
      character*(*) string
      integer       nline, ierr

      test = 'ab'
      do 10 i=1,nwdx
        words(i) = ' '
 10   continue
      nwds = nwdx
      call bwords(string,nwds,words)

      if (nwds.lt.4) goto 666
      check = words(1)
      call lower(check)
      if (check(1:2) .eq. 'nu') goto 20
      if (is2z(words(1)).eq.0) goto 666
      if (istrln(words(1)).gt.2) goto 666
 20   continue
c      if ( (nwds.lt.4).or.(is2z(words(1)).eq.0).or.
c     $            (istrln(words(1)).gt.2) ) goto 666


c  get element symbol, x, y, z, tag in that order
      elem = words(1)
      call case(test,elem)
      call getrea('atomic position', words(2), x, nline, ierr)
      call getrea('atomic position', words(3), y, nline, ierr)
      call getrea('atomic position', words(4), z, nline, ierr)
      if ((nwds.gt.4) .and. ((words(5).ne.'!').and.
     $            (words(5).ne.'%').and.(words(5).ne.'*')) ) then
          tag = words(5)
      else
          tag = words(1)
          call fixsym(tag(1:2))
      endif

      return

 666  continue
      call messag(' ')
      call messag(' ')
 400  format(' *** Warning at line ',i3, ' while reading the atom',
     $            ' or basis list')
      write(messg,400)nline
      call messag(messg)
      call messag(' ')
      call messag('     The atom list is a formatted list, and must '//
     $            'be at the end of')
      call messag('        the input file.')
      call messag('     The first column is a two character '//
     $            'elemental symbol.')
      call messag('     The next three are coordinates in the unit '//
     $            'cell.')
      call messag('     The fifth column is the optional site tag.')
      call messag('     Atoms should be able to continue, but may '//
     $            'not produce')
      call messag('        the correct output.')
      call messag('     You should edit atoms.inp and try again.-')
      call messag(' ')
      ierr = 2

c  end subroutine getatm.f
      end
