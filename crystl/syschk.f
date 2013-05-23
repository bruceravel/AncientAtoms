      subroutine syschk(isyst, isystm, spcgrp, sysmes)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-----------------------------------------------------------------
c  write a useful error message if the lattice parameters do not
c  match the space group.
c-----------------------------------------------------------------
c  isyst:     system determination from space group, from equipt
c  isystm:    system determination from lattice constants, from systm
c  spcgrp*10: space group from input file
c  sysmes*77: (4) message about lattice/space group mismatch
c-----------------------------------------------------------------
c     isystm : 1 = monoclinic     2 = orthorhombic
c              3 = <not used>     4 = tetragonal
c              5 = cubic          6 = hexagonal
c              7 = triclinic      8 = <not used>
c  who came up with this notation convention? it sucks! -br
c-----------------------------------------------------------------
      character*10 spcgrp,sg
      character*12 ctype(8)
      character*50 prmtrs(8)

      parameter(nsysm=4)
      character*77 sysmes(nsysm)

      data (ctype(i),i=1,8) / 'monoclinic', 'orthorhombic', ' ',
     $                        'tetragonal', 'cubic', 'hexagonal',
     $                        'triclinic', ' '/
      data prmtrs(1)/'a, b, c unequal; alpha = gamma = 90; beta <> 90'/
      data prmtrs(2)/'a, b, c unequal; alpha = beta = gamma = 90'/
      data prmtrs(3), prmtrs(8) /' ', ' '/
      data prmtrs(4)/'a = b <> c ; alpha = beta = gamma = 90'/
      data prmtrs(5)/'a = b = c ; alpha = beta = gamma = 90'/
      data prmtrs(6)/'a = b <> c; alpha = beta = 90; gamma = 120'/
      data prmtrs(7)/'a, b, c unequal; alpha <> beta <> gamma <> 90'/


      ii = istrln(ctype(isystm))
      sysmes(1) = 'Your lattice constants and angles are '//
     $        'appropriate for a "'//ctype(isystm)(:ii)//'" crystal.'

      ii = istrln(ctype(isyst))
      sg = spcgrp
      call upper(sg)
      is = istrln(sg)
      sysmes(2) = 'Your input space group is "'//ctype(isyst)(:ii)//'".'

      ip = istrln(prmtrs(isyst))
      sysmes(3) = ctype(isyst)(:ii)//' groups require the '//
     $            'following parameters:'

      sysmes(4) = '       '//prmtrs(isyst)(:ip)

      return
c  end subroutine syschk
      end
