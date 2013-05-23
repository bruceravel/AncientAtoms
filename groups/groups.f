      subroutine groups
c----------------------------------------------------------------------
c  identify space group from supplied symbol and set approriate unit
c  cell parameters, as needed perform the following:
c      1.  convert rhombohedral axes to hexagonal
c      2.  convert schoenflies notation to the standard hermann-maguin symbol
c      3.  permute non-standard settings to the standard
c      4.  convert special symbols (fcc, bcc, hcp, hex, diamond, cubic,
c            salt, nacl, cscl, perovskite, zincblende, zns, graphite)
c            to corresponding h-m symbols
c
c  some more chores for this:
c      1.  handle monoclinic cells properly
c      2.  accept a number between 1 and 230 for spcgrp
c
c  this subroutine calls:
c        atspec origin rh2hex schfix settng spcchk systm
c                + string manipulation routines case and messag
c----------------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)

      include 'atparm.h'
      include 'crystl.h'
      character*2  test
      logical      lrh, ltri, lhex, lalph, lbet, lgam
      parameter(epsi=1.e-5)

      test  = 'ab'
      lrh   = .false.
      ltri  = .false.
      lhex  = .false.
      lalph = .false.
      if ( abs(cell(4)-90.e0).gt.epsi ) lalph = .true.
      lbet  = .false.
      if ( abs(cell(5)-90.e0).gt.epsi ) lbet  = .true.
      lgam  = .false.
      if ( abs(cell(6)-90.e0).gt.epsi ) lgam  = .true.
      iall  = max(iatom, ibasis)
      call case(test,spcgrp)

c  equate a, b, and c if b and c are not given.
      if (cell(2).lt.0.1) cell(2) = cell(1)
      if (cell(3).lt.0.1) cell(3) = cell(1)

c  check for special lattice types, ie diamond,fcc,hcp,etc.
      call atspec(spcgrp,cell)

c  check the spcgrp against a list of the 230 actual groups and various
c  alternate settings
      call spcchk(spcgrp,ispa)
      if (ispa.eq.0) call settng(spcgrp, iperm, ispa)

c  flag rhombohedral, trigonal, hexagonal groups
      if (spcgrp(1:1).eq.'r') lrh  = .true.
      if (spcgrp(3:3).eq.'3') ltri = .true.
      if ((spcgrp(3:3).eq.'-').and.(spcgrp(4:4).eq.'3')) ltri = .true.
      if (spcgrp(3:3).eq.'6') lhex = .true.
      if ((spcgrp(3:3).eq.'-').and.(spcgrp(4:4).eq.'6')) lhex = .true.

c  check for rhombohedral space groups, convert data if rhombohedral
c  axes are given.
      if (lrh.and.lalph) then
          call rh2hex(iat,iall,cell,x,y,z)
          goto 100
      endif

c  set angles in hexagonal or trigonal group
      if ((lhex).or.(lrh).or.(ltri)) then
          cell(4) = 90.e0
          cell(5) = 90.e0
          cell(6) = 120.e0
      endif
 100  continue

c  obtain the system from the cell parameters for later comparison
c  with the results of equipt.
      call systm(cell,isystm)

c  check to see if space group is one that might need a shift
      call origin(spcgrp, shift, shwarn)

c      print*,cell
c      print*,'atom 1',x(1),y(1),z(1)
c      print*,'atom 2',x(2),y(2),z(2)

      return
c  end subroutine groups
      end
