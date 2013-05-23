      subroutine atchck(iat,ndopx,core,dopant,edge,
     $                  iatom,ibasis,ispa,idop,
     $                  x,y,z,cell,dmax,gasses,qvect)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c--------------------------------------------------------------------
c input:
c   dopant: matrix of element symbols
c   x,y,z : arrays of unique atom fractional coordinates
c   cell:   array of cell constants
c   iatom:  number of unique atoms
c   ibasis: number of atoms in basis
c   dmax:   radius of desired cluster
c   ispa:   number between 1 and 230 denoting space group
c   pargon: percent argon in the i0 chamber
c   pnitro: percent nitrogen in the i0 chamber
c   pkrypt: percent krypton in the i0 chamber
c--------------------------------------------------------------------
c  check the consistancy of all the input parameters.
c  if any are funny write a run-time error message and
c  die gracefully.
c--------------------------------------------------------------------
c      parameter(iat=50,ndopx=4)
      parameter(epsi=0.001)
      character*2  el,dopant(iat,ndopx),dp,test,edge
      character*10 core
c      character*77 messg
      dimension    x(iat), y(iat), z(iat), cell(6), qvect(3)
      dimension    idop(iat), gasses(3)
      logical      ldie

      ldie  = .false.
      icent = 0
      inull = 0
      icore = 0
      test  = 'ab'
      iall  = iatom
      if (ibasis.gt.0) iall=ibasis

      do 100 i=1,iall
        el = dopant(i,1)
        call case(test,el)
        if ((el.eq.'nu').and.(i.ne.1)) inull = inull+1
        if (core.eq.el)                icore = icore+1

        if ((x(i).gt.1.0).or.(x(i).lt.-1.0)) then
            call messag(' ')
            call messag('Atom positions must be real numbers '//
     $                  'between -1 and 1.')
            ldie=.true.
        endif
        if ((y(i).gt.1.0).or.(y(i).lt.-1.0)) then
            call messag(' ')
            call messag('Atom positions must be real numbers '//
     $                  'between -1 and 1.')
            ldie=.true.
        endif
        if ((z(i).gt.1.0).or.(z(i).lt.-1.0)) then
            call messag(' ')
            call messag('Atom positions must be real numbers '//
     $                  'between -1 and 1.')
            ldie=.true.
        endif

        if ((is2z(el).eq.0).and.(el.ne.'nu')) then
            call messag(' ')
            call messag(dopant(i,1)//'???')
            call messag('One of your elements is not in the '//
     $                  'periodic table.')
            ldie=.true.
        endif
        if (core.eq.'nu') then
            call messag(' ')
            call messag('The core atom cannot be a null site.')
            ldie=.true.
        endif
 100  continue

c%%%        if (icent.ne.1) then
C%%%              call messag(' ')
c%%%            call messag('Feff requires one and only one absorption '//
c%%%       $                'site.')
c%%%            ldie=.true.
c%%%        endif

      if (inull.ge.1) then
          call messag(' ')
          call messag('Only one empty site is allowed and it '//
     $                'must be the first site listed.')
          ldie=.true.
      endif

c  check if core is a dopant
      do 104 i=1,iat
        do 102 j=2,idop(i)
          dp=dopant(i,j)
          call case(test,dp)
          if (core.eq.dp) icore=icore+1
 102    continue
 104  continue

      if (icore.eq.0) then
          call messag(' ')
          call messag('The central atom specified by the keyword '//
     $                '"core" is not')
          call messag('in the atom or dopant lists.')
          ldie=.true.
      endif

C%%%        if (icore.ge.2) then
C%%%            call messag(' ')
C%%%            call messag('Error reading the core atom.')
C%%%            ii=istrln(core)
C%%%            messg = core(:ii)//' appears more than once in the atom list.')
C%%%            call messag(messg)
C%%%            ldie=.true.
C%%%        endif

      if ((iatom.eq.0).and.(ibasis.eq.0)) then
          call messag(' ')
          call messag('You included no atoms in your input file.')
          ldie=.true.
      endif

      if (ispa.eq.0) then
          call messag(' ')
          call messag('Your space group does not exist.  Check '//
     $                'the international tables.')
          call messag('The ATOMS document explains adapting '//
     $                'notation for the keyboard.')
          ldie=.true.
      endif

      do 110 i=1,3
        if (cell(i).le.0) then
            call messag(' ')
            call messag('Cell constants cannot be negative or zero.')
            ldie=.true.
        endif
        if ((cell(i+3).lt.0).or.(cell(i+3)-180.0.gt.epsi)) then
            call messag(' ')
            call messag('Cell angles must be stated between '//
     $                  '0 and 180 degrees.')
            ldie=.true.
        endif
c%%%          if (qvect(i).lt.0) then
c%%%              call messag(' ')
c%%%              call messag('Values for the components of the q vector '//
c%%%       $              'in dafs must be non-negative.')
c%%%              ldie=.true.
c%%%          endif
110   continue

      if (dmax.le.0) then
          call messag(' ')
          call messag('Rmax is a radial distance.  It must be '//
     $                'positive.')
          ldie=.true.
      endif

      if ( (iatom.gt.1).and.(ibasis.gt.1) ) then
          call messag(' ')
          call messag('You may not specify both an atom list and a '//
     $                'basis list.')
          ldie=.true.
      endif

c       if ( (gasses(1)+gasses(3)+gasses(2))-1.0.gt.epsi ) then
c           call messag(' ')
c           call messag('The sum of the percentages of nitrogen, '//
c      $                'argon and krypton in the ')
c           call messag('I0 chamber cannot exceed 1.0')
c           ldie=.true.
c       endif
c
c       if ( (gasses(1).lt.-epsi).or.(gasses(3).lt.-epsi).or.
c      $            (gasses(3).lt.-epsi) ) then
c           call messag(' ')
c           call messag('The percentage of nitrogen, argon or '//
c      $         'krypton in the i0 chamber')
c           call messag('cannot be less than 0.')
c           ldie=.true.
c       endif

      call case(test,edge)
      if ( .not.((edge.eq.'k').or.(edge.eq.'l3').or.(edge.eq.'l2')
     $            .or.(edge.eq.'l1')) ) then
          call messag(' ')
          call messag('ATOMS only recognizes K and L edges.')
          ldie=.true.
      endif

      if (ldie) goto 666


      return

666   continue
      call messag('ATOMS cannot continue.  Please edit atoms.inp '//
     $            'and try again.')
      call messag(' ')
      stop

c end subroutine atchck
      end
