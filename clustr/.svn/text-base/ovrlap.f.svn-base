      subroutine ovrlap(natx,itot,atlis,lover)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision(a-h,o-z)
c-----------------------------------------------------------------
c  check if any overlapping atoms have been generated in cluster by
c  comparing coordinates.
c  remove overlapping atoms from list 
c  write run-time messages if overlapped atoms found.
c  this has to be done after the heap sort and before the hash sort
c-----------------------------------------------------------------
c  input:
c    natx: dimension parameter from calling program
c  input/output:
c    itot:  number of atoms in cluster, on output: number of 
c           nonredundant atoms
c    atlis(natx,8): 1-3 -> pos. in cell    4 -> dist to origin
c                     5 -> atom type     6-8 -> cartesian coords
c                   on output has redundant atoms removed 
c  logical work array:  lover(natx)
c-----------------------------------------------------------------
c      parameter (natx=800)
      parameter(eps=.001)
c      parameter (iou=21)
c  eps is used in floating point logicals
      dimension    atlis(natx,8)
      logical      lover(natx),lmessg
c      character*80 messg, ovrout*10

      lmessg   = .false.
      lover(1) = .false.
c  loop through atoms, get dist x y z, then loop through previous atoms
c  and compare.  mark all but one overlapping atoms
      do 30 i=2,itot
        lover(i) = .false.
        dd = atlis(i,4)
        xi = atlis(i,6)
        yi = atlis(i,7)
        zi = atlis(i,8)
        do 20 j=1,i-1
          if ( abs(atlis(j,4)-dd).le.(eps*10) ) then
              xdiff = abs(xi-atlis(j,6))
              ydiff = abs(yi-atlis(j,7))
              zdiff = abs(zi-atlis(j,8))
              if ((xdiff.le.eps).and.(ydiff.le.eps).and.(zdiff.le.eps))
     $                    then
                  lover(i) = .true.
                  lmessg   = .true.
              endif
          endif
 20     continue
 30   continue

c  remove overlapping atoms from the list
c  every time one atom is removed, bubble the remaining atoms up
c  write overlapped positions to overlp.err
      if (lmessg) then
C%%%            ovrout= 'overlp.err'
C%%%            call lower(ovrout)
C%%%            open(unit=iou,file=ovrout,status='unknown')
C%%%            write(iou,400)' Positions of overlapping atoms in '//
C%%%       $                  'fractional cell coordinates.'
C%%%   400      format(1x,a)
          i = 1
c  do..while..
 40       continue 
          i = i+1
          if (lover(i)) then
C%%%                write (iou,410) atlis(i,1), atlis(i,2), atlis(i,3)
C%%%   410          format (3(3x,f8.5))
              do 70 j=i+1,itot
                lover(j-1) = lover(j)
                do 60 k=1,8
                  atlis(j-1,k) = atlis(j,k)
 60             continue 
 70           continue 
              i    = i-1
              itot = itot - 1
          endif
          if (i.lt.itot) goto 40

c  write run-time warning if overlapped atoms found
          call messag(' ')
          call messag('* * * * WARNING * * * *')
          call messag('Your input file has generated atoms '//
     $                'overlapping in space.')
          call messag('All redundant atoms have been removed from '//
     $                'the atom list.')
          call messag('Feff.inp has been written and the atom '//
     $                'list may be correct')
          call messag('but the atom labels and McMaster '//
     $                'calculations are incorrect.')
          call messag(' ')
          call messag('The most likely causes of this are:')
          call messag('  1: Specifying a unique crystallographic '//
     $                'site more than once. (See')
          call messag('      section 3.1 of the ATOMS document.)')
          call messag('  2: Constructing a basis list that overlaps '//
     $                'itself when translated.')
          call messag('  3: Specifying incorrect unique '//
     $                'crystallographic positions.')
C%%%            call messag(' ')
C%%%            ii = istrln(ovrout)
C%%%            call messag('Positions of atoms found to overlap have '//
C%%%       $                'been written')
C%%%            call messag('to a file called '//ovrout(:ii))
          call messag(' ')
C%%%            close(iou)
      endif

      return
c  end subroutine ovrlap
      end

