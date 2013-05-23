      subroutine cellex(iat,natx,
     $                  iabs,iatom,iptful,fullcl,cell,dmax,
     $                  trmtx,itot,atlis)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c------------------------------------------------------------------------
c
c  !!!!!! presently passing ipt and st in fractional coords rather than
c         iptful and fullcl
c  !!!!!!!!! does not take overfull or cartesian info !!!!!!!!!!
c
c  given an overfull cell in the not-necessarily-orthogonal basis, 
c  expand into a cluster by tacking on adjacent cells.
c
c  check distance of each new atom, if less than dmax, convert to 
c  proper cartesian coordinates and add to atlis
c  cosum is used to check orthogonality of cell axes
c------------------------------------------------------------------------
c  input:
c    iat, natx: dimension parameters from calling program
c    iabs: index of central atom in arrays dimensioned iat
c    iatom: number of unique types
c    iptful: (iat) number of each type in overfull cell
c    fullcl(iat,192,3):  first index marks atom type
c                        second index marks each occurence in overfull
c                        third index marks x,y,z (in cell axis basis) 
c    cell: (6) array of a,b,c,alpha,beta,gamma
c    dmax: radius of cluster
c    trmtx: (3,3) transformation matrix cell axis & cartesian bases, from
c           subroutine metric
c  output:
c    itot: number of atoms in cluster
c    atlis(natx,8): 1-3 -> pos. in cell    4 -> dist to origin
c                     5 -> atom type     6,8 -> cartesian coords
c------------------------------------------------------------------

c      parameter (iat=50, natx=800)
c                natx: max # of atoms to keep
c                iat  : max # of atom types
      parameter (one=1., epsi=0.0001, zero=0)
      parameter (pi=3.141592653589793238462643)
      parameter (big=100000)

      dimension    cell(6)
      dimension    fullcl(iat,192,3), atlis(natx,8), trmtx(3,3)
      dimension    iptful(iat)
      dimension    anewpt(3), centpt(3), d(3)
      character*6  cnum
      character*7  creal
c      character*80 messg

4000  format(f7.4)
4010  format(i6)

c------------------------------------------------------------------
c  define the central atom, make cosum = sum of cosines of cell angles
      cosum = zero
      do 2 i=1,3
c        centpt(i) = fullcl(iabs,icnt,i) * cell(i)
        centpt(i) = fullcl(iabs,1,i) * cell(i)
        cosum     = cosum + abs(cos( cell(i+3)*pi/180 ))
 2    continue
 5    continue
      itot=0

c------------------------------------------------------------------------
c  determine max number of adjacent cells to fully enclose desired 
c  bubble.  this is a safe overkill
      na = int( dmax/cell(1) + 1 )
      nb = int( dmax/cell(2) + 1 )
      nc = int( dmax/cell(3) + 1 )
      dmin = big

c------------------------------------------------------------------------
c  loop through all adjacent cells and through all atoms in the overfull
c  cell.  then begin adding atoms to list
      do 80 ia=-na,na
        do 70 ib=-nb,nb
          do 60 ic=-nc,nc
            do 50 i=1,iatom
              do 40 j=1,iptful(i)
                anewpt(1) = (fullcl(i,j,1) + ia)*cell(1)
                anewpt(2) = (fullcl(i,j,2) + ib)*cell(2)
                anewpt(3) = (fullcl(i,j,3) + ic)*cell(3)
c                           distance to central atom
                do 10 icheck=1,3
                  d(icheck) = anewpt(icheck) - centpt(icheck)
10              continue
c                           too far away
                dd = dist(d,cell)
                if (dd.gt.epsi) dmin = min(dmin, dd)
                if ( dd.gt.dmax ) goto 30

c                           add to atlis
                itot          = itot + 1
                atlis(itot,1) = anewpt(1)/cell(1)
                atlis(itot,2) = anewpt(2)/cell(2)
                atlis(itot,3) = anewpt(3)/cell(3)
                atlis(itot,4) = dd
                atlis(itot,5) = i * one
c                           calculate cartesian coordinates of atom 
                if (cosum.gt.epsi) call trans(anewpt,trmtx)
                do 20 k=1,3
                        atlis(itot,k+5) = anewpt(k)
20              continue

c               try again with smaller dmax if storage exceeded 
                if (itot.ge.natx) then
                    call messag(' ')
                    call messag('* * * WARNING * * *')
                    write(cnum,4010)natx
                    call messag('You have exceeded '//cnum//' atoms.')
                    call messag('Rmax reduced by 1 Angstrom '//
     $                          'to accommodate.')
                    dmax = dmax - 1
                    write(creal,4000)dmax
                    call messag('New rmax = '//creal)
                    call messag(' ')
                    goto 5
                endif
30              continue
40            continue
50          continue
60        continue
70      continue
80    continue

c  if dmax chosen too small, there will be only one atom in the atom list
c  (the central atom) and the sort routine will barf
      if (dmin.gt.dmax) then
          call messag(' ')
          call messag('* * * WARNING * * *')
          call messag('Rmax was smaller than the nearest neighbor '//
     $                'distance.')
          call messag('Rmax doubled and cluster expansion run again.')
          call messag(' ')
          dmax = 2*dmax
          goto 5
      endif

      return
c  end subroutine cellex
      end

