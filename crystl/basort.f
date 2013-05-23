      subroutine basort(iat,ndopx,ibasis,iabs,x,y,z,dopant,tag)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c------------------------------------------------------------------
c  want the absorber at the top of the unique atom list in the
c  basis calculation.  this routine puts the absorber on the point
c  of translation symmetry and shifts the coordinates of the rest
c  of the basis accordingly.  it also removes a null site from the 
c  list.
c------------------------------------------------------------------
c  input/output:
c    iat,ndopx: parameters set in calling program
c    ibasis: number of sites in basis on input, null site removed
c            on output 
c    iabs:   index of absorber in basis list on input, 1 on output
c    x,y,z:  (iat) on input: fractional coordinates of basis as listed 
c            in input file, on output:  fractional coordinates
c            shifted to put absorber at point of symmetry
c    dopant*2: (iat,ndopx) list of atomic symbols, reordered on output 
c            so that absorber is the first element of the array
c    tag*10: (iat) list of site tags, reordered on output so that 
c            absorber is the first element of the array
c------------------------------------------------------------------
c      parameter (iat=50, ndopx=4)
      parameter(ibig=10)
c  this number needs to be bigger than ndopx, 10 should suffice

      dimension    x(iat),y(iat),z(iat)
      character*2  dopant(iat,ndopx), etoss(ibig), el, test
      character*10 tag(iat),ttoss

      test = 'ab'
      if (iabs.eq.0) goto 999

c  on input, first site is the point of symmetry and may be null
c  iabs may be any other point in the list.  want iabs listed first
c  and the null site at the end for easy removal.

c  get vector to shift absorber onto basis point for point transform
      xshift = x(iabs) - x(1)
      yshift = y(iabs) - y(1)
      zshift = z(iabs) - z(1)

c  swap first and absorber and last
      do 10 i=1,ndopx
        etoss(i) = dopant(1,i)
 10   continue 
      ttoss = tag(1)
      xtoss = x(1)
      ytoss = y(1)
      ztoss = z(1)

      do 20 i=1,ndopx
        dopant(1,i) = dopant(iabs,i)
 20   continue 
      tag(1)    = tag(iabs)
      x(1)      = x(iabs)
      y(1)      = y(iabs)
      z(1)      = z(iabs)

      do 30 i=1,ndopx
        dopant(iabs,i) = dopant(ibasis,i)
 30   continue 
      tag(iabs)    = tag(ibasis)
      x(iabs)      = x(ibasis)
      y(iabs)      = y(ibasis)
      z(iabs)      = z(ibasis)

      do 40 i=1,ndopx
        dopant(ibasis,i) = etoss(i)
 40   continue 
      tag(ibasis)    = ttoss
      x(ibasis)      = xtoss
      y(ibasis)      = ytoss
      z(ibasis)      = ztoss

c  remove null atom site from list
      el = dopant(ibasis,1)
      call case(test,el)
      if (el.eq.'nu') ibasis = ibasis-1

c  shift coordinates
      do 100 i=1,ibasis
        x(i) = x(i) - xshift
        y(i) = y(i) - yshift
        z(i) = z(i) - zshift
 100  continue

c  reset iabs to reflect the switch just made
      iabs = 1

 999  continue
      return
c  end subroutine basort
      end

