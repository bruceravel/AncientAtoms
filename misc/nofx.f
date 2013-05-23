      function  nofx(x,array,npts)

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

c
c   function nofx
c
c   purpose
c     given a value x and an array of values, find the index
c     corresponding to the array element closest to x
c
c   usage
c     n = nofx(x,array,npts)
c
c   parameters
c     x     - a given value
c     array - array of values, assumed to be stored in order of
c             increasing value
c     npts  - number of elements in array
c
c   subroutines and function subprograms required
c     none
c
c   written  8/11/81 by j.m. tranquada
c
      dimension  array(npts)
      imin = 1
      imax = npts
      inc = ( imax - imin ) / 2
   10 continue
      it  = imin + inc
      xit = array(it)
      if ( x .lt. xit ) then
         imax = it
      else if ( x .gt. xit ) then
         imin = it
      else
         nofx = it
         return
      endif
      inc = ( imax - imin ) / 2
      if ( inc .gt. 0 ) go to 10
      xave = ( array(imin) + array(imin+1) ) / 2.
      if ( x .lt. xave ) then
         nofx = imin
      else
         nofx = imin + 1
      endif
      return
c end function nofx
      end
