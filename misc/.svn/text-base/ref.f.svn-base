      function ref(x,y,z)
c      double precision function ref(x,y,z)
c====================================================================
c  hash cartesian coordinates, function returns a hash key based on 
c  ordered coordinate absolute values,  three digits of precision
c  (-3,0,4) gives same hash as (4,3,0) etc. but different from (0,5,0)
c====================================================================
c  input:
c    x,y,z: cartesian coordinates
c  output:
c    hash value fo cartesian coordinates
c====================================================================

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      parameter (bb=7.238e0, cc=12.092e0)
      parameter (mildix = 10000, icent = 100)

c  sort coordinates by size
      sum = abs(x) + abs(y) + abs(z)
      a   = min( abs(x), abs(y), abs(z) )
      c   = max( abs(x), abs(y), abs(z) )
      b   = sum - a - c

c  integer part of real number, nint might get rounded differently for 
c  numbers that are only different due to numerical precision.
      ia  = int(a*mildix)
      ib  = int(b*mildix)
      ic  = int(c*mildix)

      ia  = int(ia/icent)
      ib  = int(ib/icent)
      ic  = int(ic/icent)

      ref = ia + bb*ib + cc*ic
      return
c  end function ref
      end
