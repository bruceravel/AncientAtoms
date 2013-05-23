      subroutine rh2hex(iat,iatom,cell,x,y,z)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision(a-h,o-z)
c----------------------------------------------------------------------
c  converts from a rhombohedral basis to a hexagonal basis
c  the hexagonal a and c are calculated from the rhombohedral a and 
c  alpha, and the hexagonal gamma is 120.
c  the new cell constants and angles are calculated from the old
c  and the unique atom positions are changed as indicated by the itxc
c----------------------------------------------------------------------
c  input:
c    iat:   dimension of x,y,z
c    iatom: number of unique positions
c  input/output:
c    cell:  array containing a,b,c,alpha,beta, and gamma
c    x,y,z: fractional coordinates of unique positions
c----------------------------------------------------------------------

c      parameter (iat=50)
      parameter (one = 1.e0, two = 2.e0, three = 3.e0, one80 = 180.e0)
      parameter (third=one/three, twoth=two/three)
      parameter (thirdm=-third, twothm=-twoth)
      parameter (pi=3.141592653589793238462643)

      dimension cell(6),x(iat),y(iat),z(iat),convrt(3,3)

      data ((convrt(i,j),j=1,3),i=1,3) /twoth,thirdm,thirdm,
     $                                  third, third,twothm,
     $                                  third, third, third/

      a      = cell(1)
      alpha  = cell(4)*pi/one80

c----------------------------------------------------------------------
c  the hexagonal a is the third side of the triangle formed by two of 
c  the rhombohedral cell axes, a, aprime/2, and alpha/2 form a right 
c  triangle
      aprime = 2 * a * sin(alpha/2)

c----------------------------------------------------------------------
c  the hexagonal c is the long diagonal of the rhombohedron.  the 
c  following is repeated uses of the law of cosines.  csqr is the 
c  square of the short face diagonal.  bsqr is the square of the long
c  face diagonal.  gsqr is the square of the short body diagonal.  cosg
c  is the cosine of the angle between the long face diagonal and the 
c  opposing axis.  this is an ugly mess.
      bsqr = a**2 * (2 + 2*cos(alpha))
      csqr = a**2 * (2 - 2*cos(alpha))
      gsqr = csqr + a**2
      cosg = (bsqr+a**2-gsqr) / (2*a*sqrt(bsqr))

      cprime = sqrt( bsqr + a**2 + 2*a*sqrt(bsqr)*cosg )

c      print*,'aprime,cprime: ',aprime,cprime

c----------------------------------------------------------------------
c  repack cell with new a,b,c,alpha,beta,gamma
      cell(1) = aprime
      cell(2) = aprime
      cell(3) = cprime
      cell(4) = 90.
      cell(5) = 90.
      cell(6) = 120.

c----------------------------------------------------------------------
c  now fix the unique atom fractional coordinates
c  convrt is the metric tensor connecting a rhomohedron with its 
c  equivalent heagonal prism

      do 20 i=1,iatom
        xnew = x(i)*convrt(1,1)+y(i)*convrt(1,2)+z(i)*convrt(1,3)
        ynew = x(i)*convrt(2,1)+y(i)*convrt(2,2)+z(i)*convrt(2,3)
        znew = x(i)*convrt(3,1)+y(i)*convrt(3,2)+z(i)*convrt(3,3)
c        print 400,'rh2hex: old x,y,z= ',x(i),y(i),z(i)
c        print 400,'rh2hex: new x,y,z= ',xnew,ynew,znew
c400     format(a,3(2x,f8.4))
        x(i) = xnew
        y(i) = ynew
        z(i) = znew
20    continue

30    continue
      return 
c end subroutine rh2hex
      end

