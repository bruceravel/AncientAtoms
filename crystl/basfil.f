      subroutine basfil(iat,ibasis,x,y,z,ipt,st)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------
c  after the first atom in a basis list has been expanded according
c  to the group summetries, fill in the remainder of the atoms in 
c  the basis at each point in the unit cell
c----------------------------------------------------------------
c  input:
c    iat:     parameter set in calling program
c    ibasis:  number of atoms in basis
c    x,y,z:   (iat) fractional coordinates of representative atoms
c  i/o:
c    ipt:     (iat) multiplicity of each site (all equal for a basis)
c    st:      (iat,192,3) fractional coordinates of all atoms in unit cell
c----------------------------------------------------------------
c      parameter (iat=50)

      dimension x(iat), y(iat), z(iat), st(iat,192,3), ipt(iat)

      do 50 i=2,ibasis
        do 40 j=1,ipt(1)
          st(i,j,1) = st(1,j,1) + x(i)
          st(i,j,2) = st(1,j,2) + y(i)
          st(i,j,3) = st(1,j,3) + z(i)
 40     continue 
        ipt(i) = ipt(1)
 50   continue 

      return
c  end subroutine basfil
      end

