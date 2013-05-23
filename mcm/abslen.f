      subroutine abslen(iat,ndopx,
     $                  iatom,ipt,iedge,idop,core,dopant,
     $                  percnt, v, amu, delmu, spgrav)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c  calculate the absorption of the crystal.
c  recall that barn=10e-24 cm^2/atom and [v]=10e-24 cm^3, thus the
c  numerical factor cancels exactly.  yeah!
c  also the numerical factor cancels in the specific gravity calculation.
c----------------------------------------------------------------------
c  input
c     iat, ndopx: parameters set in calling program
c     iatom:  number of unique atoms
c     ipt:    number of repititions of each unique atom
c     iedge:  desired edge, defaulted by z
c     idop:   (iat) number of species at each site, 1=no dopants, 2+=dopants
c     core*2: core atom symbol
c     dopant*2: (iat,ndopx) symbol of doping element
c     percnt: (iat,ndopx) percnt of dopant, real
c     v:      volume of unit cell, real
c  output
c     amu:    mu above the edge, real
c     delmu:  delta mu of core, real
c     spgrav: the specific gravity of the material assuming the density
c             of the bulk is the same as the density of the cell, real
c----------------------------------------------------------------------

c      parameter(iat=50,ndopx=4)
      parameter(factor=1.66053)
c  atomic mass unit = 1.66053e-24 gram

      logical      lerr
      character*2  units*1, dopant(iat,ndopx), dp, test
      character*10 core,cr
      dimension    ipt(iat), idop(iat)
      dimension    energy(9), xsec(10)
      dimension    percnt(iat,ndopx)

c----------------------------------------------------------------------
c  get edge energy of core atom, w/ units=b mucal returns barns
      test  = 'ab'
      ener  = -1000
      units = 'b'
      lerr  = .false.
      cr    = core
      call case(test, cr)
      iz    = is2z(core)
      call mucal(ener,core,iz,units,xsec,energy,lerr,ier)
      ener  = energy(iedge)

c----------------------------------------------------------------------
c  calculate total mu above and below the edge and calculate the density
c  of the material.  mu calculated at +/- estep=50 ev.
c  outer loop over sites, inner loop over principle and doping atoms
      estep  = 0.05
      ier    = 0
      amu    = 0
      delmu  = 0
      camu   = 0
      cbmu   = 0
      spgrav = 0
      do 20 i=1,iatom
        do 15 j=1,idop(i)
          dp = dopant(i,j)
          call case(test,dp)
          if (dp.eq.'nu') then
              goto 15
          else
              iz = is2z(dopant(i,j))
          endif
          call mucal(ener+estep,dopant(i,j),iz,units,
     $               xsec,energy,lerr,ier)
          amu    = amu    + ipt(i)*percnt(i,j)*xsec(4)
          if (dp.eq.cr) then
              camu = camu + ipt(i)*percnt(i,j)*xsec(4)
              call mucal(ener-estep,cr,iz,units,xsec,energy,lerr,ier)
              cbmu = cbmu + ipt(i)*percnt(i,j)*xsec(4)
          endif
          spgrav = spgrav + ipt(i)*percnt(i,j)*xsec(7)
 15     continue
 20   continue

      amu    = amu / v
      delmu  = (camu - cbmu) / v
      spgrav = factor*spgrav / v

      call mucerr(ier)

      return
c end subroutine abslen
      end
