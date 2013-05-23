      subroutine dopfix(iat,ndopx,ndpmax,
     $                  iatom,ndop,tag,doplis,replcd,pclis,
     $                  elemnt,dopant,percnt,idop)

c-----------------------------------------------------------------------
c  inputs:
c    iat:    number of unique crystallographic sites
c    ndop:   total number of dopants
c    tag:    array of site tags
c    doplis: array of all doping atoms
c    replcd: array of all replaced atoms
c    pclis:  array of doping percentages
c    elemnt: array of atoms from atom list
c  output:
c    dopant: iat x ndopx matrix of doping atoms indexed by site
c    percnt: iat x ndopx matrix of doping percentages indexed by site
c    idop:   array indicating number of dopants at each site
c-----------------------------------------------------------------------
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

c      parameter(iat=50, ndopx=4, ndpmax=10)

      character*2  test
      character*2  doplis(ndpmax), dopant(iat,ndopx), elemnt(iat)
      character*10 tag(iat), tagup, replcd(ndpmax), repup
      dimension    pclis(ndpmax),  percnt(iat,ndopx)
      dimension    idop(iat)

      test = 'ab'
      
      do 100 i=1,ndop
        repup = replcd(i)
        call case(test,repup)
        do 50 j=1,iatom
          tagup = tag(j)
          call case(test,tagup)
          if (repup.eq.tagup) then
              idop(j) = idop(j) + 1
              dopant(j,idop(j)) = doplis(i)
              percnt(j,idop(j)) = pclis(i)
          endif
 50     continue
 100  continue 

      do 200 i=1,iatom
        dopant(i,1) = elemnt(i)
        sumdop = 0
        do 150 j=2,idop(i)
          sumdop = sumdop + percnt(i,j)
 150    continue 
        percnt(i,1) = 1. - sumdop 
 200  continue 

      return
c  end subroutine dopfix
      end

