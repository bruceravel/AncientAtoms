      subroutine genmul(ng,isymce,ibravl,igen)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------
c  calculates general multiplicity (igen) of unit cell, i.e. the 
c  number of times all of the symmetries of a cell generate a new 
c  point from any arbitrary point
c-------------------------------------------------------------------
c  input:
c    ng     = multiplicity before bravais translations and 
c             centrosymmetry operation
c    isymce = centrosymmetry flag
c    ibravl = bravais lattice type, 1..7 = p i r f a b c
c  output:
c    igen   = general multiplicity, product of simple multiplicity, 
c             mult. of bravais lattice, and mult. of centrosymmetry
c-------------------------------------------------------------------
c    icmf   = mult. factor due to centrosymmety
c    ibrmf  = mult. factor due to bravais lattices
c    igen   = total multiplicity
c-------------------------------------------------------------------

      icmf  = 1
      igen  = ng
      ibrmf = 1
c                                          if centrosymmetric...
      if (isymce.eq.1) then
          icmf = 2
          igen  = ng*icmf
      endif
c                                          2=i, body center
c                                          3=r, rhombohedral
c                                          4=f, face center
      if ((ibravl.ge.2).and.(ibravl.le.4)) then
          ibrmf = ibravl
          igen  = ng*icmf*ibrmf
c                                          5..7=a,b,c, single face center
      elseif ((ibravl.ge.5).and.(ibravl.le.7)) then
          ibrmf = 2
          igen   = ng*icmf*ibrmf
      endif

      return
c end subroutine genmul
      end
