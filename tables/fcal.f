      subroutine fcal(elem,e,f1,f2)

c================================================================
c  get a values for f' and f" given an element and an energy
c================================================================
c  input:
c    elem: atomic symbol of desired element, character*2
c    e:    energy in eV, real
c  output:
c    f1:   value of f' at energy e for element elem, real
c    f2:   value of f" at energy e for element elem, real
c================================================================
c  the values are obtained by doing a four term interpolation
c  from tabulated data (either sasaki_block.f or chantler_block.f)
c  care is taken to avoid interpolating through steps and cusps
c================================================================

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      parameter (zero=0)
      character*12 cnum
      character*2  elem, el, edge, test
      logical      there
c------------------------------------------------------------------------
c  from f'/f" block data 
      parameter(nelem=92, ndatx=233)
c      ndatx is 86 for chant., ndatx is 233 for sas.)
      common /fdat/ nfdata(nelem), engrd(nelem,ndatx),
     $            fp(nelem,ndatx), fpp(nelem,ndatx)
c------------------------------------------------------------------------
      dimension    fpin(ndatx),  fppin(ndatx),  ein(ndatx)
      dimension    egrid(ndatx)

 4000 format(f12.5)

      el = elem
      test = 'ab'
c  return 0 for null site
      call case(test, el)
      if (el.eq.'nu') then
          f1 = 0
          f2 = 0
          goto 999
      endif
c  be sure the element called is available
      call avlble(el, there)
      if ( (.not.there).and.(el.ne.'nu') ) goto 666
c  values are tabulated between 3 and 60 keV
      if ( (e.lt.3000) .or. (e.gt.60000) ) goto 667

      iz   = is2z(el) 
      ndat = nfdata(iz)

c  transfer energy grid into an array for convenience below
      do 10 i=1,ndat
        egrid(i) = engrd(iz,i)
 10   continue 

c-----------------------------------------------------------------------
c  want to interpolate but avoid going through resonant energies in the
c  interpolation.  find proper boundries for interpolation by comparing 
c  input energy to edge energies
c  careful that closest energy is not on wrong side
c  s2e is tabulated in keV so mult by 1000

      edge= 'k'
      ek  = s2e(el, edge) * 1000
      nk  = nofx(ek, egrid, ndat)
      if (e.ge.ek) then
          nfirst = nk
          if (ek.gt.egrid(nk))  nfirst = nfirst+1
          nlast  = ndat
          goto 40
      endif

      edge= 'l1'
      el1 = s2e(el, edge) * 1000
      nl1 = nofx(el1, egrid, ndat)
      if ( (e.lt.ek).and.(e.ge.el1) ) then
          nfirst = nl1
          if (el1.gt.egrid(nl1)) nfirst = nfirst+1
          nlast  = nk
          if (ek .lt.egrid(nk))  nlast  = nlast-1
          goto 40
      endif

      edge= 'l2'
      el2 = s2e(el, edge) * 1000
      nl2 = nofx(el2, egrid, ndatx)
      if ( (e.lt.el1).and.(e.ge.el2) ) then
          nfirst = nl2
          if (el2.gt.egrid(nl2)) nfirst = nfirst+1
          nlast  = nl1
          if (el1.lt.egrid(nl1)) nlast  = nlast-1
          goto 40
      endif

      edge= 'l3'
      el3 = s2e(el, edge) * 1000
      nl3 = nofx(el3, egrid, ndatx)
      if ( (e.lt.el2).and.(e.ge.el3) ) then
          nfirst = nl3
          if (el3.gt.egrid(nl3)) nfirst = nfirst+1
          nlast  = nl2
          if (el2.lt.egrid(nl2)) nlast  = nlast-1
      else
          nfirst = 1
          nlast  = nl3
          if (el3.lt.egrid(nl3)) nlast  = nlast-1
      endif

 40   continue 
      nval = nlast - nfirst + 1
      do 50 i=1, nval
        il       = i+nfirst-1
        ein(i)   = egrid(il)
        fpin(i)  = fp(iz,il)
        fppin(i) = fpp(iz,il)
 50   continue 

c-----------------------------------------------------------------------
c  interpolate for values, (META: want an interpolation with binary search.)
      nterms = 4
      call interp(ein, fpin,  nval, nterms, e, f1)
      call interp(ein, fppin, nval, nterms, e, f2)

 999  continue 
      return

 666  continue
c  error message for unsupported elements 
      call messag(el//' is unavailable for a0 calculation.')
      call messag('Returning 0 for fp and fpp.')
      f1 = 0
      f2 = 0
      goto 999

 667  continue 
c  error message for out-of-range energies
      call messag('Fcal only calculates for energies between '//
     $            '3000 and 60000 eV.')
      write(cnum,4000)e
      call messag(cnum//' is out of range')
      call messag('Returning 0 for fp and fpp.')
      f1 = 0
      f2 = 0
      goto 999

c  end subroutine fcal
      end

c*********************************************************************

      subroutine avlble(elem, there)

c=======================================================================
c  elements 3-83 and 92 are tabulated in Sas. 
c  all elements <= 92 are in Chant.
c  there = true if elem in table
c  there = false if elem not in table
c=======================================================================

      character*2 elem
      logical     there

      there = .false.
      ielem = is2z(elem)

c  sasaki      
      if ( ((ielem.ge.3).and.(ielem.le.83))
     $            .or. (ielem.eq.92) ) there = .true.

c  chantler
c      if (ielem.le.92) there = .true.

      return
c  end subroutine avlble
      end
