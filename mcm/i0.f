      subroutine i0(gasses,core,iedge,li0,sigi0,qrti0)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c  calculate the i0 corrections for a measurement taken in
c  in fluorescence.
c  this correction is needed to compensate for the energy dependence of
c  the gasses in the i0 chamber.  since fluorescent photons are a single
c  energy, this energy dependence is not divided out when mu is constructed.
c  it is well approximated by a cubic polynomial regression in ln(E) to
c  the absorption of the gasses in the i0 chamber
c  corrections are given as sigma^2 and fourth cumulant.
c----------------------------------------------------------------------
c  input
c     pargon: percent by pressure of argon in i0 chamber, real
c     pnitro: percent by pressure of nitrogen in i0 chamber, real
c     pnitro: percent by pressure of krypton in i0 chamber, real
c     core*2: symbol of central atom
c     iedge:  desired edge, defaulted by z, currently k or l3
c     li0:    undocumented flag to write out i0 diagnostic file, logical
c  output
c     sigi0:  sigma2 for i0 correction, real
c     qrti0:  sigma4 (quartic term) for i0 correction, real
c----------------------------------------------------------------------
c      implicit double precision (a-h,o-z)

      parameter ( iou = 10 )
      parameter ( nvals = 50, nfit = 5 )
      parameter ( etok=0.26246 82917 )
c  nvals: # of points in regression
c  nfit:  max order of polynomial
c  etok:  conversion btwn ev and invang

      logical      lerr,li0
c  lerr:  used for error checking in mucal
c  lself: true means to write out the numbers in the self absorption
c         correction calculation.  this is undocumented and will
c         always be that way.
      character*2  units*1,core,el
      character*6  fname
      dimension    xlnxmu(nvals), qsqr(nvals)
      dimension    energy(9), xsec(10), gasses(3)
      dimension    afit(nfit)

      pargon = gasses(1)
      pkrypt = gasses(2)
      pnitro = gasses(3)
c----------------------------------------------------------------------
c  get edge energy of core atom,  e0 = edge energy
      e0     = -1000
      units  = 'b'
      lerr   = .false.
      iz     = is2z(core)
      call mucal(e0,core,iz,units,xsec,energy,lerr,ier)
      e0     = energy(iedge)

c     estep: energy step in kev, elimit: end of post-edge region = .5kev
      estep  = .010
      elimit = nvals * estep
c      print*,'elimit = ',elimit

c----------------------------------------------------------------------
c  get pressure percentage of helium from argon and nitrogen
c  percentages.  then convert percentages to reflect the number of
c  absorbers.  recall that these gases behave like ideal gases at room
c  temp. and pressure -- nitrogen is diatomic, the other two are
c  monoatomic, thus the percentages of absorbers is different from the
c  pressure percentages.
      phelm  = 1 - pargon - pnitro - pkrypt
      pnorm  =      phelm + pargon + 2*pnitro + pkrypt
      phelm  =      phelm / pnorm
      pargon =     pargon / pnorm
      pnitro =   2*pnitro / pnorm
      pkrypt =     pkrypt / pnorm

c----------------------------------------------------------------------
c  calculate total mu above the edge at several energy values
c  and keep log(mu - preedge extrapolation) v. k^2 (k in inverse angs)
      ier    = 0
      do 300 i = 1, nvals
         ener    = i * estep
         e       = e0 + ener

         el = 'h'
         iz = is2z(el)
         call mucal(e,el,iz,units,xsec,energy,lerr,ier)
         xmuhe   = xsec(4)

         el = 'ar'
         iz = is2z(el)
         call mucal(e,el,iz,units,xsec,energy,lerr,ier)
         xmuar   = xsec(4)

         el = 'n'
         iz = is2z(el)
         call mucal(e,el,iz,units,xsec,energy,lerr,ier)
         xmun    = xsec(4)

         el = 'kr'
         iz = is2z(el)
         call mucal(e,el,iz,units,xsec,energy,lerr,ier)
         xmukr   = xsec(4)

         xmu = phelm*xmuhe + pargon*xmuar + pnitro*xmun + pkrypt*xmukr

         xlnxmu(i) = log( xmu )
         qsqr(i)   = etok * ener * 1000.000
 300  continue

c--------------------------------------------------------------------
c  fit log(xmu) v. energy with a quadratic
c  the linear coef is the i0 sigma2, the quad. coef. is the i0 sigma4
      nterms = 3
      do 310 i=1,nfit
        afit(i) = 0.0
310   continue
      call polyft(qsqr(1),qsqr(nvals),qsqr,xlnxmu,nvals,nterms,afit)
      sigi0 = - afit(2) / 2.0
      qrti0 =   afit(3) * 3.0 / 2.0

      call mucerr(ier)

c  undocumented diagnostic file
c  if you do not already know what this file is for, then you don't want
c  to use it.
      if (li0) then
          fname = 'i0.dat'
          call lower(fname)
          open (unit=iou,file=fname,status='unknown')
          write(iou,*)'afit[2,3]',afit(2),afit(3)
          write(iou,*)'converted',sigi0,qrti0
          do 1000 i=1,nvals
            ener = i * estep
            e    = e0 + ener
            val  = afit(1) + qsqr(i)*(afit(2) + afit(3)*qsqr(i))
            write(iou,400)e,xlnxmu(i),val
400         format(3(2x,f10.6))
1000      continue
          close(iou)
      endif

      return
c end subroutine i0
      end
