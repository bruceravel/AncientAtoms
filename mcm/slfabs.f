      subroutine slfabs(iat, ndopx, iatom, ipt, iedge, idop, core,
     $            dopant, percnt, lself,
     $            ampslf, sigslf, qrtslf, xmuf, xmub)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c  calculate the self absorption corrections for a measurement taken in
c  in fluorescence.  this correction is calculated assuming the crystal
c  is fully concentrated, i.e. not diluted in any medium.  the correction
c  is needed to compensate for the variable absorption depth due to the
c  chi(E) wiggles above the edge.  it is well approximated by regressing
c  a cubic polynomial in ln(E) to McMaster data for the full sample.
c  see Chapter 10 for details
c  corrections are given as a constant, sigma^2 and fourth cumulant.
c----------------------------------------------------------------------
c  input
c     iat, ndopx: parameters set in calling program
c     iatom:  number of unique atoms
c     ipt:    (iat) number of repititions of each unique atom,
c             array of length iat
c     iedge:  desired edge, defaulted by z, currently k or l3
c     core*2: core atom symbol
c     dopant*2: (iat, ndopx) symbol of doping element
c     percnt: (iat, ndopx) real, percentage of dopant
c     lself:  undocumented diagnostic file flag, logical
c  output
c     ampslf: amplitude factor from self absorption, real
c     sigslf: sigma2 for self absorption correction, real
c     qrtslf: sigma4 (quartic term) for self absorption correction, real
c     muf:    absoprtion of whole sample at fluorescent energy
c     mub:    absorption of rest of sample + non-res electrons at E0
c----------------------------------------------------------------------
c      parameter (iat=50, ndopx=4)
      parameter (iou = 10)
      parameter (nvals = 50, nfit = 5)
      parameter (etok=0.26246 82917, eps = 0.0001)
c  nvals: # of points in regression
c  npre:  # of points in preedge regression
c  nfit:  max order of polynomial
c  etok:  conversion btwn ev and invang
c  eps:   small number for floating point logic

      logical      lerr,lself
c  lerr:  used for error checking in mucal
c  lself: true means to write out the numbers in the self absorption
c         correction calculation.  this is undocumented and will
c         always be that way.
      character*2  units*1,dopant(iat,ndopx), dp, test
      character*8  fname
      character*10 core
      dimension    xlnxmu(nvals), qsqr(nvals), afit(nfit), xmurat(nvals)
      dimension    ipt(iat), idop(iat)
      dimension    energy(9), xsec(10), percnt(iat,ndopx)

c----------------------------------------------------------------------
c  get edge energy of core atom,  e0 = edge energy
      test   = 'ab'
      e0     = -1000
      units  = 'b'
      lerr   = .false.
      iz     = is2z(core)
      call mucal(e0,core,iz,units,xsec,energy,lerr,ier)
      e0     = energy(iedge)
      e0m10  = e0 - 0.01e0
      efluor = energy(6)
      if (iedge.ge.2) efluor = energy(8)
c      print*,'in slfabs: core,iedge=',core,iedge

c     estep: energy step in kev, elimit: end of post-edge region = .5kev
      estep  = .010
      elimit = nvals * estep
c      print*,'elimit = ',elimit

c----------------------------------------------------------------------
c  check to see if post-edge region will run into another
c  absorption edge in the problem.  if it does, stop post-edge
c  region 10 volts before that edge by resetting estep
      do 50 i = 1, iatom
        do 40 j=1,idop(i)
          iz = is2z(dopant(i,j))
          call mucal(e, dopant(i,j), iz,units,xsec,energy,lerr,ier)
          do 30 k=1,5
            if ( (energy(k) - e0).gt.eps ) then
                echeck = (energy(k) - e0) - .010
                elimit = min(elimit, echeck)
            endif
 30       continue
 40     continue
 50   continue

      estep = elimit/nvals

c----------------------------------------------------------------------
c  calculate total absorption at the flourescent energy and 10 volts
c  below the edge --> xmuf and xmub
      ier  = 0
      xmuf = 0
      do 150 j=1,iatom
        do 140 k=1,idop(j)
          dp = dopant(j,k)
          call case(test,dp)
          if (dp.eq.'nu') then
              goto 140
          else
              iz = is2z(dopant(j,k))
          endif
          call mucal(efluor,dopant(j,k),iz,units,xsec,energy,lerr,ier)
          xmuf  = xmuf + ipt(j)*xsec(4)*percnt(j,k)
          call mucal(e0m10, dopant(j,k),iz,units,xsec,energy,lerr,ier)
          xmub  = xmub + ipt(j)*xsec(4)*percnt(j,k)
 140    continue
 150  continue

c----------------------------------------------------------------------
c  calculate total mu above the edge at several energy values
c  and keep log(mu - preedge extrapolation) v. k^2 (k in inverse angs)
      ier    = 0
      do 300 i = 1, nvals
         ener    = i * estep
         e       = e0 + ener
c --- xmu is the absorption of the rest of the sample
c --- xmucor is the absorption of the resonant atom
         xmu     = 0.0
         xmucor  = 0.0
         do 250 j = 1, iatom
           do 240 k=1,idop(j)
             if (dopant(j,k).eq.'nu') then
                 iz = 1
             else
                 iz = is2z(dopant(j,k))
             endif
             call mucal(e, dopant(j,k),iz,units,xsec,energy,lerr,ier)
c%%%%%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%
c%%%  the old, bad code.  I leave it here in case I need to know how
c%%%  big a mistake i made before
c%%%               if (dopant(j,k).eq.core) then
c%%%                   xmucor = xmucor + ipt(j)*xsec(4)*percnt(j,k)
c%%%               endif
c%%%               xmu  = xmu  + ipt(j)*xsec(4)*percnt(j,k)
c%%%  to test it, uncomment these lines and comment out the next 5
c%%%%%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%wrong%%%%
             if (dopant(j,k).eq.core) then
                 xmucor = xmucor + ipt(j)*xsec(4)*percnt(j,k)
             else
                 xmu  = xmu  + ipt(j)*xsec(4)*percnt(j,k)
             endif
 240       continue
 250     continue
c         if (i.eq.1) then
c             print*,'xmu=',xmu/ipt(1)
c             print*,'xmucor=',xmucor/ipt(1)
c         endif

         xmurat(i) = (xmuf + xmu + xmucor) / (xmuf + xmu)
c         if (xmurat(i).le.0) xmurat(i) = 0.0001
         xlnxmu(i) = log( xmurat(i) )
         qsqr(i)   = etok * ener * 1000.e0
 300  continue

c--------------------------------------------------------------------
c  fit log(xmucore / xmutotal) v. energy with a quadratic
c  the linear coef is the mcmaster sigma2, the quad. coef. is the
c  mcmaster sigma4
      nterms = 3
      do 310 i=1,nfit
        afit(i) = 0.0
310   continue
      call polyft(qsqr(1),qsqr(nvals),qsqr,xlnxmu,nvals,nterms,afit)
      ampslf =   exp(afit(1))
      sigslf = - afit(2) / 2.0
      qrtslf =   afit(3) * 3.0 / 2.0

      call mucerr(ier)

c  undocumented diagnostic file
c  if you do not already know what this file is for, then you don't want
c  to use it.
      if (lself) then
          fname = 'self.dat'
          call lower(fname)
          open (unit=iou,file=fname,status='unknown')
          write(iou,*)'afit(1,2,3)',afit(1),afit(2),afit(3)
          write(iou,*)'converted',ampslf,sigslf,qrtslf
          do 1000 i=1,nvals
            ener = i * estep
            e    = e0 + ener
            val  = exp( afit(1) + qsqr(i)*(afit(2) + afit(3)*qsqr(i)) )
            write(iou,400)e,xmurat(i),val
400         format(3(2x,f10.6))
1000      continue
          close(iou)
      endif

      return
c end subroutine slfabs
      end
