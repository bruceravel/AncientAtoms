      subroutine mcmast(core,iedge,sigmm,qrtmm,lmm)
c---------------------------------------------------------------------------
c  copyright 1993 university of washington     matt newville and bruce ravel
c---------------------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c  calculate the mcmaster correction for the central atom.  e and e^2
c  correction terms are given as sigma^2 and fourth cumulant.
c----------------------------------------------------------------------
c  input
c     core*2: central atom in problem
c     iedge:  desired edge, defaulted by z, currently k or l3
c     lmm:    undocumented logical keyword, write mcmast.dat, logical
c  output
c     sigmm:  sigma2 for mcmaster correction, real
c     qrtmm:  sigma4 (quartic term) for mcmaster correction, real
c----------------------------------------------------------------------
      parameter ( iou = 10 )
      parameter ( nvals = 50, npre = 10, nfit = 5 )
      parameter ( etok=0.26246 82917 )
c  nvals: # of points in postedge regression
c  npre:  # of points in preedge regression
c  nfit:  max order of polynomial
c  etok:  conversion btwn ev and invang

      logical      lerr,lmm
c  lerr: used for error checking in mucal
c  lmm : true means to write out the numbers in the mcmaster correction
c        calculation.  this is undocumented and will always be that way.
      character*2  units*1,core
      character*10 fname
      dimension    xlnxmu(nvals), qsqr(nvals), xmupre(npre), epre(npre)
      dimension    energy(9), xsec(10), afit(nfit)

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

c----------------------------------------------------------------------
c  calculate total mu below the edge (200ev to 100ev),

c  undocumented diagnostic file
c  if you do not already know what this file is for, then you don't want
c  to use it.
      if (lmm) then
          fname = 'mcmast.dat'
          call lower(fname)
          open(unit = iou, file=fname, status= 'unknown')
          write(iou,*) ' iz = ', iz
          write(iou,*) ' epre(i), xmupre(i) '
      endif

      do 100 i = 1, npre
         epre(i)   = e0 - 0.200 + (i-1) * estep
         xmu       = 0.0
         iz        = is2z(core)
         call mucal(epre(i),core,iz,units,xsec,energy,lerr,ier)
         xmupre(i) = xsec(4)

         if (lmm) write(iou,*) epre(i), xmupre(i)
100   continue

c----------------------------------------------------------------------
c  fit pre-edge with a straight line.
      nterms = 2
      do 110 i=1,nfit
        afit(i) = 0.0
110   continue
      call polyft(epre(1),epre(npre),epre,xmupre,npre,nterms,afit)
      bpre  = afit(1)
      slope = afit(2)

      if (lmm) then
          write(iou,*) 'slopre, bpre = ', afit(2), afit(1)
          write(iou,*) ' '
      endif

c----------------------------------------------------------------------
c  calculate total mu above the edge at several energy values
c  and keep log(mu - preedge extrapolation) v. k^2 (k in inverse angs)
      ier    = 0
      if (lmm) write(iou,*) ' qsqr(i), delxmu(i) , xlnxmu(i)'
      do 300 i = 1, nvals
         ener = i * estep
         e    = e0 + ener
         iz   = is2z(core)
         call mucal(e,core,iz,units,xsec,energy,lerr,ier)
         xmu  = xsec(4)

         delxmu    = xmu - (bpre + e*slope)
         if (delxmu.le.0) delxmu = 0.0001
         xlnxmu(i) = log( delxmu )
         qsqr(i)   = etok * ener * 1000.000
         if (lmm) write(iou,*)  qsqr(i), delxmu  , xlnxmu(i)
 300  continue

c--------------------------------------------------------------------
c  fit log(mu - preedge extrapolation) v. energy with a quadratic
c  the linear coef is the mcmaster sigma2, the quad. coef. is the
c  mcmaster sigma4
      nterms = 3
      do 310 i=1,nfit
        afit(i) = 0.0
310   continue
      call polyft(qsqr(1),qsqr(nvals),qsqr,xlnxmu,nvals,nterms,afit)
      sigmm = - afit(2) / 2.0
      qrtmm =   afit(3) * 3.0 / 2.0

      if (lmm) then
          write(iou,*) 'qrtmm, sigmm, b = ', qrtmm, sigmm, afit(1)
          write(iou,*) ' '
          write(iou,*) ' fit results:'
          write(iou,*) ' qsqr , delxmu, log(xmu) '
          do 555 i = 1, nvals
             ener    = i * estep
             qsqr(i) = etok * ener * 1000.000
             temp    = afit(1) + qsqr(i) * (afit(2) + afit(3)*qsqr(i))
             fit     = exp( temp )
             write(iou,*)  qsqr(i), fit, temp
555       continue
          close(iou)
      endif

      call mucerr(ier)

      return
c end subroutine mcmast
      end
