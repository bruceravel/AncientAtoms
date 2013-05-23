      subroutine metric (cell,trmtx)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c------------------------------------------------------------------
c  calculate metric tensor (trmtx) used in transforming from
c  a cell axis basis to an orthogonal basis
c------------------------------------------------------------------
c  input:
c    cell:  (6) array containing a,b,c,alpha,beta,gamma
c  output:
c    trmtx: (3,3) metric tensor
c------------------------------------------------------------------
      dimension cell(6),co(3),si(3),trmtx(3,3),cosqr(3)
      parameter (pi = 3.14159265358979323844e0)
      parameter (radian = 180.e0/pi)
      parameter (zero = 0.e0, one = 1.e0)
      parameter (eps=1.e-6)
c------------------------------------------------------------------
c  calculate and store sines and cosines of the cell angles
      do 10 i=1,3
        co(i)    = cos(cell(i+3)/radian)
        si(i)    = sin(cell(i+3)/radian)
        cosqr(i) = co(i)**2
 10   continue

c------------------------------------------------------------------
c  calculate various trigonometric quantities for use in the three 
c  dimensional transformation
      cosxx = (co(1)*co(3) - co(2)) / (si(1)*si(3))
      cosyy = (co(1)*co(2) - co(3)) / (si(1)*si(2))
      if ((one-cosxx**2).lt.eps) then
          sinxx = 0.e0
      else
          sinxx = sqrt(one - cosxx**2)
      endif
      if ((one-cosyy**2).lt.eps) then
          sinyy = 0.e0
      else
          sinyy = sqrt(one - cosyy**2)
      endif

c------------------------------------------------------------------
c  evaluate the transformation matrix elements
      trmtx(1,1) = sinyy*si(2)
      trmtx(1,2) = zero
      trmtx(1,3) = zero

      trmtx(2,1) = -((cosyy/(sinyy*si(1)))+(co(1)*cosxx)/(sinxx*si(1)))
     $              *(sinyy*si(2))
      trmtx(2,2) = one
      trmtx(2,3) = co(1)

      trmtx(3,1) = -(cosxx*sinyy*si(2))/sinxx
      trmtx(3,2) = zero
      trmtx(3,3) = si(1)

      return
c end subroutine metric
      end

