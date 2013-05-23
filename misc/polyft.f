       subroutine polyft(xfit1,xfit2,xdata,ydata,ndata,nterms,aout)
c
c  get coefficients for polynomial fit :
c      ydata = aout(1) + aout(2)*xdata  + aout(3) *xdata^2 + ...
c  the fit is done between xdata = xfit1 and xfit2
c
c  inputs :
c    xfit1    lower bound of fitting range
c    xfit2    upper bound of fitting range
c    xdata    array of abscissa values for data
c    ydata    array of ordinate values for data
c    ndata    length of data arrays
c    nterms   number of terms in polynomial
c
c  outputs :
c    aout     coefficients of fitted polynomial to data
c
c   copyright 1992  university of washington :          matt newville
c
c  requires function nofx
c
c  see bevington pg 104 for expalanation of these variables
c
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

       parameter (max= 5, max2m1 = 2*max-1, zero = 0.)
       dimension         xdata(ndata), ydata(ndata), aout(nterms)
       double precision  sumx(max2m1), sumy(max)
       double precision  array(max,max), ain(max), delta, determ
       external          determ
c
c find points closest to endpoints of fitting range
       nfit1 = nofx(xfit1,xdata,ndata)
       nfit2 = nofx(xfit2,xdata,ndata)
       if (nfit1.gt.nfit2) then
            ntemp = nfit1
            nfit1 = nfit2
            nfit2 = ntemp
       elseif(nfit1.eq.nfit2) then
            go to 300
       end if
c
c   initialize internal arrays
       nmax   = 2 * nterms - 1
       do 100 i=1, nmax
          sumx(i) = zero
 100   continue
       do 110 i = 1, nterms
          ain(i) = zero
          sumy(i) = zero
          do 110 j = 1,  nterms
            array(i,j) = zero
  110  continue
c
c  collect sums of data, sum of squares of data, etc.
       do 200 i = nfit1, nfit2
          xi = xdata(i)
          yi = ydata(i)
          xterm = 1.0
          do 180 n=1, nmax
             sumx(n) = sumx(n) + xterm
             xterm   = xterm * xi
  180     continue
          yterm = yi
          do 190 n=1,nterms
             sumy(n) = sumy(n) + yterm
             yterm   = yterm * xi
  190     continue
  200  continue
c
c construct matrices and evaluate coefficients
c
       do 210 j=1,nterms
         do 210 k=1,nterms
            array(j,k) = sumx(j + k - 1)
  210  continue
       delta = determ(array,nterms,max)
       if (delta.ne.zero) then
           do 260 l=1,nterms
              do 250 j=1,nterms
                 do 240 k=1,nterms
                    array(j,k) = sumx(j+k-1)
  240            continue
                 array(j,l) = sumy(j)
  250         continue
              ain(l) = determ(array,nterms,max)/delta
  260      continue
       end if
  300  continue
       do 400 i = 1, nterms
          aout(i) = sngl(ain(i))
  400  continue
       return
c end  subroutine polyft
       end
