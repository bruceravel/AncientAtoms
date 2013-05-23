
      function volume(cell)
c      double precision function volume(cell)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)

c--------------------------------------------------------------
c  calculate the volume of any primitive cell
c
c  cell:   real array (6) containing a,b,c,alpha,beta,gamma
c--------------------------------------------------------------
      parameter (pi=3.141592653589793238462643e0)
      parameter (radian = pi/180.e0)

      dimension cell(6)
      
      cosal = cos( cell(4)*radian )
      cosbe = cos( cell(5)*radian )
      cosga = cos( cell(6)*radian )

c  term = 1 for orthogonal cells     
      term = 1 - cosal**2 - cosbe**2 - cosga**2 + 
     $       2*cosal*cosbe*cosga

      volume = cell(1)*cell(2)*cell(3)*sqrt(term)

      return
c end function volume
      end
