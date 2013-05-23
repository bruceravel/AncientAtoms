      subroutine spcing(qvect,cell,d)
c----------------------------------------------------------------
c  calculate interplanar spacing for any reflection in any cell
c
c  input:
c    qvect:  (3) reflection, real
c    cell:   (6) array of a,b,c,alpha,beta,gamma; real
c  output:
c    d:      interplanar spacing, real
c----------------------------------------------------------------
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      parameter (one80 = 180., epsi=0.0001)
      parameter (pi=3.141592653589793238462643)

      dimension qvect(3), cell(6)

      v     = volume(cell)
      alpha = cell(4)*pi/one80
      beta  = cell(5)*pi/one80
      gamma = cell(6)*pi/one80
      cosum = abs(cos(alpha)) + abs(cos(beta)) + abs(cos(gamma))

c  these are all the yucky angular terms that are needed handle a 
c  general unit cell.
c  notice that the cross terms all vanish with right angles
      s11 = (cell(2)*cell(3)*sin(alpha))**2
      s22 = (cell(1)*cell(3)*sin(beta))**2
      s33 = (cell(1)*cell(2)*sin(gamma))**2
      s12 = 0.
      s23 = 0.
      s13 = 0.
      if (cosum.gt.epsi) then
          s12 =  cell(1)*cell(2)*(cell(3)**2)*
     $                (cos(alpha)*cos(beta)-cos(gamma))
          s23 =  cell(2)*cell(3)*(cell(1)**2)*
     $                (cos(beta)*cos(gamma)-cos(alpha))
          s13 =  cell(3)*cell(1)*(cell(2)**2)*
     $                (cos(gamma)*cos(alpha)-cos(beta))
      endif

      sumsqr = s11*(qvect(1)**2)        + s22*(qvect(2)**2) +
     $         s33*(qvect(3)**2)        + 2.*s12*qvect(1)*qvect(2) +
     $         2.*s23*qvect(2)*qvect(3) + 2.*s13*qvect(1)*qvect(3)

      d = v / sqrt(sumsqr)

      return 
c  end subroutine spcing
      end

