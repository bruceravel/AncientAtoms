      function dist(vector,cell)
c      double precision function dist (vector,cell)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------
c   calculate distance between (0,0,0) and a position
c   in the cell axis basis, cosum checks for orthogonal axes
c
c input:
c   vector: (3) position in cell axis basis
c   cell:   (6) array of a,b,c,alpha,beta,gamma
c output:
c   dist:   distance from origin of the position
c-------------------------------------------------------------------
      parameter (zero=0.e0, epsi=0.0001e0)
      parameter (pi=3.141592653589793238462643e0)
      parameter (radian = pi/180.e0)
      dimension vector(3), cell(6)

      dist  = zero
      cosum = zero
      do 10 i=1,3
        dist  = dist + vector(i)**2
        cosum = cosum + abs(cos( cell(i+3)*radian ))
10    continue

c     correct for non-orthogonal cell axis basis
      if (cosum.gt.epsi) then
         dist = dist + 2*vector(1)*vector(2)*cos(cell(6)*radian)
     $               + 2*vector(1)*vector(3)*cos(cell(5)*radian)
     $               + 2*vector(2)*vector(3)*cos(cell(4)*radian)
      endif
      dist = sqrt(dist)

      return
c end function dist
      end
