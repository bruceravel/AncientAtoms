      subroutine fperm(iat, nat, iperm, cell, x, y, z)
      
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      dimension cell(6), x(iat), y(iat), z(iat)
      dimension abc(3), perm(3,3)

      do 20 j=1,3
        do 10 i=1,3
          perm(i,j) = 0
 10     continue 
        abc(j) = 0.e0
 20   continue 

c  set elements fo permutation matrix
c             orthorhombic, cab
      if (iperm.eq.2) then
          perm(1,2) = 1.e0
          perm(2,3) = 1.e0
          perm(3,1) = 1.e0
c             orthorhombic, bca
      elseif (iperm.eq.3) then
          perm(1,3) = 1.e0
          perm(2,1) = 1.e0
          perm(3,2) = 1.e0
c             orthorhombic, a-cb
      elseif (iperm.eq.4) then
          perm(1,1) = 1.e0
          perm(2,3) = 1.e0
          perm(3,2)= -1.e0
c             orthorhombic, ba-c
      elseif (iperm.eq.5) then
          perm(1,2) = 1.e0
          perm(2,1) = 1.e0
          perm(3,3) = -1.e0
c             orthorhombic, -cba
      elseif (iperm.eq.6) then
          perm(1,3) = 1.e0
          perm(2,2) = 1.e0
          perm(3,1) = -1.e0
c             monoclinic, acb(2nd) to abc(1st)
      elseif (iperm.eq.11) then
          perm(1,1) = 1.e0
          perm(2,3) = 1.e0
          perm(3,2) = 1.e0
c             tetragonal, rotated --> standard
      elseif (iperm.eq.22) then
          perm(1,1) = 1.e0
          perm(1,2) = 1.e0
          perm(2,1) = -1.e0
          perm(2,2) = 1.e0
          perm(3,3) = 1.e0
      endif

c     --- orthorhombic or monoclinic
      if (iperm.lt.20) then
          do 40 i=1,3
            do 30 j=1,3
              abc(i) = abc(i) + cell(j) * abs(perm(i,j))
 30         continue 
 40       continue 
c     --- tetragonal
      else
          abc(1) = cell(1) / sqrt(2.e0) 
          abc(2) = abc(1)
          abc(3) = cell(3)
      endif
c       print*,'iperm=',iperm
c       print*,'before, after'
c       print*,cell(1),abc(1)
c       print*,cell(2),abc(2)
c       print*,cell(3),abc(3)
c       print*,'angles:',cell(4),cell(5),cell(6)

      do 45 i=1,3
        cell(i) = abc(i)
 45   continue 

      do 50 i=1,nat
        xx = x(i)*perm(1,1) + y(i)*perm(1,2) + z(i)*perm(1,3)
        yy = x(i)*perm(2,1) + y(i)*perm(2,2) + z(i)*perm(2,3)
        zz = x(i)*perm(3,1) + y(i)*perm(3,2) + z(i)*perm(3,3)
c         print*,i,x(i),xx
c         print*,i,y(i),yy
c         print*,i,z(i),zz
        x(i) = xx
        y(i) = yy
        z(i) = zz
 50   continue 

      return 
c  end subroutine fperm
      end
