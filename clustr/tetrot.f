      subroutine tetrot(natx, itot, cell, atlis)
      
c  permute tetragonal cystals back to the non-standard setting

c    natx:  dimension set in celling routine
c    itot:  number of atoms in cluster
c    cell:  (6) array of a,b,c,alpha,beta,gamma
c    atlis: (natx, 8) array of all atoms in cluster
c           1-3 -> pos. in cell     4 -> dist to origin
c              5 -> atom type     6-8 -> cartesian coords

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      parameter (root2=1.41421356237309504880e0)
      dimension atlis(natx,8), cell(6), perm(2,2)

      perm(1,1) = 1.e0
      perm(1,2) = -1.e0
      perm(2,1) = 1.e0
      perm(2,2) = 1.e0

      ab = cell(1)
      cell(1) = ab * sqrt(2.e0)
      cell(2) = cell(1)

      do 70 i=1,itot

c  permute cartesian coordinates
        xx = atlis(i,6)*perm(1,1) + atlis(i,7)*perm(1,2) 
        yy = atlis(i,6)*perm(2,1) + atlis(i,7)*perm(2,2) 
        atlis(i,6) = xx / root2
        atlis(i,7) = yy / root2

c  permute fractional coordinates
        xx = atlis(i,1)*perm(1,1) + atlis(i,2)*perm(1,2) 
        yy = atlis(i,1)*perm(2,1) + atlis(i,2)*perm(2,2) 
        atlis(i,1) = xx 
        atlis(i,2) = yy 

 70   continue 

      return 
c  end subroutine tetrot
      end

