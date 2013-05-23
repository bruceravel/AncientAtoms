      subroutine bperm(iat, nat, iperm, cell, ipt, st)
      
c  permute monoclinic and orthorhombic back to the non-standard setting

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      dimension st(iat,192,3), ipt(iat), cell(6)
      dimension abc(3), perm(3,3)

c       if (iperm.eq.22) then
c           print*,iperm
c           return 
c       endif

      do 20 j=1,3
        do 10 i=1,3
          perm(i,j) = 0
 10     continue 
        abc(j) = 0.e0
 20   continue 

c  set elements fo permutation matrix
c             orthorhombic, cab
      if (iperm.eq.2) then
          perm(2,1) = 1.e0
          perm(3,2) = 1.e0
          perm(1,3) = 1.e0
c             orthorhombic, bca
      elseif (iperm.eq.3) then
          perm(3,1) = 1.e0
          perm(1,2) = 1.e0
          perm(2,3) = 1.e0
c             orthorhombic, a-cb
      elseif (iperm.eq.4) then
          perm(1,1) = 1.e0
          perm(3,2) = 1.e0
          perm(2,3)= -1.e0
c             orthorhombic, ba-c
      elseif (iperm.eq.5) then
          perm(2,1) = 1.e0
          perm(1,2) = 1.e0
          perm(3,3) = -1.e0
c             orthorhombic, -cba
      elseif (iperm.eq.6) then
          perm(3,1) = 1.e0
          perm(2,2) = 1.e0
          perm(1,3) = -1.e0
c             monoclinic, acb(2nd) to abc(1st)
      elseif (iperm.eq.12) then
          perm(1,1) = 1.e0
          perm(2,3) = 1.e0
          perm(3,2) = 1.e0
      endif

c     --- orthorhombic or monoclinic
      if (iperm.lt.20) then
          do 40 i=1,3
            do 30 j=1,3
              abc(i) = abc(i) + cell(j) * abs(perm(i,j))
 30         continue 
 40       continue 
          do 45 i=1,3
            cell(i) = abc(i)
 45       continue 
      endif

c     --- tetragonal
c       else
c           abc(1) = cell(1) * sqrt(2.e0) 
c           abc(2) = abc(1)
c           abc(3) = cell(3)

c       print*,'iperm=',iperm
c       print*,'before, after'
c       print*,cell(1),abc(1)
c       print*,cell(2),abc(2)
c       print*,cell(3),abc(3)
c       print*,'angles:',cell(4),cell(5),cell(6)


      do 70 i=1,nat
        do 60 j=1,ipt(i)
          xx = 0
          yy = 0
          zz = 0
c  permute fractional coordinates
          xx = st(i,j,1)*perm(1,1) + st(i,j,2)*perm(1,2) +
     $                st(i,j,3)*perm(1,3)
          yy = st(i,j,1)*perm(2,1) + st(i,j,2)*perm(2,2) +
     $                st(i,j,3)*perm(2,3)
          zz = st(i,j,1)*perm(3,1) + st(i,j,2)*perm(3,2) +
     $                st(i,j,3)*perm(3,3)
c           print*,i,st(i,j,1),xx
c           print*,i,st(i,j,2),yy
c           print*,i,st(i,j,3),zz
          st(i,j,1) = xx
          st(i,j,2) = yy
          st(i,j,3) = zz
 60     continue 
 70   continue 

      return 
c  end subroutine bperm
      end
