      subroutine trans(vector,trmtx)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------
c  transform real lengths along cell axes into cartesian coordinates.
c  trmtx(3,3) is the metric tensor calculated in the subroutine 
c  `metric'.  vector(3) is input (cell axis basis) and output 
c  (cartesian)
c-------------------------------------------------------------------
      dimension trmtx(3,3),vector(3),toss(3)

      do 10 i=1,3
         toss(i)   = vector(i)
         vector(i) = 0
10    continue 

      do 30 i=1,3
         do 20 j=1,3
            vector(i) = vector(i) + trmtx(i,j)*toss(j)
20       continue
30    continue

      return
c end subroutine trans
      end
