      subroutine systm(cell,isystm)
c--------------------------------------------------------------
c  copyright (c) 1998 Bruce Ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c     determine system from cell constants for comparison
c     to determination for space group.
c     isystm : 1 = monoclinic     2 = orthorhombic
c              3 = <not used>     4 = tetragonal
c              5 = cubic          6 = hexagonal
c              7 = triclinic      8 = <not used>
c  who came up with this notation convention? why is cubic 5 and not 3?
c----------------------------------------------------------------------
      parameter (eps=.0001)
      dimension cell(6)

      if (abs(cell(6)-120.0).lt.eps) then
         isystm=6
      elseif ( (abs(cell(4)-90).gt.eps) .and.
     $               (abs(cell(5)-90).lt.eps) .and.
     $               (abs(cell(6)-90).lt.eps) .and.
     $               (abs(cell(3)-cell(1)).lt.eps) .and.
     $               (abs(cell(2)-cell(1)).lt.eps) ) then
         isystm=6
      elseif ( ( (abs(cell(5)-90).gt.eps) .and.
     $                   (abs(cell(6)-90).gt.eps) )
     $      .or. (abs(cell(4)-90).gt.eps)
     $      .or. (abs(cell(6)-90).gt.eps)) then
         isystm=7
      elseif (abs(cell(5)-90).gt.eps) then
         isystm=1
      elseif (abs(cell(3)-cell(1)).lt.eps) then
         isystm=5
      elseif (abs(cell(1)-cell(2)).lt.eps) then
         isystm=4
      else
         isystm=2
      endif

      return
c end subroutine systm
      end
