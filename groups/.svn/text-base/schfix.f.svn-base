      subroutine schfix(spcgrp)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision(a-h,o-z)
c  make sure shoenflies designation has the subscript before the 
c  superscript.  this way the user can enter d_3^3 or d^3_3
c  and d^3_3 will be converted to d_3^3.
c  'v' groups will be changed to 'd' groups
      character*10 spcgrp, second, third, first*1, test*2

      test   = 'ab'
      call case(test,spcgrp)
      iunder = index(spcgrp,'_')
      iover  = index(spcgrp,'^')
      first  = spcgrp(1:1)

c  change V groups to more modern D notation
      if (first.eq.'v') then
          spcgrp(1:1) = 'd'
          first = 'd'
      endif

c  certain T or O groups
      if (istrln(spcgrp).eq.3) then
          if ( (first.eq.'t').and.(iover.eq.2) ) then
              return
          elseif ( (first.eq.'o').and.(iover.eq.2) ) then
              return
          endif
      endif

c  underscore preceeds carot, as it should
      if ((iunder.lt.iover).and.(iunder.ne.0)) return

c  else switch underscore part and carot part
      if (iunder.eq.0) iunder=9
      second = spcgrp(2:iunder-1)
      third  = spcgrp(iunder:)
      if (first.eq.'d') then
c         --- groups 16-24
          if (third.eq.' ') then
              third = '_2'
c         --- groups 47-74
          elseif (third.eq.'_h') then
              third = '_2h'
          endif
      endif

      i2 = istrln(second)
      i3 = istrln(third)
      spcgrp = first//third(1:i3)//second(1:i2)

      return
c end subroutine schfix
      end
