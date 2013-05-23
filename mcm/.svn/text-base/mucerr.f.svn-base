      subroutine mucerr(ier)

c  generate a warning message from the mucal error code

      character*2  ending
      character*22 begin
      character*78 messg

      if (ier.eq.0) return

      begin = ' *** Warning: (mucal) '
      ending = '.-'

      if (ier.eq.1) then
          messg = begin//'Energy input to mucal is zero'//ending
      elseif (ier.eq.2) then
          messg = begin//'Element name does not match Z'//ending
      elseif ((ier.eq.3).or.(ier.eq.4)) then
          messg = begin//'No data for Po, At, Fr, Ra, Ac, Pa, '//
     $                'Np, or above Pu'//ending
      elseif (ier.eq.5) then
          messg = begin//'McMaster only provides l1 data for Z<30'//
     $                ending
      elseif (ier.eq.6) then
          messg = begin//'Energy input to mucal is on an edge'//ending
      elseif (ier.eq.7) then
          messg = begin//'No element name of Z supplied'//ending
      endif

      call messag(messg)
      return
c  end subroutine mucerr
      end
