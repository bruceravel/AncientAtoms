      integer function nxtunt(iunit)

c  this function returns the value of the next unopened unit number equal
c  to or larger than iunit.  it will return neither unit numbers 0, 5, 
c  or 6 nor a negative unit number 
      
      integer iunit, iun
      logical open

      iun = iunit
      if (iun.le.0) iun = 1
 10   continue 
      if ((iun.eq.5).or.(iun.eq.6)) then
          iun = 7
          goto 10
      endif
      inquire (unit=iun, opened=open)
      if (open) then
          iun = iun + 1
          goto 10
      endif
      
      nxtunt = iun
      return 
c  end integer function nxtunt
      end

c=======================================================================
