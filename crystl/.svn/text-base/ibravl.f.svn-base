      integer function ibrav(csymbr)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit integer(i-n)
c-------------------------------------------------------------------
c  this translates between the integer convention for bravais lattice
c  type used in equipt and that used elsewhere in the program
c
c  csymbr is the bravais lattice type, ibrav is a number 
c  corresponding to it
c-------------------------------------------------------------------
      character bra(7),csymbr
      data (bra(i),i=1,7)/'p','i','r','f','a','b','c'/

      do 10 i=1,7
        if (csymbr.eq.bra(i)) then
            ibrav = i
            goto 99
        endif
 10   continue
      ibrav = 1
 99   return
c end integer function ibravl
      end

