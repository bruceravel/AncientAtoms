      subroutine card(word,feffcd,ii)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c------------------------------------------------------------------
c  feff requires the cards to be upper case.  this routine assures 
c  that the card will be regardless of the case of this source code
c
c  input:  
c    word:   feff card for writing to feff.inp
c  output:
c    feffcd: word in all upper case characters
c    ii:     length of word
c------------------------------------------------------------------
      character*(*) word, feffcd

      feffcd = word
      call upper(feffcd)
      ii = istrln(feffcd)

      return
c  end subroutine card
      end

