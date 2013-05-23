      subroutine spcchk (spcgrp,ispace)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------
c  make sure spcgrp is a valid hermann-mauguin symbol
c  ispace is the index of the space group as listed below
c  also check against the schoenflies notation and substitute the
c  corresponding h-m symbol.
c  these are listed such that space(i) = schoen(i) for i=1,230
c  isn't this cool?
c----------------------------------------------------------------
      character*10 spcgrp,space(230),schoen(230),test

c====================================================================
c  hermann-maguin notation
c  triclinic and monoclinic
      data (space(i),i=1,15)/
     $'p 1','p -1','p 2','p 21','c 2','p m','p c','c m','c c','p 2/m',
     $'p 21/m','c 2/m','p 2/c','p 21/c','c 2/c'/
c  orthorhombic
       data(space(i),i=16,74)/
     $'p 2 2 2','p 2 2 21','p 21 21 2','p 21 21 21','c 2 2 21',
     $'c 2 2 2','f 2 2 2','i 2 2 2','i 21 21 21','p m m 2','p m c 21',
     $'p c c 2','p m a 2','p c a 21','p n c 2','p m n 21','p b a 2',
     $'p n a 21','p n n 2','c m m 2','c m c 21','c c c 2','a m m 2',
     $'a b m 2','a m a 2','a b a 2','f m m 2','f d d 2','i m m 2',
     $'i b a 2','i m a 2','p m m m','p n n n','p c c m','p b a n',
     $'p m m a','p n n a','p m n a','p c c a','p b a m','p c c n',
     $'p b c m','p n n m','p m m n','p b c n','p b c a','p n m a',
     $'c m c m','c m c a','c m m m','c c c m','c m m a','c c c a',
     $'f m m m','f d d d','i m m m','i b a m','i b c a','i m m a'/
c  tetragonal
       data(space(i),i=75,142)/
     $'p 4','p 41','p 42','p 43','i 4','i 41','p -4','i -4','p 4/m',
     $'p 42/m','p 4/n','p 42/n','i 4/m','i 41/a','p 4 2 2','p 4 21 2',
     $'p 41 2 2','p 41 21 2','p 42 2 2','p 42 21 2','p 43 2 2',
     $'p 43 21 2','i 4 2 2','i 41 2 2','p 4 m m','p 4 b m','p 42 c m',
     $'p 42 n m','p 4 c c','p 4 n c','p 42 m c','p 42 b c','i 4 m m',
     $'i 4 c m','i 41 m d','i 41 c d','p -4 2 m','p -4 2 c','p -4 21 m',
     $'p -4 21 c','p -4 m 2','p -4 c 2','p -4 b 2','p -4 n 2',
     $'i -4 m 2','i -4 c 2','i -4 2 m','i -4 2 d','p 4/m m m',
     $'p 4/m c c','p 4/n b m','p 4/n n c','p 4/m b m','p 4/m n c',
     $'p 4/n m m','p 4/n c c','p 42/m m c','p 42/m c m','p 42/n b c',
     $'p 42/n n m','p 42/m b c','p 42/m n m','p 42/n m c','p 42/n c m',
     $'i 4/m m m','i 4/m c m','i 41/a m d','i 41/a c d'/
c  trigonal
       data(space(i),i=143,167)/
     $'p 3','p 31','p 32','r 3','p -3','r -3','p 3 1 2','p 3 2 1',
     $'p 31 1 2','p 31 2 1','p 32 1 2','p 32 2 1','r 3 2','p 3 m 1',
     $'p 3 1 m','p 3 c 1','p 3 1 c','r 3 m','r 3 c','p -3 1 m',
     $'p -3 1 c','p -3 m 1','p -3 c 1','r -3 m','r -3 c'/
c  hexagonal
       data(space(i),i=168,194)/
     $'p 6','p 61','p 65','p 62','p 64','p 63','p -6','p 6/m','p 63/m',
     $'p 6 2 2','p 61 2 2','p 65 2 2','p 62 2 2','p 64 2 2','p 63 2 2',
     $'p 6 m m','p 6 c c','p 63 c m','p 63 m c','p -6 m 2','p -6 c 2',
     $'p -6 2 m','p -6 2 c','p 6/m m m','p 6/m c c','p 63/m c m',
     $'p 63/m m c'/
c  cubic
       data(space(i),i=195,230)/
     $'p 2 3','f 2 3','i 2 3','p 21 3','i 21 3','p m 3','p n 3','f m 3',
     $'f d 3','i m 3','p a 3','i a 3','p 4 3 2','p 42 3 2','f 4 3 2',
     $'f 41 3 2','i 4 3 2','p 43 3 2','p 41 3 2','i 41 3 2','p -4 3 m',
     $'f -4 3 m','i -4 3 m','p -4 3 n','f -4 3 c','i -4 3 d','p m 3 m',
     $'p n 3 n','p m 3 n','p n 3 m','f m 3 m','f m 3 c','f d 3 m',
     $'f d 3 c','i m 3 m','i a 3 d'/

c====================================================================
c  schoenflies notation
c  triclinic and monoclinic
      data (schoen(i),i=1,15)/
     $'c_1^1','c_i^1','c_2^1','c_2^2','c_2^3','c_s^1','c_s^2','c_s^3',
     $'c_s^4','c_2h^1','c_2h^2','c_2h^3','c_2h^4','c_2h^5','c_2h^6'/
c  orthorhombic
       data(schoen(i),i=16,74)/
     $'d_2^1','d_2^2','d_2^3','d_2^4','d_2^5',
     $'d_2^6','d_2^7','d_2^8','d_2^9','c_2v^1','c_2v^2',
     $'c_2v^3','c_2v^4','c_2v^5','c_2v^6','c_2v^7','c_2v^8',
     $'c_2v^9','c_2v^10','c_2v^11','c_2v^12','c_2v^13','c_2v^14',
     $'c_2v^15','c_2v^16','c_2v^17','c_2v^18','c_2v^19','c_2v^20',
     $'c_2v^21','c_2v^22','d_2h^1','d_2h^2','d_2h^3','d_2h^4',
     $'d_2h^5','d_2h^6','d_2h^7','d_2h^8','d_2h^9','d_2h^10',
     $'d_2h^11','d_2h^12','d_2h^13','d_2h^14','d_2h^15','d_2h^16',
     $'d_2h^17','d_2h^18','d_2h^19','d_2h^20','d_2h^21','d_2h^22',
     $'d_2h^23','d_2h^24','d_2h^25','d_2h^26','d_2h^27','d_2h^28'/
c  tetragonal
       data(schoen(i),i=75,142)/
     $'c_4^1','c_4^2','c_4^3','c_4^4','c_4^5','c_4^6','s_4^1','s_4^2',
     $'c_4h^1','c_4h^2','c_4h^3','c_4h^4','c_4h^5','c_4h^6',
     $'d_4^1','d_4^2','d_4^3','d_4^4','d_4^5','d_4^6','d_4^7',
     $'d_4^8','d_4^9','d_4^10','c_4v^1','c_4v^2','c_4v^3',
     $'c_4v^4','c_4v^5','c_4v^6','c_4v^7','c_4v^8','c_4v^9',
     $'c_4v^10','c_4v^11','c_4v^12','d_2d^1','d_2d^2','d_2d^3',
     $'d_2d^4','d_2d^5','d_2d^6','d_2d^7','d_2d^8',
     $'d_2d^9','d_2d^10','d_2d^11','d_2d^12','d_4h^1',
     $'d_4h^2','d_4h^3','d_4h^4','d_4h^5','d_4h^6',
     $'d_4h^7','d_4h^8','d_4h^9','d_4h^10','d_4h^11',
     $'d_4h^12','d_4h^13','d_4h^14','d_4h^15','d_4h^16',
     $'d_4h^17','d_4h^18','d_4h^19','d_4h^20'/
c  trigonal
       data(schoen(i),i=143,167)/
     $'c_3^1','c_3^2','c_3^3','c_3^4','c_3i^1','c_3i^2',
     $'d_3^1','d_3^2','d_3^3','d_3^4','d_3^5','d_3^6','d_3^7',
     $'c_3v^1','c_3v^2','c_3v^3','c_3v^4','c_3v^5','c_3v^6',
     $'d_3d^1','d_3d^2','d_3d^3','d_3d^4','d_3d^5','d_3d^6'/
c  hexagonal
       data(schoen(i),i=168,194)/
     $'c_6^1','c_6^2','c_6^3','c_6^4','c_6^5','c_6^6',
     $'c_3h^1','c_6h^1','c_6h^2',
     $'d_6^1','d_6^2','d_6^3','d_6^4','d_6^5','d_6^6',
     $'c_6v^1','c_6v^2','c_6v^3','c_6v^4','d_3h^1','d_3h^2',
     $'d_3h^3','d_3h^4','d_6h^1','d_6h^2','d_6h^3','d_6h^4'/
c  cubic
       data(schoen(i),i=195,230)/
     $'t^1','t^2','t^3','t^4','t^5','t_h^1','t_h^2','t_h^3',
     $'t_h^4','t_h^5','t_h^6','t_h^7','o^1','o^2','o^3',
     $'o^4','o^5','o^6','o^7','o^8','t_d^1',
     $'t_d^2','t_d^3','t_d^4','t_d^5','t_d^6','o_h^1',
     $'o_h^2','o_h^3','o_h^4','o_h^5','o_h^6','o_h^7',
     $'o_h^8','o_h^9','o_h^10'/
c%%%
c====================================================================
c  the second character of the schoenflies notation must be _ or ^,
c  these symbols are not used in the hermann-maguin notation.
c  change hm hexagonal 'c' groups to 'p'
c  the value of test *must* be in the same case as the data above!!
      test = 'ab'
      call case(test,spcgrp)
      ispace=0

      if ((spcgrp(2:2).ne.'_').and.(spcgrp(2:2).ne.'^')) then
c                                              hermann-maguin
          do 10 i=1,230
            if ((spcgrp(3:3).eq.'6').and.(spcgrp(1:1).eq.'c'))
     $                  spcgrp(1:1) = 'p'
            if ((spcgrp(3:3).eq.'-').and.(spcgrp(4:4).eq.'6').and.
     $                  (spcgrp(1:1).eq.'c'))   spcgrp(1:1) = 'p'
            if (spcgrp.eq.space(i)) then
                ispace=i
                goto 30
            endif
 10       continue
      else
c                                              schoenflies
          call schfix(spcgrp)
          do 20 i=1,230
            if (spcgrp.eq.schoen(i)) then
                ispace = i
                spcgrp = space(i)
                goto 30
            endif
 20       continue
      endif

 30   continue

      return
c  end subroutine spcchk
      end
