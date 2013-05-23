      subroutine settng(spcgrp, iperm, ispa)

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)
c---------------------------------------------------------------------------
c recognize symbols for alternate settings of monoclinic, orthorhombic,
c and tetragonal space groups
c
c requires subroutines case, messag
c
c input:
c   spcgrp:  symbol for group as recognized in spcchk
c output:
c   iperm:   index of permutation matrix
c            1-6:   6 orthorhombic settings (abc, cab, bca, a-cb, ba-c, -cba)
c            11-12: 2 monoclinic settings, z-axis unique, y-axis unique
c            21-22: 2 tetragonal settings, standard and rotated
c   ispa:    index of space group
c---------------------------------------------------------------------------
      character*10 spcgrp, oldgrp, altern(6,3:142), test*2

c  herman maguinn symbols for alternate settings of monoclinic groups
c  atoms uses y-axis unique (2nd setting)
      data (altern(i,3),  i=1,2) / 'p 2', 'p 2'/
      data (altern(i,4),  i=1,2) / 'p 21', 'p 21'/
      data (altern(i,5),  i=1,2) / 'b 2', 'c 2'/
      data (altern(i,6),  i=1,2) / 'p m', 'p m'/
      data (altern(i,7),  i=1,2) / 'p b', 'p c'/
      data (altern(i,8),  i=1,2) / 'b m', 'c m'/
      data (altern(i,9),  i=1,2) / 'b b', 'c c'/
      data (altern(i,10), i=1,2) / 'p 2/m', 'p 2/m'/
      data (altern(i,11), i=1,2) / 'p 21/m', 'p 21/m'/
      data (altern(i,12), i=1,2) / 'b 2/m', 'c 2/m'/
      data (altern(i,13), i=1,2) / 'p 2/b', 'p 2/c'/
      data (altern(i,14), i=1,2) / 'p 21/b', 'p 21/c'/
      data (altern(i,15), i=1,2) / 'b 2/b', 'c 2/c'/

c  herman maguinn symbols for alternate settings of orthorhombic groups
      data (altern(i,16), i=1,6) / 'p 2 2 2', 'p 2 2 2', 'p 2 2 2',
     $            'p 2 2 2', 'p 2 2 2', 'p 2 2 2'/
      data (altern(i,17), i=1,6) / 'p 2 2 21', 'p 21 2 2', 'p 2 21 2',
     $            'p 2 21 2', 'p 2 2 21', 'p 21 2 2'/
      data (altern(i,18), i=1,6) / 'p 21 21 2', 'p 2 21 21',
     $            'p 21 2 21', 'p 21 2 21', 'p 21 21 2', 'p 2 21 21'/
      data (altern(i,19), i=1,6) / 'p 21 21 21', 'p 21 21 21',
     $        'p 21 21 21', 'p 21 21 21', 'p 21 21 21', 'p 21 21 21'/
      data (altern(i,20), i=1,6) / 'c 2 2 21', 'a 21 2 2', 'b 2 21 2',
     $            'b 2 21 2', 'c 2 2 21', 'a 21 2 2'/
      data (altern(i,21), i=1,6) / 'c 2 2 2', 'a 2 2 2', 'b 2 2 2',
     $            'b 2 2 2', 'c 2 2 2', 'a 2 2 2'/
      data (altern(i,22), i=1,6) / 'f 2 2 2', 'f 2 2 2', 'f 2 2 2',
     $            'f 2 2 2', 'f 2 2 2', 'f 2 2 2'/
      data (altern(i,23), i=1,6) / 'i 2 2 2', 'i 2 2 2', 'i 2 2 2',
     $            'i 2 2 2', 'i 2 2 2', 'i 2 2 2'/
      data (altern(i,24), i=1,6) / 'i 21 21 21', 'i 21 21 21',
     $        'i 21 21 21', 'i 21 21 21', 'i 21 21 21', 'i 21 21 21'/
      data (altern(i,25), i=1,6) / 'p m m 2', 'p 2 m m', 'p m 2 m',
     $            'p m 2 m', 'p m m 2', 'p 2 m m'/
      data (altern(i,26), i=1,6) / 'p m c 21', 'p 21 m a', 'p b 21 m',
     $            'p m 21 b', 'p c m 21', 'p 21 a m'/
      data (altern(i,27), i=1,6) / 'p c c 2', 'p 2 a a', 'p b 2 b',
     $            'p b 2 b', 'p c c 2', 'p 2 a a'/
      data (altern(i,28), i=1,6) / 'p m a 2', 'p 2 m b', 'p c 2 m',
     $            'p m 2 a', 'p b m 2', 'p 2 c m'/
      data (altern(i,29), i=1,6) / 'p c a 21', 'p 21 a b', 'p c 21 b',
     $            'p b 21 a', 'p b c 21', 'p 21 c a'/
      data (altern(i,30), i=1,6) / 'p n c 2', 'p 2 n a', 'p b 2 n',
     $            'p n 2 b', 'p c n 2', 'p 2 a n'/
      data (altern(i,31), i=1,6) / 'p m n 21', 'p 21 m n', 'p n 21 m',
     $            'p m 21 n', 'p n m 21', 'p 2 n m'/
      data (altern(i,32), i=1,6) / 'p b a 2', 'p 2 c b', 'p c 2 a',
     $            'p c 2 a', 'p b a 2', 'p 2 c b'/
      data (altern(i,33), i=1,6) / 'p n a 21', 'p 21 n b', 'p c 21 n',
     $            'p n 21 a', 'p b n 21', 'p 2 c n'/
      data (altern(i,34), i=1,6) / 'p n n 2', 'p 2 n n', 'p n 2 n',
     $            'p n 2 n', 'p n n 2', 'p 2 n n'/
      data (altern(i,35), i=1,6) / 'c m m 2', 'a 2 m m', 'b m 2 m',
     $            'b m 2 m', 'c m m 2', 'a 2 m m'/
      data (altern(i,36), i=1,6) / 'c m c 21', 'a 21 m a', 'b b 21 m',
     $            'b m 21 b', 'c c m 21', 'a 21 a m'/
      data (altern(i,37), i=1,6) / 'c c c 2', 'a 2 c a', 'b b 2 c',
     $            'b b 2 b', 'c c c 2', 'a 2 a a'/
      data (altern(i,38), i=1,6) / 'a m m 2', 'b 2 m m', 'c m 2 m',
     $            'a m 2 m', 'b m m 2', 'c 2 m m'/
      data (altern(i,39), i=1,6) / 'a b m 2', 'b 2 c m', 'c m 2 a',
     $            'a c 2 m', 'b m a 2', 'c 2 m b'/
      data (altern(i,40), i=1,6) / 'a m a 2', 'b 2 m b', 'c c 2 m',
     $            'a m 2 a', 'b b m 2', 'c 2 c m'/
      data (altern(i,41), i=1,6) / 'a b a 2', 'b 2 c b', 'c c 2 a',
     $            'a c 2 a', 'b b a 2', 'c 2 c b'/
      data (altern(i,42), i=1,6) / 'f m m 2', 'f 2 m m', 'f m 2 m',
     $            'f m 2 m', 'f m m 2', 'f 2 m m'/
      data (altern(i,43), i=1,6) / 'f d d 2', 'f 2 d d', 'f d 2 d',
     $            'f d 2 d', 'f d d 2', 'f 2 d d'/
      data (altern(i,44), i=1,6) / 'i m m 2', 'i 2 m m', 'i m 2 m',
     $            'i m 2 m', 'i m m 2', 'i 2 m m'/
      data (altern(i,45), i=1,6) / 'i b a 2', 'i 2 c b', 'i c 2 a',
     $            'i c 2 a', 'i b a 2', 'i 2 c b'/
      data (altern(i,46), i=1,6) / 'i m a 2', 'i 2 m b', 'i c 2 m',
     $            'i m 2 a', 'i b m 2', 'i 2 c m'/
      data (altern(i,47), i=1,6) / 'p m m m', 'p m m m', 'p m m m',
     $            'p m m m', 'p m m m', 'p m m m'/
      data (altern(i,48), i=1,6) / 'p n n n', 'p n n n', 'p n n n',
     $            'p n n n', 'p n n n', 'p n n n'/
      data (altern(i,49), i=1,6) / 'p c c m', 'p m a a', 'p b m b',
     $            'p b m b', 'p c c m', 'p m a a'/
      data (altern(i,50), i=1,6) / 'p b a n', 'p n c b', 'p c n a',
     $            'p c n a', 'p b a n', 'p n c b'/
      data (altern(i,51), i=1,6) / 'p m m a', 'p b m m', 'p m c m',
     $            'p m a m', 'p m m b', 'p c m m'/
      data (altern(i,52), i=1,6) / 'p n n a', 'p b n n', 'p n c n',
     $            'p n a n', 'p n n b', 'p c n n'/
      data (altern(i,53), i=1,6) / 'p m n a', 'p b m n', 'p n c m',
     $            'p m a n', 'p n m b', 'p c n m'/
      data (altern(i,54), i=1,6) / 'p c c a', 'p b a a', 'p b c b',
     $            'p b a b', 'p c c b', 'p c a a'/
      data (altern(i,55), i=1,6) / 'p b a m', 'p m c b', 'p c m a',
     $            'p c m a', 'p b a m', 'p m c b'/
      data (altern(i,56), i=1,6) / 'p c c n', 'p n a a', 'p b n b',
     $            'p b n b', 'p c c n', 'p n a a'/
      data (altern(i,57), i=1,6) / 'p b c m', 'p m c a', 'p b m a',
     $            'p c m b', 'p c a m', 'p m a b'/
      data (altern(i,58), i=1,6) / 'p n n m', 'p m n n', 'p n m n',
     $            'p n m n', 'p n n m', 'p m n n'/
      data (altern(i,59), i=1,6) / 'p m m n', 'p n m m', 'p m n m',
     $            'p m n m', 'p m m n', 'p n m m'/
      data (altern(i,60), i=1,6) / 'p b c n', 'p n c a', 'p b n a',
     $            'p c n b', 'p c a n', 'p n a b'/
      data (altern(i,61), i=1,6) / 'p b c a', 'p b c a', 'p b c a',
     $            'p c a b', 'p c a b', 'p c a b'/
      data (altern(i,62), i=1,6) / 'p n m a', 'p b n m', 'p m c n',
     $            'p n a m', 'p m n b', 'p c m n'/
      data (altern(i,63), i=1,6) / 'c m c m', 'a m m a', 'b b m m',
     $            'b m m b', 'c c m m', 'a m a m'/
      data (altern(i,64), i=1,6) / 'c m c a', 'a b m a', 'b b c m',
     $            'b m a b', 'c c m b', 'a c a m'/
      data (altern(i,65), i=1,6) / 'c m m m', 'a m m m', 'b m m m',
     $            'b m m m', 'c m m m', 'a m m m'/
      data (altern(i,66), i=1,6) / 'c c c m', 'a m a a', 'b b m b',
     $            'b b m b', 'c c c m', 'a m a a'/
      data (altern(i,67), i=1,6) / 'c m m a', 'a b m m', 'b m c m',
     $            'b m a m', 'c m m b', 'a c m m'/
      data (altern(i,68), i=1,6) / 'c c c a', 'a b a a', 'b b c b',
     $            'b b a b', 'c c c b', 'a c a a'/
      data (altern(i,69), i=1,6) / 'f m m m', 'f m m m', 'f m m m',
     $            'f m m m', 'f m m m', 'f m m m'/
      data (altern(i,70), i=1,6) / 'f d d d', 'f d d d', 'f d d d',
     $            'f d d d', 'f d d d', 'f d d d'/
      data (altern(i,71), i=1,6) / 'i m m m', 'i m m m', 'i m m m',
     $            'i m m m', 'i m m m', 'i m m m'/
      data (altern(i,72), i=1,6) / 'i b a m', 'i m c b', 'i c m a',
     $            'i c m a', 'i b a m', 'i m c b'/
      data (altern(i,73), i=1,6) / 'i b c a', 'i b c a', 'i b c a',
     $            'i c a b', 'i c a b', 'i c a b'/
      data (altern(i,74), i=1,6) / 'i m m a', 'i b m m', 'i m c m',
     $            'i m a m', 'i m m b', 'i c m m'/

c  herman maguinn symbols for alternate settings of tetragonal groups
      data (altern(i,75),  i=1,2) / 'p 4', 'c 4'/
      data (altern(i,76),  i=1,2) / 'p 41', 'c 41'/
      data (altern(i,77),  i=1,2) / 'p 42', 'c 42'/
      data (altern(i,78),  i=1,2) / 'p 43', 'c 43'/
      data (altern(i,79),  i=1,2) / 'i 4', 'f 4'/
      data (altern(i,80),  i=1,2) / 'i 41', 'f 41'/
      data (altern(i,81),  i=1,2) / 'p -4', 'c -4'/
      data (altern(i,82),  i=1,2) / 'i -4', 'f -4'/
      data (altern(i,83),  i=1,2) / 'p 4/m', 'c 4/m'/
      data (altern(i,84),  i=1,2) / 'p 42/m', 'c 42/m'/
      data (altern(i,85),  i=1,2) / 'p 4/n', 'c 4/a'/
      data (altern(i,86),  i=1,2) / 'p 42/m', 'c 42/a'/
      data (altern(i,87),  i=1,2) / 'i 4/m', 'f 4/m'/
      data (altern(i,88),  i=1,2) / 'i 41/a', 'f 41/d'/
      data (altern(i,89),  i=1,2) / 'p 4 2 2', 'c 4 2 2'/
      data (altern(i,90),  i=1,2) / 'p 4 2 21', 'c 4 2 21'/
      data (altern(i,91),  i=1,2) / 'p 41 2 2', 'c 41 2 2'/
      data (altern(i,92),  i=1,2) / 'p 41 2 21', 'c 41 2 21'/
      data (altern(i,93),  i=1,2) / 'p 42 2 2', 'c 42 2 2'/
      data (altern(i,94),  i=1,2) / 'p 42 2 21', 'c 42 2 21'/
      data (altern(i,95),  i=1,2) / 'p 43 2 2', 'c 43 2 2'/
      data (altern(i,96),  i=1,2) / 'p 43 2 21', 'c 43 2 21'/
      data (altern(i,97),  i=1,2) / 'i 4 2 2', 'f 4 2 2'/
      data (altern(i,98),  i=1,2) / 'i 41 2 2', 'f 41 2 2'/
      data (altern(i,99),  i=1,2) / 'p 4 m m', 'c 4 m m'/
      data (altern(i,100), i=1,2) / 'p 4 b m', 'c 4 m b'/
      data (altern(i,101), i=1,2) / 'p 42 c m', 'c 42 m c'/
      data (altern(i,102), i=1,2) / 'p 42 n m', 'c 42 m n'/
      data (altern(i,103), i=1,2) / 'p 4 c c', 'c 4 c c'/
      data (altern(i,104), i=1,2) / 'p 4 n c', 'c 4 c n'/
      data (altern(i,105), i=1,2) / 'p 42 m c', 'c 42 c m'/
      data (altern(i,106), i=1,2) / 'p 42 b c', 'c 42 c b'/
      data (altern(i,107), i=1,2) / 'i 4 m m', 'f 4 m m'/
      data (altern(i,108), i=1,2) / 'i 4 c m', 'f 4 m c'/
      data (altern(i,109), i=1,2) / 'i 41 m d', 'f 41 d m'/
      data (altern(i,110), i=1,2) / 'i 41 c d', 'f 41 d c'/
      data (altern(i,111), i=1,2) / 'p -4 2 m', 'c -4 m 2'/
      data (altern(i,112), i=1,2) / 'p -4 2 c', 'c -4 c 2'/
      data (altern(i,113), i=1,2) / 'p -4 21 m', 'c -4 m 21'/
      data (altern(i,114), i=1,2) / 'p -4 21 c', 'c -4 c 21'/
      data (altern(i,115), i=1,2) / 'p -4 m 2', 'c -4 2 m'/
      data (altern(i,116), i=1,2) / 'p -4 c 2', 'c -4 2 c'/
      data (altern(i,117), i=1,2) / 'p -4 b 2', 'c -4 2 b'/
      data (altern(i,118), i=1,2) / 'p -4 n 2', 'c -4 2 n'/
      data (altern(i,119), i=1,2) / 'i -4 m 2', 'f -4 2 m'/
      data (altern(i,120), i=1,2) / 'i -4 c 2', 'f -4 2 c'/
      data (altern(i,121), i=1,2) / 'i -4 2 m', 'f -4 m 2'/
      data (altern(i,122), i=1,2) / 'i -4 2 d', 'f -4 d 2'/
      data (altern(i,123), i=1,2) / 'p 4/m m m', 'c 4/m m m'/
      data (altern(i,124), i=1,2) / 'p 4/m c c', 'c 4/m c c'/
      data (altern(i,125), i=1,2) / 'p 4/n b m', 'c 4/a m b'/
      data (altern(i,126), i=1,2) / 'p 4/n n c', 'c 4/a c n'/
      data (altern(i,127), i=1,2) / 'p 4/m b m', 'c 4/m m b'/
      data (altern(i,128), i=1,2) / 'p 4/m n c', 'c 4/m c n'/
      data (altern(i,129), i=1,2) / 'p 4/n m m', 'c 4/a m m'/
      data (altern(i,130), i=1,2) / 'p 4/n c c', 'c 4/a c c'/
      data (altern(i,131), i=1,2) / 'p 42/m m c', 'c 42/m c m'/
      data (altern(i,132), i=1,2) / 'p 42/m c m', 'c 42/m m c'/
      data (altern(i,133), i=1,2) / 'p 42/n b c', 'c 42/a c b'/
      data (altern(i,134), i=1,2) / 'p 42/n n m', 'c 42/a m n'/
      data (altern(i,135), i=1,2) / 'p 42/m b c', 'c 42/m c b'/
      data (altern(i,136), i=1,2) / 'p 42/m n m', 'c 42/m m n'/
      data (altern(i,137), i=1,2) / 'p 42/n m c', 'c 42/a c m'/
      data (altern(i,138), i=1,2) / 'p 42/n c m', 'c 42/a m c'/
      data (altern(i,139), i=1,2) / 'i 4/m m m', 'f 4/m m m'/
      data (altern(i,140), i=1,2) / 'i 4/m c m', 'f 4/m m c'/
      data (altern(i,141), i=1,2) / 'i 41/a m d', 'f 41/d d m'/
      data (altern(i,142), i=1,2) / 'i 41/a c d', 'f 41/d d c'/


      test = 'ab'
      oldgrp = spcgrp
      call case(test, oldgrp)
      iperm = 1

c  check monoclinic settings
      do 20 i=3,15
        do 10 j=2,1,-1
          if (oldgrp.eq.altern(j,i)) then
              iperm  = j+10
              ispa   = i
              spcgrp = altern(2,i)
              call messag('* * WARNING!')
              call messag('Your crystal is monoclinic.   Atoms '//
     $                    'cannot definitively resolve')
              call messag('the ambiguities in crystal setting using '//
     $                    'the standard short symbols.')
              call messag('Consult The International Tables of '//
     $                    'X-Ray Crystallography for details on')
              call messag('permuting monoclinic axes if the output '//
     $                    'is unsatisfactory.')
              goto 999
          endif
 10     continue
 20   continue

c  check orthorhombic settings
      do 120 i=16,74
        do 110 j=1,6
          if (oldgrp.eq.altern(j,i)) then
              iperm  = j
              ispa   = i
              spcgrp = altern(1,i)
              goto 999
          endif
 110    continue
 120  continue

c  check tetragonal settings
      do 220 i=75,142
        do 210 j=1,2
          if (oldgrp.eq.altern(j,i)) then
              iperm  = j+20
              ispa   = i
              spcgrp = altern(1,i)
              goto 999
          endif
 210    continue
 220  continue


 999  continue
c       print*,'standard space group,   input space group'
c       print*,spcgrp, oldgrp, iperm

      return
c  end subroutine setort
      end
