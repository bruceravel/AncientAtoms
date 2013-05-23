      subroutine origin(spcgrp, warn, wrning)
c--------------------------------------------------------------
c  copyright (c) 1998 Bruce Ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)

c  input:  spcgrp:  character*10 containing Hermann-Maguin symbol
c  output: warn:    logical, true if need to write message
c          wrning:  (4) of character*74 containing warning message

c  search a list to see if the chosen space group is among those with
c  two origins reported in itxc.  if it is, give a run-time message
c  telling  the user to check the output and, if it is yucky, to shift
c  atom coords by a proscribed amount.

      parameter (norg=22)
c  norg: # of space groups with two origins reported in itxc
      character*10 sp(norg), spcgrp, space, test*2
      character*23 shift(norg)
      parameter(nwarn = 4)
      character*74 wrning(nwarn)
      logical      warn

      data (sp(i),i=1,norg)/
     $            'p n n n'   , 'p b a n',
     $            'p m m n'   , 'c c c a',
     $            'f d d d'   , 'p 4/n b m',
     $            'p 4/n n c' , 'p 4/n m m',
     $            'p 4/n c c' , 'p 4/n b c',
     $            'p 42/n n m', 'p 42/n m c',
     $            'p 42/n c m', 'i 41/a m d',
     $            'i 41/a c d', 'p n 3',
     $            'f d 3'     , 'p n 3 n',
     $            'p n 3 m'   , 'f d 3 m',
     $            'f d 3 c'   , ' '/

      data (shift(i),i=1,norg)/
     $            ' 0.25   0.25   0.25' , ' 0.25   0.25   0',
     $            ' 0.25   0.25   0'    , ' 0      0.25   0.25',
     $            '-0.125 -0.125 -0.125', '-0.25  -0.25   0',
     $            '-0.25  -0.25  -0.25' , '-0.25   0.25   0',
     $            '-0.25   0.25   0'    , '-0.25   0.25  -0.25',
     $            '-0.25   0.25  -0.25' , '-0.25   0.25  -0.25',
     $            '-0.25   0.25  -0.25' , ' 0      0.25  -0.125',
     $            ' 0      0.25  -0.125', '-0.25  -0.25  -0.25',
     $            '-0.125 -0.125 -0.125', '-0.25  -0.25  -0.25',
     $            '-0.25  -0.25  -0.25' , '-0.125 -0.125 -0.125',
     $            '-0.375 -0.375 -0.375', ' '/

      test = 'ab'
      call case(test,spcgrp)

      warn=.false.
      ishift=norg
      do 10 i=1,norg
        if (spcgrp.eq.sp(i)) then
            ishift = i
            space = sp(i)
            warn=.true.
        endif
 10   continue

      if (warn) then
          m1=istrln(space)
          call upper(space)
          m2=istrln(shift(ishift))

          wrning(1) = '     Space group "'//space(:m1)//'" may be '//
     $                'referenced to a different origin.'

          wrning(2) = '     If the atom list seems incorrect, '//
     $                'put this line in your input file'

          wrning(3) = '            shift  '//shift(ishift)(:m2)

          wrning(4) = '     and run atoms again.-'

      endif
      return
c end subroutine origin
      end
