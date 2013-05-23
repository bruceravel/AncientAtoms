      subroutine atspec(spcgrp,cell)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c  spcgrp: character string with space group designation
c  cell:   lattice parameters, a,b,c,alpha,beta & gamma
c----------------------------------------------------------------------
c  interpret words corresponding to very common space groups
c  set gamma=120 for hcp and graphite
c----------------------------------------------------------------------

      character*10 spcgrp, test*2
      dimension cell(6)

c  test must be in the same case as the options listed below
      test = 'ab'
      call case(test,spcgrp)

      if ((spcgrp(1:3).eq.'hex').or.(spcgrp(1:3).eq.'hcp')) then
          spcgrp='p 63/m m c'
          cell(6)=120
      elseif (spcgrp.eq.'fcc') then
          spcgrp='f m 3 m'
      elseif (spcgrp.eq.'bcc') then
          spcgrp='i m 3 m'
      elseif (spcgrp(1:3).eq.'cub') then
          spcgrp='p m 3 m'
      elseif ((spcgrp.eq.'salt').or.(spcgrp.eq.'nacl')) then 
          spcgrp='f m 3 m'
      elseif ((spcgrp.eq.'cscl').or.(spcgrp(1:3).eq.'ces')) then
          spcgrp='p m 3 m'
      elseif (spcgrp(1:5).eq.'perov') then
          spcgrp='p m 3 m'
      elseif ((spcgrp(1:5).eq.'zincb').or.(spcgrp.eq.'zns')) then
          spcgrp='f -4 3 m'
      elseif (spcgrp(1:3).eq.'dia') then
          spcgrp='f d 3 m'
      elseif (spcgrp(1:3).eq.'gra') then
          spcgrp='p 63 m c'
          cell(6)=120
      endif

      return
c end subroutine atspec
      end

