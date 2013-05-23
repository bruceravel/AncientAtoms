      subroutine feffpr(iat,natx,ntitx,ndopx,maxln,nexafs,nlogx,
     $                ifeff,iunit,itot,ntit,imult,
     $                edge,core,tag,titl,tagcen,exafs,atlis,
     $                logic,
     $                idop,dopant,percnt,elemnt,
     $                tglist,xwrite,ywrite,zwrite,rwrite,index,npot)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c  this produces a feff.inp file that will run to completion from a
c  list of atomic coordinates.  the feff.inp file will make resonable
c  guesses at atom lables, ipots, hole type, and rmax.  it will write
c  all ones for the control card and all zeros for the print card.
c  it will also write several useful keywords but comment them with an
c  asterisk.  the absorption and the density of the material and the
c  cluster size will be written as title lines 2 and 3.
c----------------------------------------------------------------------
c  input:
c    iunit:  scratch file unit number
c    edge:   k, l1, l2, or l3
c    exafs:  array containing exafs calculations from mcm module
c    dopant: symbol of doping atom
c    replcd: symbol of atom replaced by dopant
c    percnt: occupation of dopant
c    itot:   # of atoms in cluster
c----------------------------------------------------------------------
      parameter (nw=6, maxty=7, nall=2*maxty)
c maxln: max # of atom lines for feff
c maxty: max # of unique potentials for feff
c      parameter (iat=50, ntitx=9, ndopx=4, maxln=400)

      character*2  edge, dopant(iat,ndopx), cnum, edgeup, elemnt(iat)
      character*2  at, attype, atlist(0:maxty), test
      character*10 dwarf(0:nall), core, toss
      character*10 tglist(maxln), tglast, tag(iat), tagcen
      character*13 tgword
c       character*20 words(nw)
      character*20 feffcd
      character*23 holewd
      character*72 titl(ntitx)
c       character*80 str
      dimension    atlis(natx,8)
      dimension    xwrite(maxln), ywrite(maxln), zwrite(maxln),
     $             rwrite(maxln)
      dimension    percnt(iat,ndopx), exafs(nexafs)
      dimension    npot(maxln), index(maxln), nz(0:maxty), idop(iat),
     $             imult(iat), istoi(0:maxty)
      logical      logic(nlogx), ldop, lgcd
      external     s2e

      data(dwarf(i),i=0,7) /'snow white','dopey','grumpy','bashful',
     $                      'sleepy','doc','happy','sneezy'/
      data(dwarf(i),i=8,nall) /maxty*'ed'/

 4000 format(1x,a)
 4002 format(a)
 4005 format(i2.2)
 4010 format(1x,'*       ',a2,' is the central atom.'/,1x,'*       ',
     $'the atom list describes the undoped structure.')
 4100 format(1x,14('* -- '),'*')
 4200 format(1x,a,2x,f9.5)
 4400 format(bn,f10.0)
 4500 format(6x,i2,3x,i2,3x,a)
 4505 format(6x,i2,3x,i2,4x,a,5x,i2,6x,i2,6x,i2)
 4600 format(1x,'*       total mu = ',f10.1,' cm^-1, delta mu = ',f10.1,
     $' cm^-1')
 4650 format(1x,'*       specific gravity = ',f6.3,', cluster',
     $' contains ',i4,' atoms.')
 4655 format(1x,'*       mcmaster corrections: ',f8.5,' ang^2 and ',
     $e10.3,' ang^4')
 4660 format(1x,'*       self-abs. corrections: amplitude factor = ',
     $f6.3)
 4662 format(1x,'*                             ',f8.5,' ang^2 and ',
     $e10.3,' ang^4')
 4670 format(1x,'*       i0 corrections:       ',f8.5,' ang^2 and ',
     $e10.3,' ang^4')
 4680 format(1x,'*       sum of corrections:   ',f8.5,' ang^2 and ',
     $e10.3,' ang^4')
 4685 format(1x,'*       xanes corrections, [barns/cell]:   mu(EF) = ',
     $          g11.6,/,1x,'*              mu_back = ',g11.6,
     $          '       mu_cen = ',g11.6)
 4690 format(1x,'*       ',a2,' substituted ',f6.1,'% for ',a10)

 4700 format(1x,3(f9.5,3x),i2,3x,a13,3x,f8.5)
 4800 format(1x,a,3x,a)
 4850 format(1x,a6,'   geom.dat file generated with one',
     $       ' bounce flags.')

c----------------------------------------------------------
c  initialize things
      ldop   = .false.
      ipot   = 0
      test   = 'ab'
      attype = ' '
      nwds   = nw
      do 100 i=1,maxty
        atlist(i) = ' '
        istoi(i)  = 0
 100  continue

c----------------------------------------------------------
c  read data from scratch file
c       rewind(iunit)
c  get the titles
c       do 105 i=1,ntit
c         read(iunit,4002,end=140)titl(i)
c         call triml(titl(i))
c 105   continue
c  read the line containing the central atom
c      read(iunit,4000,end=140)str
c      call bwords(str,nwds,words)
      tglist(1) = tagcen
      atlist(0) = core
      istoi(0)  = 0
c       atlist(0) = words(1)
      nz(0)     = is2z(core)
      npot(1)   = 0
      xwrite(1) = atlis(1,6)
      ywrite(1) = atlis(1,7)
      zwrite(1) = atlis(1,8)
      rwrite(1) = atlis(1,4)
c      read(words(2),4400)atlis(1,6)
c      read(words(3),4400)atlis(1,7)
c      read(words(4),4400)atlis(1,8)
c      read(words(6),4400)atlis(1,4)
c  read in all the rest, keep track of ipots by crude criteria
c  keep track of atom types after the central atom.  if an atom type is
c  repeated do not give it a new ipot.
      nat = min(itot,maxln)
      do 130 i=2,nat
c        nwds=nw
c        read(iunit,4000,end=140)str
c        call triml(str)
c        call bwords(str,nwds,words)
        kat = nint(atlis(i,5))
        at = dopant(kat,1)
        call case(test,at)
        xwrite(i) = atlis(i,6)
        ywrite(i) = atlis(i,7)
        zwrite(i) = atlis(i,8)
        rwrite(i) = atlis(i,4)
c         read(words(2),4400)atlis(i,6)
c         read(words(3),4400)atlis(i,7)
c         read(words(4),4400)atlis(i,8)
c         read(words(6),4400)atlis(i,4)
c       check the just read atom symbol against the last one.  if it is
c         different then check against the entire list.  if it is new then
c         increment ipot.  if it has been seen before then mark jpot with
c         its previous ipot value and store it in npot(i).
c       simplest possible assignment of ipots: one unique potential
c         for every unique atomic species (not crystallographic site!)
c       at:  the present atomic symbol
c       attype:  the last symbol
c       npot, tglist:  lists of ipots and tags for the atom list
c       nz, atlist:  list of z's and symbols for the potentials list
c        jpot = 0
        if (at.ne.attype) then
            do 110 k=1,ipot
              toss=atlist(k)
              call case(test,toss)
              if (at.eq.toss) then
                  jpot=k
                  goto 120
              endif
 110        continue
            ipot         = ipot+1
            jpot         = ipot
            atlist(ipot) = at
            nz(ipot)     = is2z(at)
        endif
 120    continue
        attype    = at
        tglist(i) = tag(kat)
        npot(i)   = jpot
 130  continue
 140  continue
c  nat=i-1

      if (.not.logic(19)) goto 1045
      do 1010 ia=0,iat
        do 1000 ip=1,ipot
          if (elemnt(ia).eq.atlist(ip)) then
              istoi(ip) = istoi(ip)+imult(ia)
          endif
 1000   continue
 1010 continue
      minsto=10000
      isfact = 1
      do 1020 ip=1,ipot
        minsto = min(minsto, istoi(ip))
 1020 continue
c      print *, minsto
      if (minsto.ne.1) then
         lgcd = .false.
         do 1040 is=minsto,1,-1
           do 1030 ip=1,ipot
             anum = real(istoi(ip))/real(is)
             if ((anum - int(anum)) .lt. 0.00001 ) then
                isfact = is
                lgcd = .true.
             else
                lgcd = .true.
                goto 1035
             endif
 1030      continue
           do 1032 ip=1,ipot
             istoi(ip) = istoi(ip) / isfact
 1032      continue
 1035      continue
 1040    continue
      endif
 1045 continue

c--------------------------------------------------------------------
c  interpret the atom indexing if lindex=true
c  search through the list of atoms.  for each atom, compare the tag
c  to the previous tag (tglast).  if it is the same (ie same crys.
c  site) then check distance.  if the distance is different then
c  increment the index.  if distance is the same then use the same
c  index.  if the present tag is different from tglast then search
c  backwards through all previous atoms to see if that tag (site)
c  has been used before.  if it has, then continue indexing from the
c  index last used for that site.  if it has not been seen before,
c  then begin indexing from 1.
      if (logic(4)) then
          index(1) = 0
          tglast   = tglist(1)
          do 148 j=2,nat
            if (tglist(j).eq.tglast) then
                jr     = int(rwrite(j)   * 10000)
                jrprev = int(rwrite(j-1) * 10000)
                if (jr.gt.jrprev) then
                    index(j) = index(j-1) + 1
                else
                    index(j) = index(j-1)
                endif
            else
                tglast = tglist(j)
                index(j) = 1
                do 143 k=j-2,1,-1
                  if (tglist(j).eq.tglist(k) ) then
                      index(j) = index(k) + 1
                      goto 146
                  endif
143             continue
146             continue
            endif
148       continue
      endif

c----------------------------------------------------------------
c  results of mcmaster module calculations
      if (logic(28)) then
          write(ifeff,4100)
          write(ifeff,4600)exafs(1),exafs(2)
          write(ifeff,4650)exafs(3),itot
          write(ifeff,4100)
          write(ifeff,4655)exafs(4),exafs(5)
          if (logic(3)) then
              write(ifeff,4660)exafs(6)
              write(ifeff,4662)exafs(7),exafs(8)
              write(ifeff,4670)exafs(9),exafs(10)
              sumsig = exafs(4) + exafs(7) + exafs(9)
              sumqrt = exafs(5) + exafs(8) + exafs(10)
              write(ifeff,4100)
              write(ifeff,4680)sumsig,sumqrt
c           write(ifeff,4100)
c           write(ifeff,4685)exafs(11),exafs(12),exafs(13)
          endif
          write(ifeff,4100)
      endif

c----------------------------------------------------------------
c  dopant information
      do 152 i=1,iat
        do 150 j=2,idop(i)
          call fixsym(dopant(i,j))
          call fixsym(dopant(i,1))
          write(ifeff,4690)dopant(i,j),percnt(i,j)*100.,tag(i)
          ldop = .true.
 150    continue
 152  continue
      if (ldop) then
          call fixsym(core)
          write(ifeff,4010)core
          write(ifeff,4100)
      endif

c----------------------------------------------------------------
c  titles
      write(ifeff,4000)' '
      call card('title',feffcd,ii)
      do 155 i=1,ntit
        ij = istrln(titl(i))
        ij = min(ij, 70)
        write(ifeff,4800)feffcd(:ii),titl(i)(:ij)
155   continue
      write(ifeff,4000)' '

c  -- **  -- **  -- **  -- **  -- **  -- **  -- **  -- **  -- **  --
c     logic(19) is the feff8 flag, control and print take six args
c     logic(18) is the xanes flag
c     logic(6)  is the nogeom flag
c  -- **  -- **  -- **  -- **  -- **  -- **  -- **  -- **  -- **  --

c----------------------------------------------------------------
c  hole
      if (logic(19)) then
          call card('edge',feffcd,ii)
          edgeup = edge
          call upper(edgeup)
          write(ifeff,4000) feffcd(:ii) // '      ' // edgeup
          call card('s02',feffcd,ii)
          write(ifeff,4000) feffcd(:ii) // '       1.0'
      else
          call card('hole',feffcd,ii)
          call fixsym(atlist(0))
          if (edge.eq.'l3') then
              holewd = ' 4   1.0     '//atlist(0)//' L3 edge'
              edgen = s2e(atlist(0), 'l3')
          elseif (edge.eq.'l2') then
              holewd = ' 3   1.0     '//atlist(0)//' L2 edge'
              edgen = s2e(atlist(0), 'l2')
          elseif (edge.eq.'l1') then
              holewd = ' 2   1.0     '//atlist(0)//' L1 edge'
              edgen = s2e(atlist(0), 'l1')
          elseif (edge.eq.'k') then
              holewd = ' 1   1.0     '//atlist(0)//' K edge'
              edgen = s2e(atlist(0), 'k')
          else
              if (nz(0).gt.57) then
                  holewd = ' 4   1.0     '//atlist(0)//' L3 edge'
                  edgen = s2e(atlist(0), 'l3')
              else
                  holewd = ' 1   1.0     '//atlist(0)//' K edge'
                  edgen = s2e(atlist(0), 'k')
              endif
          endif
          write(ifeff,4900)feffcd(:ii), holewd, edgen
 4900     format (1x, a4, a23, 1x, '(', f7.3,
     $                ' keV), second number is S0^2')
      endif

c----------------------------------------------------------------
c  control, print, rmax cards
      call card('control',feffcd,ii)
      write(ifeff,4000)' '
      if (logic(19)) then
          write(ifeff,4000)
     $                '*         pot    xsph  fms   paths genfmt ff2chi'
      else
          write(ifeff,4000)     '*         mphase,mpath,mfeff,mchi'
      endif
      if (logic(18)) then
          write(ifeff,4000)feffcd(:ii)//'   1      0     0     0'
      elseif (logic(19)) then
          write(ifeff,4000)feffcd(:ii)//
     $                '   1      1     1     1     1      1'
      else
          write(ifeff,4000)feffcd(:ii)//'   1      1     1     1'
      endif

      call card('print',feffcd,ii)
      if (logic(6)) then
          write(ifeff,4000)feffcd(:ii)//'     1      2     0     3'
      elseif (logic(19)) then
          write(ifeff,4000)feffcd(:ii)//
     $                '     1      0     0     0     0      0'
      else
          write(ifeff,4000)feffcd(:ii)//'     1      0     0     3'
      endif
      write(ifeff,4000)' '

      wrmin = 2.2*rwrite(2)
      wrmax = rwrite(nat)+0.00001
      if (.not.logic(19)) then
          call card('rmax',feffcd,ii)
          write(ifeff,4200)feffcd(:ii)//'   ',wrmax
      else
          call card('scf',feffcd,ii)
          write(ifeff,4000) '*         r_scf   [ l_scf  n_scf  ca ]'
          write(ifeff,4910)feffcd(:ii)//
     $                '       ', wrmin, '   0      15     0.1'
 4910     format(1x,a,f7.5,a)
      endif
      if (logic(6)) then
          call card('nogeom',feffcd,ii)
          write(ifeff,4850)feffcd(:ii)
      endif
      write(ifeff,4000)' '

c-----------------------------------------------------------------
c  various feff cards
      if (logic(18)) then
          call card('xanes',feffcd,ii)
          write(ifeff,4200)feffcd(:ii)//'    ',rwrite(nat)+0.00001
          call card('vintfix',feffcd,ii)
          write(ifeff,4000)'*'//feffcd(:ii)//'  10.0'
          call card('egrid',feffcd,ii)
          write(ifeff,4000)'*'//feffcd(:ii)//'    1  80'
          call card('emesh',feffcd,ii)
          write(ifeff,4000)'*'//feffcd(:ii)//'    1'
          call card('exchange',feffcd,ii)
          write(ifeff,4000)feffcd(:ii)//'  2  0  0'
      elseif (logic(19)) then
          call card('exchange',feffcd,ii)
          write(ifeff,4000)    '*         ixc  [ Vr  Vi ]'
          write(ifeff,4000)feffcd(:ii)//'  0      0   0'
          write(ifeff,4000)' '
          call card('exafs',feffcd,ii)
          write(ifeff,4000)feffcd(:ii)
          call card('rpath',feffcd,ii)
          write(ifeff,4200)feffcd(:ii)//' ', 2*wrmin
          write(ifeff,4000)' '
          call card('xanes',feffcd,ii)
          write(ifeff,4000)'*         kmax  [ delta_k  delta_e ]'
          write(ifeff,4000)
     $                '*'//feffcd(:ii)//'     4.0     0.07     0.5'
          call card('fms',feffcd,ii)
          write(ifeff,4000)'*         r_fms     [ l_fms ]'
          write(ifeff,4207)'*'//feffcd(:ii), wrmin, 0
 4207     format(1x,a,4x,2x,f8.5,4x,i2)
          write(ifeff,4000)'*'
          call card('rpath',feffcd,ii)
          write(ifeff,4200)'*'//feffcd(:ii)//' ', 0.1
          call card('ldos',feffcd,ii)
          write(ifeff,4000)'*         emin  emax  resolution'
          write(ifeff,4000)'*'//feffcd(:ii)//'      -20    20   0.1'
      else
          call card('criteria',feffcd,ii)
          write(ifeff,4000)'*'//feffcd(:ii)//'     curved   plane'
          call card('debye',feffcd,ii)
          write(ifeff,4000)'*'//feffcd(:ii)//
     $                '        temp     debye-temp'
          call card('nleg',feffcd,ii)
          write(ifeff,4000)'*'//feffcd(:ii)//'         8'
      endif
      write(ifeff,4000)' '

c----------------------------------------------------------
c  write out the potential list, list is longer for feff8
      call card('potentials',feffcd,ii)
      write(ifeff,4000)feffcd(:ii)
      if (logic(19)) then
          write(ifeff,4000)
     $                '*   ipot   z [ label   l_scmt  l_fms  '//
     $                'stoichiometry ]'
      else
          write(ifeff,4000)'*   ipot   z  label'
      endif
      do 160 i=0,ipot
c        istoi = 0
c        if (i .ne. 0) istoi = imult(i)
        if (logic(5)) then
            if (logic(19)) then
                llmm = 3
                if (nz(i).le.36) llmm = 2
                if (nz(i).le.10) llmm = 1
c                write(ifeff,4505)i,nz(i),dwarf(i),llmm,llmm
                write(ifeff,4505)i,nz(i),dwarf(i),-1,-1,istoi(i)
            else
                write(ifeff,4500)i,nz(i),dwarf(i)
            endif
        else
          call fixsym(atlist(i))
          if (logic(19)) then
              llmm = 3
              if (nz(i).le.36) llmm = 2
              if (nz(i).le.10) llmm = 1
              write(ifeff,4505)i,nz(i),atlist(i),-1,-1,istoi(i)
          else
              write(ifeff,4500)i,nz(i),atlist(i)
          endif
        endif
160   continue

c----------------------------------------------------------
c  write out the atom list
      write(ifeff,4000)' '
      call card('atoms',feffcd,ii)
      write(ifeff,4000)feffcd(:ii)
      do 170 i=1,nat
        if (logic(4)) then
            itgwd  = istrln(tglist(i))
            write(cnum,4005)index(i)
            tgword = tglist(i)(1:itgwd)//'_'//cnum
        else
            tgword = tglist(i)
        endif
        write(ifeff,4700)xwrite(i),ywrite(i),zwrite(i),npot(i),
     $              tgword,rwrite(i)
170   continue

      call card('end',feffcd,ii)
      write(ifeff,4000)feffcd(:ii)

      return
c end subroutine feffpr
      end
