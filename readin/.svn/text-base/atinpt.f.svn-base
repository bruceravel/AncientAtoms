      subroutine atinpt(iat,ntitx,ndopx,nlogx,
     $            title,ntit,iatom,ibasis,iabs,dmax,ispa,iperm,
     $            idop,iedge,nepts,nrefl,isystm,nnoan,
     $            elemnt,tag,noantg,edge,core,spcgrp,inpgrp,
     $            outfil,shwarn,afname,refile,dopant,
     $            x,y,z,cell,percnt,gasses,qvect,egr,
     $            logic, stdout, expnd, shift)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------
c  iatom   : number of unique atoms in cell
c  ibasis  : number of atoms in basis
c  ispa    : number between 1 and 230, index of s.g. in int'l tables
c  iperm   : index of permutation matrix for non-standard setting
c  title   : user supplied comment line
c  elemnt  : character array of atom types
c  tag     : character array of site tags
c  spcgrp  : chararcter string with hermann-maguin space group designation
c  edge    : absorption edge of core atom, k or l3
c  outfil  : output file name
c  refile  : reflection amplitude file name
c  dopant  : dopant element symbol
c  x,y,z   : positions of atoms in cell-axis coordinates
c  cell    : lattice constants, a,b,c,alpha,beta & gamma
c  dmax    : maximum distance in cluster
c  percent : percent of dopant
c  logic   : various flags, see arrays.map
c  noantg  : tags for which anomalous correction is to be neglected
c---------------------------------------------------------------------
c  this parses the lines of the command file looking for keywords then
c  reads in the value for that keyword from the next word.  the
c  structure is nearly free format except for the atom list which must
c  come at the end of the command file,  the reason for this is that
c  "b" and "c" are keywords and atomic symbols thus enforcing some
c  structure in the command file is the easiest way to distinguish them.
c---------------------------------------------------------------------

      parameter(ndpmax=10, nwdx=20)
c      parameter(iat=50, ntitx=9, ndopx=4 )
      parameter(zero=0)
      character*2  elemnt(iat),edge,dopant(iat,ndopx),
     $             doplis(ndpmax),test
      character*10 spcgrp,inpgrp,tag(iat),replcd(ndpmax),core,co,tagup
      character*10 noantg(iat)
      character*20 words(nwdx)
      character*72 title(ntitx), titln, outfil, fname, afname, refile
      character*80 string,toss,messg*78
      logical      logic(nlogx), stdout, expnd
      logical      inpnul, lcore, lshift, there
      dimension    x(iat), y(iat), z(iat), cell(6), qvect(3)
      dimension    idop(iat), nrefl(3), gasses(3)
      dimension    percnt(iat,ndopx), pclis(ndpmax)
      parameter(nshwrn=4)
      character*74 shwarn(nshwrn)
      logical      shift

 1410 format(bn,f10.0)
 1440 format(a)
 1450 format(a,i2)
 4000 format(i2)
 4010 format(' *** Notice at line ',i3)
 4020 format(' *** Warning at line ',i3)

c---------------------------------------------------------------------
c  initialize some things used only in this routine
      inpnul = .true.
      lcore  = .false.
      lshift = .false.
      icore  = 0
      ndop   = 0
      ntit   = 0
      xshift = zero
      yshift = zero
      zshift = zero
      nnoan  = 0
      test   = 'ab'
      nline  = 0
      ierr   = 0
      iertot = 0
c  the value of the variable test must be in the same case as the
c  keyword names in the long block of elseif's

c---------------------------------------------------------------------
c  open atom.inp, look for upper and lower case file names for both
c  atom.inp and atoms.inp
      if (.not.stdout) then
          fname = 'atoms.inp'
          call lower(fname)
          inquire(file=fname,exist=there)
          if (.not.there) then
              call upper(fname)
              inquire(file=fname,exist=there)
          endif
          if (.not.there) then
              fname = 'atom.inp'
              call lower(fname)
              inquire(file=fname,exist=there)
              if (.not.there) then
                  call upper(fname)
                  inquire(file=fname,exist=there)
              endif
          endif
          if (.not.there) then
              call messag('Input file for ATOMS is not found. '//
     $                    'Hasta luego.')
              stop
          endif
          open(unit=1,file=fname,status='old')
      endif

c---------------------------------------------------------------------
c  begin reading the input file, words must be cleared each time to
c  avoid unintentionally labling more than one atom type as the
c  central atom.
 101  continue
      if (stdout) then
          read (*,1440,end=191,err=191)string
      else
          read (1,1440,end=191,err=191)string
      endif
      nline = nline+1
      call untab(string)
c       call uncomm(string)

 120  continue
      nwds = nwdx
      do 122 iw=1,nwds
        words(iw)=' '
 122  continue
      call triml(string)
c                                           - denotes end of job
      if  (string(1:1).eq.'-') goto 191
c                                           skip a comment line
      if  ((string(1:1).eq.'!').or.(string(1:1).eq.'*')
     $ .or.(string(1:1).eq.'%').or.(string(1:1).eq.'#')) goto 101
c                                           skip a blank line
      if  (string.eq.' ') goto 101
c                                           begin reading line
      i=1
      call bwords(string,nwds,words)
      inpnul = .false.

c  ******************** input file parsing *************************
c  read a word, identify it, assign the value from the following word(s)
c  increment i and come back.  i points to position in string, when i
c  exceeds nwds go read a new line.
130   continue
      call case(test,words(i))
c                                           skip a blank line
      if     (words(i).eq.' ') then
          goto 101
c                                           ignore everything after !,*,%
      elseif ((words(i)(1:1).eq.'!').or.(words(i)(1:1).eq.'*').or.
     $        (words(i)(1:1).eq.'%').or.(words(i)(1:1).eq.'#')) then
          goto 101
      elseif (words(i).eq.'end') then
          goto 191
c                                           title and comment are synonyms
      elseif ((words(i).eq.'comment').or.(words(i).eq.'title')) then
          if (ntit.lt.ntitx) then
              call gettit(words(i), string, titln, ntit, stdout)
              title(ntit) = titln
          endif
          goto 101
c                                           read the next ten characters
c                                           into space, continue reading
c                                           rest of the line.
c                                           this one is perverse and must
c                                           be handled specially
      elseif (words(i).eq.'space') then
          toss = string
          call case(test,toss)
          m = index(toss, 'space')
          toss = string(m+6:)
          call triml(toss)
          string = toss
          if (string(:1).eq.'=') toss = string(2:)
          call triml(toss)
          spcgrp = toss(1:10)
          call triml(spcgrp)
          call case(test,spcgrp)
          string = toss(11:)
          inpgrp = spcgrp
          goto 120
c                                           outfile, default=feff.inp
      elseif (words(i)(1:3).eq.'out') then
          outfil=words(i+1)
          i=i+2
c                                           specified edge, default by z
      elseif ((words(i).eq.'edge').or.(words(i).eq.'hole'))  then
          edge=words(i+1)
          call case(test,edge)
          i=i+2
c                                           specified core
      elseif ((words(i).eq.'core').or.(words(i)(1:4).eq.'cent')) then
          core=words(i+1)
          call case(test,core)
          lcore = .true.
          i=i+2
c                                           dopants
      elseif (words(i)(1:3).eq.'dop') then
          ndop = ndop + 1
          if (ndop.gt.ndopx) then
              write(messg, 4020)nline
              call messag(messg)
              call messag('     You have exceeded the '//
     $                    'maximum number of dopants.')
              call messag('     ATOMS will ignore this and all '//
     $                    'further dopants.-')
              ierr = 2
              goto 137
          endif
          doplis(ndop) = words(i+1)
          call case(test,doplis(ndop))
          replcd(ndop) = words(i+2)
          call case(test,replcd(ndop))
          call getrea(words(i), words(i+3), pclis(ndop),
     $                nline, ierr)
 137      continue
          i=i+4
c                                           argon, fluorescence
      elseif (words(i)(1:3).eq.'arg') then
          call getrea(words(i), words(i+1), gasses(1), nline, ierr)
          logic(3) = .true.
          i=i+2
c                                           krypton, fluorescence
      elseif (words(i)(1:3).eq.'kry') then
          call getrea(words(i), words(i+1), gasses(2), nline, ierr)
          logic(3) = .true.
          i=i+2
c                                           nitrogen, fluorescence
      elseif (words(i)(1:3).eq.'nit') then
          call getrea(words(i), words(i+1), gasses(3), nline, ierr)
          logic(3) = .true.
          i=i+2
c                                           flag for indexing
      elseif (words(i)(1:3).eq.'ind') then
          call getlgc(words(i), words(i+1), logic(4), nline, ierr)
          i=i+2
c                                           geom.dat
      elseif (words(i)(1:3).eq.'geo') then
          call getlgc(words(i), words(i+1), logic(6), nline, ierr)
          i=i+2
c                                           xanes keywords (& not feff8)
c       elseif (words(i)(1:3).eq.'xan') then
c           call getlgc(words(i), words(i+1), logic(18), nline, ierr)
c           logic(19) = .false.
c           i=i+2
c                                           feff8 keywords (& not xanes)
      elseif (words(i)(1:5).eq.'feff8') then
          call getlgc(words(i), words(i+1), logic(19), nline, ierr)
          logic(18) = .false.
          i=i+2
c                                           run clustr & output (note order)
      elseif (words(i)(1:4).eq.'feff') then
          call getlgc(words(i), words(i+1), logic(7), nline, ierr)
          i=i+2
c                                           calc and write McM corrs.
      elseif (words(i)(1:4).eq.'corr') then
          call getlgc(words(i), words(i+1), logic(28), nline, ierr)
          i=i+2
c                                           heh, heh, heh!
      elseif (words(i)(1:3).eq.'dwa') then
          call getlgc(words(i), words(i+1), logic(5), nline, ierr)
          i=i+2
c * * * * * * * * * * diagnostic functions * * * * * * * * * * * *
c                                           mcmast.dat
      elseif (words(i)(1:3).eq.'mcm') then
          call getlgc(words(i), words(i+1), logic(14), nline, ierr)
          i=i+2
c                                           self.dat
      elseif (words(i).eq.'self') then
          call getlgc(words(i), words(i+1), logic(15), nline, ierr)
          i=i+2
c                                           i0.dat
      elseif (words(i).eq.'i0') then
          call getlgc(words(i), words(i+1), logic(16), nline, ierr)
          i=i+2
c                                           unit.dat
      elseif (words(i)(1:3).eq.'uni') then
          call getlgc(words(i), words(i+1), logic(9), nline, ierr)
          i=i+2
c                                           p1.inp
      elseif (words(i).eq.'p1') then
          call getlgc(words(i), words(i+1), logic(8), nline, ierr)
          i=i+2
c                                           f.dat, diagnostic for f'/"
      elseif (words(i)(1:4).eq.'fdat') then
          call getlgc(words(i), words(i+1), logic(17), nline, ierr)
          i=i+2
c                                           print module messages
      elseif (words(i)(1:3).eq.'mod') then
          call getlgc(words(i), words(i+1), logic(20), nline, ierr)
          i=i+2
c                                           print location messages
      elseif (words(i).eq.'message') then
          call getint(words(1),words(i+1), imess, nline, ierr)
          if ((imess.eq.0).or.(imess.eq.2)) logic(21) = .true.
          if ((imess.eq.0).or.(imess.eq.3)) logic(22) = .true.
          if ((imess.eq.0).or.(imess.eq.4)) logic(23) = .true.
          if ((imess.eq.0).or.(imess.eq.5)) logic(24) = .true.
          if ((imess.eq.0).or.(imess.eq.6)) logic(25) = .true.
          if ((imess.eq.0).or.(imess.eq.7)) logic(26) = .true.
          if (imess.eq.0) logic(20) = .true.
          i=i+2
c * * * * * * * * end diagnostic functions * * * * * * * * * * * *
c                                           rmax, default=5.0
      elseif (words(i)(1:3).eq.'rma') then
          call getrea(words(i), words(i+1), dmax, nline, ierr)
          i=i+2
c                                           the lattice constants
      elseif (words(i).eq.'a') then
          call getrea(words(i), words(i+1), cell(1), nline, ierr)
          i=i+2
      elseif (words(i).eq.'b') then
          call getrea(words(i), words(i+1), cell(2), nline, ierr)
          i=i+2
      elseif (words(i).eq.'c') then
          call getrea(words(i), words(i+1), cell(3), nline, ierr)
          i=i+2
c                                           the latice angles
      elseif (words(i)(1:3).eq.'alp') then
          call getrea(words(i), words(i+1), cell(4), nline, ierr)
          i=i+2
      elseif (words(i)(1:3).eq.'bet') then
          call getrea(words(i), words(i+1), cell(5), nline, ierr)
          i=i+2
      elseif (words(i)(1:3).eq.'gam') then
          call getrea(words(i), words(i+1), cell(6), nline, ierr)
          i=i+2
c                                           shift vector
      elseif (words(i).eq.'shift') then
          call getrea(words(i), words(i+1), xshift, nline, ierr)
          call getrea(words(i), words(i+2), yshift, nline, ierr)
          call getrea(words(i), words(i+3), zshift, nline, ierr)
          lshift = .true.
          i=i+4

c     ************ DAFS STUFF ****************

c                                           q vector for dafs
      elseif ((words(i)(1:4).eq.'qvec').or.(words(i).eq.'dafs')) then
          call getrea(words(i), words(i+1), qvect(1), nline, ierr)
          call getrea(words(i), words(i+2), qvect(2), nline, ierr)
          call getrea(words(i), words(i+3), qvect(3), nline, ierr)
          logic(10) = .true.
          i=i+4
      elseif (words(i).eq.'feout') then
          afname = words(i+1)
          i=i+2
c                                           reflection amplitudes
      elseif (words(i)(1:4).eq.'refl') then
          call getint(words(i), words(i+1), nrefl(1), nline, ierr)
          call getint(words(i), words(i+2), nrefl(2), nline, ierr)
          call getint(words(i), words(i+3), nrefl(3), nline, ierr)
          logic(11) = .true.
          i=i+4
c                                           reflection amplitude file name
      elseif (words(i).eq.'refile') then
          refile = words(i+1)
          i=i+2
c                                           # of grid points for dafs
      elseif (words(i)(1:3).eq.'nep') then
          call getint(words(i), words(i+1), nepts, nline, ierr)
          i=i+2
c                                           grid spacing for dafs
      elseif (words(i)(1:3).eq.'egr') then
          call getrea(words(i), words(i+1), egr, nline, ierr)
          i=i+2
c                                           neglect anomalous correction
      elseif (words(i)(1:4).eq.'noan') then
          nnoan = nnoan+1
          noantg(nnoan) =  words(i+1)
          i=i+2

c                                           beneath the word atom is a
c                                           5 col. list of unique atoms
c                                           in this order:
c                                              sym  x  y  z   center?
c                                           atom info is read in until eof
c                                           or - is found
      elseif (words(i)(1:3).eq.'ato') then
 140      continue
          if (stdout) then
              read (*,1440,end=191,err=191)string
          else
              read (1,1440,end=191,err=191)string
          endif
          nline = nline+1
          call untab(string)
          call triml(string)
          if (string(1:1).eq.'-') goto 191
          if ( (string(1:1).eq.' ').or.(string(1:1).eq.'!').or.
     $         (string(1:1).eq.'*').or.(string(1:1).eq.'%').or.
     $         (string(1:1).eq.'#'))
     $                goto 140
          iatom=iatom+1
          if (iatom.gt.iat) then
              ierr = 3
 400          format(' *** Error: Your atoms list exceeds ',i3,
     $                    ', the hardwired limit.')
              write(messg, 400)iat
              call messag(' ')
              call messag(messg)
              call messag('     Reset iat in the source code then ')
              call messag('     recompile ATOMS to handle a larger '//
     $                    'atom list.-')
              goto 191
          endif
          if (expnd) then
              call messag(' *** Error: No expanded atoms list yet.-')
              stop
          else
              call getatm(nline, ierr, string, elemnt(iatom),
     $                    x(iatom), y(iatom), z(iatom),
     $                    tag(iatom))
              if (ierr.eq.2) iatom=iatom-1
          endif
          if (ierr.eq.3) iertot = 1
          ierr = 0
          goto 140
c                                           handle basis
c                                           iatom set to one, one atom
c                                           will be expanded as a
c                                           point and rest will be added
      elseif (words(i)(1:3).eq.'bas') then
          logic(2) = .true.
          iatom  = iatom+1
 155      continue
          if (stdout) then
              read (*,1440,end=191,err=191)string
          else
              read (1,1440,end=191,err=191)string
          endif
          nline = nline+1
          call untab(string)
          call triml(string)
          if (string(1:1).eq.'-') goto 191
          if ( (string(1:1).eq.' ').or.(string(1:1).eq.'!').or.
     $         (string(1:1).eq.'*').or.(string(1:1).eq.'%').or.
     $         (string(1:1).eq.'#'))
     $                goto 155
          ibasis=ibasis+1
          if (ibasis.gt.iat) then
              ierr = 3
 410          format(' *** Error: Your basis list exceeds ',i3,
     $                    ', the hardwired limit.')
              write(messg, 400)iat
              call messag(' ')
              call messag(messg)
              call messag('     Reset iat in the source code then ')
              call messag('     recompile ATOMS to handle a larger '//
     $                    'atom list.-')
              goto 191
          endif
          if (expnd) then
              call messag(' *** Error: No expanded atoms list yet.-')
              stop
          else
              call getatm(nline, ierr, string, elemnt(ibasis),
     $                    x(ibasis), y(ibasis), z(ibasis),
     $                    tag(ibasis))
              if (ierr.eq.2) iatom=iatom-1
          endif
          if (ierr.eq.3) iertot = 1
          ierr = 0
          goto 155
      else
          write(messg, 4020)nline
          call messag(messg)
          iunk  = istrln(words(i))
          messg = '     "'//words(i)(:iunk)//'" is an unknown keyword.-'
          call messag(messg)
          ierr = 2
          i=i+2
      endif
c     if read entire line then read next line else read next word in line
      if (ierr.eq.3) iertot = 1
      ierr = 0
      if (i.ge.nwds) goto 101
      goto 130

c     done reading lines
191   continue
      if (ierr.eq.3) iertot = 1
      if (iertot.ne.0) then
          call messag(' ')
          call messag(' *** Ending early due to faulty input file.-')
          call messag(' ')
          ierr = 1
          stop
      endif
      if (inpnul) then
          logic(1)=.true.
          goto 300
      endif

c--------------------------------------------------------------------
c----- do a few things with the keyword values before leaving -------
c--------------------------------------------------------------------

c  turn off all messages if program compiled for standard in/output
c  also disable dafs, geom.dat, unit.dat, p1.dat outputs
      if (stdout) then
          logic(6)  = .false.
          logic(8)  = .false.
          logic(9)  = .false.
          logic(10) = .false.
          logic(11) = .false.
          logic(14) = .false.
          logic(15) = .false.
          logic(16) = .false.
          logic(17) = .false.
          logic(20) = .false.
          logic(21) = .false.
          logic(22) = .false.
          logic(23) = .false.
          logic(24) = .false.
          logic(25) = .false.
          logic(26) = .false.
      endif

c  iall is the number of atom in the atom or basis list
      iall = iatom
      if (logic(2)) iall = ibasis

c  add shift vector to unique atom coordinates
      do 250 i=1,iall
        x(i) = x(i) + xshift
        y(i) = y(i) + yshift
        z(i) = z(i) + zshift
 250  continue

c  identify space group from supplied symbol
      call groups

c  organize matrices containing site contents
      call dopfix(iat,ndopx,ndpmax,
     $            iall,ndop,tag,doplis,replcd,pclis,
     $            elemnt,dopant,percnt,idop)

c  try to determine the atomic symbol of the core atom
c  it need not be specified if the atom list has one site in it...
      if ((iatom.eq.1).and.(ibasis.le.1).and.(.not.lcore))  then
          core  = elemnt(1)
          iabs  = 1
          call case(test,core)
          lcore = .true.
      endif

c      print*,'initial value of core: ', core, ' iabs=', iabs
c  search through tag list for the specified core...
      co   = core
      call case(test,co)
      do 200 i=1,iall
        tagup=tag(i)
        call case(test,tagup)
        if (tagup.eq.co) then
            icore = icore + 1
            core  = elemnt(i)
            iabs  = i
            lcore = .true.
        endif
 200  continue
 210  continue
c      print*,'after searching tag list: ', core, ' iabs=', iabs

c  if the specified core was not found in the tag list, search the dopants
c  (use the variable name tagup to avoid defining another variable)
      if (iabs.eq.0) then
          do 280 i=1,iall
            do 270 j=2,idop(i)
              co   = core
              call case(test,co)
              tagup=dopant(i,j)
              call case(test,tagup)
              if (tagup.eq.co) then
                  iabs  = i
                  icore = icore+1
              endif
 270        continue
 280      continue
      endif

c----------------------------------------------------------------------
c  choose l3 (=4) or k (=1) edge, the numbers are chosen to suit mucal
      if (edge.eq.'k') then
          iedge = 1
      elseif (edge.eq.'l1') then
          iedge = 2
      elseif (edge.eq.'l2') then
          iedge = 3
      elseif (edge.eq.'l3') then
          iedge = 4
      elseif (edge.eq.' ') then
          if (is2z(core).gt.57) then
              edge  = 'l3'
              iedge = 4
          else
              edge  = 'k'
              iedge = 1
          endif
      endif

c  check if core has been specified, this is to satisfy backwards
c  incompatibility concerns
      if (.not.lcore) then
          call messag(' ')
          call messag(' *** Error:  while reading atom or basis list')
          call messag('     In this version of ATOMS, the central '//
     $                'atom is specified by the keyword "core"')
          call messag('     The fifth column of the atom list is '//
     $                'reserved for the site tag.')
          call messag('Please edit the input file and run atoms '//
     $                'again.-')
          call messag(' ')
          stop
      endif

c  check if 0 sites match core
      if (icore.eq.0) then
          ii = istrln(co)
          call messag(' ')
          call messag(' *** Error:  while reading central atom.')
          call messag('     '//co(:ii)//' matches no sites.-')
          call messag(' ')
          stop
      endif

c  check if 2 or more sites match core
      if (icore.ge.2) then
          ii = istrln(co)
          call messag(' ')
          call messag(' *** Error:  while reading central atom.')
          call messag('     '//co(:ii)//' matches more than one site.-')
          call messag(' ')
          stop
      endif

c  check if outfil or afname will overwrite input file
      if ((outfil.eq.fname).or.(afname.eq.fname).or.
     $            (refile.eq.fname)) then
          call messag(' *** Error: ')
          call messag('     One of your specified output files names '//
     $                'will overwrite the input file.')
          call messag('     This is not allowed.-')
          call messag(' ')
          stop
      endif

 300  continue
      return

c 666  continue
c      stop

c end subroutine atinpt
      end
