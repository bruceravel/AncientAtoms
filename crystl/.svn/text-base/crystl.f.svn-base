      subroutine crystl(idebug, ierr)
c=====================================================================
c  atoms module 2:  decode space group, determine unit cell contents
c=====================================================================
c  this module consists of the following subroutines and functions:
c     crystl basfil basort equipt genmul ibravl metric multip syschk
c  and requires
c     case dbglvl istrln messag positn upper
c
c  most variables are passed in common in crystl.h
c  idebug: an integer denoting the debug level, this is interpreted
c          into an array of binary bits which are used as logical flags.
c          multiple debugging features can be enables by specifying a
c          sum of bits
c     0:  disable all debuging function
c     1:  enable positional run-time messages
c  ierr:   output error code (0=no problem, 1=info, 2=warning, 3=error)
c------------------------------------------------------------------------
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      parameter(zero=0.e0, one=1.e0)
      include 'atparm.h'
      include 'crystl.h'
      character*1  csymbr
      logical      lbasis
      dimension    ts(3,24), fs(3,3,24)
      integer      idebug, idbg(0:ndbgx)
      character*78 module, string

      lbasis=.false.
      if (ibasis.gt.0) lbasis=.true.
      call dbglvl(idebug,idbg)
c       print*,'idebug=',idebug,' idbg=(',idbg,')'
      module = 'crystl'
      ierr = 0

c  initialize centrosymmetry flag
c  isymce=1 marks a centrosymmetric atom, this is determined in
c  equipt and used elsewhere
      isymce = 0
      ns     = 0
c  initialize the transformation matrices
      do 30 i=1,3
        do 20 j=1,24
          ts(i,j) = zero
          do 10 k=1,3
            fs(i,k,j) = zero
 10       continue
 20     continue
 30   continue
      fs(1,1,1) = one
      fs(2,2,1) = one
      fs(3,3,1) = one

c------------------------------------------------------------
c  reorganize basis list so absorber is at symmetry point in cell
      nat = iatom
      if (lbasis.and.(iabs.ne.1)) then
          string = 'sorting basis'
          if (idbg(0).gt.0) call positn(module, string)
          call basort(iat,ndopx,ibasis,iabs,x,y,z,dopant,tag)
          nat = ibasis
      endif

c------------------------------------------------------------
c  permute atoms to standard setting
      if (iperm.gt.1) call fperm(iat, nat, iperm, cell, x, y, z)

c------------------------------------------------------------
c  calculate the transformation matrix
      string = 'getting matrix'
      if (idbg(0).gt.0) call positn(module, string)
      call metric(cell, trmtx)

c------------------------------------------------------------
c  call the routine that does all the group theory to find
c  individual atom positions from the symmetry properties of
c  the space group and the positions of the unique atoms within
c  the primitive cell
      string = 'entering equipt'
      if (idbg(0).gt.0) call positn(module, string)
      call equipt(isyst, isymce, csymbr, ns, ts, fs, spcgrp)

c------------------------------------------------------------
c  two independent checks on system of crystal
      string = 'checking system'
      if (idbg(0).gt.0) call positn(module, string)
      syserr = .false.
      if (isyst.ne.isystm) then
          call syschk(isyst, isystm, spcgrp, sysmes)
          syserr = .true.
          ierr = 1
      endif

c------------------------------------------------------------
c  get number of bravais lattice for use in multip
      ibravl=ibrav(csymbr)

c------------------------------------------------------------
c  get igen, the multiplicity of the general position then get
c  multiplicities of atom positions,
c  for bases iatom = 1, want to run multip only on first atom in basis
      string = 'calling genmul'
      if (idbg(0).gt.0) call positn(module, string)
      call genmul(ns,isymce,ibravl,igen)
      string = 'callinp multip'
      if (idbg(0).gt.0) call positn(module, string)
      call multip(iatom, ibravl, x, y, z, tag, fs, ts, isymce,
     $            ns, igen, cell, st, ipt, imult, ierr)

c------------------------------------------------------------
c  fill in basis at each point
      if (lbasis) then
          string = 'filling basis'
          if (idbg(0).gt.0) call positn(module, string)
          call basfil(iat, ibasis, x, y, z, ipt, st)
      endif

c------------------------------------------------------------
c  permute atoms back to non-standard setting, but not tetragonal
      if ((iperm.ne.22).and.(iperm.gt.1))
     $            call bperm(iat, nat, iperm, cell, ipt, st)

      return
c  end of module crystl
      end
