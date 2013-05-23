      subroutine mcm(idebug)
      implicit integer(i-n)
      implicit real(a-h,o-z)
c=====================================================================
c  Atoms Module 3:  perform various calculations using mcmaster tables
c=====================================================================
c  this module consists of the following subroutines and functions:
c     mcm abslen i0 mcmast slfabs
c  this module also requires: mucal volume positn messag lower dbglvl
c                             polyft case
c=====================================================================
c  input describing the unit cell comes via crystl.h
c  input controlling exafs function comes via exafs.h
c
c  idebug: an integer denoting the debug level, this is interpreted
c          into an array of binary bits which are used as logical flags.
c          multiple debugging features can be enables by specifying a
c          sum of bits
c     0      :  disable all debuging function
c     1 (+1) :  enable positional run-time messages
c     2 (+2) :  write mcmaster correction diagnostic file
c     3 (+4) :  write self-absorption correction diagnostic file
c     4 (+8) :  write i0 correction diagnostic file
c------------------------------------------------------------------------

      include 'atparm.h'
      include 'crystl.h'
      include 'exafs.h'

      parameter(epsi = 1.e-3)
      integer idebug, idbg(0:ndbgx)
      logical lmcm, lself, li0
      character*78 module, string, messg

      call dbglvl(idebug,idbg)
c       print*,'idebug=',idebug,' idbg=(',idbg,')'
      lmcm  = .false.
      if (idbg(1).gt.0) lmcm  = .true.
      lself = .false.
      if (idbg(2).gt.0) lself = .true.
      li0   =.false.
      if (idbg(3).gt.0) li0   = .true.
      module = 'mcm'

 400  format(' *** Warning: ', a, ' set to a negative number.')
 410  format(' *** Warning: ', a, ' set to larger than 1.')
c  check nitrogen
      if (gasses(1).lt.0.) then
          write(messg,400)'argon'
          call messag(messg)
          call messag('     The value was reset to 0.-')
          gasses(1) = 0
          iexerr = 2
      endif
      if (gasses(1).gt.1.) then
          write(messg,410)'argon'
          call messag(messg)
          call messag('     The value was reset to 1.-')
          gasses(1) = 1
          iexerr = 2
      endif
c  check argon
      if (gasses(2).lt.0.) then
          write(messg,400)'krypton'
          call messag(messg)
          call messag('     The value was reset to 0.-')
          gasses(2) = 0
          iexerr = 2
      endif
      if (gasses(2).gt.1.) then
          write(messg,410)'krypton'
          call messag(messg)
          call messag('     The value was reset to 1.-')
          gasses(2) = 1
          iexerr = 2
      endif
c  check krypton
      if (gasses(3).lt.0.) then
          write(messg,400)'nitrgen'
          call messag(messg)
          call messag('     The value was reset to 0.-')
          gasses(3) = 0
          iexerr = 2
      endif
      if (gasses(3).gt.1.) then
          write(messg,410)'nitrogen'
          call messag(messg)
          call messag('     The value was reset to 1.-')
          gasses(3) = 1
          iexerr = 2
      endif
c  check sum of gasses
      lfluo = .false.
      sum = gasses(1) + gasses(2) + gasses(3)
      if (sum.gt.1.+epsi) then
          call messag(' *** Warning: the sum off gasses exceeds 1.')
          call messag('     Turning off fluorescence corrections.-')
          iexerr = 2
      elseif (sum.gt.epsi) then
          lfluo = .true.
      endif

c------------------------------------------------------------
c  calculate absorption and density of the crystal
c  calculate mcmaster correction for central atom
c  calculate self absorption correction for fluorescence
c  calculate i0 correction for fluorescence

      if ( (iedge.gt.1) .and. (is2z(core).lt.30) ) then
          call messag(' ')
          call messag(' *** Warning: McMaster calculations are '//
     $                'unreliable for L edges of Z<30.-')
          call messag(' ')
          iexerr = 2
      endif

      string = 'computing absorption length'
      if (idbg(0).gt.0) call positn(module, string)
      v = volume(cell)
      call abslen(iat, ndopx, nsites, ipt, iedge, idop, core, dopant,
     $            percnt, v, amu, delmu, spgrav)

      string = 'computing McMaster correction'
      if (idbg(0).gt.0) call positn(module, string)
      call mcmast(core,iedge,sigmm,qrtmm,lmcm)

      if (lfluo) then
          string = 'computing self absorption correction'
          if (idbg(0).gt.0) call positn(module, string)
          call slfabs(iat, ndopx, nsites, ipt, iedge, idop, core,
     $                dopant, percnt, lself,
     $                ampslf, sigslf, qrtslf, xmuf, xmub)

          string = 'computing I0 correction'
          if (idbg(0).gt.0) call positn(module, string)
          call i0(gasses, core, iedge, li0, sigi0, qrti0)
      endif

      exafs(1)  = amu
      exafs(2)  = delmu
      exafs(3)  = spgrav
      exafs(4)  = sigmm
      exafs(5)  = qrtmm
      exafs(6)  = ampslf
      exafs(7)  = sigslf
      exafs(8)  = qrtslf
      exafs(9)  = sigi0
      exafs(10) = qrti0
      exafs(11) = xmuf
      exafs(12) = xmub
      exafs(13) = delmu*v

      return
c  end of module mcm
      end
