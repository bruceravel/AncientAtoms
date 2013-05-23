      subroutine unit(idebug, ntit, title, vaxflg, ierr)
c=====================================================================
c  atom module 4:  construct unit.dat and/or p1.inp
c=====================================================================
c  this module consists of the following subroutines and functions:
c     unit p1out realcl unitcl writcl
c=====================================================================
c  n.b  this module is a skeleton right now -- eventually it will be
c       used for distortions
c=====================================================================
c * denotes output
c
c passed via unit.h:
c  * iptful: (iat) number of positions of unique atom in overfull unit cell
c  * fullcl: (iat,192,3) fractional coords of all atoms in overfull unit cell
c
c  idebug: an integer denoting the debug level, this is interpreted
c          into an array of binary bits which are used as logical flags.
c          multiple debugging features can be enables by specifying a
c          sum of bits
c     0:  disable all debuging function
c     1:  enable positional run-time messages (2**0)
c     2:  enable writing p1.inp               (2**1)
c     3:  enable writing unit.inp             (2**2)
c  title*72: (ntitx) array of title lines
c  vaxflg:   set to true for opening files on a vax
c  ierr:   output error code (0=no problem, 1=info, 2=warning, 3=error)
c----------------------------------------------------------------------

      implicit integer(i-n)
      implicit real(a-h,o-z)

      include 'atparm.h'
      include 'crystl.h'
      include 'exafs.h'
      include 'unit.h'
      include 'version.h'

      integer      idebug, idbg(0:ndbgx)
      character*72 title(ntitx)
      character*78 module, string
      logical      vaxflg

      call dbglvl(idebug,idbg)
c       print*,'idebug=',idebug,' idbg=(',idbg,')'
      module = 'unit'
      ierr = 0

c------------------------------------------------------------
c  if space group P1 input file desired, write it out
      if (idbg(1).gt.0) then
          string = 'calling p1out'
          if (idbg(0).gt.0) call positn(module, string)
          call p1out(iat,ndopx,ntitx,nsites,ntit,iabs,dopant,edge,
     $               cell,rmax,st,ipt,title,vrsion,vaxflg)
      endif

c------------------------------------------------------------
c  if unit cell file desired, calculate overfull cell and write
c  out unit.dat
      if (idbg(2).gt.0) then
          string = 'calling unitcl'
          if (idbg(0).gt.0) call positn(module, string)
          call unitcl(iat,nsites,st,ipt,fullcl,iptful)

          string = 'calling realcl'
          if (idbg(0).gt.0) call positn(module, string)
          call realcl(iat,nsites,iptful,cell,fullcl)

          string = 'calling writcl'
          if (idbg(0).gt.0) call positn(module, string)
          call writcl(iat,ndopx,nsites,dopant,tag,iptful,fullcl,
     $                vrsion, vaxflg)
      endif

      return
c  end of module unit
      end
