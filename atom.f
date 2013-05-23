      program atoms
      implicit real(a-h,o-z)
      implicit integer(i-n)
      include 'version.h'
c----------------------------------------------------------------------
c             copyright 1998  Bruce Ravel
c        copyright 1993-1997  University of Washington
c             written by      Bruce Ravel
c                             Ceramics Division
c                             National Institute of Standards and Technology
c                             Bldg 223, Room 329
c                             Gaithersburg, MD 20899
c
c                  phone      (301) 975-5759
c                  fax        (301) 975-5334
c                  e-mail     bruce.ravel@nist.gov
c                          or ravel@phys.washington.edu
c           please use email for communication with the author
c
c  insert permissions notice here
c----------------------------------------------------------------------
c      brief description of the code:
c
c  atoms writes a list of atomic coordinates for any crystal given its
c  crystallographic information. the list will be sorted by radial
c  distance from an atom chosen as the central atom. atoms also estimates
c  the bulk absorption and density of the material and various corrections
c  to xafs data due to experimental effects.
c----------------------------------------------------------------------
c  comments blocks that follow:
c       glossary of variables
c       descriptions of runtime error codes
c       sample input file
c       version history
c
c  the code begins after the first occurance of this string:  %%%%
c  the first executable statement begins after its second occurance
c----------------------------------------------------------------------
c >>> glossary of variables
      include 'glossary'
c=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
c >>> runtime error codes
      include 'runtime'
c=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
c >>> sample input file:
      include 'sample'
c=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
c >>> history
      include 'history'
c----------------------------------------------------------------------
c%%%%
      include 'atparm.h'

c  stdout=.true. for reading from standard input and writing to standard
c  output
      logical stdout
      parameter (stdout=.false.)

c  expnd=.true. for reading expanded atoms list, false for normal atoms list
      logical expnd
      parameter (expnd=.false.)

c  vaxflg=.true. for compilation on a VMS machine, used for opening
c  files in such a way that VMS version numbers are used
      logical vaxflg
      parameter (vaxflg=.false.)

c  crystl.h contains all information relevant to a unit cell
      include 'crystl.h'
c  exafs.h contains information relevant to exafs calculations
      include 'exafs.h'
c  unit.h contains information over-full unit cell information
      include 'unit.h'

c  see glossary for descriptions of variables, see arrays.map for
c  description of logic, gasses, exafs
      parameter(m=8)
      character*2  elemnt(iat), test, vrsn*9
      character*10 inpgrp, noantg(iat)
      character*72 title(ntitx), outfil, afname, refile
      character*78 messg
      character*78 module, string
      logical      logic(nlogx)
      complex      anot(neptx)
      dimension    ngeom(ngeomx)
      dimension    atlis(natx,8), qvect(3)
      dimension    nrefl(3)


c  reserve some array space for use in modules cluster and output
      dimension    cltmp1(natx),  cltmp2(natx),  cltmp3(natx,m)
      logical      cltmp4(natx)
      character*10 optmp1(maxln)
      dimension    optmp2(maxln), optmp3(maxln), optmp4(maxln),
     $             optmp5(maxln), ioptp6(maxln), ioptp7(maxln)

c------------------------------------------------------------------------
c-- dafs stuff
c       dimension    f0(iat,ndopx), usqr(iat)
c       complex      fcore(neptx)
c c  from block data sasaki
c       parameter(nelem=92, ndatx=233)
c       common /fdat/ nfdata(nelem), engrd(nelem,ndatx),
c      $            fp(nelem,ndatx), fpp(nelem,ndatx)
c c  this is the crommer-mann common block
c       common /sk/ aa(1926)
c------------------------------------------------------------------------

c------------------------------------------------------------------------
c  formats needed at top level: version # (4000), line of '=' (4010),
c                 ... aN,(55-N)x ...
 4000 format('ATOMS ',a9,46x,'by Bruce Ravel')
 4010 format(75('='))

c%%%%
c---------------------------------------------------------------------
c  run time messages: lines of '=' and version number
c  set unit number of feff.inp
      if (.not.stdout) then
          ifeff = 2
          write(messg,4010)
          call messag(messg)
          write(messg,4000)vrsion
          call messag(messg)
          write(messg,4010)
          call messag(messg)
      else
c  this clause may trigger a compiler warning regarding there being no
c  possible path to this spot in the code.  this is an unavoidable
c  result of using F77 and no preprocessor and should pose no problem
c  during execution of the code.
          ifeff = 6
      endif

c---------------------------------------------------------------------
c  initialize some top-level variables
      vrsn  = vrsion
      test  = 'ab'
      logic(1) = .false.
      module = 'Atoms'

c 10    continue

c=====================================================================
c  atoms module 1:  initialize, read input file, error checking
c=====================================================================
c      print*,'beginning module 1'
      call readin(iat, natx, ntitx, ndopx, ngeomx, neptx, nlogx,
     $            ifeff, ntit, iatom, ibasis, iabs, isystm,
     $            iperm, ipt, iptful, idop, ngeom, iedge,
     $            nepts, nsites, nrefl, nnoan,
     $            vrsn, spcgrp, inpgrp, title, outfil, afname,
     $            refile, tag, noantg, edge, core, elemnt, dopant,
     $            x, y, z, cell, rmax, st, fullcl, atlis,
     $            percnt, gasses, qvect, egr, anot,
     $            logic, stdout, vaxflg, expnd)

      if (logic(1)) goto 99

c=====================================================================
c  atoms module 2:  decode space group, determine unit cell contents
c=====================================================================
      string = 'the crystl module'
      if (logic(20)) call positn(module, string)
      iunidb = 0
      if (logic(21)) icrydb = icrydb + 2**0
      call crystl(icrydb, iercry)

      if ( (.not.stdout).and.syserr ) then
          call messag(' ')
          call messag(' *** Warning ')
          do 20 i=1,nsysm
            call messag(sysmes(i))
 20       continue
          call messag('The calculation will be finished, but you '//
     $                'might want to edit your')
          call messag('crystallographic input data and try again.-')
          call messag(' ')
      endif

c=====================================================================
c  atoms module 3:  perform various calculations using mcmaster tables
c=====================================================================
      if (logic(28)) then
          string = 'the mcm module'
          if (logic(7)) then
              if (logic(20)) call positn(module, string)
              iexadb = 0
              if (logic(22)) iexadb = iexadb + 2**0
              if (logic(14)) iexadb = iexadb + 2**1
              if (logic(15)) iexadb = iexadb + 2**2
              if (logic(16)) iexadb = iexadb + 2**3
              call mcm(iexadb)
          endif
          if (.not.lfluo) logic(3)=.false.
      endif

c=====================================================================
c  atoms module 4:  construct unit.dat
c=====================================================================
      string = 'the unit module'
      if (logic(20)) call positn(module, string)
      iunidb = 0
      if (logic(23)) iunidb = iunidb + 2**0
      if (logic(8))  iunidb = iunidb + 2**1
      if (logic(9))  iunidb = iunidb + 2**2
      call unit(iunidb, ntit, title, vaxflg, ieruni)

c=====================================================================
c  atoms module 5:  construct structure factor for dafs applications
c=====================================================================
c       if (logic(10).or.logic(11)) then
c           string = 'the ascat module'
c           if (logic(20)) call positn(module, string)
c           call ascat(iat, ntitx, ndopx, neptx, nlogx,
c      $             nsites, ipt, idop, ntit, nepts, nrefl, nnoan,
c      $             edge, core, dopant, title, afname, refile, vrsn,
c      $             tag, noantg,
c      $             st, usqr, cell, percnt, qvect, egr, fcore, f0,
c      $             anot, logic, vaxflg)
c       endif

c=====================================================================
c  atoms module 6:  expand cluster around central atom
c=====================================================================
      if (logic(7)) then
          string = 'the clustr module'
          if (logic(20)) call positn(module, string)
          call clustr(iat,natx,ngeomx,nlogx,
     $            iabs,nsites,iperm,ipt,itot,ngeom,
     $            cell,trmtx,rmax,st,atlis,
     $            cltmp1,cltmp2,cltmp3,cltmp4,
     $            logic)

c=====================================================================
c  atoms module 7:  write feff.inp, geom.dat
c=====================================================================
          string = 'the output module'
          if (logic(20)) call positn(module, string)
          call output(iat, natx, ntitx, ndopx, ngeomx, maxln, nlogx,
     $            nexafs,
     $            ifeff, iabs, itot, ntit, idop, ngeom, imult,
     $            title, tag, edge, core, dopant, outfil, elemnt,
     $            vrsn, percnt, exafs, atlis,
     $            logic, vaxflg,
     $            optmp1,optmp2,optmp3,optmp4,optmp5,ioptp6,ioptp7)

c#####################################################################

c---------------------------------------------------------------------
c  run time message: output file and line of '='
          if ((logic(6)).and.(outfil.ne.'list')) then
              if (logic(13)) then
                  call messag('geom.dat overwritten.')
              else
                  call messag('geom.dat written.')
              endif
          endif
          if (.not.stdout) then
              ii = istrln(outfil)
              if (logic(12)) then
                  call messag(outfil(:ii)//' overwritten.')
              else
                  call messag('Output written to '//outfil(:ii))
              endif
          endif
      else
          call messag('not writing feff.inp.')
      endif
      if (.not.stdout) then
          write(messg,4010)
          call messag(messg)
      endif
c---------------------------------------------------------------------

c      if (.not.logic(1)) goto 10
99    continue
      stop
c  end main program atom
      end
