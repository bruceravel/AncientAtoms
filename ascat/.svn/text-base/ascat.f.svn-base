      subroutine ascat(iat, ntitx, ndopx, neptx, nlogx,
     $             nsites, ipt, idop, ntit, nepts, nrefl, nnoan,
     $             edge, core, dopant, title, afname, refile, vrsn,
     $             tag, noantg,
     $             st, usqr, cell, percnt, qvect, egr, fcore, f0,
     $             anot, logic, vaxflg)
c=====================================================================
c  atom module 5:  construct a0 for a given reflection
c=====================================================================
c  this module consists of the following subroutines and functions:
c     a0 makea0 spcing wrta0
c  this module also requires the cromann and fcal packages
c  this module requires function volume
c=====================================================================
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)
c---------------------------------------------------------------------
c calling the form factor a0 is an historical thing.  Hans Stragier
c adopted that convention in his thesis and it has been carried
c around the Stern and Sorensen labs ever since.
c---------------------------------------------------------------------
c makea0 calls getf0 (tabulation of Cromer-Mann aikman coefficients)
c and fcal (tabulation of calculations of f' anf f" for interpolation)
c---------------------------------------------------------------------
c integers:
c   iatx, ntitx, ndopx, neptx, nlogx: parameters set in main program
c   nsites: number of atoms in atom or basis list
c   ipt:    (iat) number of positions of unique atom in unit cell
c   idop:   (iat) number of species at each site, 1=no dopants, 2+=dopants
c   ntit:   number of title lines
c   nepts:  number of energy points in grid for a0
c   nnoan:  number of tags for which anomalous correction should be neglected
c
c characters:
c   edge*2:    excitation edge measured, k or l3
c   core*10:   atomic symbol of resonant atom
c   noantg*10: (iat) tags of atoms for which an. corr. should be neglected
c   dopant*2:  (iat,ndopx) matrix with all host and dopant atomic symbols
c   title*72:  (ntitx) user supplied title lines
c   afname*72: output file name for structure factor as a function of energy
c   refile*72: output file name for reflection amplitudes
c   vrsn*5:    atoms version number
c
c reals:
c   st:     (iat,192,3) fractional coords of all atoms in unit cell
c   usqr:   (iat) crystallographic thermal factors
c   cell:   (6) a,b,c,alpha,beta,gamma of cell
c   percnt: (iat,ndopx) percent occupancies of hosts and dopants
c   qvect:  (3) q vector of reflection
c   egr:    energy grid spacing for a0
c   fcore:  (neptx) f" for resonant atom (output)
c   f0:     (iat,ndopx) work array for regular part of scattering factor
c
c complex (and output):
c   anot:   (neptx) complex array containing a0
c
c logicals:
c   logic:  array of flags, see arrays.map
c---------------------------------------------------------------------
c      parameter (iat=50, ndopx=4, neptx=2**11, ntitx=9)
      parameter (eps=0.0001, iou=17)
      parameter (pi = 3.1415 92653 58979 32384 62643e0)

      character*2  edge, dopant(iat,ndopx)
      character*5  vrsn
      character*10 core, whose*8, tag(iat), noantg(iat)
      character*72 title(ntitx), afname, refile
      dimension    ipt(iat), idop(iat), nrefl(3)
      dimension    percnt(iat,ndopx), qvect(3)
      dimension    st(iat,192,3), usqr(iat)
      complex      fcore(neptx)
      dimension    f0(iat,ndopx)
      logical      logic(nlogx), vaxflg
      complex      anot(neptx)

c  which f'/f" data set is used?
      whose = 'Sasaki'
c  ignore crystallographic dwf for now
      do 10 i=1,iat
        usqr(i) = 0
 10   continue
c  set defaults for nepts (50 above and below) and egrid (20 eV)
      if (nepts.eq.0) nepts = 50
      if (egr.lt.eps) egr   = 0.02

      if (logic(10)) then
          ii = istrln(whose)
          call messag('Computing anomalous scattering factor using '//
     $                whose(:ii)//' data.')
          if (logic(24)) call messag('  calling makea0...')
          call makea0(iat,ndopx,neptx,
     $                nsites, ipt, idop, nepts, nnoan,
     $                edge, core, dopant, tag, noantg, egr, ecusp,
     $                st, usqr, cell, percnt, qvect, anot, fcore, f0,
     $                logic(17))

          if (logic(24)) call messag('  calling wrta0...')
          call wrta0(iat, ntitx, neptx,
     $                nepts, ntit, nnoan ,qvect, egr, ecusp, anot,
     $                fcore, afname, noantg, edge, core, title,
     $                vrsn, vaxflg)
      endif

      if (logic(11)) then
          call refls(iat, ndopx, ntitx, nsites, nrefl, ntit,
     $            refile, core, dopant, tag, title, vaxflg, st, ipt)
      endif

      return
c  end of module ascat
      end
