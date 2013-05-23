      subroutine atinit(iat,natx,ntitx,ndopx,ngeomx,neptx,nlogx,
     $            title,elemnt,tag,noantg,spcgrp,edge,core,
     $            dopant,outfil,afname,geodat,refile,
     $            iatom,ibasis,dmax,ispa,iperm,idop,ngeom,iedge,iabs,
     $            ipt,iptful,isyst,nepts,nrefl,nnoan,
     $            x,y,z,cell,st,fullcl,atlis,
     $            percnt,gasses,qvect,egr,anot,
     $            logic,stdout)

c      $            lbasis,lindex,lfluor,ldwarf,lmm,lunit,lp1,
c      $            lself,li0,lgeom,ldafs,lf,lmod,lfeff,lfex,lgex,
c      $            lcrys, lmcm, lun, la0, lclus, lout, lrefl)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------
c  title   : user supplied comment line
c  elemnt  : character array of atom types
c  tag     : character array of site tags
c  spcgrp  : chararcter string with hermann-maguin space group designation
c  edge    : absorption edge of core atom, k or l3
c  x,y,z   : positions of atoms in cell-axis coordinates
c  cell    : lattice constants, a,b,c,alpha,beta & gamma
c  iatom   : number of unique atoms in cell
c  ibasis  : number of atoms in basis
c  dmax    : maximum distance in cluster
c  ispa    : number between 1 and 230, index of s.g. in int'l tables
c  iperm   : index of permutation matrix for non-standard setting
c  outfil  : output file name
c  refile  : reflection amplitudes file name
c  st:     : arrays containing positions in unit cell
c  isystm  : bravais lattice type
c  dopant  : dopant element symbol
c  percnt  : percent of dopant
c  gasses  : percent of argon, krypton, nitrogen in i0 chamber
c  logic   : logic flags, see arrays map
c-------------------------------------------------------
c  initialize everything there is to initialize
c-------------------------------------------------------
      parameter (zero=0.000000000, one=1.000000000)
c      parameter (iat=50,natx=800,neptx=2**11)
c      parameter (ntitx=9,ndopx=4,ngeomx=800)
      character*2  elemnt(iat),edge,dopant(iat,ndopx)
      character*10 spcgrp,tag(iat),core,geodat,noantg(iat)
      character*72 title(ntitx),outfil, afname, refile
      complex      anot(neptx)
      dimension    x(iat), y(iat), z(iat), cell(6), qvect(3)
      dimension    st(iat,192,3),
     $             fullcl(iat,192,3), atlis(natx,8)
      dimension    ipt(iat), iptful(iat), idop(iat),
     $             ngeom(ngeomx), nrefl(3)
      dimension    percnt(iat,ndopx), gasses(3)
      logical      logic(nlogx), stdout

c--------- characters -------------------------------------------------
      do 10 i=1,ntitx
        title(i) = ' '
10    continue
      spcgrp = ' '
      edge   = ' '
      core   = ' '
      outfil = 'feff.inp'
      call lower(outfil)
      afname = ' '
      geodat = 'geom.dat'
      call lower(geodat)
      refile = 'reflect.dat'
      call lower(refile)

c--------- logicals ---------------------------------------------------
      do 15 i=1,nlogx
        logic(i) = .false.
 15   continue
c     7: feff   19: feff8   28: mcm
      logic(7)  = .true.
      logic(19) = .true.
      logic(28) = .true.
c      stdout = .true.

c--------- reals ------------------------------------------------------
      dmax   = 5.
      do 17 i=1,3
        gasses(i) = zero
 17   continue
      egr    = zero

c--------- integers ---------------------------------------------------
      ispa   = 0
      iperm  = 1
      iatom  = 0
      ibasis = 0
      isyst  = 0
      iabs   = 0
      iedge  = 0
      nepts  = 0
      nnoan  = 0
      do 20 i=1,3
        nrefl(i) = 0
 20   continue

c--------- arrays ------------------------------------------------------
      do 100 i=1,iat
        elemnt(i) = ' '
        tag(i)    = ' '
        noantg(i) = ' '
        x(i)      = zero
        y(i)      = zero
        z(i)      = zero
        ipt(i)    = 0
        iptful(i) = 0
        idop(i)   = 1
        do 90 j=1,ndopx
          dopant(i,j) = ' '
          percnt(i,j) = zero
 90     continue
100   continue
      do 110 i=1,3
        cell(i)   = zero
        cell(i+3) = 90.
        qvect(i)  = zero
110   continue
      do 120 i=1,ngeomx
        ngeom(i) = 0
 120  continue

c  unit transformation for trivial position k=1
      do 180 i=1,3
        do 170 j=1,iat
          do 160 k=1,192
            st(j,k,i)     = zero
            fullcl(j,k,i) = zero
160       continue
170     continue
180   continue

      do 200 i=1,natx
        do 190 j=1,8
          atlis(i,j) = zero
190     continue
200   continue

      do 210 i=1,neptx
        anot(i) = cmplx(zero,zero)
 210  continue

      return
c end subroutine atinit
      end
