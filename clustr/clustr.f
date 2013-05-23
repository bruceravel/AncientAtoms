      subroutine clustr(iat,natx,ngeomx,nlogx,
     $         iabs,nsites,iperm,ipt,itot,ngeom,
     $         cell,trmtx,dmax,st,atlis,
     $         radii,reftmp,temp,lover,
     $         logic)
c=====================================================================
c  atom module 6:  expand cluster around central atom
c=====================================================================
c  this module consists of the following subroutines and functions:
c     clustr atheap cellex dist ovrlap subshl tetrot trans
c  this module calls function ref
c=====================================================================
c input integers:
c    iat,natx,ngeomx,nlogx: parameters set in calling program
c    iabs:   index of central atom in arrays dimensioned (iat)
c    nsites: number of unique positions
c    ipt:    (iat) number of positions of unique atom in unit cell  
c
c output integers:
c    itot:   total number of atoms in cluster
c    ngeom:  (ngeomx) single bounce geometry of each atom in cluster
c
c input reals:
c    cell:   (6) array of a,b,c,alpha,beta,gamma
c    trmtx:  (3,3) transformation matrix between cell axis and cartesian 
c            coordinates, from subroutine metric
c    dmax:   max cluster size in angstroms
c    st:     (iat,192,3) array of all atomic coordinates in unit cell, 
c            output of module crystl
c    radii:  (natx) workspace
c    reftmp: (natx) workspace
c    temp    (natx,8) workspace
c
c output reals:
c    atlis: (natx, 8) array of all atoms in cluster
c           1-3 -> pos. in cell     4 -> dist to origin
c              5 -> atom type     6-8 -> cartesian coords
c
c input logicals
c    lover:  (natx) workspace
c    logic : array of flags, see arrays.map
c---------------------------------------------------------------------
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

c      parameter (iat=50, natx=800, ngeomx=800)
      parameter(m=8)

      dimension ipt(iat), ngeom(ngeomx)
      dimension cell(6), trmtx(3,3)
      dimension st(iat,192,3), atlis(natx,8)
      dimension radii(natx),reftmp(natx),temp(natx,m)
      logical   logic(nlogx), lover(natx)

4000  format(a6)

c------------------------------------------------------------
c  expand the cluster, sort it by radial distance, check for overlap

      if (logic(25)) call messag('  calling cellex...')
      call cellex(iat,natx,iabs,nsites,ipt,st,cell,dmax,
     $            trmtx,itot,atlis)

c expand overfull cell, but not right now...
c      call cellex(iat,natx,iabs,nsites,iptful,fullcl,cell,dmax,
c     $            trmtx,itot,atlis)

c  load distances to atoms, (mode=1 in atheap means sort by distance)
      mode   = 1
      ifirst = 1
      do 50 i=1,itot
        radii(i) = atlis(i,4)
 50   continue 
      if (logic(25)) call messag('  calling atheap...')
      call atheap(natx,mode,ifirst,itot,atlis,radii,reftmp,temp)

c  translate all the atoms so that the first entry in atlis is at (0,0,0)
      xcent = atlis(1,6)
      ycent = atlis(1,7)
      zcent = atlis(1,8)
      do 70 i=1,itot
        atlis(i,6) = atlis(i,6) - xcent
        atlis(i,7) = atlis(i,7) - ycent
        atlis(i,8) = atlis(i,8) - zcent
 70   continue 

c  flag down and remove overlapping atoms
      if (logic(25)) call messag('  calling ovrlap...')
      call ovrlap(natx,itot,atlis,lover)
c  sort by subshells
      if (logic(25)) call messag('  calling subshl...')
      call subshl(natx,ngeomx,itot,atlis,ngeom,logic(6),
     $            radii,reftmp,temp)

c  permute tetragonal crystal back to original setting
      if (iperm.eq.22) call tetrot(natx, itot, cell, atlis)

      return
c  end of module clustr
      end
