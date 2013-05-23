      subroutine subshl(natx,ngeomx,itot,atlis,ngeom,lgeom,
     $            refer,reftmp,temp)

c==================================================================
c integers:
c    natx, ngeomx:  dimension parameters from calling program
c    itot:  number of atoms in cluster
c
c reals:
c    atlis: (natx, 8) array of all atoms in cluster
c           1-3 -> pos. in cell     4 -> dist to origin
c              5 -> atom type     6-8 -> cartesian coords
c    
c output:
c    ngeom:  (ngeomx) one bounce flags for geom.dat; integer array
c
c work space, real arrays:
c    refer(natx),  reftmp(natx), temp(natx,8)
c==================================================================
c  this routine sorts equidistant atoms by a hash value
c  of the cartesian coordinates.  thus atoms at (0,0,5) etc will
c  be grouped together and different from those at (0,3,-4) etc.
c  also fills ngeom for use as the one bounce flag in geom.dat
c==================================================================

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

c      parameter (natx=800, ngeomx=800)
      parameter (eps=0.001e0, m=8)

      logical   lgeom
      dimension atlis(natx,8), ngeom(ngeomx)
      dimension refer(natx), reftmp(natx), temp(natx,m)

      mode   = 2
      ifirst = 0
      ig     = 0
      do 10 i=1,itot
        refer(i) = ref( atlis(i,6), atlis(i,7), atlis(i,8) )
c        print*,i,atlis(i,6), atlis(i,7), atlis(i,8), refer(i)
 10   continue 

      rlast  = atlis(1,4)
c  step through atom list, when distance changes sort the preceding block
c  then reset all the markers and continue
      do 20 i=2,itot
        if (abs(atlis(i,4)-rlast).gt.eps) then
            ilast     = i-1
            if ((ilast-ifirst).gt.2) then
                call atheap(natx,mode,ifirst,ilast,atlis,
     $                      refer,reftmp,temp)
            endif
            ifirst = i
            rlast  = atlis(i,4)
        endif
 20   continue 

c     sort the last group!
      if ((itot-ifirst).gt.2) call atheap(natx,mode,ifirst,itot,atlis,
     $            refer,reftmp,temp)

c  fill ngeom for use in geom.dat
      if (lgeom) then
          ifirst = 1
          ig     = 0
          do 50 i=2,itot
            if (abs(refer(i)-refer(i-1)).gt.eps)  then
                ig        = ig + 1
                ngeom(ig) = i - ifirst
                ifirst    = i
            endif
 50       continue
          ig        = ig +1 
          ngeom(ig) = itot-ifirst+1
      endif

      return
c  end subroutine subshl
      end
