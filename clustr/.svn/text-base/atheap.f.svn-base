      subroutine atheap(natx,mode,nfirst,nlast,atlis,
     $            refer,reftmp,temp)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------
c  heapsort adapted from numerical recipes.  sort rows of atlis 
c  comparison to be done with refer, pertinent rows of atlis
c  first transfered into temp.  refer and temp are sorted together.
c  all the pesky little do loops are for transferring rows
c  of temp into toss.
c
c  requires function ref
c-------------------------------------------------------------------
c  natx:   dimension parameter from calling program
c  mode:   contents of refer, mode=1: sort by radial distance
c                             mode=2: sort by hash of cart. coords.
c  nfirst: first element of array to sort
c  nlast:  last element to sort
c  atlis(natx,8):  array of (position,distance,atom type,cart. coords)
c                  sorted on output
c  refer(natx): array of sorting criteria values, sorted on output
c  reftmp(natx), temp(natx,8): work space
c-------------------------------------------------------------------
c  atlis(natx,8): 1-3 -> pos. in cell    4 -> dist to origin
c                   5 -> atom type     6,8 -> cartesian coords
c-------------------------------------------------------------------
      parameter (m=8)
c      parameter (natx=800)
      dimension atlis(natx,m),toss(m)
      dimension refer(natx),reftmp(natx),temp(natx,m)

c  switch atoms to be sorted into temp
      natom = nlast-nfirst+1
      do 40 i=1,natom
        reftmp(i) = refer(nfirst+i-1)
        do 20 j=1,m
          temp(i,j) = atlis(nfirst+i-1,j)
 20     continue 
 40   continue 

      l  = natom/2+1
      ir = natom
 110  continue
        if (l.gt.1) then
            l = l-1
            do 120 index=1,m
              toss(index)=temp(l,index)
 120        continue
            rref = reftmp(l)
        else
            do 130 index=1,m  
              toss(index) = temp(ir,index)
 130        continue
            rref = reftmp(ir)
            do 140 index=1,m
              temp(ir,index)=temp(1,index)
 140        continue
            reftmp(ir)=reftmp(1)
            ir=ir-1
            if(ir.eq.1)then
               do 150 index=1,m
                 temp(1,index)=toss(index)
 150           continue
               reftmp(1)=rref
c              sort is finished
               goto 300
            endif
        endif
        i=l
        j=l+l

 160    if (j.le.ir) then
          if (j.lt.ir) then
            if (reftmp(j).lt.reftmp(j+1)) j=j+1
          endif
c                                             * * choose sorting mode * *
          
          reftos = toss(4)
          if (mode.eq.1) then
              reftos = toss(4)
          elseif (mode.eq.2) then
              reftos = ref( toss(6), toss(7), toss(8) )
          endif
          if(reftos.lt.reftmp(j))then
             do 170 index=1,m
               temp(i,index)=temp(j,index)
 170         continue
             reftmp(i)=reftmp(j)
             i=j
             j=j+j
          else
             j=ir+1
          endif
        goto 160     
        endif
        do 180 index=1,m
          temp(i,index)=toss(index)
 180    continue
        reftmp(i)=rref
      goto 110
    
c  switch atoms back into atlis
 300  continue 
      do 320 i=1,natom
        refer(nfirst+i-1) = reftmp(i)
        do 310 j=1,m
          atlis(nfirst+i-1,j) = temp(i,j)
 310    continue 
 320  continue 

      return 
c end subroutine atheap
      end
