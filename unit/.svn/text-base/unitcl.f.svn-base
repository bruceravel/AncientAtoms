      subroutine unitcl(iat,iatom,st,ipt,fullcl,iptful)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c--------------------------------------------------------------------
c  create the overfull unit cell.  this is defined as the true unit 
c  cell plus additional atoms on all the walls.  thus any atom that 
c  is fully or fractionally within the central cell will have its 
c  coordinates in fullcl.
c--------------------------------------------------------------------
c  input:
c    iatom:  number of unique atoms
c    st:     fractional coords of each atom in the true unit cell
c    ipt:    number of occurrences of each unique atom in the true cell
c  output:
c    fullcl: fractional coords of each atom in the overfull unit cell
c    iptful: number of occurences of each unique atom in the overfull
c            cell
c--------------------------------------------------------------------
c  because 0 is a special position the largest number of atoms on 
c  true unit cell walls is 96 -- doubled is 192, thus ipt and fullcl
c  can have the same dimensions
c--------------------------------------------------------------------

c      parameter (iat=50)
      parameter(eps=0.001)

      dimension st(iat,192,3),fullcl(iat,192,3)
      dimension ipt(iat),iptful(iat)

c--------------------------------------------------------------------
c  copy true cell contents into overfull cell arrays
      do 30 i=1,iatom
        do 20 j=1,ipt(i)
          do 10 k=1,3
            fullcl(i,j,k) = st(i,j,k)
10        continue
20      continue
        iptful(i) = ipt(i)
30    continue

c--------------------------------------------------------------------
c  search for atoms on a cell edge, wall, or corner.  increment iptful 
c  for that atom type and put an atom on the opposite cell edge, wall
c  or corner.
c  need to check each position for which wall it is at, i.e. i want
c  to translate atoms at coordinate 0 to coordinate 1 and those at 
c  coordinate 1 back to 0.  ix/y/z =0 when the atom is not near a wall.
c  ix/y/z =1 when atom is at 0 because i want to add 1 to the coordinate
c  ix/y/z =-1 when atom is at 1 because i want to sub 1 from the coord.
c  up to seven new atoms can be generated from a single position for 
c  the overfull cell
c                           ----------
c  so 0 means don't translate (not at wall), 1 means translate from 0 
c  wall to 1 wall, -1 means translate from 1 wall to 0 wall
c                           ----------
      do 50 i=1,iatom
        do 40 j=1,ipt(i)
          ix = 0
          iy = 0
          iz = 0
c                                              load ix/y/z for each atom
          if ( abs(st(i,j,1))  .lt.eps ) ix =  1
          if ( abs(st(i,j,1)-1).lt.eps ) ix = -1
          if ( abs(st(i,j,2))  .lt.eps ) iy =  1
          if ( abs(st(i,j,2)-1).lt.eps ) iy = -1
          if ( abs(st(i,j,3))  .lt.eps ) iz =  1
          if ( abs(st(i,j,3)-1).lt.eps ) iz = -1
c                                              x at wall
          if (ix.ne.0) then                
                  iptful(i)      = iptful(i) + 1
                  jj             = iptful(i)
c                  print*,'unitcl: jj=',jj
                  fullcl(i,jj,1) = st(i,j,1) + ix
                  fullcl(i,jj,2) = st(i,j,2)
                  fullcl(i,jj,3) = st(i,j,3)
c                                              x and y at wall
                  if (iy.ne.0) then
                          iptful(i)      = iptful(i) + 1
                          jj             = iptful(i)
c                          print*,'unitcl: jj=',jj
                          fullcl(i,jj,1) = st(i,j,1) + ix
                          fullcl(i,jj,2) = st(i,j,2) + iy
                          fullcl(i,jj,3) = st(i,j,3)
c                                              x, y and z at wall
                          if (iz.ne.0) then
                                  iptful(i)      = iptful(i) + 1
                                  jj             = iptful(i)
c                                  print*,'unitcl: jj=',jj
                                  fullcl(i,jj,1) = st(i,j,1) + ix
                                  fullcl(i,jj,2) = st(i,j,2) + iy
                                  fullcl(i,jj,3) = st(i,j,3) + iz
                          endif
                  endif
c                                              x and z at wall
                  if (iz.ne.0) then
                          iptful(i)      = iptful(i) + 1
                          jj             = iptful(i)
c                          print*,'unitcl: jj=',jj
                          fullcl(i,jj,1) = st(i,j,1) + ix
                          fullcl(i,jj,2) = st(i,j,2) 
                          fullcl(i,jj,3) = st(i,j,3) + iz
                  endif
          endif
c                                              y at wall
          if (iy.ne.0) then
                  iptful(i)      = iptful(i) + 1
                  jj             = iptful(i)
c                  print*,'unitcl: jj=',jj
                  fullcl(i,jj,1) = st(i,j,1)
                  fullcl(i,jj,2) = st(i,j,2) + iy
                  fullcl(i,jj,3) = st(i,j,3)
c                                              y and z at wall
                  if (iz.ne.0) then
                          iptful(i)      = iptful(i) + 1
                          jj             = iptful(i)
c                          print*,'unitcl: jj=',jj
                          fullcl(i,jj,1) = st(i,j,1)
                          fullcl(i,jj,2) = st(i,j,2) + iy
                          fullcl(i,jj,3) = st(i,j,3) + iz
                   endif
          endif
c                                              z at wall
          if (iz.ne.0) then
                  iptful(i)      = iptful(i) + 1
                  jj             = iptful(i)
c                  print*,'unitcl: jj=',jj
                  fullcl(i,jj,1) = st(i,j,1)
                  fullcl(i,jj,2) = st(i,j,2) 
                  fullcl(i,jj,3) = st(i,j,3) + iz
          endif
40      continue
50    continue

      return
c  end subroutine unitcl
      end

