      subroutine realcl(iat,iatom,iptful,cell,fullcl)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision(a-h,o-z)
c---------------------------------------------------------------------
c  convert coordinates of atoms in overfull cell from fractional to 
c  angstroms along cell axis basis
c  do not transform yet!
c---------------------------------------------------------------------
c  input:
c    iatom:  number of unique atoms
c    iptful: number of each unique atom in overfull cell
c    cell:   array containing lattice params a,b,c,alpha,beta,gamma
c  i/o:
c    fullcl: fractional positions of atoms in overfull cell on input,
c            cartesian positions on output
c---------------------------------------------------------------------
c      parameter (iat=50)

      dimension iptful(iat),cell(6),fullcl(iat,192,3)

      do 40 i=1,iatom
        do 30 j=1,iptful(i)
c                                        read in position of an atom 
c                                        and multiply by axis lengths
          do 10 icoord=1,3
            fullcl(i,j,icoord)  = fullcl(i,j,icoord)  * cell(icoord)
10        continue
30      continue
40    continue

      return
c  end subroutine realcl
      end

