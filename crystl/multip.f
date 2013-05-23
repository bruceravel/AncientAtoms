      subroutine multip (iatom,ibravl,x,y,z,tag,fs,ts,isymce,ns,
     $                   igen,cell,st,ipt,imult,ierr)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)
c-----------------------------------------------------------------
c     finds multiplicities of unique positions in primitive cell
c-----------------------------------------------------------------
c   iat = max # of unique atom types
c   ia  = index used below for counting through all atom types in problem
c
c  input:
c    iat:    dimension parameter set in calling program
c    iatom:  # of atom types in problem (iatom<=iat)
c    ibravl: number denoting bravais lattice type
c    x,y,z:  (iat) fractional coordinates of unique atom positions
c    tag*10: (iat) tags of each unique site
c    fs:     (3,3,24) rotational symmetries from equipt
c    ts:     (3,24) fractional coords of simply multiple positions from equipt
c    isymce: flag for centrosymmetry, (1=yes 0=no)
c    ns:     simple multiplicity (ng<=24)
c    igen:   general multiplicity of cell
c    cell:   (6) array containing a,b,c,alpha,beta, and gamma
c  output:
c    st(iat,192,3):  first index flags atom type
c                    second index flags each occurence of that atom type
c                      there are never more that 192 of an atom type in a cell
c                    third index contains x,y,z (in cell axis basis) and dist.
c    ipt:    (iat) array telling how many of each unique position in cell
c    ierr:   error code, 0 if ok, 2 if decoding problem
c------------------------------------------------------------------------
c      parameter (iat=50)
      include 'atparm.h'
      parameter (zero  = 0.e0, one   = 1.e0, eps = 0.0001e0)
      parameter (three = 3.e0, third = one/three)
      parameter (two   = 2.e0, half  = one/two)
      parameter (twoth = two/three)

      dimension    x(iat),y(iat),z(iat)
      dimension    tsb(3),xyz(3),xyz1(3),cell(6)
      dimension    brvais(3,9),st(iat,192,3),ts(3,24),fs(3,3,24)
      dimension    ipt(iat), imult(iat)
      character*10 tag(iat)

c----------------------------------------------------------------------
c  b1 contains numbers appropriate to calculating translation symmetries
c  of i, r, f, and abc type cells.  i type cells require all positions to
c  be translated by (1/2,1/2,1/2).  r type cells require all positions to
c  be translated by (2/3,1/3,1/3) and (1/3,2/3,2/3).  f type cells require
c  all positions to be translated by (1/2,1/2,0),(1/2,0,1/2),& (0,1/2,1/2).
c  abc type cells require (1/2,1/2,0) where the 0 is in the abc position.
      data ((brvais(i,j),i=1,3),j=1,9)/
     $        zero,zero,zero,  half,half,half,    twoth,third,third,
     $        half,half,zero,  zero,half,half,    half,zero,half,
     $        half,half,zero,  third,twoth,twoth, zero,zero,zero/

c  the above is more precise but watcom compiler bug barks
c      data ((brvais(i,j),i=1,3),j=1,9)/
c     $        0.0,0.0,0.0,     0.5,0.5,0.5,     0.6667,0.3333,0.3333,
c     $        0.5,0.5,0.0,     0.0,0.5,0.5,     0.5,0.0,0.5,
c     $        0.5,0.5,0.0,     0.3333,0.6667,0.6667,    0.0,0.0,0.0/

c----------------------------------------------------------------------
c  symmetry associated with bravais lattice types (ie p,i,abc,r,f)
      nbra=1
      if (ibravl.gt.1) nbra=2
      if (ibravl.eq.3) nbra=3
      if (ibravl.eq.4) nbra=4

c----------------------------------------------------------------------
c     loop over all atoms
      do 110 ia=1,iatom

c       --- initialize flags and counters and load cell-axis-basis position
c       --- vector each time through
        xyz1(1)   = x(ia)
        xyz1(2)   = y(ia)
        xyz1(3)   = z(ia)
        imult(ia) = 0
        ipt(ia)   = 1

c       --- put coordinates in cell in first octant
        do 10 i=1,3
          if ((xyz1(i)+eps).lt.0.e0) xyz1(i) = xyz1(i)+one
          if (xyz1(i).ge.1.e0) xyz1(i) = xyz1(i) - 1.e0
          st(ia,ipt(ia),i) = xyz1(i)
 10     continue

c       --- loop over equivalent positions
        do 100 j=1,ns
c         --- loop over bravais lattice points
          do 90 nb=1,nbra
            ianti=1
c           --- return here to perform centrosymmetry operation
 20         continue

c           --- calculate atom coordinates in the cell axis basis,
c               perform rotation symmetries
            do 30 i=1,3
              xyz(i) = ianti * (ts(i,j) + fs(1,i,j)*xyz1(1) +
     $                 fs(2,i,j)*xyz1(2) + fs(3,i,j)*xyz1(3))

c             --- load appropriate values of bravais translation vector
              if (nb.eq.1)                     tsb(i) = brvais(i,1)
              if (nb.eq.2)                tsb(i) = brvais(i,ibravl)
              if (nb.eq.4)                     tsb(i) = brvais(i,6)
              if ((nb.eq.3).and.(ibravl.eq.3)) tsb(i) = brvais(i,8)
              if ((nb.eq.3).and.(ibravl.eq.4)) tsb(i) = brvais(i,5)

c             --- translate by the bravais lattice vector
              xyz(i) = xyz(i)+tsb(i)

c             --- put position back into first octant
              if ((xyz(i)+eps).ge.2.0)  xyz(i) = xyz(i)-two
              if ((xyz(i)+eps).ge.1.0)  xyz(i) = xyz(i)-one
              if ((xyz(i)+eps).le.-1.0) xyz(i) = xyz(i)+two
              if ((xyz(i)+eps).lt.0.0)  xyz(i) = xyz(i)+one
 30         continue

c  check if we had that position already in memory
c  nb: need to recognize that 0 and 1 are the same while considering
c      floating point precision problems
c  deincrement istred each time a coordinate is found
c  if all three are found skip over the storage block
            do 50 ii=1,ipt(ia)

              istred=3
              do 40 ichk=1,3
                posnew = xyz(ichk)
                posold = st(ia,ii,ichk)
                if ( ( abs( posnew-posold ).lt.eps)
     $            .or.(abs( abs(posnew-posold)-1 ).lt.eps) )
     $                     istred=istred-1
 40           continue
              if (istred.eq.0) goto 60
 50         continue

c           --- increment pointer and store new position
            ipt(ia) = ipt(ia)+1
            do 70 inew=1,3
              st(ia,ipt(ia),inew) = xyz(inew)
 70         continue

 60         continue
c           --- calculate number imult of coinciding atoms
c               imult counts how many coincidences, that is the
c               number of times the symmetry operations reproduced
c               the same point.
            dif=0.0
            do 80 k=1,3
              dif  = dif + abs( xyz(k)-xyz1(k) )*cell(k)
 80         continue
            if (abs(dif).lt.eps) then
                imult(ia) = imult(ia)+1
            endif

c           --- check for centrosymmetry
            if (isymce.ne.1) goto 90
            if (ianti.eq.-1) goto 90
            ianti=-1
            goto 20
 90       continue
 100    continue

c       --- error check, if ok then convert imult to a true multiplicity
c           and check consistency
        ierrl = 0
        if (imult(ia).eq.0) then
            call messag(' *** Warning: Multiplicity=0 for atom '//
     $                  tag(ia)//'.-')
            ierrl = 2
        else
            imult(ia) = igen/imult(ia)
            if (imult(ia).ne.ipt(ia)) then
                call messag(' *** Warning: Decoding error for atom '//
     $                      tag(ia)//'.-')
                ierrl = 2
            endif
c            imult(ia) = igen/imult(ia)
        endif
        if (ierrl.ne.0) ierr = ierrl
        if (ierrl.eq.2) then
            call messag(' *** Caution: Check your input file to make '//
     $                  'sure your')
            call messag('crystallographic information is correct.-')
        endif

c     --- go to next atom type
 110  continue
c       do 990 ii=1,iatom
c         print*, imult(ii)
c  990  continue

      return
c end subroutine multip
      end
