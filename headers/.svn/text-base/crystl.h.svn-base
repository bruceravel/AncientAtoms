c-*-fortran-*-

c  various parameters used by module crystl
      common /cryint/ iabs, iatom, ibasis, isystm, ispa, iperm, nsites,
     $            ipt(iat), idop(iat), imult(iat)
      save /cryint/

      parameter(nsysm=4, nshwrn=4)
      character*2  dopant(iat,ndopx)
      character*10 spcgrp, tag(iat)
      character*74 shwarn(nshwrn)
      character*77 sysmes(nsysm)
      common /crystr/ shwarn, sysmes, dopant, tag, spcgrp
      save /crystr/

      logical syserr, shift
      common /crylog/ syserr, shift
      save /crylog/

      dimension trmtx(3,3), st(iat,192,3)
      dimension cell(6), x(iat), y(iat), z(iat)
      dimension percnt(iat,ndopx)
      common /cryflt/ trmtx, st, cell, x, y, z, percnt
      save /cryflt/
c -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c  GLOSSARY:  (* = user input, % = error handling, ! = output needed to
c              construct cluster, the rest are used internally)
c
c * iabs:   index of absorber in unique coordinate list
c * iatom:  >=1 if atoms list is used, else =1
c * ibasis: =1 if basis list is used, else =0
c   isystm: index of crystal system (1..7)=(mono,orth,<not used>,tetr,
c           cubic,hex,triclinic)
c   ispa:   space group index, 1-230 from IXTC
c             0:       not recognized, error in input symbol
c             1-2:     triclinic
c             3-15:    monoclinic
c             16-74:   orthorhombic
c             75-142:  tetragonal
c             143-167: trigonal
c             168-194: hexagonal
c             195-230: cubic
c   iperm:  permutation index for non-standard settings, used in crystl
c             1:     default value -- no permutation necessary
c             1-6:   6 orthorhombic settings (abc, cab, bca, a-cb, ba-c, -cba)
c             11-12: 2 monoclinic settings, z-axis unique, y-axis unique
c             21-22: 2 tetragonal settings, standard and rotated
c ! ipt:    (iat) number of positions of unique atom in unit cell (1..192)
c   imult:  (iat) workspace for calculating multiplicities
c
c % sysmes: 4 line message if space group and axes/angles don't match
c % syserr: true if space group and axes/angles don't match
c % shwarn: 3 line message if space group may require a shift vector
c % shift:  true if space group may require a shift vector
c
c * dopant*2:  (iat,ndopx) matrix with all host and dopant atomic symbols
c * percnt:    (iat,ndopx) matrix with occupancies of hosts and dopants
c * tag*10:    (iat) character tag for each unique site in cell
c * spcgrp*10: space group symbol.  On output it is the short
c              Hermann-Maguin symbol in standard setting.  On input
c              spcgrp can be any short HM, Schoenflies, a number
c              between 1 and 230, or one of a small set of special
c              words (fcc, bcc, etc.).  Other symbol conventions
c              (full HM symbol, Shubnikov, 1935 ITXC, etc.) are not
c              and never will be used.
c
c * x,y,z:  (iat) arrays of fractional coordinates of unique positions
c                 in unit cell
c * cell:   (6) array of a,b,c,alpha,beta,gamma
c
c ! trmtx:  (3,3) transformation matrix between cell-axis and cartesian
c                 bases, see subroutine trans in clustr
c ! st:     (iat,192,3) fractional coordinates of all atoms in unit cell,
c                       first arg refers to unique atom list, second
c                       to position in cell, third is xyz.
c
c           fyi: 192 is the largest possible number of equivalent
c                positions in a cell of any symmetry. see, for example,
c                cubic f m 3 c.
c----------------------------------------------------------------------
