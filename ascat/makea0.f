      subroutine makea0(iat, ndopx, neptx,
     $            nsites, ipt, idop, nepts, nnoan,
     $            edge, core, dopant, tag, noantg, egr, ecusp,
     $            st, usqr, cell, percnt, qvect, anot, fcore,f0,
     $            lf)

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)
c-----------------------------------------------------------------
c  construct anomoulous scattering factor from full atomic scattering
c  factors and crystallographis information.
c-----------------------------------------------------------------
c  makea0 calls getf0 (tabulation of Cromer-Mann aikman coefficients)
c  and fcal (tabulation of calculations of f' and f" for
c  interpolation)
c-----------------------------------------------------------------
c  integers:
c    iat,ndopx,neptx:  dimension parameters set in main program
c    nsites: number of unique sites in input file
c    ipt:    number of positions of a unique atom in unit cell
c    idop:   number of dopants at each site
c    nepts:  number of points in energy for calculation
c    nnoan:  number of tags for which anomalous correction should be neglected
c  characters:
c    edge*2:   edge of central atom
c    core*2:   symbol of central atom
c    dopant*2: (iat, ndopx) atoms at each site
c    tag*10:   (iat) cite tags
c    noantg*10: (iat) tags of atoms for which an. corr. should be neglected
c  reals
c    egr:    grid size in energy
c    ecusp:  resonant energy in problem
c    st:     (iat,192,3) fractional coords of all atoms in unit cell
c    usqr:   (iat) crystallographic DW factors
c    cell:   (6) array of a,b,c,alpha,beta,gamma
c    percnt: (iat,ndopx) occupancies of dopants
c    qvect:  (3) reflection in problem
c  logicals:
c    lf:     true --> write diagnostic file
c  output:
c    anot:   (neptx) complex anomalous scattering factor
c    fcore:  (neptx) f" for central atom
c-----------------------------------------------------------------
c      parameter (iat=50, ndopx=4, neptx=2**11)
      parameter (iouf=17)
      parameter (pi = 3.1415 92653 58979 32384 62643)
      parameter (one = 1, mone = -1, zero = 0, twopi = pi+pi)

      character*2  edge,dopant(iat,ndopx), el, test
      character*10 core, fname*8, tag(iat), noantg(iat),
     $             tagup, noan, dp
      dimension    ipt(iat),idop(iat)
      dimension    percnt(iat,ndopx),qvect(3)
      dimension    st(iat,192,3), usqr(iat)
      complex      fcore(neptx)
      dimension    f0(iat,ndopx)
      complex      anot(neptx), phfact
      logical      lf, there

      test   = 'ab'
      ecusp  = s2e(core,edge)
      call spcing(qvect,cell,d)
      s      = one / (2*d)
      qsqr   = qvect(1)**2 + qvect(2)**2 + qvect(3)**2
      dwf    = one

c  fill matrix of f0 values used below, and check that all the elements
c  are available for f' and f"
      do 20 i=1,nsites
        do 10 j=1,idop(i)
          el = dopant(i,j)
          f0(i,j) = getf0(el, s)
          call avlble(el, there)
          call case(test,el)
          if ( (.not.there).and.(el.ne.'nu') ) then
              call fixsym(el)
              call messag(' ')
              call messag('* * * CAUTION * * *')
              call messag(el//' is not available in the tables of '
     $                    //'anomalous scattering.')
              call messag('0 was used for f'' and f".')
              call messag(' ')
          endif
 10     continue
 20   continue

c  ie loops thru energy values, i loops thru sites, j loops thru
c  multiplicity of sites, k loops thru dopants at a site
      do 200 ie=1,2*nepts
        energy = (ecusp + (ie-nepts+0.5) * egr) * 1000
c        print*,ie,energy

        do 100 i=1,nsites
c         crystallographic debye-waller factor, certainly wrong for dopants
c          dwf  = exp( mone * qsqr * usqr(i) )

          do 70 k=1,idop(i)
c            print*,dopant(i,k),energy

c           neglect anomalous correction for some atoms

c           --- case 1 is when a site tag is specified
            if (k.eq.1) then
                tagup=tag(i)
                call case(test,tagup)
                do 30 nn=1,nnoan
                  noan=noantg(nn)
                  call case(test,noan)
                  if (noan.eq.tagup) then
                      fp  = zero
                      fpp = zero
                      goto 40
                  endif
 30             continue
            endif

c           --- case 2 is when an element symbol is specified
            dp   = dopant(i,k)
            call case(test,dp)
            do 35 nn=1,nnoan
              noan = noantg(nn)
              call case(test,noan)
              if (noan.eq.dp) then
                  fp  = zero
                  fpp = zero
                  goto 40
              endif
 35         continue
c           --- get here if both of the tests above are failed
            call fcal(dopant(i,k), energy, fp, fpp)
 40         continue
c            print*,dopant(i,k),fp,fpp
            if (dopant(i,k).eq.core) then
c               --- want this even if noantg is set for central atom
                call fcal(dopant(i,k), energy, ffp, ffpp)
                fcore(ie) = cmplx(ffp, ffpp)
            endif

            do 50 j=1,ipt(i)
c              print*,dopant(i,k),st(i,j,1),st(i,j,2),st(i,j,3)
              xq = st(i,j,1) * qvect(1)
              yq = st(i,j,2) * qvect(2)
              zq = st(i,j,3) * qvect(3)
c             crystallographic phase factor, 2*pi*(q dot r)
              phase  = twopi*(xq + yq + zq)
              phfact = cmplx( cos(phase), sin(phase) )

              anot(ie) = anot(ie) + cmplx( (f0(i,k)+fp), fpp ) *
     $                   phfact * percnt(i,k) * dwf

 50         continue
 70       continue
 100    continue

 200  continue

      if (lf) then
 400      format(a)
 410      format(a,a,a,f10.6)
 420      format(a,f8.4)
 430      format(a,3(f5.1))
          ii = istrln(core)
          fname = 'f-'//core(:ii)//'.dat'
          call lower(fname)
          open (unit=iouf, file=fname, status='unknown')
          write(iouf,400)'# f'' and f" interpolated from data'
          write(iouf,410)'# core,edge,ecusp: ',core,edge,ecusp
          write(iouf,420)'# plane spacing in angstroms: ',d
          write(iouf,430)'# qvect: ',qvect
          write(iouf,400)'#-------------------------------------------'
          write(iouf,400)'#        e              f''             f"'
          do 1000 i=1,2*nepts
            energy = (ecusp + (i-nepts+0.5) * egr) * 1000
            call fcal(core(:2), energy, fp, fpp)
            write(iouf,*)energy, fp, fpp
 1000     continue
          close(iouf)
          call messag('  wrote '//fname//'...')
      endif

      return
c  end subroutine makea0
      end
