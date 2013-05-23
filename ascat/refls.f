      subroutine refls(iat, ndopx, ntitx, nsites, nrefl, ntit,
     $            refile, core, dopant, tag, title, vaxflg, st, ipt)

      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

c------------------------------------------------------------------
c  calculate exp[ 2pi Q.R ] for each resonant site in the structure
c  for different reflections
c  this can be used to determine weight from each site in a dafs
c  spectrum
c------------------------------------------------------------------
c  input:
c     iat, ndopx, ntitx:  parameters set in calling routine
c     nsites:  number of unique crystallographic sites
c     nrefl:   (3) largest h k and l of calculation
c     ntit:    number of title lines to write to output file
c     refile:  output file name
c     core:    atomic symbol of central atom
c     dopant:  (iat,ndopx) list of atoms in cell
c     tag:     (iat) identifying tag of each site
c     title:   (ntitx) title lines
c     vaxflg:  vax flag for how to open output file
c     st:      (iat,192,3) array of coordinates
c     ipt:     (iat) number of positions of each site in cell
c------------------------------------------------------------------

      parameter (pi = 3.1415 92653 58979 32384 62643e0)
      parameter (twopi = pi+pi)
      parameter (nresx=5)
      character*2  core, dopant(iat,ndopx)
      character*10 tag(iat)
      character*72 title(ntitx), refile
      complex      phfact(nresx)
      dimension    st(iat,192,3), ipt(iat), nrefl(3)
      dimension    ires(100)
      integer      h
      logical      vaxflg, there

 400  format(a)

      if ((nrefl(1).eq.0).and.(nrefl(2).eq.0).and.(nrefl(3).eq.0)) then
          call messag(' ')
          call messag('* * * ERROR in reflection amplitude '//
     $                'calculation!')
          call messag('Number of requested reflections must be '//
     $                'positive.')
          call messag('Reflection amplitude calculation skipped.')
          call messag(' ')
          return 
      endif

c       if ((nrefl(1).gt.10).or.(nrefl(2).gt.10).or.(nrefl(3).gt.10)) then
c           call messag(' ')
c           call messag('* * * WARNING in reflection amplitude '//
c      $                'calculation!')
c           call messag('You asked for a very large number of '//
c      $                'reflections.')
c           call messag('This will produce a very large output file.')
c           call messag(' ')
c       endif

      irefl = nxtunt(12)
      inquire(file=refile, exist=there)
      if (vaxflg) then
          open(unit=irefl, file=refile, status='new')
      else
          open(unit=irefl, file=refile, status='unknown')
      endif
      
      nres = 0
      do 10 i=1,nsites
        if (dopant(i,1).eq.core) then
            nres = nres+1
            ires(nres) = i
        endif
 10   continue 
      
      do 20 i=1,ntit
        write(irefl,400)'# '//title(i)
 20   continue 
      write(irefl,400)'# exp[2pi Q dot R] for each resonant site at '//
     $            'various reflections'
      write(irefl,410)
 410  format('#',71('-'))
      write(irefl,420)(tag(ires(i)), i=1,nres)
 420  format('# h  k  l  ',5(a10,6x))

      do 130 l=0,nrefl(3)
        do 120 k=0,nrefl(2)
          do 110 h=0,nrefl(1)

            if ((h+k+l).eq.0) goto 110
            do 100 isite=1,nres
              phfact(isite) = cmplx(0.e0, 0.e0)
              do 50 jpt=1,ipt(ires(isite))

                xq = st(ires(isite),jpt,1) * h
                yq = st(ires(isite),jpt,2) * k
                zq = st(ires(isite),jpt,3) * l
c               crystallographic phase factor, 2*pi*(q dot r)
                phase  = twopi*(xq + yq + zq)
                phfact(isite) = phfact(isite) +
     $                      cmplx( cos(phase), sin(phase) )

 50           continue 
 100        continue 

            write(irefl,430)h,k,l,(phfact(i),i=1,nres)
 430        format(3(1x,i2),5(1x, '(', f6.2, ',', f6.2, ')'))

 110      continue 
 120    continue 
 130  continue 

      close(irefl)
      ii = istrln(refile)
      if (there) then
          call messag(refile(:ii)//' overwritten.')
      else
          call messag(refile(:ii)//' written.')
      endif

      return 
c  end subroutine refls
      end
