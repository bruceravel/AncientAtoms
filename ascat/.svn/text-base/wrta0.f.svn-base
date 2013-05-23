      subroutine wrta0(iat, ntitx, neptx,
     $                 nepts, ntit, nnoan, qvect, egrid, ecusp, anot,
     $                 fcore, afname, noantg, edge, core, title,
     $                 vrsn, vaxflg)
c===========================================================================
c write results of complex scattering factor calculation.
c output file is a five column file of energy in keV, energy with respect
c to the edge, real(A0), imag(A0), f" of resonant atom
c=====================================================================
c integers:
c    ntitx, neptx:  dimension parameters set in calling program
c    nepts:  number of energy points used in calculation
c    ntit:   number of title lines
c    nnoan:  number of tags or atoms for which an. corr. was neglected
c
c  reals:
c    qvect:  (3) momentum transfer vector
c    egrid:  energy grid for calculation
c    ecusp:  resonant energy of resonant atom
c    fcore:  (neptx) f" for resonant atom
c
c complex:
c    anot:   (neptx) complex scattering factor for crystal
c
c characters:
c    afname*72: name of output file
c    noantg*10: (iat) tags or element symbols for which an. corr. was neglected
c    edge*2:    resonant edge, K, L3
c    core*10:   site tag of resoanant atom
c    title*72:  (ntitx) title lines
c    vrsn*5:    version number as a character string
c
c  logical:
c     vaxflg:   vax flag for opening output file
c===========================================================================
      implicit integer(i-n)
      implicit real(a-h,o-z)
c      implicit double precision(a-h,o-z)

      parameter (iou=13)
c      parameter(neptx=2**11, ntitx=9)

      dimension    qvect(3)
      complex      anot(neptx), fcore(neptx)
      character*1  refl1
      character*2  edge, refl2
      character*5  vrsn
      character*10 core, noantg(iat)
      character*72 title(ntitx), afname, foo
      logical      there, vaxflg

 4000 format(2x,f12.4,2x,f12.4,4(2x,f10.4))
 4050 format('#',5x,'e',4x,'e-e0',4x,'real(a0)',4x,'imag(a0)',
     $            4x,'f''_core',4x,'f"_core')
 4100 format('#',1x,a)
 4200 format('#',60('-'))
 4300 format('# reflection: (', i2, ',', i2, ',', i2, ')' )
 4350 format('# core: ',a10,' edge: ',a2)
 4360 format(i1)
 4370 format(i2)
 4400 format('# ATOMS ',a5,3x,'by Bruce Ravel')
 4450 format('# anomalous correction neglected for ', a)

      if (afname.eq.' ') then
          afname = 'f'
          foo = afname
          do 10 i=1,3
            if ( (qvect(i).lt.9.5).and.(qvect(i).gt.-0.5) ) then
                write(refl1,4360)nint(qvect(i))
                ii = istrln(foo)
                afname = foo(:ii)//refl1
                foo = afname
            else
                write(refl2,4370)nint(qvect(i))
                ii = istrln(foo)
                afname = foo(:ii)//refl2
                foo = afname
            endif
 10       continue
          ii = istrln(foo)
          afname = foo(:ii)//'.dat'
          call lower(afname)
      endif

      there = .false.
      inquire(file=afname, exist=there)
      ii = istrln(afname)
      if (there) then
          call messag(afname(:ii)//' overwritten.')
      else
          call messag('Anomalous scattering written to '//afname(:ii))
      endif

      if (vaxflg) then
          open (unit=iou, file=afname, status='new')
      else
          open (unit=iou, file=afname, status='unknown')
      endif

      write(iou,4400)vrsn
      do 50 i=1,ntit
        write(iou,4100)title(i)
 50   continue
      write(iou,4300)nint(qvect(1)),nint(qvect(2)),nint(qvect(3))
      write(iou,4350)core, edge
      do 70 i=1,nnoan
        ii = istrln(noantg(i))
        write(iou,4450)noantg(i)(:ii)
 70   continue
      write(iou,4200)
      write(iou,4050)

      do 100 i=1,2*nepts
        e   = ecusp + (i-nepts+0.5) * egrid
        val = abs(anot(i))**2
        write(iou,4000)e*1000,(e-ecusp)*1000,
     $                 real(anot(i)),  aimag(anot(i)),
     $                 real(fcore(i)), aimag(fcore(i))
 100  continue

      return
c  end subroutine wrta0
      end
