      subroutine geout(ngeomx, ntitx, natx,
     $            iug, itot, ntit, ngeom, maxln,
     $            atlis, vrsion, title, vaxflg)

c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision (a-h,o-z)

c      parameter(ngeomx=800)
      parameter (nwdx=20,maxty=7)
      parameter (eps=0.0001)

      dimension    atlis(natx,8)
      character*2  elist(maxty),test,vrsion*5
      character*20 fname,feffcd
c      character*20 words(nwdx)
      character*72 title(ntitx)
      dimension    ngeom(ngeomx)
      logical      vaxflg

 4000 format (a)
 4010 format (1x,a)
 4020 format (bn,f10.0)
 4030 format (1x,70('-'))
 4040 format (1x,a,2x,a)
 4100 format (1x,i3,3(2x,f10.6),2x,i2,1x,i3)

      test  = 'ab'
      fname = 'geom.dat'
      call lower(fname)
      if (.not.vaxflg) then
          open (unit=iug,file=fname,status='unknown')
      else
          open (unit=iug,file=fname,status='new')
      endif

      write(iug,4010)'**this geom.dat file was made by atoms '//vrsion
      write(iug,4010)'**use this file with the feff.inp file '//
     $               'that was created at the same time.'
c----------------------------------------------------------
c  write the titles
      call card('title',feffcd,ii)
      do 105 i=1,ntit
        write(iug,4040)feffcd(:ii),title(i)(:75-ii)
105   continue
      write(iug,4030)

      rlast = -1.
      icnt  = 0
      ne    = 0
      do 10 i=1,maxty
        elist(i) = ' '
 10   continue
      do 100 i=1,itot
c       want geom.dat to be the same length as feff.inp
        if (itot.gt.maxln) goto 200
        igeo = 0
        x = atlis(i, 6)
        y = atlis(i, 7)
        z = atlis(i, 8)
        r = ref(x,y,z)
        if (abs(r-rlast).gt.eps) then
            ipot = 0
            do 50 j=1,ne
c              if ( (words(1)(1:2).eq.elist(j)).and.(i.ne.1) ) ipot = j
              if (i.ne.1) ipot = j
 50         continue
            if ( (ipot.eq.0).and.(i.ne.1)) then
                ne        = ne + 1
                ipot      = ne
c                elist(ne) = words(1)(1:2)
            endif
            icnt  = icnt + 1
            igeo  = ngeom(icnt)
            rlast = r
        endif
        write(iug,4100) (i-1), x, y, z, ipot, igeo
 100  continue

 200  continue
      close (iug)

      return

c  end subroutine geout
      end
