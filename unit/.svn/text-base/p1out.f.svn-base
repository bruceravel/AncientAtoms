      subroutine p1out(iat,ndopx,ntitx,nsites,ntit,iabs,dopant,edge,
     $                 cell,rmax,st,ipt,title,vrsn,vaxflg)
c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c--------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer(i-n)
c      implicit double precision(a-h,o-z)

      parameter(ip1=19)

      character*2  dopant(iat,ndopx), edge, fname*6, vrsn*9, elmn
      character*3  ci, cj, tag*15, cabs
      character*72 title(ntitx)
      dimension    ipt(iat), cell(6), st(iat,192,3)
      logical      there, vaxflg

 400  format(a)
 405  format(a5,3x,a4)
 410  format(3(a5,2x,f9.4,3x))
 420  format('core ',3x,a2,a,a2,6x,'edge ',3x,a2,9x,'rmax ',3x,f7.4)
 430  format(2x,a2,3(3x,f8.4),3x,a)
 440  format(2x,a2,3(3x,f8.4),3x,a2,i1,'_',i2.2)
 450  format(i3)

      fname = 'p1.inp'
      call lower(fname)

      inquire(file=fname, exist=there)
      ii = istrln(fname)
      if (there) then
          call messag(fname(:ii)//' overwritten.')
      else
          call messag('Entire unit cell written as an input file to '//
     $                fname(:ii))
      endif

      if (.not.vaxflg) then
          open(unit=ip1,file=fname,status='unknown')
      else
          open(unit=ip1,file=fname,status='new')
      endif

c  write titles
      write(ip1,400)'! ATOMS '//vrsn//' by Bruce Ravel'
      do 10 i=1,ntit
        ii = istrln(title(i))
        ii = min(ii,71)
        write(ip1,400)'title  '//title(i)(:ii)
 10   continue
      write(ip1,400)'title  unit cell written as space group p 1'

c  write keywords
      write(ip1,405)'space', 'p 1'
      write(ip1,410)'a    ', cell(1), 'b    ', cell(2), 'c    ',cell(3)
      write(ip1,410)'alpha', cell(4), 'beta ', cell(5), 'gamma',cell(6)

c  construct tag for core, write more keywords
      ii=istrln(dopant(iabs,1))
      write(cabs,450)iabs
      call triml(cabs)
      ia=istrln(cabs)
      write(ip1,420)dopant(iabs,1)(:ii),cabs(:ia),'_1', edge, rmax
      write(ip1,400)'atoms'

c  write atom positions
      do 30 i=1,nsites
        do 20 j=1,ipt(i)

c         construct tag for each atom:  Abi_j
          id  = istrln(dopant(i,1))
          write(ci,450)i
          call triml(ci)
          ii  = istrln(ci)
          write(cj,450)j
          call triml(cj)
          ij  = istrln(cj)
          tag = dopant(i,1)(:id)//ci(:ii)//'_'//cj(:ij)
          elmn = dopant(i,1)
          call fixsym(elmn)

          write(ip1,430)elmn,st(i,j,1),st(i,j,2),st(i,j,3),tag

 20     continue
 30   continue

      close(ip1)

      return
c  end subroutine p1out
      end
