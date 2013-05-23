      subroutine equipt (isystm,isymce,cbravl,ng,tx,fx,spcgrp)
      implicit real (a-h,o-z)
c      implicit double precision (a-h,o-z)
c----------------------------------------------------------------------
c     interprets the  herman-mauguin symbol for space group             
c     modified from a program by h.burzlaff, j.appl.cryst, v.15, 
c     p.464 (1982)
c  double precision and character variable manipulations: br 8/92
c----------------------------------------------------------------------
c  input:
c    spcgrp*10:  hermann-maguin space group notation
c  output:
c    isystm:     integer denoting bravais lattice, (1..7)=(a,b,r,c,i,f,p)
c    isymce:     centrosymmetry flag, 1=centrosymmetric
c    cbravl*1:   character denoting bravais type, 1st letter in spcgrp
c    ng:         simple multiplicity of lattice type
c    tx(3,24):   fractional coordinates of simply multiple atoms
c    fx(3,3,24): rotational symmetries
c----------------------------------------------------------------------
      parameter (eps=.001e0)
      parameter (one = 1.e0, two = 2.e0, four = 4.e0)
      parameter (half = one/two, quart = one/four)
      integer      e(3,3), sys, isys(7), nul(3)
      dimension    sh(3), te(3)
      dimension    ss(24,3,3), tx(3,24), ts(24,3), fx(3,3,24)
      character*10 spcgrp

      character*1  cspace(10), cbra(7), cbravl, csym(3,4), cbr
      character*1  cm, cn
      character*2  test
      logical      blank

      data (sh(i),  i=1,3)  /0.25e0, 0.25e0, 0.25e0/
      data (nul(i), i=1,3)  /0, 0, 0/
      data (isys(i),i=1,7)  /7, 1, 2, 4, 6, 5, 3/
      data (cbra(i),i=1,7)  /'p','a','b','c','f','i','r'/

c%%%%%
      test  = 'ab'
c --- this stuff was added after running the code through g77 --------
      blank = .false.
      ichar = 0
      nbr   = 0
c --------------------------------------------------------------------
c  the value of test *must* be of the same case as the values of cbra
c-----------------------------------------------------
c  break spcgrp up into its 10 components
c  cspace contains the (up to) 10 characters in the space group designation
      do 9001 i=1,10
        cspace(i) = ' '
9001  continue
      do 9011 i=1,10
        read (spcgrp(i:i),4000) cspace(i)
        call case(test,cspace(i))
9011  continue
4000  format(a1)

c-----------------------------------------------------
c     initialize some things
      nsym = -1
      ng   =  1
      ns   =  1
      do 24 i=1,24
         do 22 j=1,3
           ts(i,j) = 0.e0
           tx(j,i) = 0.e0
           do 20 k=1,3
              e(j,k)    = 0
              ss(i,j,k) = 0.e0
              fx(j,k,i) = 0.e0
 20        continue
 22      continue
 24   continue
      ss(1,1,1) = 1.e0
      ss(1,2,2) = 1.e0
      ss(1,3,3) = 1.e0
      do 32 i=1,3 
         te(i) = 0.e0
         do 30 j=1,4
            csym(i,j) = ' '
 30      continue
 32   continue

c-----------------------------------------------------
c  begin interpreting the space group

      do 80 i=1,10
c        reset blank (flag) and ichar each time a blank is found
c        then go on to the next character
         if (cspace(i).eq.' ') then
             blank  = .true.
             ichar  = 0
             goto 80
         endif

c  n.b. this could be cleaned up since spcgrp is trimled before this

c  if bravais lattice type has been found...
         if (nsym.ge.0) goto 70

c  look for bravais lattice type    1..7 = p,a,b,c,f,i,r
c  cbr for internal use, i/cbravl for external use
         do 50 j=1,7
            cbr = cbra(j)
            nbr = j
            if (cspace(i).eq.cbra(j)) goto 60
50       continue

60       continue
         nsym    = 0
         cbravl  = cspace(i)
         goto 80 

c  nsym counts the three symbols in the space group after the lattice type.
c  ichar counts the characters in that symbol.  nsym=[1..3],ichar=[1..4].
c  blank=t flags the beginning of a new symbol.
70       continue
         if (blank) nsym  = nsym+1
         ichar            = ichar+1
         csym(nsym,ichar) = cspace(i)
         blank            = .false.
         if (cspace(i).eq.'/') ns = 0
 80   continue

c      print*,'spcgrp=',spcgrp
c      print*,'cspace=',cspace

c-----------------------------------------------------
c  determine the bravais system
c  sys=1...6 ==> (tric,mono,ortho,tetra,hex,cubic)

c  cubes have triad axes.
      do 90 i=1,4
         if (csym(2,i).ne.'3') goto 90
         sys = 6
         goto 145
 90   continue  

c  hexagonal cells have hexads, trigonal cells have triads
c  tetragonal cells have tetrads.
      do 110 i=1,4
         if ( (csym(1,i).eq.'3').or.(csym(1,i).eq.'6') ) then
               sys = 5
               goto 145
         endif
         if (csym(1,i).eq.'4') then
             sys = 4
             goto 145
         endif
110   continue

c     choose between triclinic and monoclinic
      if (nsym.eq.1) then
c         monoclinic might only have one symbol, but it's not "1"
          if ( (csym(1,1).eq.'1').or.(csym(1,1).eq.'-') ) then
                sys = 1 
                goto 145
          endif
          sys = 2
c         fix csym for point groups 2 and m, but not 2/m
          do 130 i=1,4
             csym(2,i) = csym(1,i)
130       continue
          csym(1,1) = '1'
          csym(3,1) = '1'
      endif

c%#@%#@  this seems to be wrong for group 2/m
c  orthorhombic by default unless monoclinic pt. grp. 2 or m
      sys=3
      if ( (csym(1,1).eq.'1').or.(csym(2,1).eq.'1') ) sys = 2

c---------------------------------------------------------------------
c   leap ahead to the appropriate section of code for the bravais 
c   type in question 
c   jump to these line numbers: (1000,1100,1200,1300,1400,1500)
c   for appropriate sys, sys=[1..6] according to lattice type:
c   (tric,mono,ortho,tetra,hex,cubic)
145   continue

      if     (sys.eq.1) then
              goto 1000
      elseif (sys.eq.2) then
              goto 1100
      elseif (sys.eq.3) then
              goto 1200
      elseif (sys.eq.4) then
              goto 1300
      elseif (sys.eq.5) then
              goto 1400
      elseif (sys.eq.6) then
              goto 1500
      endif

c---------------------------------------------------------------------
c  triclinic,  described by group containing only the identity (p1)
c              or the identity and a parity operation (p-1)
1000  if (csym(1,1).eq.'-') ns = 0
      goto 760

c---------------------------------------------------------------------
c  monoclinic, primitive (p) or one-face-centered (c).
c              point gps: diad (p2,p21,c2),mirror (pm,cm),glide (pc,cc), 
c                         combined (p2/m,p21/m,p2/c,p21/c,c2/m,c2/c)
1100  ng = 2
      ind = 0
      do 180 i=1,3
        if (csym(i,1).ne.'1') ind = i
180   continue
      id = 1
      if (csym(ind,1).eq.'2') id = -1
      do 190 i=1,3
        ss(2,i,i) =  ss(1,i,i)*id
190   continue
      ss(2,ind,ind) = -ss(2,ind,ind)
      do 220 i=1,3
        if ( (csym(i,1).eq.'2').and.(csym(i,2).eq.'1') ) ts(2,i)=0.5
        do 210 j=1,4
          if (csym(i,j).eq.'a') ts(2,1) = 0.5
          if (csym(i,j).eq.'b') ts(2,2) = 0.5
          if (csym(i,j).eq.'c') ts(2,3) = 0.5
          if (csym(i,j).eq.'n') goto 200
          goto 210
 200      continue 
          k=i+1
          if (k.gt.3) k = k-3
          ts(2,k)     = 0.5
          ts(2,6-k-i) = 0.5
 210    continue
 220  continue 
      goto 760

c---------------------------------------------------------------------
c --- orthorhombic
1200  continue
      ng  = 4
      ic  = 0
      ind = 1

c  if not point group 222...
      if ( (csym(1,1).ne.'2').and.(csym(2,1).ne.'2').and.
     $     (csym(3,1).ne.'2') )  then
            ind = -1
            ns = 0
      endif

      do 240 i=1,3
         id = 1
c        if any diads exist...
         if (csym(i,1).eq.'2') id = -1
         do 230 j=1,3
           ss(1+i,j,j) = ss(1,j,j)*id*ind
230      continue
         ss(1+i,i,i) = -ss(1+i,i,i)
240   continue

      do 282 i=1,3
c        if screw diads exist...
         if ( (csym(i,1).eq.'2').and.(csym(i,2).eq.'1') )
     $          ts(1+i,i) = 0.5
         do 280 j=1,4
c           if a glide plane...; if b glide plane...; if c glide plane...;
c           if diagonal glide plane or 1/4 glide plane
            if  (csym(i,j).eq.'a') ts(1+i,1) = 0.5
            if  (csym(i,j).eq.'b') ts(1+i,2) = 0.5
            if  (csym(i,j).eq.'c') ts(1+i,3) = 0.5
            if ((csym(i,j).eq.'n').or.(csym(i,j).eq.'d')) goto 250
            goto 280

250         continue
            k = i+1
            if (k.gt.3) k = k-3
            if (csym(i,j).eq.'d') goto 260 

c           case of not diagonal glide plane 
            ts(1+i,k)     = 0.5
            ts(1+i,6-k-i) = 0.5
            goto 280

c           case of 1/4 glide plane 
260         ic = 1
c           if no centrosymmetry...
            if (ns.eq.1) goto 270
            ts(1+i,k)     = 0.25
            ts(1+i,6-k-i) = 0.25
            goto 280
 270        ts(1+i,1)     = 0.25
            ts(1+i,2)     = 0.25
            ts(1+i,3)     = 0.25
 280     continue
 282  continue

c
      if (ic.eq.1) goto 760 
c     if no centrocymmetry...
      if (ns.eq.1) goto 360  

      do 290 i=1,3
         k = 1+i
         if (k.gt.3) k = k-3
         tc = ts(1+k,i)+ts(1+6-i-k,i)
         if (tc.eq.1.0) tc = 0.0
         ts(1+i,i) = tc
290   continue

c----------------------------------------------------------------
c  special treatment of c m m a, c m c a, i m m a
c  ma counts the number of "m's" (ie 2, 1, or 2)
      if ((cbr.eq.'p').or.(cbr.eq.'f')) goto 760
      ma = 0
      do 310 i=1,3 
         do 300 j=1,4
           if (csym(i,j).eq.'m') nul(i)=1
300      continue
         ma = ma+nul(i)
310   continue
c     if i m m a...do origin shift
      if ( (cbr.eq.'i').and.(ma.eq.2) ) goto 330
c     weed out remaining orthorhombic c and i space groups
      if ( (ma.eq.0).or.(ma.eq.3).or.(cbr.eq.'i') ) goto 760
c%#%#%#%#%#%#%#%#
      do 320 i=1,3
c  is this a typo?!   what is the purpose of the do loop -- this is equally
c                     mysterious in sexie's equipt.  will deal with this 
c                     when time comes.
         if (nul(nbr-1).eq.1) goto 760
         sh(nbr-1) = 0
320   continue

c --- origin shift
 330  do 352 i=1,ng
         do 350 j=1,3
            do 340 k=1,3
               id = 1
               if (j.ne.k) id = 0
               ts(i,j) = ts(i,j) + (id-ss(i,j,k))*sh(k)
340         continue
            if (ts(i,j).gt.1.)  ts(i,j) = ts(i,j)-1
            if (ts(i,j).lt.1.0) ts(i,j) = ts(i,j)+1
350      continue
352   continue
      goto 760

c...where to come for no centrocymmetry
360   ic = 0
      do 370 i=1,3
c                                                  .eq.-1
        if ( abs((ss(1+i,1,1)*ss(1+i,2,2)*ss(1+i,3,3)) +1).lt.eps) ic=1
370   continue
      if (ic.eq.1) goto 410
      tc = ts(2,1)+ts(3,2)+ts(4,3)
      if (abs(tc).lt.eps) goto 760
      do 400 i=1,3
         k = i+1
         if (k.gt.3) k = k-3 
         if (tc.gt.0.5) goto 380
         if (abs(ts(1+i,i)).lt.eps) goto 400 
         m = i-1
         if (m.eq.0) m = m+3
         ts(1+m,i) = 0.5
         goto 400
 380     if (tc.gt.1.0) goto 390
         if (abs(ts(1+i,i)).gt.eps) goto 400
         l = k+1
         if (l.gt.3) l = l-3  
         ts(1+k,l) = 0.5
         ts(1+l,k) = 0.5
         goto 400      
 390     ts(1+i,k) = 0.5
 400  continue 
      goto 760

410   do 420 i=1,3
        if (abs(ss(1+i,1,1)*ss(1+i,2,2)*ss(1+i,3,3)-1).lt.eps) id=i
420   continue
      do 450 i=1,3
         tc = ts(2,i)+ts(3,i)+ts(4,i)
         if ( (abs(tc).lt.eps).or.(abs(tc-1.0).lt.eps) ) goto 450
c        if mirror and diagonal glide planes exist... 
         if ( ((csym(1,1).eq.'m').and.(csym(2,1).eq.'n')) .or.
     $        ((csym(2,1).eq.'m').and.(csym(3,1).eq.'n')) .or.
     $        ((csym(3,1).eq.'m').and.(csym(1,1).eq.'n')) ) goto 440
         do 430 j=1,3
            if (id.eq.j) goto 430
            if (abs(ts(1+j,i)-0.5).lt.eps) goto 430
            ts(1+j,i) = 0.5
 430     continue
         goto 450
 440     k=i-1
         if (k.eq.0) k = k+3
         ts(1+k,i) = 0.5
 450  continue 
      goto 760

c---------------------------------------------------------------------
c --- tetragonal


1300  ng = 4
      if (nsym.eq.3) ng=8
      ss(2,1,2) = -1
      ss(2,2,1) =  1
      ss(2,3,3) =  1
      cm        = csym(1,1)
      cn        = csym(1,2)

      do 472 i=1,3  
         do 470 j=1,3
c           if enantiomorphous...
            if (cm.eq.'-') ss(2,i,j) = -ss(2,i,j)
 470     continue
 472  continue

      if (cm.eq.'-') goto 500
c     if screw axes of 1/4, 1/2. or 3/4 rotation exist...
      if (cn.eq.'1') ts(2,3) = 0.25
      if (cn.eq.'2') ts(2,3) = 0.5 
      if (cn.eq.'3') ts(2,3) = 0.75

c   if horizontal diagonal glide plane or hdgp & # of symbols=3...
      if ( (csym(1,3).eq.'n') .or. ((csym(1,4).eq.'n').and.(nsym.eq.3)))
     $                 ts(2,1) = 0.5

c   if hdgp & # of symbols=1 (p42/n) or 1/4 screw axis & no cent. & i...
      if ( ((csym(1,4).eq.'n').and.(nsym.eq.1)) .or. 
     $     ((cn.eq.'1').and.(ns.eq.1).and.(cbr.eq.'i')) ) ts(2,2) = 0.5

c   if 1/4 screw axis & no cent. & i...
      if (  (cn.eq.'1').and.(ns.eq.0).and.(cbr.eq.'i') ) goto 480

c   if secondary 1/4 screw axis or no hdgp & vert. diag. glide plane &
c                                  third mirror plane...
      if ( (csym(2,2).eq.'1') .or. 
     $    ((csym(1,4).ne.'n').and.(csym(2,1).eq.'n').and.
     $     (csym(3,1).eq.'m')) ) goto 490

      goto 500

480   continue
      ts(2,1) = 0.25
      ts(2,2) = 0.75
c     if point group 4, -4, or 4/m...
      if (nsym.eq.1) ts(2,1) = 0.75
      if (nsym.eq.1) ts(2,2) = 0.25
      goto 500

490   continue
      ts(2,1)   = 0.5
      ts(2,2)   = 0.5

500   continue
      ss(3,1,1) = -1
      ss(3,2,2) = -1
      ss(3,3,3) =  1
      ts(3,1)   = ss(2,1,2)*ts(2,2)+ts(2,1)
      ts(3,2)   = ss(2,2,1)*ts(2,1)+ts(2,2)
      ts(3,3)   = ss(2,3,3)*ts(2,3)+ts(2,3)

      do 510 i=1,3  
        if ( (cbr.eq.'i').and.(abs(ts(3,1)-0.5).lt.eps).and.
     $       (abs(ts(3,2)-0.5).lt.eps).and.(abs(ts(3,3)-0.5).lt.eps) )  
     $                   ts(3,i) = 0.0
510   continue

      do 524 i=1,3  
         ts(4,i) = ts(2,i)
         do 522 j=1,3
            ts(4,i) = ts(4,i)+ss(2,i,j)*ts(3,j)
            do 520 k=1,3
               ss(4,i,j) = ss(4,i,j)+ss(2,i,k)*ss(3,k,j)
 520        continue
 522     continue
 524  continue

c     if point group 4, -4, or 4/m then done...
      if (nsym.eq.1) goto 760
c     if centrosym...
      if (ns.eq.0) goto 560

      cm = csym(2,1)
      cn = csym(3,1)
c     if no secondary diad exists...
      if ( (cm.ne.'2').and.(cn.ne.'2') ) goto 550
c     if 2 secondary diads exist...
      if ( (cm.eq.'2').and.(cn.eq.'2') ) goto 540
c     if 2nd symbol 1st char = 2 skip to 530..
      if (cm.eq.'2') goto 530

c     if c glide plane or diag. glide plane exists...
      if ( (cm.eq.'c').or.(cm.eq.'n') )  te(3) = 0.5
      e(1,1) = -1
      e(2,2) =  1
      e(3,3) =  1
c     skip to 570 if no diag or b glide plaes exist...
      if ( (cm.ne.'n').and.(cm.ne.'b') ) goto 570
      te(1)  = 0.5
      te(2)  = 0.5
      goto 570

530   continue
      e(1,1) =  1
      e(2,2) = -1
      e(3,3) = -1
c     if c or 1/4 glide plane...
      if (cn.eq.'c') te(3) = 0.5
      if (cn.eq.'d') te(3) = 0.25
      if (cn.eq.'d') te(2) = 0.5
c     if secondary diad is not screw... 
      if (csym(2,2).ne.'1') goto 570
      te(1)  = 0.5
      te(2)  = 0.5
      goto 570

540   continue
      e(1,2) =  1
      e(2,1) =  1
      e(3,3) = -1
c     if no secondary screw or i... 
      if ( (csym(2,2).ne.' ').or.(cbr.eq.'i').or.(csym(1,2).eq.' ') ) 
     $     goto 570
c     if principle screw axis is 1/4, 1/2, 3/4...
      if (csym(1,2).eq.'1') te(3) = 0.75
      if (csym(1,2).eq.'2') te(3) = 0.5 
      if (csym(1,2).eq.'3') te(3) = 0.25
      goto 570

550   continue
      cm     = csym(2,1)
      e(1,1) = -1
      e(2,2) =  1
      e(3,3) =  1
c     if parallel diag(n), b, or c mirror planes...
      if ( (cm.eq.'c').or.(cm.eq.'n') ) te(3) = 0.5
      if ( (cm.eq.'n').or.(cm.eq.'b') ) te(1) = 0.5
      if ( (cm.eq.'n').or.(cm.eq.'b') ) te(2) = 0.5
      goto 570

560   continue
      e(1,1) = -1
      e(2,2) =  1
      e(3,3) =  1
c     if parallel glide plane...
      if ( (csym(2,1).eq.'c').or.(csym(2,1).eq.'n') ) te(3) = 0.5
      if ( (csym(2,1).eq.'b').or.(csym(2,1).eq.'n') ) te(2) = 0.5
      if ( (csym(2,1).eq.'b').or.(csym(2,1).eq.'n') ) te(1) = 0.5
c     if parallel glide plane...
      if ( (csym(1,3).eq.'n').or.(csym(1,4).eq.'n') ) te(1)= te(1)+0.5
 570  ne = 4
      goto 740

c---------------------------------------------------------------------
c --- hexagonal
1400  continue
c      print*,'starting hex'
      ng = 3
      cm  = csym(1,1)
      cn  = csym(1,2)

c     if centrosymmetric...
      if ( (cm.eq.'-').and.(cn.eq.'3') ) ns = 0
c     if a hexad exists...
      if (cm.eq.'6') goto 610
      ss(2,1,2) = -1
      ss(2,2,1) =  1
      ss(2,2,2) = -1
      ss(2,3,3) =  1

c     if a 1/3, 2/3 screw triad exists...
      if (cn.eq.'1') ts(2,3) = 1.0/3.0 
      if (cn.eq.'2') ts(2,3) = 2.0/3.0
      ss(3,1,1) = -1
      ss(3,2,1) = -1
      ss(3,1,2) =  1
      ss(3,3,3) =  1
      ts(3,3)   =  2*ts(2,3) 
      if (ts(3,3).ge.1.0) ts(3,3) = ts(3,3)-1

c     if point group -6...
      if ( (nsym.eq.1).and.(cn.ne.'6') ) goto 760
c     if not point group -6m2...
      if (cn.ne.'6') goto 600
      ng = ng+ng
      do 594 i=1,3
         do 592 j=1,3
            do 590 k=1,3
               ss(3+i,j,k) = ss(i,j,k)
               ss(3+i,3,3) = -1
 590        continue
 592     continue
 594  continue

600   continue
      if (nsym.eq.1) goto 760
      if ( (csym(2,1).ne.'c').and.(csym(3,1).ne.'c') ) goto 630
      ts(4,3)   = 0.5
      ts(5,3)   = 0.5
      ts(6,3)   = 0.5
      goto 630

610   continue
      ng        = ng+ng
      ss(2,1,1) =  1
      ss(2,1,2) = -1
      ss(2,2,1) =  1
      ss(2,3,3) =  1
c     if 1/6..5/6 screw hexad...
      if (cn.eq.'1') ts(2,3) = 1.0/6.0
      if (cn.eq.'2') ts(2,3) = 2.0/6.0
      if (cn.eq.'3') ts(2,3) = 3.0/6.0
      if (cn.eq.'4') ts(2,3) = 4.0/6.0
      if (cn.eq.'5') ts(2,3) = 5.0/6.0

      do 626 i=1,4
         do 624 j=1,3 
            ts(2+i,j)=ts(2,j)
            do 622 k=1,3
               ts(2+i,j) = ts(2+i,j)+ss(2,j,k)*ts(1+i,k)
               if (ts(2+i,j).gt.1.0) ts(2+i,j) = ts(2+i,j)-1.0
               do 620 l=1,3
                  ss(2+i,j,k) = ss(2+i,j,k)+ss(2,j,l)*ss(1+i,l,k)
 620           continue
 622        continue
 624     continue
 626  continue
      if (nsym.eq.1) goto 760

630   continue
      ng  = ng+ng
      cm  = csym(2,1)
      cn  = csym(3,1)
c     if x,y directions are identity axes...
      if (cm.eq.'1') goto 650
c     if secondary diad exists...
      if (cm.eq.'2') goto 640
      e(1,2) = -1
      e(2,1) = -1
      e(3,3) =  1
c     if c glide plane exists...
      if (cm.eq.'c') te(3)=0.5
      goto 670

640   continue
      e(1,2) =  1
      e(2,1) =  1
      e(3,3) = -1
      te(3)  = 2.0*ts(2,3)

c   group p 31 1 2 and p 32 1 2
      if ( (csym(1,1).eq.'3') .and. 
     $    ((csym(1,2).eq.'1').or.(csym(1,2).eq.'2')) ) te(3) = 0
      goto 670

650   continue
c     if secondary diad exists...
      if (cn.eq.'2') goto 660
      e(1,2) = 1
      e(2,1) = 1
      e(3,3) = 1
c     if c glide plane exists...
      if (cn.eq.'c') te(3) = 0.5
      goto 670

660   continue
      e(1,2) = -1
      e(2,1) = -1
      e(3,3) = -1
      te(3)  = 2.0*ts(2,3)
      if (te(3).gt.1.0) te(3) = te(3)-1.0

670   continue
      ne = 6
c     if trigonal (3**) or (-3**)...
      if ( (csym(1,1).eq.'3') .or.
     $    ((csym(1,2).eq.'3').and.(csym(1,1).eq.'-')) ) ne = 3 
      goto 740

c---------------------------------------------------------------------
c   cubic
1500  ng = 12
      if (nsym.eq.3) ng = 24
c     if no primary diad or tetrad exists...
      if ( (csym(1,1).ne.'2').and.(csym(1,1).ne.'4').and.
     $     (csym(1,1).ne.'-') ) ns = 0

      do 692 i=1,3
         do 690 j=1,3
            ss(1+i,j,j) =  1
            if (i.eq.j) goto 690
            ss(1+i,j,j) = -1
            if (csym(1,1).eq.'n') ts(1+i,j) = half
            if (csym(1,1).eq.'d') ts(1+i,j) = quart
 690      continue
 692   continue

c     if (no a and no 1/4 glide plane and no 1/4 and no 3/4 screw) or
c     if (face centered)...
      if ( ((csym(1,1).ne.'a').and.(csym(3,1).ne.'d').and.
     $      (csym(1,2).ne.'3').and.(csym(1,2).ne.'1')) .or. 
     $      (cbr.eq.'f') )   goto 710

      do 700 i=1,3
         ts(1+i,i) = half
         k = i+1
         if (k.eq.4) k = 1
         ts(1+i,k) = half
700   continue

710   continue
      do 724 i=1,4
         do 722 j=1,3
            do 720 k=1,3
               l = j+1
               if (l.eq.4) l = 1
               m = j-1
               if (m.eq.0) m = 3
               ss(4+i,j,k) = ss(i,l,k)
               ss(8+i,j,k) = ss(i,m,k)
               ts(4+i,j)   = ts(i,l)
               ts(8+i,j)   = ts(i,m)
720         continue
722      continue
724   continue

      if (ng.eq.12) goto 760
      ne     = 12
      e(1,2) = 1 
      e(2,1) = 1 
      e(3,3) = 1 
c     if secondary diad exists...
      if (csym(3,1).eq.'2') e(3,3) = -1
c     if c glide plane exists...
      if (csym(3,1).eq.'c') te(3)  = half

c        if diag glide plane or 2/4 screw diad..
c        if 1/4 glide or 1/4 screw or 3/4 screw..
      do 730 i=1,3
         if ( (csym(3,1).eq.'n').or.(csym(1,2).eq.'2') ) te(i) = half
         if ( (csym(3,1).eq.'d').or.(csym(1,2).eq.'1').or.
     $        (csym(1,2).eq.'3') )                       te(i) = quart
730   continue

c     if 1/4 screw or p...
      if (  (csym(1,2).eq.'1').and.(cbr.eq.'p') )      te(1) = quart * 3
c     if (not 1/4 screw or not i) and (not 3/4 screw or not p)... 
      if ( ((csym(1,2).ne.'1').or.(cbr.ne.'i')) .and. 
     $     ((csym(1,2).ne.'3').or.(cbr.ne.'p')) )        goto 740

      te(2) = quart * 3
      te(3) = quart * 3

c---------------------------------------------------------------------
 740  do 756 i=1,ne
         do 754 j=1,3
            ts(ne+i,j) = te(j)
            do 752 k=1,3
               ts(ne+i,j) = ts(ne+i,j)+e(j,k)*ts(i,k)
               do 750 l=1,3
                  ss(ne+i,j,k) = ss(ne+i,j,k)+e(j,l)*ss(i,l,k)
750            continue
752         continue
754      continue
756   continue

c---------------------------------------------------------------------
c  just about done.  put all atoms in central cell, finish up fx, mark
c  system and centrosymmetry, and make sure ng>=1.
760   continue  

c++++++++++++++++
c  put located atoms back into central cell
      do 812 i=1,ng
         do 810 j=1,3
 770        if (ts(i,j).ge.1.0) then
                ts(i,j) = ts(i,j)-1
                goto 770
            endif
 790        if (ts(i,j).lt.0.0) then
                ts(i,j) = ts(i,j)+1
                goto 790
            endif
 810     continue
 812  continue
 
c================
c---------------------------------------------------------------------
c  tell the rest of the program what system the lattice is and whether
c  is has centrosymmetry
      isystm = isys(sys)
      if (ns.eq.0) isymce = 1
      if (ns.eq.1) isymce = 0

c---------------------------------------------------------------------
c --- swap ss(i,j,k) to fx(k,i,j)
      do 840 k=1,ng
         do 830 j=1,3
            do 820 i=1,3
              fx(i,j,k) = ss(k,j,i)
820         continue
            tx(j,k)   = ts(k,j)
830      continue
840   continue

c---------------------------------------------------------------------
c  there is always at least one element in a group 
c  (identity transformation)
      if (ng.eq.0) ng = 1

      return
c end subroutine equipt
      end

