c-*-fortran-*-
      common /exaflt/ gasses(3), exafs(nexafs), rmax
      save /exaflt/

      common /exaint/ iedge, iexerr
      save /exaint/

      character*10 core, edge*2
      common /exastr/ core, edge
      save /exastr/

      logical lfluo
      common /exalog/ lfluo
      save /exalog/
c -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c  GLOSSARY:  (* = user input, ! = output, % = error handling)
c
c  * gasses:  (3) percent pressure of argon, krypton, nitrogen in i0 chamber
c         gasses(1)       percentage of argon
c         gasses(2)       percentage of krypton
c         gasses(3)       percentage of nitrogen
c  ! exafs:   (13) amu,delmu,spgrav,sigmm,qrtmm,ampslf,sigslf,qrtslf,
c                  sigi0,qrti0,muf,mub,mue
c         exafs(1):       total mu
c         exafs(2):       delta mu
c         exafs(3):       speciffic gravity
c         exafs(4):       mcmaster sigma^2
c         exafs(5):       mcmaster C4
c         exafs(6):       self absorption amplitude correction
c         exafs(7):       self absorption sigma^2
c         exafs(8):       self absorption C4
c         exafs(9):       i0 corrcetion sigma^2
c         exafs(10):      i0 correction C4
c  * iedge:   edge for calulation, 1=K 2=L1 3=L2 4=L3
c  * iexerr:  exit error code, 0=no prob, 1=info, 2=warning, 3=error
c  * core*10: tag of absorbing atom
c    lfluo:   true if fluorescence corrections were calculated
c----------------------------------------------------------------------
