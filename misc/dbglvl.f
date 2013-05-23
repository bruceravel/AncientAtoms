      subroutine dbglvl(idebug, idbg)

c  This converts an integer into a binary representation

      include 'atparm.h'
      integer idebug, idbg(0:ndbgx)
      integer isave, level

      do 5 i=0,ndbgx
        idbg(i) = 0
 5    continue
      if (idebug.ge.1) then
          isave = idebug
          do 10 i=ndbgx,0,-1
            level = 2**i
            idbg(i) = idebug/level
            idebug  = mod(idebug,level)
 10       continue
          idebug = isave
      endif

      return
c  end subroutine dbglvl
      end
