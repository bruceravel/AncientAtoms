      integer function is2z(sym)
c---------------------------------------------------------------------------
c  copyright 1993 university of washington     matt newville and bruce ravel
c---------------------------------------------------------------------------

c     returns z number given atomic symbol: default is 0
      character*2 symbol(103), sym, sy, test
c
      data (symbol(i),i=1,103) /'h','he','li','be','b','c','n','o',
     $'f','ne','na','mg','al','si','p','s','cl','ar','k','ca','sc',
     $'ti','v','cr','mn','fe','co','ni','cu','zn','ga','ge','as','se',
     $'br','kr','rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag',
     $'cd','in','sn','sb','te','i','xe','cs','ba','la','ce','pr','nd',
     $'pm','sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta',
     $'w','re','os','ir','pt','au','hg','tl','pb','bi','po','at','rn',
     $'fr','ra','ac','th','pa','u','np','pu','am','cm','bk','cf','es',
     $'fm','md','no','lr'/

c test case of this routine : note that 'ab' must be the same
c                             case as all of the symbols above
      sy   = sym
      test = 'ab'
      call case(test,sy)

      is2z = 0
      do 110 i=1,103
        if (sy.eq.symbol(i)) then
          is2z = i
          goto 120
        end if
110   continue
120   return
c end integer function is2z
      end
