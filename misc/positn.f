      subroutine positn(module, string)

      character*78 module, string, messg
      im = istrln(module)
      is = istrln(string)
      messg = ' *** Position in ' // module(:im) // ': ' //
     $            string(:is) // '.-'
      call messag(messg)
      return
c end subroutine crypos
      end
