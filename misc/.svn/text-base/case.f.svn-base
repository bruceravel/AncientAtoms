      subroutine smcase (str, contrl)
c  convert case of string *str*to be the same case
c  as the first letter of string *contrl*
c  if contrl(1:1) is not a letter, *str* will be made lower case.
      character*(*) str, contrl, s1*1, t1*1
      s1 = contrl(1:1)
      t1 = s1
      call lower(t1)
      if (t1.eq.s1)  call lower(str)
      if (t1.ne.s1)  call upper(str)
      return
c end subroutine smcase
      end
c----------------------------------------------------------------
      subroutine case(test,word)
c  returns *word* in the same case as *test*
c  note that this is just the reverse of smcase !
      character*(*) test, word
      call smcase (word, test)
      return
c  end subroutine case
      end
