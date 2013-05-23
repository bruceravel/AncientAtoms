      subroutine bwords (s, nwords, words)
c
c     breaks string into words.  words are seperated by one or more
c     blanks, or a comma or equal sign and zero or more blanks.
c
c     args        i/o      description
c     ----        ---      -----------
c     s            i       char*(*)  string to be broken up
c     nwords      i/o      input:  maximum number of words to get
c                          output: number of words found
c     words(nwords) o      char*(*) words(nwords)
c                          contains words found.  words(j), where j is
c                          greater then nwords found, are undefined on
c                          output.
c
c      written by:  steven zabinsky, september 1984
c
c**************************  deo soli gloria  **************************
c-- no floating point numbers in this routine.
      implicit integer (a-z)
      character*(*) s, words(nwords)
      character blank, comma, equal
      parameter (blank = ' ', comma = ',', equal = '=')

c-- betw    .true. if between words
c   comfnd  .true. if between words and a comma or equal has
c                                         already been found
      logical betw, comfnd
c-- maximum number of words allowed
      wordsx = nwords

c-- slen is last non-blank character in string
      slen = istrln (s)

c-- all blank string is special case
      if (slen .eq. 0)  then
         nwords = 0
         return
      endif

c-- begc is beginning character of a word
      begc = 1
      nwords = 0
      betw   = .true.
      comfnd = .true.
      do 10  i = 1, slen
         if (s(i:i) .eq. blank)  then
            if (.not. betw)  then
               nwords = nwords + 1
               words (nwords) = s (begc : i-1)
               betw = .true.
               comfnd = .false.
            endif
         elseif ((s(i:i).eq.comma).or.(s(i:i).eq.equal))  then
            if (.not. betw)  then
               nwords = nwords + 1
               words (nwords) = s(begc : i-1)
               betw = .true.
            elseif (comfnd)  then
               nwords = nwords + 1
               words (nwords) = blank
            endif
            comfnd = .true.
         else
            if (betw)  then
               betw = .false.
               begc = i
            endif
         endif
         if (nwords .ge. wordsx)  return
   10 continue
c
      if (.not. betw  .and.  nwords .lt. wordsx)  then
         nwords = nwords + 1
         words (nwords) = s (begc :slen)
      endif
      return
c end subroutine bwords
      end
