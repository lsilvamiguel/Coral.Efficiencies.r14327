
      INTEGER  ict,ixt,iyt,MAXPAD
c0104      INTEGER  ihcc                      !   1900
      INTEGER*8  ihcc                      !   1900
      REAL  pht,padpht
      PARAMETER (MAXPAD=5000)
c
      COMMON / CLU / ict(MAXPAD),ixt(MAXPAD),iyt(MAXPAD),
c1900     +               pht(MAXPAD),padpht(16,72,72)
     +               pht(MAXPAD),ihcc(MAXPAD),padpht(16,72,72)
