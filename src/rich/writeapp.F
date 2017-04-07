      SUBROUTINE writeapp( filename, n1, n2, f1, f2 )
c
c---- date  18/05/04
c     --------------
c
c
      IMPLICIT NONE
c
      INTEGER  n1, n2
      REAL   f1, f2
c
      CHARACTER*(*) filename
      INTEGER  icww
c
c
c---- write unformatted data on g-file (for linux) :
c
      icww = 50
c
ccc      print *, n1, n2, f1, f2
c
      write (icww) n1, n2, f1, f2
c
c
      RETURN
      END

