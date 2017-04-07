
      SUBROUTINE openfile( filename, icww )
c     -------------------------------------
c
c---- date  22/02/10
c
c
      IMPLICIT NONE
c
      CHARACTER*(*) filename
      INTEGER  icww
c
      open( icww , file=filename, form='unformatted' )
c     ------------------------------------------------
c
      RETURN
      END
