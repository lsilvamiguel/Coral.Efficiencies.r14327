      integer function  inithbook( zebra )

      implicit none

      integer        zebra
      integer        memory
      parameter    ( memory = 10000000 ) 
      integer        hmem( memory ) 
      common /PAWC/  hmem
      integer        iquest( 100 )
      common /QUEST/ iquest

      if( zebra .gt. 0 ) then
         call hlimit( -memory )
      else
         call hlimit( memory )
      endif
      iquest(10) = 65000

      write( 6, 20 ) memory
 20   format( ' ---------------------------------------', /,
     /        ' HBOOK Initialization ', /,
     .        ' HLIMIT set to ', I10, /, 
     .        ' IQUEST(10) = 65000' )

      inithbook = 1

      return
      end
