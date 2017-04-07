
      SUBROUTINE readbuff( icww,
c     --------------------------
     +                     nparev,
     +                     xpav, ypav, zpav, alpav, ampav, anpav,
     +                     xpov, ypov, zpov,
     +		           xpemv, ypemv, zpemv, alpemv, ampemv, anpemv,
     +		           xpewv, ypewv, zpewv, alpewv, ampewv, anpewv,
     +                     ppav,
     +                     nphoc, thetac,
     +                     zhelix0, zhelix1, nclus, chisq,
     +                     itravv, ipartvv, iowinvv,
     +                     npatt,
     +                     ict, ixt, iyt, pht )
c
c---- date  22/02/10
c
c
      IMPLICIT NONE
c
      INTEGER  NPAEV
      PARAMETER  (NPAEV = 1)
      INTEGER  nparev
      REAL  xpav(NPAEV), ypav(NPAEV), zpav(NPAEV)
      REAL  alpav(NPAEV), ampav(NPAEV), anpav(NPAEV)
      REAL  xpov(NPAEV), ypov(NPAEV), zpov(NPAEV)
      REAL  xpemv(NPAEV), ypemv(NPAEV), zpemv(NPAEV)
      REAL  alpemv(NPAEV), ampemv(NPAEV), anpemv(NPAEV)
      REAL  xpewv(NPAEV), ypewv(NPAEV), zpewv(NPAEV)
      REAL  alpewv(NPAEV), ampewv(NPAEV), anpewv(NPAEV)
      REAL  ppav(NPAEV)
      INTEGER  nphoc(NPAEV)
      REAL  thetac(NPAEV)
      INTEGER  itravv(NPAEV), ipartvv(NPAEV), iowinvv(NPAEV)
      REAL  zhelix0(NPAEV), zhelix1(NPAEV), chisq(NPAEV)
      INTEGER  nclus(NPAEV)
c
      INTEGER  MAXPAD
      PARAMETER (MAXPAD = 1)
      INTEGER  npatt
      INTEGER  ict(MAXPAD),ixt(MAXPAD),iyt(MAXPAD)
      REAL  pht(MAXPAD)
c
      INTEGER  k
      INTEGER  icww
      INTEGER  kcall
      LOGICAL  print
c
      print = .true.
ccc      print = .false.
c
      data kcall / 0 /
c
c
c---- read unformatted gfile for linux :
c
      print *, "readbuff"
c
      if( kcall .eq. 0 ) then
        kcall = 1
      endif
c
c-----------------------------------------------------------------
      read (icww) nparev
      print *, nparev
c     +            ( xpav(k), ypav(k), zpav(k),
c     +            alpav(k), ampav(k), anpav(k),
c     +            xpov(k), ypov(k), zpov(k),
c     +            xpemv(k), ypemv(k), zpemv(k),
c     +            alpemv(k), ampemv(k), anpemv(k),
c     +            xpewv(k), ypewv(k), zpewv(k),
c     +            alpewv(k), ampewv(k), anpewv(k),
c     +            ppav(k),
c     +            nphoc(k), thetac(k),
c     +            zhelix0(k), zhelix1(k), nclus(k), chisq(k),
c     +            itravv(k), ipartvv(k), iowinvv(k), k=1,nparev ),
c     +            npatt,
c     +            ( ict(k), ixt(k), iyt(k), pht(k), k=1,npatt )
c-----------------------------------------------------------------
      if( print ) then
        print *, nparev
        do   k=1,nparev
          print *, xpav(k), ypav(k), zpav(k),
     +             alpav(k), ampav(k), anpav(k),
     +             xpov(k), ypov(k), zpov(k),
     +             ppav(k),
     +             nphoc(k), thetac(k),
     +             zhelix0(k), zhelix1(k), nclus(k), chisq(k),
     +             itravv(k), ipartvv(k), iowinvv(k)
        enddo
        print *, npatt
        do   k=1,npatt
          print *, ict(k), ixt(k), iyt(k), pht(k)
        enddo
      endif
c
c
      RETURN
      END
