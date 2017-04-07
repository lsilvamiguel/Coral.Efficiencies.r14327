      SUBROUTINE writebuff(filename,
     +                     nparev,
     +                     xpav,ypav,zpav,alpav,ampav,anpav,
     +                     xpov,ypov,zpov,
     +		           xpemv, ypemv, zpemv, alpemv, ampemv, anpemv,
     +		           xpewv, ypewv, zpewv, alpewv, ampewv, anpewv,
     +                     ppav,
     +                     nphoc,thetac,
     +                     zhelix0,zhelix1,nclus,chisq,     !  020322
     +                     itravv,ipartvv,iowinvv,
     +                     npatt,
     +                     ict,ixt,iyt,pht)
c
c---- date  14/02/01
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
      CHARACTER*(*) filename
      INTEGER  icww
      INTEGER  kcall
c
      data kcall / 0 /
c
c
c---- write unformatted gfile for linux :
c
      icww = 50
c
      if(kcall.eq.0) then
        kcall = 1
c
        open(icww,file=filename,form='unformatted')
c-----------------------------------------------------------
        print *,"RICHONE, writebuff: Write File  ",filename,"  opened"
c
        write (icww) icww
c
      endif
c
c-----------------------------------------------------------------
      write (icww) nparev,
     +             (xpav(k),ypav(k),zpav(k),
     +             alpav(k),ampav(k),anpav(k),
     +             xpov(k),ypov(k),zpov(k),
     +             xpemv(k), ypemv(k), zpemv(k),
     +             alpemv(k), ampemv(k), anpemv(k),
     +             xpewv(k), ypewv(k), zpewv(k),
     +             alpewv(k), ampewv(k), anpewv(k),
     +             ppav(k),
     +             nphoc(k),thetac(k),
     +             zhelix0(k),zhelix1(k),nclus(k),chisq(k),   !   020322
     +             itravv(k),ipartvv(k),iowinvv(k),k=1,nparev),
     +             npatt,
     +             (ict(k),ixt(k),iyt(k),pht(k),k=1,npatt)
c-----------------------------------------------------------------
ccc      print *,     nparev,
ccc     +             (xpav(k),ypav(k),zpav(k),
ccc     +             alpav(k),ampav(k),anpav(k),
ccc     +             xpov(k),ypov(k),zpov(k),
ccc     +             ppav(k),
ccc     +             nphoc(k),thetac(k),
ccc     +             zhelix0(k),zhelix1(k),nclus(k),chisq(k),   !   020322
ccc     +             itravv(k),ipartvv(k),iowinvv(k),k=1,nparev),
ccc     +             npatt,
ccc     +             (ict(k),ixt(k),iyt(k),pht(k),k=1,npatt)
c
c
      RETURN
      END
c
c
      SUBROUTINE closewrtfile
c
c---- date  22/03/02
c
c
      IMPLICIT NONE
c
      INTEGER  icww
c
      icww = 50
c
      close(icww)
c
      RETURN
      END
