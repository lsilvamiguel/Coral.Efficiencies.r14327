c       Set of programs for shower separation in ECAL2.  A.Lednev 3.11.2008.
c       Complete set includes 6 files:
c       separa.f, badcheck.f, ecal2reco.f, Ecell.f, fastreco.f, showerfit.f.
c
c       Subset of programs for shower separation in ECAL2.  A.Lednev 17.11.2010.
c       blockdata, separa, book, store, trans, calib, readcal, writcal,
c       varthr, readthr, writhr, thrcut, mindis, icalad, xm2g.
c       Last correction 3.08.2011. (indmc and lnclust are saved in separahist)

      blockdata
      implicit none
#include "parameters.inc"
#include "badlist.inc"
      data rnoisy/numc*1./
#include "cells.inc"
#include "chkbad.inc"
      data rmn/numc*0./,rsq/numc*0./
#include "cluster.inc"
      data idx/48/
#include "ccuts.inc"
      data egcut/.3/,time0,timew/0.,5./
#include "debug.inc"
      data ldebug/.false./,ndebug/0/,listdebug/10*0/
#include "fcuts.inc"
      data discut/25.00/,xm2cut/.04/
      data Elow/15./,Ehigh/300./
      data chicut/999.0/,coscut/.9/
#include "files.inc"
      data filehbook/'separa000000.hbook'/
      data flthresh /'thresh000000.dat'/
      data badcells /'badcells000000.dat'/
#include "geom.inc"
      data Zecal2/33252./
      data iholel/34/,iholer/35/,iholeu/24/,iholed/23/
#include "iflags.inc"
      data ibad/0/,icalib/0/,ifast/0/,ihist/0/,ivarthr/0/
#include "leds.inc"
      data spillled/numc*0./,nledinspill/numc*0/,ameanled/numc*0./,ledstat/numc*0/
#include "mincells.inc"
      data Amin/numc*5./,astat/numc*5./
#include "rawcells.inc"
#include "showerprof.inc"
cc      data a/1.191, -.191,  1.113, -.113/     ! old parameters
cc      data b/6.664, 43.777, 7.311, 85.56/     ! mm
      data a/0.586, 0.414,  0.272, 0.728/     ! new parameters
      data b/3.26,  9.05,   1.517, 7.29 /     ! mm
#include "thresholds.inc"
      data thrcommon/.2/,cellthr/numc*.001/
#include "thr_test.inc"
      data thrtest/numc*0./
      end
c
c-----------------------------------------------------------
      subroutine separa
      implicit none
#include "errcell2.inc"
#include "icalad.inc"
#include "parameters.inc"
#include "badlist.inc"
#include "ccuts.inc"
#include "cells.inc"
#include "chkbad.inc"
#include "clb.inc"
#include "cluster.inc"
#include "debug.inc"
#include "fcuts.inc"
#include "gamma.inc"
#include "geom.inc"
#include "iflags.inc"
#include "isitlg.inc"
#include "misc.inc"
#include "rawcells.inc"
#include "mincells.inc"
#include "runpar.inc"
#include "showerprof.inc"
#include "thresholds.inc"
#include "thr_test.inc"
      integer i, icnt, ig, ix, iy, j, n, nev, nburstprev, nevthr
      logical lini
      real Atot, abst, dnear, distig, Emain, Etot, scale, xx, yy, thr0
      data lini/.true./,nev/0/,nburstprev/-1/,nevthr/1/
c
      if(lini) then
        print *,'*** separa started ***  Run',NRUN,'  time0 =',time0,'  timew=',timew
        print *,'shower parameters',a,b
c        print *,'init cuts =',egcut, discut, xm2cut, Elow, Ehigh, chicut, coscut, time0, timew
        if(icalib.ne.0) then
          ihist =0
          call read_ecalib('calib.dat')
        endif
        if(ihist.ne.0) call book
        print *,'ibad=',ibad,' icalib=',icalib,' ifast=',ifast,' ihist=',ihist,' ivarthr=',ivarthr
        print *,'egcut =',egcut,'   thrcommon =',thrcommon
C initialize XY map used to map cells by their number to their position
        call inimapxy
cc        call ecal2thr
        call readthr(athr)  !//AL to measure thresholds
        do i=1,numc
          cellthr(i)=athr(i)/cal(i)
        enddo
c        print *,'separa init. thrcommon =',thrcommon,'  thresholds ='
c        print '(24f6.3)',cellthr
        lini=.false.
      endif

      ngam=0
      ncl=0
      if(ncell.le.0) return
      if(ihist.ne.0) then
        call hreset(15,' ')
        do i=1,numc
          call hf1(15,cellthr(i),1.)
          if(cellthr(i).lt.thrcommon) cellthr(i)=thrcommon
        enddo
      endif
c      print *,'thrcommon =',thrcommon,'  thresholds ='
c      print '(24f6.3)',cellthr

      if(nburst.ne.nburstprev)then
        nburstprev=nburst
c        print *,'  NBURST =',NBURST
      endif
c      ldebug=.false.
c      ldebug=.true.
      do i=1,ndebug
        if(numbev.eq.listdebug(i)) ldebug=.true.
      enddo
      if(ldebug) print '(a,2i8)','*** separa started ***  ncell =',ncell
c      print *,'** separa **  numbev',numbev,'  ncell =',ncell,'  icailb =',icalib
      nev=nev+1
      if(ncell.lt.1) return
      do i=1,ncell
        iaddr(i)=icalad(x(i),y(i))
        icnt=iaddr(i)+1
        ix=iaddr(i)/48
        iy=iaddr(i)-ix*48
        err(i)=sqrt(errcell2(e(i)))
        if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle
! thr_test
!        if(thrtest(icnt).le.0..or.e(i).le.thrtest(icnt)/2.) thrtest(icnt)=e(i)
!        if(thrtest(icnt).gt.2.) thrtest(icnt)=2.
!        thr0=thrtest(icnt)
!        if(e(i).ge.thrtest(icnt).and.e(i).lt.thrtest(icnt)+0.050) then
!          thr0=thrtest(icnt)*1.03
!        endif
!        if(e(i).ge.0..and.e(i).lt.thrtest(icnt)) then
!          thr0=thrtest(icnt)/1.03
!        endif
!        thrtest(icnt)=thr0
! end thr_test
        if(ihist.ne.0) then
          call hfill(2000+ix,iy+.5,time(i),sqrt(e(i)))
          call hf2(50,time(i),e(i),1.)
          call hf1(14,e(i),1.)
          if(isitLG(iaddr(i)).eq.1) then
            call hf1(51,e(i),1.)
            call hf2(56,time(i),e(i),1.)
          endif
          if(isitLG(iaddr(i)).eq.0) then
            call hf1(52,e(i),1.)
            call hf2(57,time(i),e(i),1.)
          endif
          if(ix.ge.23.and.ix.le.42.and.iy.ge.16.and.iy.le.31) then
c            if(abs(x(i)-dcell).le.2.*dcell.and.abs(y(i)).lt.2.*dcell) then
            call hf1(53,e(i),1.)
            call hf2(58,time(i),e(i),1.)
          endif
          if(abs(time(i)-time0).le.2.*timew) call hf1(54,e(i),1.)
        endif
        if((acell(i).lt.Amin(icnt).or.Amin(icnt).lt.0.1).and.acell(i).gt.1.) Amin(icnt)=acell(i)

        if(ihist.ne.0.and.e(i).gt.cellthr(icnt)) call hf2(44,ix+.5,iy+.5,1.)
        abst=abs(time(i)-time0)/timew
        if(abst.gt.100.) abst=100.
        if(abst.gt.1.) then
          err(i)=err(i)*abst
c          e(i)=e(i)/abst
        endif
        if(ihist.ne.0.and.e(i).gt.cellthr(icnt).and.abst.lt.2.) call hf2(45,ix+.5,iy+.5,1.)
        if(iaddr(i).lt.0) then
          print *,'**reco_ecal2**  input x,y =',x(i),y(i),' (out of ECAL2)'
          return
        endif
      enddo
      if(ldebug) then
        print *,'*separa*  event ',numbev,'  input ecell ='
        print '(10f8.2)',(e(i),i=1,ncell)
        print *,'acell ='
        print '(10f8.2)',(acell(i),i=1,ncell)
        print *,'iaddr ='
        print '(10i8)',(iaddr(i),i=1,ncell)
        print *,'xcell ='
        print '(10f8.1)',(x(i),i=1,ncell)
        print *,'ycell ='
        print '(10f8.1)',(y(i),i=1,ncell)
        print *,'err ='
        print '(10f8.1)',(err(i),i=1,ncell)
        print *,'time ='
        print '(10f8.1)',(time(i),i=1,ncell)
      endif

      Atot=0.
      do i=1,ncell
        Atot=Atot+acell(i)
      enddo

      if(icalib.ne.0) then
        call joujou
        return
      endif

c      if(icalib.eq.0) call ecalib
      ngam=0
cc      call badcheck
cc      call varthr(nevthr)
cc      call thrcut
      ngam=0
      Etot=0.
      do i=1,ncell
        Etot=Etot+e(i)
        ix=iaddr(i)/48
        iy=iaddr(i)-ix*48
        if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle
        if(ihist.ne.0) then
          call hf2(46,ix+.5,iy+.5,1.)
          if(abs(time(i)-time0).lt.2.*timew) call hf2(4,ix+.5,iy+.5,1.)
        endif
        do j=1,noisy
          if(iaddr(i).eq.inoisy(j)-1) then
            if(ibad.ne.0) then
              scale=rsq(i+1)+1.
              if(rsq(i+1).lt..02) scale=10.
            else
              scale=rnoisy(i)
            endif
            err(i)=err(i)*scale
          endif
        enddo
      enddo
c      call badcheck
      if(ihist.ne.0) then
cc        call hf1(1,Atot,1.)
        call hf1(2,Etot,1.)
        call hf1(3,float(ncell),1.)
      endif
      if(Etot.lt.Elow.or.Etot.gt.Ehigh) then
c        print *,'Etot,Elow,Ehigh =',Etot,Elow,Ehigh
        return
      endif
c
      call ecal2reco
c      print *,'event',numbev,ngam,' gamma are found'
      if(ihist.ne.0) then
        call hf1(111,float(ndead+noisy),1.)
      endif
      call killgamma(egcut,discut,xm2cut)
c        some unnecessary things...
      do ig=1,ngam
        Eerr(ig)=sqrt(errcell2(egam(ig)))
        Xerr(ig)=50.*Eerr(ig)/(egam(ig)+.2)
        Yerr(ig)=Xerr(ig)
        xx=xgam(ig)
        yy=ygam(ig)
        Emain=0.
        dnear=1.E30
        do i=1,ncell
          ix=iaddr(i)/48
          iy=iaddr(i)-ix*48
          if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle
          distig=sqrt((xx-X(i))**2+(yy-Y(i))**2)
          if(distig.lt.dnear.and.E(i).gt.0.) then
            dnear=distig
            indmc(ig)=iaddr(i)
            gtime(ig)=time(i)
            Emain=E(i)
          endif
c          if(iaddr(i).eq.indmc(ig)) Emain=E(i)
        enddo
        if(sqrt(dnear).gt.38.) then
          print *,'dnear,Emain,egam,xgam,ygam=',dnear,Emain,egam(ig),xgam(ig),ygam(ig)
          print '(10f8.2)',(e(i),i=1,ncell)
          print '(10f8.2)',(x(i),i=1,ncell)
          print '(10f8.2)',(y(i),i=1,ncell)
        endif
        if(ihist.ne.0) then
          call hf1(9,egam(ig),1.)
          call hf2(13,Emain,egam(ig),1.)
          call hf1(16,float(lnclust(ig)),1.)
        endif
      enddo

      if(ihist.ne.0) then
        call separahist
      else
        if(nev.ge.1000) then
          if(ibad.ne.0) call writbad
          if(ihist.ne.0) call store
          if(ivarthr.eq.1) call writhr(athr)
          nev=0
        endif
      endif
cc      print *,'EvInSpill',numbinburst,'  ngam =',ngam,'  Egam =',(egam(i),i=1,ngam)
      return
      end
c
c-----------------------------------------------------------
c  separahist does not change output of separa. It is for histogramming only.
      subroutine separahist
      implicit none
#include "xm2g.inc"
#include "parameters.inc"
#include "badlist.inc"
#include "cells.inc"
#include "chkbad.inc"
#include "cluster.inc"
#include "ccuts.inc"
#include "fcuts.inc"
#include "debug.inc"
#include "iflags.inc"
#include "gamma.inc"
#include "geom.inc"
#include "misc.inc"
#include "runpar.inc"
#include "thresholds.inc"
      integer i, ig, ig1, ig2, ix, iy, nev, ngamsv, numb, ih0
      integer indmcsv(numg), ndfsv(numg), lnclustsv(numg)
      logical lini
      real chimax, cost, d, d2g, dx, dy, d2gmin, e2g, Etot, Px, Py, Pz
      real xm12, xm12min, xmass, egc, disc, xm2c
      real egamsv(numg),xgamsv(numg),ygamsv(numg),chigamsv(numg),gtimesv(numg)
      data lini/.true./,nev/0/,numb/-1/

      if(lini) then
c        print '(a,2i8)','      *** separahist started ***  Run',NRUN
c        call readcal('calib.dat')
        call book
        lini=.false.
      endif
      if(numb.eq.numbev) return
      numb=numbev

      nev=nev+1
      ngamsv=ngam
      do ig=1,ngam
        egamsv(ig)  =egam(ig)
        xgamsv(ig)  =xgam(ig)
        ygamsv(ig)  =ygam(ig)
        chigamsv(ig)=chigam(ig)
        gtimesv(ig) =gtime(ig)
        indmcsv(ig) =indmc(ig)
        ndfsv(ig)   =ndf(ig)
        lnclustsv(ig)=lnclust(ig)
        call hf1(19,gtime(ig),1.)
      enddo
      egc=egcut
      disc=discut
      xm2c=xm2cut
      if(ihist.ne.0) then
        do ig=1,ngam
          if(egam(ig).gt.20..and.xgam(ig).lt.dcell) then
            ix=nint(xgam(ig)/dcell)
            dx=abs(xgam(ig)-ix*dcell)
            iy=nint(ygam(ig)/dcell)
            dy=abs(ygam(ig)-iy*dcell)
            call hf1(401,dx,1.)
            call hf1(402,dy,1.)
            call hf1(403,xgam(ig),1.)
            call hf1(404,ygam(ig),1.)
          endif
        enddo
      endif

      ih0 = 0
    1 ih0=ih0+100
      call killgamma(egc,disc,xm2c)
      call hf1(ih0+7,float(ngam),1.)
      chimax=0.
      do ig=1,ngam
        call hf1(ih0+9,egam(ig),1.)
        call hf2(ih0+14,xgam(ig)/dcell+31.,ygam(ig)/dcell+23.,1.)
        call hf2(ih0+15,xgam(ig)/dcell+31.,ygam(ig)/dcell+23.,1.)
      enddo
c
      Etot=0.
      Px=0.
      Py=0.
      Pz=0.
      chimax=0.
      do ig=1,ngam
        d=sqrt(Zecal2**2+xgam(ig)**2+ygam(ig)**2)
        Etot=Etot+egam(ig)
        Px=Px+egam(ig)*xgam(ig)/d
        Py=Py+egam(ig)*ygam(ig)/d
        Pz=Pz+egam(ig)*Zecal2/d
        if(chigam(ig)/ndf(ig).gt.chimax) chimax=chigam(ig)/ndf(ig)
        call hf2(ih0+8,chigam(ig)/ndf(ig),egam(ig),1.)
      enddo
      xmass=0.
      if(ngam.gt.1) then
        xmass=Etot**2-Px**2-Py**2-Pz**2
        if(xmass.gt.0.) xmass=sqrt(xmass)
      endif
      do ig1=1,ngam-1
        if(abs(gtime(ig1)-time0).gt.2.*timew) cycle
        do ig2=ig1+1,ngam
          if(abs(gtime(ig2)-time0).gt.2.*timew) cycle
          xm12=xm2g(ig1,ig2)
          if(egam(ig1).lt.egam(ig2)) then
            call hf2(ih0+17,xm12,egam(ig1),1.)
          else
            call hf2(ih0+17,xm12,egam(ig2),1.)
          endif
        enddo
      enddo
      d2gmin=99999.
      xm12min=99999.
      do ig1=1,ngam-1
cc        if(chigam(ig1).gt.chicut) cycle
        if(abs(gtime(ig1)-time0).gt.2.*timew) cycle
        if(egam(ig1).lt.5.) cycle
        do ig2=ig1+1,ngam
cc          if(chigam(ig2).gt.chicut) cycle
          if(abs(gtime(ig2)-time0).gt.2.*timew) cycle
          if(egam(ig2).lt.5.) cycle
          d2g=sqrt((xgam(ig1)-xgam(ig2))**2+(ygam(ig1)-ygam(ig2))**2)
          e2g=egam(ig1)+egam(ig2)
          cost=abs(egam(ig1)-egam(ig2))/e2g
          xm12=xm2g(ig1,ig2)
          call hf1(ih0+21,xm12,1.)
          if(egam(ig1).lt.egam(ig2)) then
            call hf2(ih0+17,xm12,egam(ig1),1.)
          else
            call hf2(ih0+17,xm12,egam(ig2),1.)
          endif
          if(d2g.lt.d2gmin) then
            d2gmin=d2g
            xm12min=xm12
          endif
          if(ngam.eq.2) then
            call hf1(ih0+20,e2g,1.)
            if(e2g.gt.40..and.e2g.le.100.) then
              if(cost.le.coscut)  call hf1(ih0+23,xm12,1.)
              if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+24,cost,1.)
            elseif(e2g.gt.100..and.e2g.le.150.) then
              if(cost.le.coscut) call hf1(ih0+25,xm12,1.)
              if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+26,cost,1.)
            elseif(e2g.gt.150..and.e2g.le.200.) then
              if(cost.le.coscut) call hf1(ih0+27,xm12,1.)
              if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+28,cost,1.)
            endif
          endif
          if(e2g.gt.10..and.e2g.le.40.) then
            if(cost.le.coscut) call hf1(ih0+31,xm12,1.)
            if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+32,cost,1.)
          elseif(e2g.gt.40..and.e2g.le.100.) then
            if(cost.le.coscut)  call hf1(ih0+33,xm12,1.)
            if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+34,cost,1.)
          elseif(e2g.gt.100..and.e2g.le.150.) then
            if(cost.le.coscut) call hf1(ih0+35,xm12,1.)
            if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+36,cost,1.)
          elseif(e2g.gt.150..and.e2g.le.200.) then
            if(cost.le.coscut) call hf1(ih0+37,xm12,1.)
            if(xm12.gt..12.and.xm12.lt..15) call hf1(ih0+38,cost,1.)
          endif
        enddo
      enddo
      call hf1(ih0+16,d2gmin,1.)
      call hf2(ih0+18,xm12min,d2gmin,1.)
      if(ih0.le.100) then
        egc=egc*2.
        disc=disc*2.
        xm2c=xm2c*2.
        goto 1
      endif
      if(nev.ge.1000) then
        call store
        nev=0
      endif
      ngam=ngamsv
      do ig=1,ngam
        egam(ig)  =egamsv(ig)
        xgam(ig)  =xgamsv(ig)
        ygam(ig)  =ygamsv(ig)
        chigam(ig)=chigamsv(ig)
        gtime(ig) =gtimesv(ig)
        indmc(ig) =indmcsv(ig)
        ndf(ig)   =ndfsv(ig)
        lnclust(ig)=lnclustsv(ig)
      enddo
      return
      end
c
c-----------------------------------------------------------

      subroutine book
      implicit none
#include "parameters.inc"
#include "ccuts.inc"
#include "fcuts.inc"
#include "files.inc"
#include "pawc.inc"
#include "runpar.inc"
      character*32 tit2
      character*24 tit9
      character*28 tit17
      character*20 tit
      character*20 tittime
      integer i, ih0
      logical lexist
      data tit/'Eg/Efit    Column XX'/
      data tittime/'Time for   Column XX'/
      data tit2/'ECAL2 total energy. (Elow =    )'/
      data tit9/'Gamma energy (cut =    )'/
      data tit17/'Split distance (Cut =     ) '/
      data lexist/.false./
c      print *,' *** book started ***'
      call hlimit(mem)
c      call hbook1(1,'ECAL2 energy. Init. calibr.',400,0.,200.,0.)
      write(tit2(27:30),'(f4.0)') Elow
      call hbook1(2,tit2,600,0.,60.,0.)
      call hbook1(3,' Hitted cell number',200,0.,200.,0.)
      call hbook2(4,'ECAL2 hitted cells',64,0.,64.,48,0.,48.,0.)
      CALL HBOOK1(5,'Number of clusters',50,0.,50.,0.)
      CALL HBOOK1(6,'Number of gamma',50,0.,50.,0.)
      CALL HBOOK1(8,'CHI**2 for gamma',600,0.,30.,0.)
      write(tit9(20:23),'(f4.1)') egcut
      call hbook1(9,tit9,200,0.,100.,0.)
c      call hbook2(10,'Drift of calib. coeff.',64,0.,64.,48,0.,48.,0.)
      call hbook2(11,'Individual thresholds',64,0.,64.,48,0.,48.,0.)
      call hbook1(12,'Iteration number',40,0.,40.,0.)
      call hbook2(13,'Main cell energy  vs Egamma',100,0.,5.,100,0.,5.,0.)
      call hbook1(14,'Cell energy',200,0.,2.,0.)
      call hbook1(15,'Thresholds',200,0.,2.,0.)
      CALL HBOOK1(16,'Cluster lengths for gamma',50,0.,50.,0.)
      write(tit17(22:26),'(f5.2)') discut
      call hbook1(17,tit17,100,0.,200.,0.)
      call hbook1(18,'chisq for egam gt 20',250,0.,25.,0.)
      call hbook1(19,'gamma time',200,-50.,50.,0.)
      call hbook2(44,'ECAL2 hitted cells no cut',64,0.,64.,48,0.,48.,0.)
      call hbook2(45,'ECAL2 hitted cells time cut',64,0.,64.,48,0.,48.,0.)
      call hbook2(46,'ECAL2 hitted cells thresholds cut',64,0.,64.,48,0.,48.,0.)
      call hbook2(50,'Time  vs  cell energy',100,-50.,50.,100,0.,5.,0.)
      call hbook1(51,'LG cell energy',200,0.,2.,0.)
      call hbook1(52,'SH cell energy',200,0.,2.,0.)
      call hbook1(53,'Cell energy for central zone',200,0.,2.,0.)
      call hbook1(54,'Cell energy in time',200,0.,2.,0.)
      call hbook2(55,'Minimal amplitudes',64,0.,64.,48,0.,48.,0.)
      call hbook2(56,'Time  vs  LG cell energy',100,-50.,50.,100,0.,5.,0.)
      call hbook2(57,'Time  vs  SH cell energy',100,-50.,50.,100,0.,5.,0.)
      call hbook2(58,'Time  vs  central zone cell energy',100,-50.,50.,100,0.,5.,0.)
      call hbook2(60,'Chisq  vs  E  (1g cluster)',100,0.,20.,20,0.,200.,0.)
      call hbook2(61,'Individual thresholds',64,0.,64.,48,0.,48.,0.)
      ih0=0
    1 ih0=ih0+100
      CALL HBOOK1(ih0+07,'Final number of gamma',50,0.,50.,0.)
      call hbook2(ih0+08,'chisq vs Egam',100,0.,50.,50,0.,100.,0.)
      call hbook1(ih0+09,tit9,200,0.,100.,0.)
      call hbook1(ih0+11,'Number of bad cells',100,0.,100.,0.)
      call hbook2(ih0+12,'Dead cell positions',64,0.,64.,48,0.,48.,0.)
      call hbook2(ih0+13,'Noisy cell positions',64,0.,64.,48,0.,48.,0.)
      call hbook2(ih0+14,'Gamma frequency',64,0.,64.,48,0.,48.,0.)
      call hbook2(ih0+15,'Gamma freq. at center',40,28.,38.,40,19.,29.,0.)
      call hbook1(ih0+16,'minimal distance between gammas  (mm)',100,0.,200.,0.)
      call hbook2(ih0+17,'2g mass  vs  min Eg',50,0.,.3,50,0.,5.,0.)
      call hbook2(ih0+18,'2g mass  vs 2g min dist',50,0.,0.3,50,0.,150.,0.)
      call hbook1(ih0+19,'2g min distance',200,0.,200.,0.)
      call hbook1(ih0+20,'Energy of 2g events',400,0.,200.,0.)
      call hbook1(ih0+21,'2gamma mass (any combination)',100,0.,1.,0.)
      call hbook1(ih0+23,'Mass of 2gamma events.   E .gt. 40 GeV',500,0.,1.,0.)
      call hbook1(ih0+24,'cos(theta) in pi0.   E .gt. 40 GeV',100,0.,1.,0.)
      call hbook1(ih0+25,'Mass of 2gamma events.   E .gt. 100 GeV',300,0.,.75,0.)
      call hbook1(ih0+26,'cos(theta) in pi0.   E .gt. 100 GeV',100,0.,1.,0.)
      call hbook1(ih0+27,'Mass of 2gamma events.   E .gt. 150 GeV',300,0.,.75,0.)
      call hbook1(ih0+28,'cos(theta) in pi0.   E .gt. 150 GeV',100,0.,1.,0.)
c      call hbook1(ih0+31,'2 gamma mass.   E(2g) .gt. 10 GeV',500,0.,1.,0.)
c      call hbook1(ih0+32,'pi0 cos(theta)   E(2g) .gt. 10 GeV',100,0.,1.,0.)
      call hbook1(ih0+33,'2 gamma mass.   E(2g) .gt. 40 GeV',500,0.,1.,0.)
      call hbook1(ih0+34,'pi0 cos(theta)   E(2g) .gt. 40 GeV',100,0.,1.,0.)
      call hbook1(ih0+35,'2 gamma mass.   E(2g) .gt. 100 GeV',500,0.,1.,0.)
      call hbook1(ih0+36,'pi0 cos(theta)   E(2g) .gt. 100 GeV',100,0.,1.,0.)
      call hbook1(ih0+37,'2 gamma mass.   E(2g) .gt. 150 GeV',500,0.,1.,0.)
      call hbook1(ih0+38,'pi0 cos(theta)   E(2g) .gt. 150 GeV',100,0.,1.,0.)
      if(ih0.le.100) goto 1
c      call hbook2(301,'Peak position',64,0.,64.,48,0.,48.,0.)
      call hbook1(401,'X of gamma',96,0.,dcell/2.,0.)
      call hbook1(402,'Y of gamma',96,0.,dcell/2.,0.)
      call hbook1(403,'X of gamma',192,-3.*dcell,dcell,0.)
      call hbook1(404,'Y of gamma',192,-2.*dcell,2.*dcell,0.)
      do i=1000,1063
        write(tit(19:20),'(i2)') i-1000
        call hbprof(i,tit,48,0.,48.,0.,100.,'S')
      enddo
      call hbook2(1070,'Eg/Efit mean',64,0.,64.,48,0.,48.,0.)
      call hbook2(1071,'err(Eg/Efit-mean)',64,0.,64.,48,0.,48.,0.)
      call hbook2(1072,'Eg/Efit stat',64,0.,64.,48,0.,48.,0.)
      call hbook1(1073,'Eg/Efit',200,0.,4.,0.)
      call hbook1(1074,'err(Eg/Efit-mean)',250,0.,5.,0.)
      do i=2000,2063
        write(tittime(19:20),'(i2)') i-2000
c        call hbprof(i,tittime,48,0.,48.,-200.,200.,'S')
        call hbprof(i,tittime,48,0.,48.,-200.,200.,' ')
      enddo
      write(filehbook(7:12),'(i6.6)') NRUN
      inquire(file=filehbook,exist=lexist)
      if(lexist) then
        call hrget(0,filehbook,'A')
        print *,'file ',filehbook,' has been read'
      endif
      return
      end
C
      subroutine store
      implicit none
#include "parameters.inc"
#include "chkbad.inc"
#include "files.inc"
#include "gclb.inc"
#include "rawcells.inc"
#include "mincells.inc"
#include "runpar.inc"
#include "thresholds.inc"
#include "thr_test.inc"
      integer i
      real dtr(numc)
c      logical hexist
      write(filehbook(7:12),'(i6.6)') NRUN
      call trans(gcal,dtr)
c      call hpak(10,dtr)
      call trans(cellthr,dtr)
c      print *,'dtr ='
c      print '(4(/16f6.2))',dtr
      call hpak(11,dtr)
      call trans(Amin,dtr)
      call hpak(55,dtr)
      call trans(thrtest,dtr)
      call hpak(61,dtr)
      call hreset(1073,' ')
      do i=1,numc
       call hf1(1073,rmn(i),1.)
      enddo
      call trans(rmn,dtr)
      call hpak(1070,dtr)
      call hreset(1074,' ')
      do i=1,numc
       call hf1(1074,rsq(i),1.)
      enddo
      call trans(rsq,dtr)
      call hpak(1071,dtr)
      call trans(rst,dtr)
      call hpak(1072,dtr)
      print *,' * store *  file ',filehbook
      call hrput(0,filehbook,'T')
!      write(flthresh(7:12),'(i6.6)') NRUN
!      open(unit=8,file=flthresh,status='UNKNOWN')
!      write(8,'(a20,i6,a20,f7.3)'),' Thresholds for RUN ',NRUN,'  Common threshold ',thrcommon
!      write(8,'(128(24f7.3/))') (thrtest(i),i=1,numc)
!      close(unit=8)
      return
      end
c
      subroutine trans(datain,dataout)
      implicit none
      integer i, j
      real datain(48,64),dataout(64,48)
      do i=1,64
        do j=1,48
          dataout(i,j)=datain(j,i)
        enddo
      enddo
      return
      end
c
c-----------------------------------------------------------
      subroutine calib
      implicit none
#include "icalad.inc"
#include "parameters.inc"
#include "cells.inc"
#include "clb.inc"
      integer i, iad
C
      do i=1,ncell
        iad=icalad(x(i),y(i))
        if(iad.lt.0.or.iad.ge.numc) cycle
        if(cal(iad+1).lt..01) cycle
        e(i)=e(i)/cal(iad+1)
      enddo
      return
      END
C
c
      subroutine readcal(filename)
      implicit none
#include "parameters.inc"
#include "clb.inc"
      character*(*) filename
      integer i
C
      OPEN(UNIT=8,FILE=filename,STATUS='OLD',ERR=10)
      READ(8,*,err=9,end=9) cal
      READ(8,*,err=9,end=9) eled
      READ(8,*,err=9,end=9) stat
      CLOSE(UNIT=8)
      do i=1,numc
        calsav(i)=cal(i)
      enddo
      return
    9 close(unit=8)
   10 print *,' *** READCAL ***  READ ERROR.  File ',filename,' init.'
      stop
      do i=1,numc
        cal(i)=1.
        calsav(i)=1.
      enddo
      return
      END
C
c
      subroutine writcal(filename)
      implicit none
#include "parameters.inc"
#include "clb.inc"
      character*(*) filename
      integer i
C
      OPEN(UNIT=8,FILE=filename,form='formatted',status='UNKNOWN',err=2)
      WRITE(8,'(16f8.4)') (cal(i),i=1,numc)
      WRITE(8,'(16f8.1)') (eled(i),i=1,numc)
      WRITE(8,'(16f8.1)') (stat(i),i=1,numc)
      CLOSE(UNIT=8)
      RETURN
    2 write(*,*) ' ** WRITCAL **  OPEN failure.  File "',filename,'"'
      return
      END
c
c-----------------------------------------------------------
      subroutine varthr(nevthr)
      implicit none
#include "parameters.inc"
#include "badlist.inc"
#include "cells.inc"
#include "clb.inc"
#include "cluster.inc"
#include "iflags.inc"
#include "mincells.inc"
#include "thresholds.inc"
      real buf(numc)
      integer ibuf(numc)
      integer i, icl, icnt, ini, ip, itr, ix, ixmax, ixmin, iy, j
      integer leng, LTR, nevthr
      real defthr, sc
      DATA LTR/1000/,sc/1.02/,itr/0/,ini/1/
C
      if(ivarthr.ne.0) then
        if(ini.eq.1) then
          if(ivarthr.eq.1) print *,'*** VARTHR started ***'
          ini=0
        endif
        if(ivarthr.ne.1) return
        if(ncell.lt.1) return
        nevthr=nevthr+1
        ITR=ITR+1
        IF(ITR.GT.LTR) THEN  !  Decrease all thresholds
          ITR=0
          DO I=1,numc
            cellthr(I)=cellthr(I)/sc-.002
            IF(cellthr(I).lt.0.) cellthr(I)=0.
          ENDDO
        ENDIF
        call findclust(ncell,e,iaddr,x,y) !  Search noisy channels
        ip=1
        do icl=1,ncl
          leng=0
          ixmin=99
          ixmax=-1
          do i=ip,ip+lencl(icl)-1
            j=listcells(i)
            ix=iaddr(j)/48
            iy=iaddr(j)-ix*48
            if(ix.lt.ixmin) ixmin=ix
            if(ix.gt.ixmax) ixmax=ix
            leng=leng+1
            buf(leng)=e(j)
            ibuf(leng)=iaddr(j)
          enddo
!         Increase individual thresholds
          IF(ixmin.gt.0.and.ixmax.lt.63.and.ixmax.eq.ixmin) THEN
            do i=ip,ip+lencl(icl)-1
              j=listcells(i)
              icnt=iaddr(j)+1
              if(e(j).gt.cellthr(icnt)*.95) cellthr(icnt)=cellthr(icnt)*sc+.002
            enddo
c
c            if(leng.gt.1.and.ixmin.eq.62) then
c              print *
c              call cprint(leng,buf,ibuf)
c              print '(a,200f6.2)','cellthr =',(cellthr(ibuf(i)+1),i=1,leng)
c            endif
c          ELSE
c            DO i=ip,ip+lencl(icl)-1 !  Decrease thresholds in cluster
c              j=listcells(i)
c              icnt=iaddr(j)+1
c              IF(cellthr(icnt).GT..005) cellthr(icnt)=cellthr(icnt)-.0001
c            ENDDO
          ENDIF
          ip=ip+lencl(icl)
        enddo
      endif
      return
      END
c
      subroutine readthr(athr)
      implicit none
#include "parameters.inc"
#include "files.inc"
#include "runpar.inc"
#include "thresholds.inc"
      integer i,nrunthr
      character*16 char16
      character*20 char20
      real athr(*)
C
      write(flthresh(7:12),'(a6)') '000000'
      OPEN(UNIT=8,FILE=flthresh,STATUS='OLD',ERR=100)
      print *,'File ',flthresh,' is open'
      read(8,'(a20,i6,a16,f7.2)',err=101),char20,nrunthr,char16,thrcommon
      READ(8,'(24f7.3)',err=100) (athr(i),i=1,numc)
      print *,'**readthr**  Thresholds are read for run ',nrunthr
      CLOSE(UNIT=8)
      RETURN
 100  write(flthresh(7:12),'(i6.6)') NRUN
      OPEN(UNIT=8,FILE=flthresh,STATUS='OLD',ERR=101)
      print *,'File ',flthresh,' is open'
      read(8,'(a20,i6,a16,i8)',err=101),char16,nrunthr
      READ(8,'(24f7.2)',err=102) (athr(i),i=1,numc)
      print *,'**readthr**  Thresholds are read for run ',nrunthr
      CLOSE(UNIT=8)
      RETURN
 101  print '(3a)','**readthr**  NO FILE ',flthresh,' Threshols initialized'
      goto 103
 102  print '(3a)', '**readthr** READ ERROR. File ',flthresh,' init.'
      CLOSE(UNIT=8)
 103  do i=1,numc
        athr(i)=12.
      enddo
      return
      END
C
      subroutine writhr(athr)
      implicit none
#include "parameters.inc"
#include "files.inc"
#include "runpar.inc"
#include "thresholds.inc"
      integer i
      real athr(*)
C
      write(flthresh(7:12),'(i6.6)') NRUN
      OPEN(UNIT=8,FILE=flthresh,status='UNKNOWN',err=10)
      write(8,'(a20,i6,a16,f7.2,a4)'),' Thresholds for RUN ',NRUN,'    thrcommon = ',thrcommon,' GeV'
      WRITE(8,'(128(24f7.2/))') (athr(i),i=1,numc)
      CLOSE(UNIT=8)
      print '(3a)',' ** writhr **  File "',flthresh,'" is written'
      RETURN
   10 print '(3a)',' ** writhr ** OPEN failure.  File "',flthresh,'"'
      return
      END
C
C********************************************************************
      subroutine thrcut
      implicit none
#include "icalad.inc"
#include "parameters.inc"
#include "cells.inc"
#include "geom.inc"
#include "thresholds.inc"
      integer i, ix, iy, j, ncell_thr
      ncell_thr=0
      do i=1,ncell
        j=icalad(x(i),y(i))
        if(j.le.0) cycle
        ix=j/48
        iy=j-ix*48
        if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle
        if(e(i).gt.cellthr(j+1).and.e(i).gt.thrcommon) then
          ncell_thr=ncell_thr+1
          e(ncell_thr)=e(i)
          iaddr(ncell_thr)=iaddr(i)
          x(ncell_thr)=x(i)
          y(ncell_thr)=y(i)
          err(ncell_thr)=err(i)
          time(ncell_thr)=time(i)
        endif
      enddo
      ncell=ncell_thr
      return
      end
c
c-----------------------------------------------------------
      subroutine mindis(dismin)
      implicit none
#include "parameters.inc"
#include "gamma.inc"
      integer ig1, ig2
      real dis, dismin, dist
      dismin=999999.
      do ig1=1,ngam-1
        do ig2=ig1+1,ngam
          dist=sqrt((xgam(ig1)-xgam(ig2))**2+(ygam(ig1)-ygam(ig2))**2)
          dis=2.*sqrt(egam(ig1)*egam(ig2))*dist/(egam(ig1)+egam(ig2))
          if(dis.lt.dismin) dismin=dis
        enddo
      enddo
      return
      end
c
c-----------------------------------------------------------
      integer function icalad(x,y)
      implicit none
      real x, y
#include "parameters.inc"
      integer ix, iy
      icalad=-1
      ix=nint(x/dcell-.5)+32
      iy=nint(y/dcell-.5)+24
      if(ix.lt.0.or.ix.ge.64.or.iy.lt.0.or.iy.ge.48) return
      icalad=ix*48+iy
      end
c
c-----------------------------------------------------------
      real function xm2g(ig1, ig2)
      implicit none
#include "parameters.inc"
#include "gamma.inc"
#include "geom.inc"
      integer ig1, ig2
      real d1, d2, px, py, pz, xm
      xm2g = 0.
      if(ig1.le.0.or.ig1.gt.200.or.ig2.le.0.or.ig2.gt.200) return
      d1 = sqrt(xgam(ig1)**2+ygam(ig1)**2+Zecal2**2)
      d2 = sqrt(xgam(ig2)**2+ygam(ig2)**2+Zecal2**2)
      px = egam(ig1)*xgam(ig1)/d1 + egam(ig2)*xgam(ig2)/d2
      py = egam(ig1)*ygam(ig1)/d1 + egam(ig2)*ygam(ig2)/d2
      pz = egam(ig1)*Zecal2/d1 + egam(ig2)*Zecal2/d2
      xm = (egam(ig1)+egam(ig2))**2-px**2-py**2-pz**2
      if(xm.gt.0.) xm2g = sqrt(xm)
      return
      end
c
c-----------------------------------------------------------
      subroutine ecal2thr
      implicit none
#include "parameters.inc"
#include "runpar.inc"
#include "cluster.inc"
#include "thresholds.inc"
      integer i,ix,iy
      if(NRUN.ge.81473.and.NRUN.le.82145) then
        do i=1,numc    !  put individual thresholds
          ix=(i-1)/idx
          iy=i-1-ix*idx
          cellthr(i)=0.2
          if(ix.ge.16.and.ix.le.47) then
            cellthr(i)=0.400
            if(ix.ge.23.and.ix.le.42.and.iy.ge.16.and.iy.le.31) cellthr(i)=0.600
          endif
        enddo
      endif
      return
      end
c
c-----------------------------------------------------------
      subroutine joujou
      implicit none
#include "errcell2.inc"
#include "icalad.inc"
#include "xy.inc"
#include "parameters.inc"
#include "badlist.inc"
#include "cells.inc"
#include "chkbad.inc"
#include "clb.inc"
#include "cluster.inc"
#include "ccuts.inc"
#include "fcuts.inc"
#include "debug.inc"
#include "gamma.inc"
#include "geom.inc"
#include "iflags.inc"
#include "isitlg.inc"
#include "leds.inc"
#include "misc.inc"
#include "rawcells.inc"
#include "mincells.inc"
#include "runpar.inc"
#include "showerprof.inc"
#include "thresholds.inc"
      integer iaddrmem(numc)
      integer i,j,ic,icl,iad,ibin,icalcthr,iclmax,icell,icellmax,icnt,ip,ipclmax
      integer ifit,iled
      integer ix,ixcell,iy,iycell
      integer ledstatlim,itr,n,lgmax
      integer nev,nevlong,nevmax,nspill,numspill
      integer nbeam,nbm
      logical lini,ledinf
      real abst,aclmax,athr0,athrlog,bin
      real cumLG(700),statLG(700),cumSH(700),statSH(700)
      real xwcl(numcl),ywcl(numcl),ecl(numcl),xcl(numcl),ycl(numcl),distcl(numcl),eclk(numcl)
      real distclmax,dr,dx,eclmax,eclkmax
      real emin,emax,Etot
      real Ebeam,xbeam,ybeam
      real ebmspill,xbmspill,ybmspill
      real scale,stat0,statlim,w,w_stat,sumA,sumE,sumEK
      real xc,xLG,xSH
      real xx,xcmax,xwmax,xwmin
      real yy,ycmax,ywmax,ywmin
      real SumCum(7),Xcum(7),aledmem(48,64),EK(numc)
      data Ebeam/40./,xbeam,ybeam/0.,0./,nbeam/0/
      data lini/.true./,ledinf/.false./
      data nev,nevlong/2*0/,nevmax/99999999/,nspill,numspill/2*-1/
      data emin,emax/5.,100./
      data eled/numc*0./
      data ledstatlim/1000/,iled/0/,aledmem/numc*0./,iaddrmem/numc*-1/
      data distclmax/300./,statlim/200./
      data athr0/12./

      if(lini) then
        ifit=0
        time0=0.
        timew=10.
        call read_ecalib('calib.dat')

!        icalcthr=0     ! icalcthr=0 to use constant thresholds (athr0)
!        icalcthr=1     ! icalcthr=1 to read them only.
!        icalcthr=2     ! icalcthr=2 to read and calculate thresholds
!                       ! other icalcthr to use exturnal thresholds
        icalcthr=0
        if(icalcthr.eq.1.or.icalcthr.eq.2) call readthr(athr)
        if(icalcthr.eq.0.or.icalcthr.eq.1.or.icalcthr.eq.2) then
          do i=1,numc
            if(icalcthr.eq.0) athr(i)=athr0
            cellthr(i)=athr(i)/cal(i)
            Amin(i)=athr(i)
            if(icalcthr.eq.2) Amin(i)=Amin(i)+10.
          enddo
c          print *,'athr init ='
c          print '(24f7.2)', athr
        endif

        print *,' *** joujou started *** ','  ifit =',ifit,'  icalcthr =',icalcthr
        print *,'thrcommon =',thrcommon,'  egcut =',egcut,'  emin,emax =',emin,emax
        call joujoubook(ifit)
        call readled('leds.dat')
        call hreset(14,' ')
        do i=1,numc
          call hf1(14,cal(i),1.)
        enddo
        if(icalib.ne.0) then
          ihist =0
        endif
        lini=.false.
      endif

      nevlong=nevlong+1
      nev=nev+1
      ngam=0
      if(ncell.gt.numc/2) then
        do i=1,ncell
          itr=iand(ntrig,2**i)
          if(itr.ne.0) call hf1(98,i+.5,1.) ! LED trigger bits
        enddo
      endif
      if(ncell.le.0.or.ncell.gt.numc) goto 99
      if(icalcthr.eq.2) then ! to calculate thresholds
        do i=1,ncell
          icnt=icalad(x(i),y(i))+1
          if(icnt.le.0.or.icnt.gt.numc) cycle
          if(acell(i).lt.2.*athr(icnt)) then
            astat(icnt)=astat(icnt)+1.
            athrlog=alog(acell(i))
            athrlog=athrlog+(alog(athr(icnt))-athrlog)/astat(icnt)
            athr(icnt)=exp(athrlog)
          endif
        enddo
      endif
      do i=0,15
        itr=iand(ntrig,2**i)
        if(itr.ne.0) call hf1(99,i+.5,1.) ! Trigger bits
      enddo
      if(ntrig.eq.4096) then ! LED signals
c        print *,' ** joujou **  ntrig =',ntrig,'  ncell =',ncell
        if(ncell.lt.1000) goto 99
        if(iled.lt.32767) iled=iled+1
c        print *,' iled =',iled
        if(nburst.ne.numspill) then
          numspill=nburst
          do i=1,numc
            nledinspill(i)=0
            ix=(i-1)/48
            iy=i-1-ix*48
            call hfill(2000+ix,iy+.5,spillled(i),1.)
            if(ledstat(i).gt.0) call hf2(3000+ix,iy+.5,spillled(i)/ameanled(i),1.)
          enddo
          call hf1(10082,float(nburst),spillled(83))
          call hf1(10083,float(nburst),spillled(84))
          call hf1(10084,float(nburst),spillled(85))
        endif
        ledinf=.true.
        do i=1,ncell
          iaddr(i)=icalad(x(i),y(i))
          if(iaddrmem(iaddr(i)+1).eq.iaddr(i)) then
            print *,'cell X =',iaddr(i)/48,' Y =',iaddr(i)-iaddr(i)/48*48,' is filled already !!!'
            cycle
          endif
          if(iaddr(i).lt.0) then
            print *,'icell,x,y =',i,x(i),y(i)
            cycle
          endif
          ix=iaddr(i)/48
          iy=iaddr(i)-ix*48
          bin=iled+.5
          call hf1(20000+ix*100+iy,bin,acell(i))
cc          if(ix.eq.1.and.iy.eq.35) then
cc            print *,'histo',20000+ix*100+iy,' filled with bin=',bin,'  acell =',acell(i),'  time =',time(i)
cc          endif
          call hf1(30000,abs(acell(i)-aledmem(iy+1,ix+1)),1.)
          aledmem(iy+1,ix+1)=acell(i)
          nledinspill(iaddr(i)+1)=nledinspill(iaddr(i)+1)+1
          if(ledstat(iaddr(i)+1).lt.ledstatlim) ledstat(iaddr(i)+1)=ledstat(iaddr(i)+1)+1
          spillled(iaddr(i)+1)=spillled(iaddr(i)+1)+(acell(i)-spillled(iaddr(i)+1))/nledinspill(iaddr(i)+1)
          ameanled(iaddr(i)+1)=ameanled(iaddr(i)+1)+(acell(i)-ameanled(iaddr(i)+1))/ledstat(iaddr(i)+1)
c          if(ix.eq.1.and.iy.eq.35) print *,'spillled,ameanled=',spillled(iaddr(i)+1),ameanled(iaddr(i)+1)
        enddo
        do i=1,numc
          iaddrmem(i)=-1
        enddo
        goto 99
      endif ! end of LED signals

cc      if(.not.ledinf) return
      if(ldebug) then
        print *,'ecell ='
        print '(10f8.2)',(e(i),i=1,ncell)
cc        print *,'cellthr ='
cc        print '(10f8.2)',(cellthr(iaddr(i)+1),i=1,ncell)
!        print *,'iaddr ='
!        print '(10i8)',(iaddr(i),i=1,ncell)
        print *,'xcell ='
!        print '(10f8.2)',(x(i),i=1,ncell)
        print '(10i8)',(nint(x(i)/dcell+.5)+31,i=1,ncell)
        print *,'ycell ='
!        print '(10f8.2)',(y(i),i=1,ncell)
        print '(10i8)',(nint(y(i)/dcell+.5)+23,i=1,ncell)
      endif

      if(nevlong-nevlong/1000*1000.eq.0) then
        print '(a,i9,a,f5.1,a,f5.1,a)',' ***** event ',nevlong
     ,,'  Beam at X =',xbeam/dcell+31.5,'  Y =',ybeam/dcell+23.5,' (cell)'
      endif
c      print *,' ***joujou***  run ',nrun,'  event ',numbev
c     ,  ,'  burst ',nburst,'  event in burst ',numbinburst

      sumA=0.
      sumEK=0.
      do i=1,ncell
        sumA=sumA+acell(i)
c        if(acell(i).gt.10.) print *,'i =',i,'  acell =',acell(i),'  time =',time(i)
        ek(i)=e(i) ! Kol calib
        if(ek(i).gt.thrcommon) sumEK=sumEK+ek(i)
        e(i)=acell(i)
      enddo
      call hf1(1,sumA,1.)
      call hf1(2,sumEK,1.)

      call ecalib
      sumE=0.
      do i=1,ncell
        if(e(i).gt.thrcommon) sumE=sumE+e(i)
      enddo
      call hf1(17,sumE,1.)
      do i=1,ncell
        iad=icalad(x(i),y(i))
        if(iad.lt.0.or.iad.ge.numc) cycle
        err(i)=sqrt(errcell2(e(i)))
        call hf2(50,time(i),e(i),1.)
        time(i)=time(i)-time0
        abst=abs(time(i))/timew ! time correction
        if(abst.gt.100.) abst=100.
        if(abst.gt.1.) then
          err(i)=err(i)*abst
          e(i)=e(i)/abst
          ek(i)=ek(i)/abst
        endif
      enddo
      do i=1,ncell
        iaddr(i)=icalad(x(i),y(i))
        if(iaddr(i).lt.0) then
          print *,'icell,x,y =',i,x(i),y(i)
        endif
      enddo
!      print *,'e =',(e(i),i=1,ncell)
      call hf1(3,float(ncell),1.)
      n=0
      do i=1,ncell
        icnt=iaddr(i)+1
        if(e(i).gt.thrcommon.and.e(i).gt.cellthr(icnt)) then
          n=n+1
          e(n)=e(i)
          iaddr(n)=iaddr(i)
          x(n)=x(i)
          y(n)=y(i)
          err(n)=err(i)
          time(n)=time(i)
          ek(n)=ek(i)
        endif
      enddo
      ncell=n
      if(ldebug) then
        print *,'ecell ='
        print '(10f8.3)',(e(i),i=1,ncell)
!        print *,'iaddr ='
!        print '(10i8)',(iaddr(i),i=1,ncell)
        print *,'xcell ='
        print '(10f8.2)',(x(i),i=1,ncell)
        print *,'ycell ='
        print '(10f8.2)',(y(i),i=1,ncell)
      endif

cc      call varthr(ivarthr)
cc      call thrcut
cc      call badcheck

      if(ncell.lt.6) goto 99
      Etot=0.
      do i=1,ncell
        Etot=Etot+e(i)
        scale=1.
        if(ledinf.and.spillled(iaddr(i)+1).gt.10.) scale=eled(iaddr(i)+1)/spillled(iaddr(i)+1)
        if(ledinf) e(i)=e(i)*scale
c        if(e(i).gt.50.) then
c          print *,'iaddr,e,eled,spillled =',iaddr(i),e(i),eled(iaddr(i)+1),spillled(iaddr(i)+1)
c        endif
        if(ledinf.and.e(i).gt.10..and.spillled(iaddr(i)+1).gt.10.) then
          call hf1(34,eled(iaddr(i)+1)/spillled(iaddr(i)+1),1.)
        endif
      enddo

      call findclust(ncell,e,iaddr,x,y)
cc      print *,'ncl =',ncl,'  lencl(icl)=',(lencl(icl),icl=1,ncl)
cc      print *,'listcells =',(listcells(i),i=1,ncell)
cc      print '(10f8.1)',(e(i),i=1,ncell)
      call hf1(5,float(ncl),1.)
      Etot=0.
      ip=1
      Eclmax=-9999.
      do icl=1,ncl
        ecl(icl)=0.
        xcl(icl)=0.
        ycl(icl)=0.
        eclk(icl)=0.
        xwmin=9999.
        xwmax=-9999.
        ywmin=9999.
        ywmax=-9999.
        do i=ip,ip+lencl(icl)-1
          ic=listcells(i)
cc          print *,'i,ic,e(ic),Etot =',i,ic,e(ic),Etot
          Etot=Etot+e(ic)
          ecl(icl)=ecl(icl)+e(ic)
          eclk(icl)=eclk(icl)+ek(ic)
          xcl(icl)=xcl(icl)+e(ic)*x(ic)
          ycl(icl)=ycl(icl)+e(ic)*y(ic)
          if(x(ic).lt.xwmin) xwmin=x(ic)
          if(x(ic).gt.xwmax) xwmax=x(ic)
          if(y(ic).lt.ywmin) ywmin=y(ic)
          if(y(ic).gt.ywmax) ywmax=y(ic)
        enddo
        xcl(icl)=xcl(icl)/ecl(icl)
        ycl(icl)=ycl(icl)/ecl(icl)
        xwcl(icl)=(xwmax-xwmin)/dcell+1.
        ywcl(icl)=(ywmax-ywmin)/dcell+1.
        call hf1(22,float(lencl(icl)),1.)
        if(Ecl(icl).gt.Eclmax) then
          Eclmax=Ecl(icl)
          Eclkmax=Eclk(icl)
          iclmax=icl
          ipclmax=ip
        endif
        ip=ip+lencl(icl)
      enddo
cc      print *,'Etot,delta,thrcommon =',Etot,delta(Etot,thrcommon),thrcommon
cc      Etot=Etot+delta(Etot,thrcommon)
cc      if(Etot.gt.50.) then
cc        print *,'e ='
cc        print '(10f8.1)',(e(i),i=1,ncell)
cc      endif
      call hf1(32,Etot,1.)
      do icl=1,ncl
        call hf1(20,xwcl(icl)+.5,1.)
        call hf1(21,ywcl(icl)+.5,1.)
        distcl(icl)=sqrt((xcl(icl)-xcl(iclmax))**2+(ycl(icl)-ycl(iclmax))**2)
      enddo
      if(numbinburst.lt.200) call hf1(202,Etot,1.)
      if(ntrig.eq.128) call hf1(101,Etot,1.)
      if(ntrig.eq.256) call hf1(102,Etot,1.)
      if(ntrig.eq.4096) call hf1(103,Etot,1.)

      call hf1(23,float(lencl(iclmax)),1.)
      if(lencl(iclmax).le.5.or.lencl(iclmax).gt.49) then
!        print *,'iclmax =',iclmax
        goto 99
      endif
      call hf1(24,Etot,1.)
      call hf1(27,xwcl(iclmax)+.5,1.)
      call hf1(28,ywcl(iclmax)+.5,1.)
      if(xwcl(iclmax).lt.1.5.or.xwcl(iclmax).gt.10.5) goto 99
      if(ywcl(iclmax).lt.1.5.or.ywcl(iclmax).gt.10.5) goto 99
      call hf1(25,Etot,1.)
      if(iand(ntrig,256).eq.0) goto 99
      call hf1(26,Etot,1.)
      aclmax=0.
      eclmax=0.
      eclkmax=0.
      xx=0.
      yy=0
      do i=ipclmax,ipclmax+lencl(iclmax)-1
        ic=listcells(i)
        aclmax=aclmax+acell(ic)
        eclmax=eclmax+e(ic)
        eclkmax=eclkmax+ek(ic)
        xx=xx+e(ic)*x(ic)
        yy=yy+e(ic)*y(ic)
      enddo
      xcmax=xx/eclmax
      ycmax=yy/eclmax
      icellmax=icalad(xcmax,ycmax)+1
      LGmax=IsItLG(icellmax-1)
      xcmax=xy(xcmax,LGmax)
      ycmax=xy(ycmax,LGmax)
      eclmax=eclmax+3.*cellthr(icellmax)
      call hf1(29,eclmax,1.)
      if(eclmax.lt.egcut) goto 99
      do icl=1,ncl
        if(icl.eq.iclmax) then
          if(eclmax.gt.Emin) call hfill(33,1.,eclmax,1.)
          cycle
        endif
        if(eclmax.gt.Emin.and.ecl(icl).gt.egcut) then
          call hf1(19,distcl(icl)/10.,1.)
          call hfill(33,distcl(icl)/10.,Ecl(icl),1.)
        endif
        if(distcl(icl).lt.distclmax) eclmax=eclmax+Ecl(icl)
      enddo
      call hf1(4,eclmax,1.)
      call hf1(44,Eclkmax,1.)
      call hf1(104,aclmax,1.)
      if(numbinburst.lt.100) call hf1(204,eclmax,1.)
      if(nbeam.le.99) nbeam=nbeam+1
      xbeam=xbeam+(xcmax-xbeam)/nbeam
      ybeam=ybeam+(ycmax-ybeam)/nbeam
      if(nbeam.gt.99) call hf2(9,xbeam/dcell+31.5,ybeam/dcell+23.5,1.)
      if(nspill.eq.nburst.and.eclmax.gt.20.) then
        nbm=nbm+1
        xbmspill=xbmspill+(xcmax-xbmspill)/nbm
        ybmspill=ybmspill+(ycmax-ybmspill)/nbm
        ebmspill=ebmspill+(eclmax-ebmspill)/nbm
        if(nbm.gt.100) then
          dr=sqrt((xcmax-xbmspill)**2+(ycmax-ybmspill)**2)
          call hfill(210,dr,ebmspill,1.)
          call hfill(211,xcmax-xbmspill,ebmspill,1.)
          call hfill(212,ycmax-ybmspill,ebmspill,1.)
          call hf2(213,xcmax-xbmspill,ycmax-ybmspill,1.)
        endif
      else
        nbm=0
        xbmspill=0.
        ybmspill=0.
        ebmspill=0.
        nspill=nburst
      endif

      if(abs(xcmax-xbeam).gt.200..or.abs(ycmax-ybeam).gt.200.) then
        print *,'xbeam,ybeam=',xbeam,ybeam,'  xcmax,ycmax =',xcmax,ycmax
        goto 99
      endif
cc      if(iand(ntrig,256).eq.0) goto 99
      if(eclmax.gt.emin.and.eclmax.lt.emax) then
        ip=ipclmax
cc        Eclmax=Eclmax*1.02
        do i=ip,ip+lencl(iclmax)-1
          ic=listcells(i)
          if(e(ic).gt.eclmax*0.1.and.e(ic).lt.eclmax*0.95) then
            icell=icalad(x(ic),y(ic))+1
cc            if(ledstat(icell).gt.0) then
            w=(e(ic)/Ebeam)**2
            if(w.gt.9.) w=9.
            stat(icell)=stat(icell)+w
            if(stat(icell).gt.99999.) stat(icell)=99999.
            stat0=stat(icell)
            if(stat0.gt.statlim)  stat0=statlim
            w_stat=w/stat0
            if(w_stat.ge.1.) w_stat=1.
            cal(icell)=cal(icell)*(1.+(Eclmax/Ebeam-1.)*w_stat)
            eled(icell)=eled(icell)-(eled(icell)-spillled(icell))*w_stat
c            if(y(ic).gt.0.) print *,'nev,icell,e,x,y ='
c     ,        ,nev,icell,e(ic),x(ic),y(ic)
cc            endif
          endif
        enddo
        ixcell=nint(xcmax/dcell+31.5)
        iycell=nint(ycmax/dcell+23.5)
        call hf2(16,ixcell+.5,iycell+.5,1.)
        if(ixcell.gt.16.and.ixcell.lt.47) then
          xc=xcmax/dcell-nint(xcmax/dcell-.5)
          if(xc.ge..5) xc=1.-xc
          if(iycell.gt.2.and.iycell.lt.13) then ! LG
            call hf1(30,xc,1.)
          endif
          if(iycell.gt.15.and.iycell.lt.32) then ! Sha
            call hf1(31,xc,1.)
          endif
!                       Cumulative function
          do j=1,7
            Xcum(j)=(nint(xcmax/dcell-.5)+j-4.)*dcell
            Sumcum(j)=0.
          enddo
          sumE=0.
          do i=ipclmax,ipclmax+lencl(iclmax)-1
            ic=listcells(i)
            sumE=sumE+e(ic)
            do j=1,7
              if(Xcum(j).lt.x(ic)) Sumcum(j)=Sumcum(j)+e(ic)
            enddo
          enddo
          do j=1,7
            Sumcum(j)=Sumcum(j)/sumE
          enddo
          if(iycell.gt.2.and.iycell.lt.13) then ! LG
            xLG=xy(xcmax,1)
            do j=1,7
              dx=(xLG-Xcum(j))/dcell
              ibin=int(dx*100.+350.)
              if(ibin.le.0.or.ibin.gt.700) cycle
              statLG(ibin)=statLG(ibin)+1.
              cumLG(ibin)=cumLG(ibin)+(Sumcum(j)-cumLG(ibin))/statLG(ibin)
            enddo
          endif
          if(iycell.gt.15.and.iycell.lt.32) then ! Sha
            xSH=xy(xcmax,0)
            do j=1,7
              dx=(xSH-Xcum(j))/dcell
              ibin=int(dx*100.+350.)
              if(ibin.le.0.or.ibin.gt.700) cycle
              statSH(ibin)=statSH(ibin)+1.
              cumSH(ibin)=cumSH(ibin)+(Sumcum(j)-cumSH(ibin))/statSH(ibin)
            enddo
          endif
        endif
      endif
!
   99 if(nev.ge.1000) then
        call joustore(ifit)
        nev=0
        call write_ecalib('calib.dat')
        call writled('leds.dat')
        if(ifit.ne.0) call writbad
c        if(ivarthr.ne.0) call writhr(athr)
        if(icalcthr.eq.2) call writhr(athr)
      endif
      return
      end
c
c-----------------------------------------------------------
      subroutine joujoubook(ifit)
      implicit none
      integer ifit
#include "parameters.inc"
#include "files.inc"
#include "runpar.inc"
#include "pawc.inc"
      character*20 tit
      character*24 titled,ledtit
      integer i,ix,iy,memrun
      logical hexist
      data tit/'Eg/Efit    Column XX'/
      data titled/'LED fluct.     Column XX'/
      data ledtit/'LED in time X= XX  Y= YY'/
      data memrun/-1/

      print *,' *** joujoubook started ***'
      write(filehbook(1:6),'(a6)') 'joujou'
      write(filehbook(7:12),'(i6.6)') NRUN
      call hlimit(10000000)
      call hbook1(1,' Total energy before calibration',200,0.,4000.,0.)
      call hbook1(2,' Calibrated total energy ',600,0.,60.,0.)
      call hbook1(3,' GAMS length (above threshold)',100,0.,100.,0.)
      call hbook1(4,'Energy of max. cluster',600,0.,60.,0.)
      CALL HBOOK1(5,'Number of clusters',50,0.,50.,0.)
      CALL HBOOK1(6,'Number of gamma',20,0.,20.,0.)
      CALL HBOOK1(7,'Number of gamma',20,0.,20.,0.)
      CALL HBOOK1(8,'CHI**2/ndf for gamma',200,0.,20.,0.)
      call hbook2(9,'Y/X  beam positions',64,0.,64.,48,0.,48.,0.)
      call hbook2(10,'Drift of calib. coeff.',64,0.,64.,48,0.,48.,0.)
      call hbook2(11,'Individual ADC thresholds',64,0.,64.,48,0.,48.,0.)
      call hbook2(12,'Statistics',64,0.,64.,48,0.,48.,0.)
      call hbook2(13,'Leds',64,0.,64.,48,0.,48.,0.)
      call hbook1(14,'Initial calib. coeff.',250,0.,125.,0.)
      call hbook1(15,'Current calib. coeff.',250,0.,125.,0.)
      call hbook2(16,'Max cluster position',64,0.,64.,48,0.,48.,0.)
      call hbook1(17,'Total energy calibrated',600,0.,60.,0.)
      CALL HBOOK1(18,'CHI**2/ndf for Egamma.gt.20GeV',200,0.,20.,0.)
      call hbook1(19,'Distance between clusters  (cm)',200,0.,200.,0.)
      call hbook1(20,'X cluster width',20,0.,20.,0.)
      call hbook1(21,'Y cluster width',20,0.,20.,0.)
      call hbook1(22,'Cluster lengths',50,0.,50.,0.)
      call hbook1(23,'Emax cluster length',50,0.,50.,0.)
      call hbook1(24,'Total energy cluster length rej.',300,0.,60.,0.)
      call hbook1(25,'Total energy cluster x width rej.',300,0.,60.,0.)
      call hbook1(26,'Total energy with triggermask 256',300,0.,60.,0.)
      call hbook1(27,'X width of Emax cluster',20,0.,20.,0.)
      call hbook1(28,'Y width of Emax cluster',20,0.,20.,0.)
      call hbook1(29,'Energy of Emax cluster',300,0.,60.,0.)
      call hbook1(30,'Xc LG',500,0.,.5,0.)
      call hbook1(31,'Xc Sha',500,0.,.5,0.)
      call hbook1(32,'Calibrated total energy. LED corr. ',300,0.,60.,0.)
      call hbook2(33,'Cluster energy vs distance from max cluster',
     ,  50,0.,200.,100,0.,50.,0.)
      call hbook1(34,'LED stability',500,0.75,1.25,0.)
      call hbook1(40,'LG cumulative function',700,-3.5,3.5,0.)
      call hbook1(41,'Shaslyk cumulative function',700,-3.5,3.5,0.)
      call hbook1(44,'EK of max. cluster',600,0.,60.,0.)
      call hbook2(50,'Time  vs  cell energy',100,-50.,50.,100,0.,10.,0.)
      call hbook2(55,'Minimal amplitudes',64,0.,64.,48,0.,48.,0.)
      call hbook2(56,'Minimal amplitude statistics',64,0.,64.,48,0.,48.,0.)
      call hbook2(57,'Proposed thresholds',64,0.,64.,48,0.,48.,0.)
      call hbook1(98,'LED trigger bits',16,0.,16.,0.)
      call hbook1(99,'Trigger bits',16,0.,16.,0.)
c      call hbook1(100,'alog(rled) values',100,-1.,1.,0.)
      call hbook1(101,' Calibrated total energy trig=128',200,0.,100.,0.)
      call hbook1(102,' Calibrated total energy trig=256',200,0.,100.,0.)
      call hbook1(103,' Calibrated total energy trig=4096',200,0.,100.,0.)
      call hbook1(104,'Raw energy of max. cluster',300,0.,3000.,0.)
      call hbook1(202,'Calibrated total energy 200 events in spill',300,0.,60.,0.)
      call hbook1(204,'Energy of max. cluster',600,0.,60.,0.)
      call hbprof(210,'R energy distribition in spill',96,   0.,    96., 20.,60.,'  ')
      call hbprof(211,'X energy distribition in spill',96,-192.,192.,20.,60.,'  ')
      call hbprof(212,'Y energy distribition in spill',96,-192.,192.,20.,60.,'  ')
      call hbook2(213,'X,Y distribition in spill',96,-192.,192.,96,-192.,192.,0.)

      if(ifit.eq.0) goto 99 ! if there is no fit
      do i=1000,1063
        write(tit(19:20),'(i2)') i-1000
        call hbprof(i,tit,48,0.,48.,0.,100.,'S')
      enddo
      call hbook2(1070,'Eg/Efit mean',64,0.,64.,48,0.,48.,0.)
      call hbook2(1071,'err(Eg/Efit-mean)',64,0.,64.,48,0.,48.,0.)
      call hbook2(1072,'Eg/Efit stat',64,0.,64.,48,0.,48.,0.)
      call hbook1(1073,'Eg/Efit',200,0.,4.,0.)
      call hbook1(1074,'err(Eg/Efit-mean)',200,0.,4.,0.)
      do i=2000,2063
        write(titled(23:24),'(i2)') i-2000
        call hbprof(i,titled,48,0.,48.,0.,4096.,'S')
      enddo
      do i=3000,3063
        write(titled(23:24),'(i2)') i-3000
        call hbook2(i,titled,48,0.,48.,200,.8,1.2,0.)
      enddo
      call hbook1(10082,'spill number  vs  LED(1,34)',200,1.,201.,0.)
      call hbook1(10083,'spill number  vs  LED(1,35)',200,1.,201.,0.)
      call hbook1(10084,'spill number  vs  LED(1,36)',200,1.,201.,0.)
      do ix=0,63
        do iy=0,47
          write(ledtit(23:24),'(i2)') iy
          write(ledtit(15:16),'(i2)') ix
          call hbook1(20000+ix*100+iy,ledtit,500,1.,501.,0.)
        enddo
      enddo
      call hbook1(30000,'Subsequent LED difference',2050,0.,2050.,0.)
      if(nrun.ne.memrun) then
        memrun=nrun
        do ix=0,63
          do iy=0,47
            call hreset(20000+ix*100+iy,' ')
          enddo
        enddo
      endif
c      call hreset(19,' ')

  99  inquire(file=filehbook,exist=hexist)
      if(hexist) then
        call hrget(0,filehbook,'A')
        print *,'file ',filehbook,' has been read'
      endif
      return
      end
C
c-----------------------------------------------------------------
      subroutine joustore(ifit)
      implicit none
      integer ifit
#include "parameters.inc"
#include "badlist.inc"
#include "cells.inc"
#include "chkbad.inc"
#include "clb.inc"
#include "cluster.inc"
#include "ccuts.inc"
#include "fcuts.inc"
#include "debug.inc"
#include "files.inc"
#include "gamma.inc"
#include "geom.inc"
#include "iflags.inc"
#include "isitlg.inc"
#include "leds.inc"
#include "misc.inc"
#include "rawcells.inc"
#include "mincells.inc"
#include "runpar.inc"
#include "showerprof.inc"
#include "thresholds.inc"
      integer i
      logical hexist
      real r(numc),dtr(numc)
      write(filehbook(1:6),'(a6)') 'joujou'
      write(filehbook(7:12),'(i6.6)') NRUN
      call trans(cal,dtr)
      call hpak(10,dtr)
      do i=1,numc
        r(i)=cellthr(i)*cal(i)
c        r(i)=cellthr(i)
      enddo
      call trans(r,dtr)
c      print *,'dtr ='
c      print '(4(/16f6.2))',dtr
      call hpak(11,dtr)
      call trans(stat,dtr)
      call hpak(12,dtr)
      call trans(eled,dtr)
      call hpak(13,dtr)
      call hreset(15,' ')
      do i=1,numc
        call hf1(15,cal(i),1.)
      enddo
      call trans(Amin,dtr)
      call hpak(55,dtr)
      call trans(Astat,dtr)
      call hpak(56,dtr)
      call trans(Athr,dtr)
      call hpak(57,dtr)

      if(ifit.eq.0) goto 99 ! if there is no fit
      call hreset(1073,' ')
      do i=1,numc
       call hf1(1073,rmn(i),1.)
      enddo
      call trans(rmn,dtr)
      call hpak(1070,dtr)
      call hreset(1074,' ')
      do i=1,numc
       call hf1(1074,rsq(i),1.)
      enddo
      call trans(rsq,dtr)
      call hpak(1071,dtr)
      call trans(rst,dtr)
      call hpak(1072,dtr)

  99  continue
cc      filehbook = 'test.hbook'
      print *,' * store *  file ',filehbook
      call hrput(0,filehbook,' ')
      return
      end
c
c-----------------------------------------------------------------
      subroutine ecalib
      implicit none
#include "icalad.inc"
#include "parameters.inc"
#include "cells.inc"
#include "clb.inc"
#include "rawcells.inc"
      integer i,iad
      real Ecorr
      data Ecorr/1.0/
C
      do i=1,ncell
        iad=icalad(x(i),y(i))
        if(iad.lt.0.or.iad.ge.numc) then
          print *,'***ecalib***  x,y,iad =',x(i),y(i),iad
          cycle
        endif
        if(cal(iad+1).lt..01) then
          print *,'***ecalib*** cal(',iad+1,') =',cal(iad+1)
          cycle
        endif
        e(i)=Ecorr*acell(i)/cal(iad+1)
      enddo
      return
      END
c
      subroutine read_ecalib(filename)
      implicit none
      character*(*) filename
#include "parameters.inc"
#include "clb.inc"
      integer i
      real digits_per_GeV
      data digits_per_GeV/60./
C
      OPEN(UNIT=8,FILE=filename,STATUS='OLD',ERR=10)
      READ(8,'(16f8.2)',err=9,end=9) cal
      READ(8,'(16f8.1)',err=9,end=9) eled
      READ(8,'(16f8.1)',err=9,end=9) stat
      print *,'*read_ecalib* ',filename,' has been read'
      CLOSE(UNIT=8)
      do i=1,numc
        calsav(i)=cal(i)
      enddo
      return
    9 close(unit=8)
   10 print *,' *** READ_ECALIB ***  READ ERROR.  File ',filename,' init.'
      stop
      do i=1,numc
        cal(i)=digits_per_GeV
        calsav(i)=digits_per_GeV
      enddo
      return
      END
C
c
      subroutine write_ecalib(filename)
      implicit none
      character*(*) filename
#include "parameters.inc"
#include "clb.inc"
      integer i
C
      OPEN(UNIT=8,FILE=filename,form='formatted',status='UNKNOWN',err=2)
      WRITE(8,'(16f8.2)') (cal(i),i=1,numc)
      WRITE(8,'(16f8.1)') (eled(i),i=1,numc)
      WRITE(8,'(16f8.1)') (stat(i),i=1,numc)
      print *,'*write_ecalib* file ',filename,' is written'
      CLOSE(UNIT=8)
      RETURN
    2 write(*,*) ' ** WRITCAL_ECALIB **  OPEN failure.  File "',filename,'"'
      return
      END
c
c-----------------------------------------------------------
      subroutine readled(filename)
      implicit none
      character*(*) filename
#include "parameters.inc"
#include "leds.inc"
      integer i
C
      OPEN(UNIT=8,FILE=filename,STATUS='OLD',ERR=10)
      READ(8,*,err=9,end=9) spillled
      READ(8,*,err=9,end=9) ameanled
      READ(8,*,err=9,end=9) ledstat
      close(unit=8)
      return
    9 close(unit=8)
      print *,' *** READLED ***  READ ERROR.  File ',filename,' init.'
      do i=1,numc
        spillled(i)=0.
        ameanled(i)=0.
        ledstat(i)=0
      enddo
      return
   10 print *,' *** READLED ***  OPEN ERROR.  File ',filename,' init.'
      do i=1,numc
        spillled(i)=0.
        ameanled(i)=0.
        ledstat(i)=0
      enddo
      return
      END
C
c
      subroutine writled(filename)
      implicit none
      character*(*) filename
#include "parameters.inc"
#include "leds.inc"
      integer i
C
      OPEN(UNIT=8,FILE=filename,form='formatted',status='UNKNOWN',err=2)
      WRITE(8,'(16f8.1)') (spillled(i),i=1,numc)
      WRITE(8,'(16f8.1)') (ameanled(i),i=1,numc)
      WRITE(8,'(16i8)') (ledstat(i),i=1,numc)
      CLOSE(UNIT=8)
      RETURN
    2 write(*,*) ' ** WRITLED **  OPEN failure.  File "',filename,'"'
      return
      END
