c        Subset of programs for shower separation in ECAL2.  A.Lednev 3.11.2008.
c        Subroutines and functions for FULL shower fit.
c        ecal2reco, findclust, reorder, zone, printzone, startpoint, startp2g,
c        momenta, killgamma, inimapxy, xy.
c-----------------------------------------------------------------
c        The latest modification 19.12.2010. (killgamma)
      subroutine ecal2reco
      implicit none
#include "cell.inc"
#include "chisq.inc"
#include "isitlg.inc"
#include "parameters.inc"
      real xm2(50,50),chisqmem(numgcl)
#include "badlist.inc"
#include "cells.inc"
#include "chkbad.inc"
#include "cluster.inc"
#include "ccuts.inc"
#include "fcuts.inc"
#include "debug.inc"
#include "gamma.inc"
#include "geom.inc"
#include "iflags.inc"
#include "misc.inc"
#include "runpar.inc"
#include "thresholds.inc"
#include "zn.inc"
#include "zonesize.inc"
      logical fast
      real varfit(numgcl3),var(numgcl3),varmem(numgcl3),efit(numc)
      integer i, i1, i2, ic, icell, icl, icnt, ifl, ig, ig1, ig2, ipcl, iter, ix, iy
      integer j, k, lenclust, LG, nc, nc0, ncmem, ndf0, ngam0, nonzero
      real chilim, chisqfit, d12min, d2, dchisq, dis, dismin, distmin, dx2, dy2
      real e1, e2, eclsum, ee, Eexp, r2, scale, w, xx, yy,xzmax,xzmin,yzmax,yzmin
      real ecl(numc),xcl(numc),ycl(numc) !//AL test for shower profile
      data fast/.false./,LG/0/
c
c      if(numbev.eq.152044190) then
c        ldebug=.true.
c      else
c        ldebug=.false.
c      endif
      do i=1,ndebug
        if(numbev.eq.listdebug(i)) ldebug=.true.
      enddo
      if(ldebug) then
        print '(a,2i8)','*** ecal2reco started ***  numbev,ncell =',numbev,ncell
        print *,'ecell ='
        print '(10f7.2)',(e(i),i=1,ncell)
        print *,'iaddr ='
        print '(10i7)',(iaddr(i),i=1,ncell)
        print *,'xcell ='
        print '(10f7.1)',(x(i),i=1,ncell)
        print *,'ycell ='
        print '(10f7.1)',(y(i),i=1,ncell)
        print *,'err ='
        print '(10f7.2)',(err(i),i=1,ncell)
      endif
c
      ngam=0
      if(ncell.lt.1) then
c        ldebug=.false.
        return
      endif
      ngam0=0
      call findclust(ncell,e,iaddr,x,y)
      if(ldebug) then
        print '(a,i4)','ncell =',ncell
        print '(i2,a,99i4)',ncl,' clusters found *   lencl =',(lencl(i),i=1,ncl)
        print *,'listcells='
        print '(10i8)',(listcells(i),i=1,ncell)
        print *,'e ='
        print '(10f8.2)',(e(listcells(i)),i=1,ncell)
        print *,'x ='
        print '(10f8.2)',(x(listcells(i)),i=1,ncell)
        print *,'y ='
        print '(10f8.2)',(y(listcells(i)),i=1,ncell)
        print *,'err ='
        print '(10f8.2)',(err(listcells(i)),i=1,ncell)
        print *,'iaddr ='
        print '(10i8)',(iaddr(listcells(i)),i=1,ncell)
        print *,'ndead,deadcells=',ndead,(idead(i),i=1,ndead)
        print *,'noisy,noisycells=',noisy,(inoisy(i),i=1,noisy)
      endif
c
      ipcl=1
      lenclust=0
      do icl=1,ncl
        ipcl=ipcl+lenclust
        lenclust=lencl(icl)
        eclsum=0.
        if(ifast.ne.0) fast=.true.
        do i=ipcl,ipcl+lenclust-1
          j=listcells(i)
          eclsum=eclsum+e(j)
          if(err(j).gt.10.) fast=.false.        ! dead channel
          do k=1,noisy
            if(iaddr(j).eq.inoisy(k)-1) fast=.false.
          enddo
c          if(ldebug) print *,'icl,i,j,e(j),eclsum =',icl,i,j,e(j),eclsum
        enddo
        if(ldebug) print '(a,i3,a,f8.2)','cluster',icl,'  eclsum =',eclsum
        if(eclsum.lt.egcut) cycle
        if(ldebug) print *,' ifast =',ifast
        if(fast) then
          call fastreco(ipcl,lenclust)
          do ig=ngam0+1,ngam
            egam(ig)=egam(ig)+3.*thrcommon
c            if(ihist.ne.0.and.egam(ig).gt.egcut) call hf1(108,chigam(ig)/ndf(ig),1.)
c            if(ihist.ne.0) call hf1(109,chigam(ig)/ndf(ig),1.)
          enddo
c          print *,'***fastreco found',ngam-ngam0,' gamma in cluster',icl
c     ,  ,' e,x,y,chisq =',(egam(ig),xgam(ig),ygam(ig),chigam(ig),ig=ngam0+1,ngam)
          do ig=ngam0+1,ngam
            lnclust(ig)=lencl(icl)
          enddo
          ngam0=ngam
          cycle
        endif
        call zone(icl,e,x,y,err)
        if(ldebug) then
          print '(a,i3,a)','cluster',icl,'  e ='
          print '(10f8.2)',(ez(i),i=1,ncellz)
          print *,'  x ='
          print '(10f8.2)',(xz(i),i=1,ncellz)
          print *,'  y ='
          print '(10f8.2)',(yz(i),i=1,ncellz)
        endif
        nonzero=0
        xzmax=-1225. ! 32*38.3
        xzmin= 1225.
        yzmax=-919. ! 24*38.3
        yzmin= 919.
        do icell=1,ncellz
          efit(icell)=0.
cc          if(ez(icell).gt.thrcommon/1.9) nonzero=nonzero+1
          if(ez(icell).gt.errz(icell)) nonzero=nonzero+1
          if(xz(icell).gt.xzmax) xzmax=xz(icell)
          if(xz(icell).lt.xzmin) xzmin=xz(icell)
          if(yz(icell).gt.yzmax) yzmax=yz(icell)
          if(yz(icell).lt.yzmin) yzmin=yz(icell)
        enddo
        do i=1,3
          varmem(i)=0.
        enddo
c
        nc=0
   1    continue
        nc=nc+3
         if(ldebug) then
          print *
          print '(3(a,i3))','***  ',nc/3,' gamma search in cluster',icl,'  ***'
          print '(a,i4,a)','ncellz=',ncellz,'  input e ='
          call printzone(ez)
          print *,'  err ='
          call printzone(errz)
        endif
        do i=1,nc-3
          var(i)=varmem(i)
        enddo
        call startpoint(nc,var,efit,distmin)
        if(ldebug.and.nc.gt.3) print *,'*startpoint* distmin =',distmin
        if(nc.gt.3.and.distmin.lt..8*discut) then
          if(ldebug) print *,'* startpoint * distmin =',distmin
          nc=nc-3
          goto 2
        endif
        do ic=3,nc,3
          if(var(ic).lt.egcut/10.) then
            nc=nc-3
            goto 2
          endif
        enddo
        if(ldebug) then
          print *,'nc=',nc,'  var =',(var(i),i=1,nc)
        endif
        chisqfit=chisq(nc,var,efit)
        if(ldebug) print *,'input chisq=',chisqfit
!        chilim=(ncellz-nc)*2.5
        chilim=(nonzero-nc)*2.5
        if(chilim.lt.3.) chilim=3.

        ncmem=nc
        call showerfit(nc,var,varfit,efit,chisqfit,iter,ifl)
        do ic=3,nc,3
          if(varfit(ic-2).lt.xzmin) varfit(ic-2)=xzmin
          if(varfit(ic-2).gt.xzmax) varfit(ic-2)=xzmax
          if(varfit(ic-1).lt.yzmin) varfit(ic-1)=yzmin
          if(varfit(ic-1).gt.yzmax) varfit(ic-1)=yzmax
        enddo
        chisqfit=chisq(nc,varfit,efit)
        if(nc.lt.3) cycle
        chisqmem(nc/3)=chisqfit
        if(ldebug) then
          print *,'chisqfit,chisqmem,chilim=',chisqfit,chisqmem(nc/3),chilim,
     ,    ' iter=',iter,'  efit ='
          call printzone(efit)
          print *,'varfit =',(varfit(i),i=1,nc)
        endif
        d12min=9999.
        do i1=1,nc-3,3
          e1=varfit(i1+2)
          do i2=i1+3,nc,3
          e2=varfit(i2+2)
            d2=(varfit(i1)-varfit(i2))**2+(varfit(i1+1)-varfit(i2+1))**2
            dis=2.*sqrt(e1*e2*d2)/(e1+e2)
            if(dis.lt.d12min) then
              d12min=dis
            endif
          enddo
        enddo
        if(ldebug) print *,'chisqfit,d12min =',chisqfit,d12min
        dchisq=1.e30
        if(nc.gt.3) then
          dchisq=chisqmem(nc/3-1)-chisqfit
        endif
        if((dchisq.gt.4..and.d12min.gt.discut).or.dchisq.gt.10.) then
          nc0=0
          do ic=1,nc,3
            if(varfit(ic).lt.xzone(1)-9..or.varfit(ic).gt.xzone(nxzone)+9.) cycle
            if(varfit(ic+1).lt.yzone(1)-9..or.varfit(ic+1).gt.yzone(nyzone)+9.) cycle
            varmem(nc0+1)=varfit(ic)
            varmem(nc0+2)=varfit(ic+1)
            varmem(nc0+3)=varfit(ic+2)
            nc0=nc0+3
          enddo
          nc=nc0
        else
          if(nc.eq.ncmem) nc=nc-3
          do ic=1,nc
            varfit(ic)=varmem(ic)
          enddo
          goto 2
        endif
        if(ldebug) print *,'nc,chisqfit,chilim =',nc,chisqfit,chilim
        if(chisqfit.gt.chilim.and.nc.eq.ncmem.and.nc.lt.18) go to 1

    2   continue
!//AL test for shower profile
!        if(nc.eq.3.and.varfit(3).gt.20..and.lenclust.gt.4.and.lenclust.lt.20) then
!          print '(a,i10,2(a,i3),a)','numbev',numbev,'  cluster',icl,' Size',lenclust,' E,X,Y ='
!          print '(i2)',lenclust
!          do i=ipcl,ipcl+lenclust-1
!            print '(f8.3,2f8.2)',e(listcells(i)),x(listcells(i)),y(listcells(i))
!          enddo
!          j=0
!          do i=ipcl,ipcl+lenclust-1
!            j=j+1
!            ecl(j)=e(listcells(i))
!            xcl(j)=x(listcells(i))
!            ycl(j)=y(listcells(i))
!          enddo
!          call printclust(lenclust,ecl,xcl,ycl)
!        endif !//AL end of test for shower profile

        chisqfit=chisq(nc,varfit,efit)
        ndf0=(nonzero+ncellz)/2-nc
!        ndf0=ncellz-nc
        if(ndf0.lt.1) ndf0=1
        if(ihist.ne.0.and.nc.eq.3) call hf2(60,chisqfit/ndf0,varfit(3),1.)
        ngam=ngam0
        do i=1,nc,3
          if(ngam.lt.numg) then
            ngam=ngam+1
            xgam(ngam)=varfit(i)
            ygam(ngam)=varfit(i+1)
            egam(ngam)=varfit(i+2)
            chigam(ngam)=chisqfit
            ndf(ngam)=ndf0
          endif
        enddo
        if(ngam.ge.numg) exit

        if(ldebug) then
          print *,'nc,ngam0,ngam=',nc,ngam0,ngam
          print *,'***',ngam-ngam0,' gamma found in cluster',icl,' iter=',iter,'  efit='
          call printzone(efit)
          print *,'chisq =',chisqfit,'   ndf =',ndf(ngam)
          print *,'varfit=',(varfit(i),i=1,nc)
        endif
        if(ihist.ne.0) call hf1(12,iter+.5,1.)
        do ig=ngam0+1,ngam
          lnclust(ig)=lencl(icl)
          if(ihist.ne.0) then
            call hf1(8,chigam(ig)/ndf(ig),1.)
            if(egam(ig).gt.20.) call hf1(18,chigam(ig)/ndf(ig),1.)
          endif
        enddo
        ngam0=ngam
      enddo
cc      print *,'numbev =',numbev,'  ngam =',ngam,'  e,x,y ='
cc      print '(12f8.2)',(egam(i),xgam(i),ygam(i),i=1,ngam)

      icnt=0                ! check for dead and noisy channels
      do ix=0,63
        do iy=0,47
          icnt=icnt+1
          if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle
          xx=(ix-31.5)*dcell
          yy=(iy-23.5)*dcell
          LG=IsItLG(icnt-1)
          ee=0.
          do ig=1,ngam
            ee=ee+egam(ig)*cell(xx-xgam(ig),yy-ygam(ig),LG)
          enddo
          if(ee.gt.egcut) then
            Eexp=0.
            do j=1,ncell
              if(abs(xx-x(j)).lt.1..and.abs(yy-y(j)).lt.1.) Eexp=e(j)
            enddo
            scale=Eexp/(ee+.02)
            if(ihist.ne.0) call hfill(1000+ix,iy+.5,scale,ee)
            if(ihist.ne.0) call hf1(10000+icnt-1,Eexp,1.)
            w=1.
            est(icnt)=est(icnt)*.999 + w
            emn(icnt)=emn(icnt)+(ee-emn(icnt))/est(icnt)
            if(est(icnt).gt.10.) w=ee/emn(icnt)
            rst(icnt)=rst(icnt)*.999 + w
            if(rst(icnt).lt.1.) rst(icnt)=1.
            rmn(icnt)=rmn(icnt)+(scale-rmn(icnt))/rst(icnt)
            r2=rsq(icnt)**2
            rsq(icnt)=sqrt(r2+((scale-rmn(icnt))**2-r2)/rst(icnt))
          endif
        enddo
      enddo

      if(ihist.ne.0) call hf1(5,float(ncl),1.)
      if(ihist.ne.0) call hf1(6,float(ngam),1.)
      if(ngam.gt.1.and.ngam.le.50) then
        do i1=1,ngam-1
          do i2=i1+1,ngam
            dx2=(xgam(i1)-xgam(i2))**2
            dy2=(ygam(i1)-ygam(i2))**2
            xm2(i1,i2)=sqrt(egam(i1)*egam(i2)*(dx2+dy2))/Zecal2
            if(ihist.ne.0) call hf1(121,xm2(i1,i2),1.)
            if(ihist.ne.0.and.ngam.eq.2) call hf2(122,xm2(i1,i2),egam(i1)+egam(i2),1.)
          enddo
        enddo
      endif
      if(ldebug) then
        print *,'**',ngam,' gamma found in event',numbev,' before killgamma'
        print '(a,200f8.1)','egam =',(egam(ig),ig=1,ngam)
        print '(a,200f8.2)','xgam =',(xgam(ig),ig=1,ngam)
        print '(a,200f8.2)','ygam =',(ygam(ig),ig=1,ngam)
      endif
      call mindis(dismin)
      if(ihist.ne.0) call hf1(17,dismin,1.)
      call killgamma(egcut,discut,xm2cut)
      call mindis(dismin)
      if(ihist.ne.0) call hf1(119,dismin,1.)

      if(ldebug) then
        print *,'***',ngam,' gamma found in event',numbev
        print *,'chisq =',(chigam(ig),ig=1,ngam)
        print '(a,200f8.1)','egam =',(egam(ig),ig=1,ngam)
        print '(a,200f8.2)','xgam =',(xgam(ig),ig=1,ngam)
        print '(a,200f8.2)','ygam =',(ygam(ig),ig=1,ngam)
        if(ngam.gt.1.and.ngam.le.50) then
          do i1=1,ngam-1
            do i2=i1,ngam
              dx2=(xgam(i1)-xgam(i2))**2
              dy2=(ygam(i1)-ygam(i2))**2
              xm2(i1,i2)=sqrt(egam(i1)*egam(i2)*(dx2+dy2))/Zecal2
              xm2(i2,i1)=xm2(i1,i2)
            enddo
          enddo
          print *,' xm2 ='
          do ig1=1,ngam
            print '(20f8.3)',(xm2(ig1,ig2),ig2=1,ngam)
          enddo
          print *,'varfit=',(varfit(i),i=1,nc)
        endif
      endif
c      ldebug=.false.
      return
      end
c
c-----------------------------------------------------------------
      subroutine findclust(ncell,e,iaddr,x,y)
      implicit none
#include "parameters.inc"
      integer ncell, iaddr(numc)
      real e(numc),x(numc),y(numc)
#include "cluster.inc"
#include "debug.inc"
#include "runpar.inc"
      logical localdb
      integer iwork(numc)
      integer i, icl, iend, ij, ist, istart, itest, j, k, kl, kprev, last, last_cl, leng, next
      real d2_cut, dist2, dy_cut, xj, xk, xtest, yj, yk, ytest
      data localdb/.false./
c
c      localdb=ldebug
c      localdb=.true.
      if(localdb) then
        print '(a,i4)','* findclust *   ncell =',ncell
        print *,'  e ='
        print '(10f8.2)',(e(i),i=1,ncell)
        print *,'  iaddr ='
        print '(10i8)',(iaddr(i),i=1,ncell)
        print *,'  x ='
        print '(10f8.2)',(x(i),i=1,ncell)
        print *,'  y ='
        print '(10f8.2)',(y(i),i=1,ncell)
      endif
      ncl=0
      if(ncell.le.0.or.ncell.gt.numc) return
      do i=1,ncell
        listcells(i)=i
      enddo
      if(localdb) then
        print '(a/9(20i6/))','init. iaddr =',(iaddr(listcells(i)),i=1,ncell)
      endif
      if(ncell.eq.1) then
        ncl=1
        lencl(1)=1
        return
      endif
      call reorder(ncell,iaddr,listcells)
      if(localdb) then
        print '(a/9(20i6/))','reord. iaddr =',(iaddr(listcells(i)),i=1,ncell)
        print *,'listcells ='
        print '(10i8)',(listcells(i),i=1,ncell)
        print *,'  e ='
        print '(10f8.2)',(e(listcells(i)),i=1,ncell)
        print *,'  x ='
        print '(10f8.2)',(x(listcells(i)),i=1,ncell)
        print *,'  y ='
        print '(10f8.2)',(y(listcells(i)),i=1,ncell)
      endif
c
      next=1
      do k=2,ncell+1
        if(k.le.ncell) then
          if(k.le.ncell) kl=listcells(k)
          kprev=listcells(k-1)
          xk=x(kl)
          yk=y(kl)
          d2_cut=1.1*dcell**2
c          d2_cut=2.2*dcell**2 ! diagonal cells are included
          dy_cut=1.1*dcell
        endif
        dist2=(x(listcells(k-1))-xk)**2+(y(listcells(k-1))-yk)**2
        if(localdb) then
          print '(a,i4,a,2i4,a,4f8.2)','  k =',k,' icells=',listcells(k-1),listcells(k)
     ,     ,'  x1,y1,x2,y2=',x(listcells(k-1)),y(listcells(k-1)),xk,yk
          print *,'x,dist2,d2_cut =',x(listcells(k-1)),dist2,d2_cut
        endif
        if(dist2.le.d2_cut.and.k.le.ncell) cycle
        istart=next                ! first word of curr. (sub)cluster
        iend=k-1                ! last word of curr. (sub)cluster
        lencl(ncl+1)=k-istart        ! length of the (sub)cluster
        next=k                        ! first word of the next (sub)cluster
        ncl=ncl+1
        if(localdb) then
          print '(2(a,i4))','----- subcluster',ncl, '  iend =',iend
          print '(a,9(20i4/))','----- subcluster listcells=',(listcells(i),i=istart,iend)
        endif
        if(ncl.eq.1) cycle
        if(ncl.ge.500) return
        last=istart-1
        if(last.le.0) cycle
        if(localdb) print *,'a job to glue the subcluster',ncl
        last_cl=ncl-1
        ist=istart
        do icl=last_cl,1,-1        ! loop on cluster to test
          leng=lencl(icl)
          if(localdb) then
            print '(a,i4,a)','cluster to test =',icl,'  listcells ='
            print '(20i4)',(listcells(i),i=last-leng+1,last)
          endif
          do itest=ist,iend        ! loop on new subcluster
            kl=listcells(itest)
            xtest=x(kl)
            ytest=y(kl)
            if(localdb) then
              print '(a,i4,a,2f8.2)','test listcell',listcells(itest),'  x,y =',xtest,ytest
            endif
            dy_cut=1.1*dcell
            do ij=last,last-leng+1,-1                ! loop on cluster elements
              j=listcells(ij)
              xj=x(j)
              yj=y(j)
              dist2=(xtest-xj)**2+(ytest-yj)**2
              d2_cut=1.1*dcell**2
c              d2_cut=2.2*dcell**2 ! diagonal cells are included
              if(localdb) print *,'j,dist2=',j,dist2,' d2_cut=',d2_cut
              if(dist2.le.d2_cut) then                ! subcluster to be glued
                if(localdb) print *,'***subcluster is glued***'
                if(icl.lt.ncl-1) then
                  do j=1,leng
                    iwork(j)=listcells(last-leng+j)
                  enddo
                  do j=1,ist-1-last
                    listcells(j+last-leng)=listcells(last+j)
                  enddo
                  do j=1,leng
                    listcells(ist-leng-1+j)=iwork(j)
                  enddo
                  do j=icl,ncl-2
                    lencl(j)=lencl(j+1)
                  enddo
                endif
                lencl(ncl-1)=lencl(ncl)+leng
                ncl=ncl-1
                istart=istart-leng
                exit
              endif
            enddo
            if(dist2.le.d2_cut) exit
          enddo
          last=last-leng
        enddo
      enddo
c      ip=1
c      do icl=1,ncl
c        call reorder(lencl(icl),iaddr(ip),listcells(ip))
c        ip=ip+lencl(icl)
c      enddo
c      localdb=.false.
      if(localdb) then
        print *,ncl,' clusters found.   lencl =',(lencl(i),i=1,ncl)
        print *,' listcells =',(listcells(i),i=1,ncell)
      endif
      return
      end
c
c-----------------------------------------------------------------
      subroutine reorder(ncell,iaddr,listcells)
      implicit none
#include "parameters.inc"
      integer ncell, iaddr(numc),listcells(numc)
      integer i, k, kaddr, listk
      if(ncell.le.1) return
      do 2 k=2,ncell
        kaddr=iaddr(listcells(k))
        listk=listcells(k)
        do 1 i=k-1,1,-1
          if(kaddr.lt.iaddr(listcells(i))) then
            listcells(i+1)=listcells(i)
            go to 1
          else
            listcells(i+1)=listk
            go to 2
          endif
    1   continue
c        iaddr(1)=kaddr
        listcells(1)=listk
    2 continue
      return
      end
c
c-----------------------------------------------------------------
      subroutine zone(icl,e,x,y,err)
      implicit none
      integer icl
#include "parameters.inc"
      real e(numc),x(numc),y(numc),err(numc)
#include "errcell2.inc"
#include "cluster.inc"
#include "runpar.inc"
#include "thresholds.inc"
#include "xymap.inc"
#include "zn.inc"
#include "zonesize.inc"
      logical localdb
      integer i, icell, ip, ix, iy, j, jmin, k, kmin, ncellmemo
      real d, dfar, dmin, dnear, edelta
      data localdb/.false./,dnear,dfar/1.,40./
c
c      if(numbev.eq.1) localdb=.true.

      ip=1
      do i=1,icl-1
        ip=ip+lencl(i)
      enddo
      ncellz=0
      nxzone=0
      nyzone=0
      do icell=ip,ip+lencl(icl)-1
        j=listcells(icell)
        ncellz=ncellz+1
        xz(ncellz)=x(j)
        yz(ncellz)=y(j)
        ez(ncellz)=e(j)
        errz(ncellz)=err(j)
        ix=nxzone+1
        do i=1,nxzone
          if(abs(x(j)-xzone(i)).lt.dnear) then
            ix=-1
            exit
          endif
          if(x(j).lt.xzone(i)) then
            ix=i
            exit
          endif
        enddo
        if(ix.ne.-1) then
          do i=nxzone,ix,-1
            xzone(i+1)=xzone(i)
          enddo
          xzone(ix)=x(j)
          nxzone=nxzone+1
        endif
        iy=nyzone+1
        do i=1,nyzone
          if(abs(y(j)-yzone(i)).lt.dnear) then
            iy=-1
            exit
          endif
          if(y(j).lt.yzone(i)) then
            iy=i
            exit
          endif
        enddo
        if(iy.ne.-1) then
          do i=nyzone,iy,-1
            yzone(i+1)=yzone(i)
          enddo
          yzone(iy)=y(j)
          nyzone=nyzone+1
        endif
      enddo
      if(localdb) then
        print *,'ez=',(ez(i),i=1,ncellz)
        print *,'xz=',(xz(i),i=1,ncellz)
        print *,'yz=',(yz(i),i=1,ncellz)
        print *,'xzone=',(xzone(i),i=1,nxzone)
        print *,'yzone=',(yzone(i),i=1,nyzone)
      endif
      ncellmemo=ncellz
c      edelta=0.
      do k=1,numc
        dmin=9999.
        do j=1,ncellmemo
          d=sqrt((xz(j)-xmap(k))**2+(yz(j)-ymap(k))**2)
          if(d.lt.dmin) then
            dmin=d
            if(dmin.lt.dnear) exit
            kmin=k
            jmin=j
          endif
        enddo
        if(dmin.le.dnear) cycle
        if(dmin.gt.dfar*1.5) cycle
        ncellz=ncellz+1
        xz(ncellz)=xmap(kmin)
        yz(ncellz)=ymap(kmin)
        edelta=cellthr(kmin)*.5
        if(edelta.lt.0.03) edelta=0.03
        if(dmin.gt.dfar) edelta=edelta/2.
        ez(ncellz)=edelta/2.
        ez(ncellz)=edelta
        errz(ncellz)=edelta*.99
        if(errz(ncellz).lt.0.03) errz(ncellz)=0.03
        if(edelta.le.0.) errz(ncellz)=sqrt(errcell2(thrcommon/2.))
        ix=nxzone+1
        do i=1,nxzone
          if(abs(xmap(kmin)-xzone(i)).lt.dnear) then
            ix=-1
            exit
          endif
          if(xmap(kmin).lt.xzone(i)) then
            ix=i
            exit
          endif
        enddo
        if(ix.ne.-1) then
          do i=nxzone,ix,-1
            xzone(i+1)=xzone(i)
          enddo
          xzone(ix)=xmap(kmin)
          nxzone=nxzone+1
        endif
        iy=nyzone+1
        do i=1,nyzone
          if(abs(ymap(kmin)-yzone(i)).lt.dnear) then
            iy=-1
            exit
          endif
          if(ymap(kmin).lt.yzone(i)) then
            iy=i
            exit
          endif
        enddo
        if(iy.ne.-1) then
          do i=nyzone,iy,-1
            yzone(i+1)=yzone(i)
          enddo
          yzone(iy)=ymap(kmin)
          nyzone=nyzone+1
        endif
      enddo
      if(localdb) then
        print *,'ncellz,ncellmemo',ncellz,ncellmemo
        print '(a,200(/10f7.2))','ez =',(ez(i),i=1,ncellz)
        print '(a,200(10f7.2/))','xzone=',(xzone(i),i=1,nxzone)
        print '(a,200(/10f7.2))','xz =',(xz(i),i=1,ncellz)
        print '(a,200(10f7.2/))','yzone=',(yzone(i),i=1,nyzone)
        print '(a,200(/10f7.2))','yz =',(yz(i),i=1,ncellz)
        print *,'ncellz =',ncellz,'  ez ='
        call printzone(ez)
        print *,'nxzone =',nxzone,'  xz ='
        call printzone(xz)
        print *,'nyzone =',nyzone,'  yz ='
        call printzone(yz)
      endif
      return
      end
c
c-----------------------------------------------------------------
      subroutine printzone(p)
      implicit none
#include "parameters.inc"
      real p(numc)
      character*6 blank6(40)
      character*240 line,blank
#include "zn.inc"
#include "zonesize.inc"
      equivalence(blank6,blank)
      integer i, ix, iy
      data blank6/40*'      '/
c      print '(a,200(/10f7.2))','xzone=',(xzone(i),i=1,nxzone)
c      print '(a,200(/10f7.2))','yzone=',(yzone(i),i=1,nyzone)
c      print '(a,200(/10f7.2))','p =',(p(i),i=1,ncellz)
c      print '(a,200(/10f7.2))','xz =',(xz(i),i=1,ncellz)
c      print '(a,200(/10f7.2))','yz =',(yz(i),i=1,ncellz)
      do iy=nyzone,1,-1
        line=blank
        do ix=1,nxzone
          if(ix.gt.40) exit
          do i=1,ncellz
            if(abs(xz(i)-xzone(ix)).lt..1.and.abs(yz(i)-yzone(iy)).lt..1) then
              write(line(ix*6-5:ix*6),'(f6.2)') p(i)
            endif
          enddo
        enddo
        print *,line(1:nxzone*6)
      enddo
      return
      end
c
c-----------------------------------------------------------------
      subroutine printclust(leng,e,x,y)
      implicit none
      integer leng
      real e(*),x(*),y(*)
#include "parameters.inc"
      character*6 blank6(40)
      character*240 line,blank
      equivalence(blank6,blank)
      integer i,ix,ixmin,ixmax,iy,iymin,iymax
      real xc,yc
      data blank6/40*'      '/
      if(leng.le.0.or.leng.gt.1000) return
c      print *,'e =',(e(i),i=1,leng)
c      print *,'x =',(x(i),i=1,leng)
c      print *,'y =',(y(i),i=1,leng)
      ixmin=999
      ixmax=-999
      iymin=999
      iymax=-999
      do i=1,leng
        ix=int(x(i)/dcell+32.)
        iy=int(y(i)/dcell+24.)
        if(ix.lt.ixmin) ixmin=ix
        if(ix.gt.ixmax) ixmax=ix
        if(iy.lt.iymin) iymin=iy
        if(iy.gt.iymax) iymax=iy
      enddo
c      print *,'ixmin,ixmax,iymin,iymax =',ixmin,ixmax,iymin,iymax
      print '(a6,40i6)','  Y\X ',(i,i=ixmin,ixmax)
      do iy=iymax-iymin,0,-1
        line=blank
        write(line(1:4),'(i4)') iy+iymin
        yc=(iy+iymin-23.5)*dcell
        do ix=0,ixmax-ixmin
          if(ix.ge.40) exit
          xc=(ix+ixmin-31.5)*dcell
          do i=1,leng
            if(abs(x(i)-xc).lt..1.and.abs(y(i)-yc).lt..1) then
              write(line(ix*6+7:ix*6+12),'(f6.2)') e(i)
            endif
          enddo
        enddo
        print *,line(1:(ixmax-ixmin+2)*6)
      enddo
      return
      end
c
c-----------------------------------------------------------------
      subroutine startpoint(nc,var,efit,distmin)
      implicit none
#include "parameters.inc"
      integer nc
      real var(numgcl3), efit(numc), distmin
#include "cell.inc"
#include "icalad.inc"
#include "isitlg.inc"
#include "xy.inc"
      real varmem(200)
#include "ccuts.inc"
#include "fcuts.inc"
#include "debug.inc"
#include "gamma.inc"
#include "runpar.inc"
#include "zn.inc"
      integer i, i1, i2, iad, ic, icell, icellmax, icmem, icmin, LG
      real dderr, dderrmin, de, demax
      real dimin, dist, dist1, dist2, distcut, dx1, dx2, dy1, dy2, ee1, ee2, f, fact, fcell1, fcell2
      real sum1, sum1x, sum1y, sum9, sum9x, sum9y, x1, x2, xx1, xx2, y1, y2, yy1, yy2
      data distcut/8./,LG/0/

      distmin=0.
      demax=-999999.
      sum9=0.
      sum9x=0.
      sum9y=0.
      if(nc.le.3) then                ! one gamma
        do icell=1,ncellz
          if(ez(icell)/errz(icell).gt.demax) then
            demax=ez(icell)/errz(icell)
            icellmax=icell
            iad=icalad(xz(icell),yz(icell))
            LG=IsItLG(iad)
          endif
        enddo
        do icell=1,ncellz
          dist=sqrt((xz(icell)-xz(icellmax))**2+(yz(icell)-yz(icellmax))**2)
          if(dist.lt.1.6*dcell) then
            if(errz(icell).gt.10.) cycle
            sum9=sum9+ez(icell)
            sum9x=sum9x+ez(icell)*xz(icell)
            sum9y=sum9y+ez(icell)*yz(icell)
          endif
        enddo
        if(sum9.lt..01) then
          var(1)=0.
          var(2)=0.
          var(3)=0.
          return
        endif
        var(1)=xy(sum9x/sum9,LG)
        var(2)=xy(sum9y/sum9,LG)
c        var(1)=xz(icellmax)
c        var(2)=yz(icellmax)
        var(3)=sum9
        if(ldebug) print *,'icellmax,sum9,var=',icellmax,sum9,(var(i),i=1,3)
        return
      endif

      icellmax=0                ! two gammas or more
      do ic=1,nc
        varmem(ic)=var(ic)
      enddo
      do icell=1,ncellz
cc        if(abs(xz(icell)).lt.30..and.abs(yz(icell)).lt.30.) cycle
        if(errz(icell).gt.10.) cycle
        de=ez(icell)-efit(icell)
        dderr=sqrt(ez(icell)/errz(icell))
        if(ldebug.and.ez(icell).gt.0.) then
          print *,'icell,e,de,dderr=',icell,ez(icell),de,dderr
        endif
        if(de.lt.egcut/4.+.01) cycle
        dimin=1.e15
        do ic=3,nc-3,3
          dist=sqrt((xz(icell)-var(ic-2))**2+(yz(icell)-var(ic-1))**2)
          fact=2.*sqrt(de*var(ic))/(de+var(ic))
          if(fact.lt..5) fact=.5
          dist=dist*fact
          if(ldebug) print *,'icell,ic,dist,dist*dderr=',icell,ic,dist,dist*dderr
          if(dist*dderr.lt.dimin) then
            dimin=dist*dderr
            icmem=ic
            if(ldebug) print *,'dimin,icell=',dimin,icell
          endif
        enddo
        if(dimin.gt.distmin) then
          distmin=dimin
          dderrmin=dderr
          icellmax=icell
          icmin=icmem
          if(ldebug) print *,'   distmin,icellmax=',distmin,icellmax
        endif
      enddo
      if(icellmax.eq.0) return        ! de is too small

      distmin=distmin/dderrmin
      if(ldebug) then
        print *,'**startpoint**  nc=',nc,'  icmin,distmin=',icmin,distmin
        print '(a,i4,3f8.2)','icellmax,x,y,e (max)=',icellmax,xz(icellmax)
     ,    ,yz(icellmax),ez(icellmax)
      endif
      if(distmin.lt.distcut) return
      LG=IsItLG(icalad(xz(icellmax),yz(icellmax)))
      if(distmin.gt.100.) then
        do icell=1,ncellz
          if(errz(icell).gt.10.) cycle
          dist=sqrt((xz(icell)-xz(icellmax))**2+(yz(icell)-yz(icellmax))**2)
          if(dist.lt.1.6*dcell) then
            sum9=sum9+ez(icell)
            sum9x=sum9x+ez(icell)*xz(icell)
            sum9y=sum9y+ez(icell)*yz(icell)
          endif
        enddo
        var(nc-2)=xy(sum9x/sum9,LG)
        var(nc-1)=xy(sum9y/sum9,LG)
        var(nc)=sum9
      elseif(distmin.gt.40.) then
c        demax=ez(icellmax)-efit(icellmax)
        demax=ez(icellmax)
        sum1=0.
        sum1x=0.
        sum1y=0.
        if(ldebug) print *,'icellmax,ezmax,demax=',icellmax,ez(icellmax),demax
        do icell=1,ncellz
          if(errz(icell).gt.10..or.ez(icell).lt..01) cycle
          dx1=xz(icell)-var(icmin-2)
          dy1=yz(icell)-var(icmin-1)
          dist1=sqrt(dx1**2+dy1**2)
          dx2=xz(icell)-xz(icellmax)
          dy2=yz(icell)-yz(icellmax)
          dist2=sqrt(dx2**2+dy2**2)
          f=0.
          if(dist1.le.1.6*dcell.or.dist2.le.1.6*dcell) then
c            xi=var(icmin-2)
c            yi=var(icmin-1)
c            call corrxy(xz(icell),yz(icell),xi,yi,xi,yi)
            fcell1=efit(icell)
            fcell2=1.25*demax*cell(dx2,dy2,LG)
            if(fcell2.gt.0.) f=fcell2/(fcell1+fcell2)
            sum1=sum1+fcell1*(1.-f)
            sum1x=sum1x+fcell1*(1.-f)*xz(icell)
            sum1y=sum1y+fcell1*(1.-f)*yz(icell)
            sum9=sum9+ez(icell)*f
            sum9x=sum9x+ez(icell)*f*xz(icell)
            sum9y=sum9y+ez(icell)*f*yz(icell)
          endif
          if(ldebug) print *,'*startpoint*  f,sum9,sum1 =',f,sum9,sum1
        enddo
        if(ldebug) print *,'*startpoint*  var(icmin),sum9,sum1 =',var(icmin),sum9,sum1
        if(sum9.gt.egcut/4.) then
          var(icmin)=sum1
          var(icmin-1)=sum1y/sum1
          var(icmin-2)=sum1x/sum1
c        var(icmin-1)=(var(icmin-1)*var(icmin)-sum1y*f)/(var(icmin)-sum1*f)
c        var(icmin-2)=(var(icmin-2)*var(icmin)-sum1x*f)/(var(icmin)-sum1*f)
c        var(icmin)=var(icmin)-sum1*f
          var(nc)=sum9
          var(nc-1)=xy(sum9y/sum9,LG)
          var(nc-2)=xy(sum9x/sum9,LG)
          distmin=sqrt((var(nc-1)-var(icmin-1))**2+(var(nc-2)-var(icmin-2))**2)
          distmin=distmin*2.*sqrt(var(nc)*var(icmin))/(var(nc)+var(icmin))
          if(ldebug) print *,'*startpoint*  var =',(var(i),i=1,nc)
        else
          distmin=0.
          go to 99
        endif
      else
        if(ldebug) print *,'distmin =',distmin,'  startp2g called'
        x1=var(icmin-2)
        y1=var(icmin-1)
        x2=xz(icellmax)
        y2=yz(icellmax)
        call startp2g(x1,y1,x2,y2,xx1,yy1,ee1,xx2,yy2,ee2)
        var(icmin-2)=xx1
        var(icmin-1)=yy1
        var(icmin)=ee1
        var(nc-2)=xx2
        var(nc-1)=yy2
        var(nc)=ee2
        if(ldebug) print *,'*startpoint*  var =',(var(i),i=1,nc)
      endif
      distmin=999999.
      do i1=3,nc-3,3
        do i2=i1+3,nc,3
          dist=sqrt((var(i1-1)-var(i2-1))**2+(var(i1-2)-var(i2-2))**2)
          dist=dist*2.*sqrt(var(i1)*var(i2))/(var(i1)+var(i2))
          if(dist.lt.distmin) distmin=dist
        enddo
      enddo
      if(ldebug) print *,'**startpoint**  nc=',nc,'  distmin=',distmin
   99 if(distmin.lt.distcut) then
        do ic=1,nc
          var(ic)=varmem(ic)
        enddo
        distmin=0.
      endif
      return
      end
c
c-----------------------------------------------------------------
      subroutine startp2g(x1,y1,x2,y2,xx1,yy1,ee1,xx2,yy2,ee2)
      implicit none
      real x1, y1, x2, y2, xx1, yy1, ee1, xx2, yy2, ee2
#include "parameters.inc"
#include "zn.inc"
      integer i
      real d1, d2, distcut, dx, dy, dxy, e0, eps, r, rsq, rsq2, u, x0, y0, xx, yy, yx
      data distcut/32./
c
      call momenta(x1,y1,x2,y2,e0,x0,y0,xx,yy,yx)
c      print *,'*startp2g* e0,x0,y0,xx,yy,yx =',e0,x0,y0,xx,yy,yx
      ee2=0.
      xx2=0.
      yy2=0.
      if(ncellz.lt.4) return
      dxy = xx-yy
      rsq2= dxy**2 + 4.*yx**2
      if(rsq2.lt.1.e-30) rsq2=1.e-30
      rsq = sqrt(rsq2)
      dx = sqrt((rsq+dxy)*2.)
      dy =-sqrt((rsq-dxy)*2.)
      if(yx.ge.0.) dy=-dy
      r=sqrt(dx**2+dy**2)
c      print *,'*startp2g*  dxy,rsq,dx,dy,r=',dxy,rsq,dx,dy,r
      eps=0.
      do i=1,ncellz
        if(errz(i).gt.1.) cycle
        d1 =(xz(i)-x1)**2+(yz(i)-y1)**2
        d2 =(xz(i)-x2)**2+(yz(i)-y2)**2
        if(d1.lt.distcut.or.d2.lt.distcut) then
          u = (xz(i)-x0)*dx/r+(yz(i)-y0)*dy/r
          eps=eps-ez(i)*u*abs(u)
        endif
      enddo
      eps=eps/(e0*rsq)
      if(eps.gt.0.8) eps=0.8
      if(eps.lt.-.8) eps=-.8
      dx=dx/sqrt(1.-eps**2)
      dy=dy/sqrt(1.-eps**2)
      ee1=e0*(1+eps)/2.
      ee2=e0-ee1
      xx1=x0+dx*(1.-eps)/2.
      yy1=y0+dy*(1.-eps)/2.
      xx2=x0-dx*(1.+eps)/2.
      yy2=y0-dy*(1.+eps)/2.
c      print *,'*startp2g*  dx,dy,eps=',dx,dy,eps
      return
      end
c
c-----------------------------------------------------------------
      subroutine momenta(x1,y1,x2,y2,a0,x0,y0,xx,yy,yx)
      implicit none
      real x1, y1, x2, y2, a0, x0, y0, xx, yy, yx
c     The first and second momenta calculation.
c  Input values: x1,y1,x2,y2 to select zone   and common/zn/
c  Output values: a0,x0,y0 - Sum of energies and center-of-gravity coordinates
c                 xx,yy,yx - Second momenta
#include "parameters.inc"
#include "zn.inc"
      integer i
      real a, d1, d2, dist2cut
      data dist2cut/3200./
c
      a0=0.
      x0=0.
      y0=0.
      xx=0.
      yy=0.
      yx=0.
      if(ncellz.le.0) return
      do i=1,ncellz                !  first momenta
        if(errz(i).gt.10.) cycle
        d1 =(xz(i)-x1)**2+(yz(i)-y1)**2
        d2 =(xz(i)-x2)**2+(yz(i)-y2)**2
        if(d1.lt.dist2cut.or.d2.lt.dist2cut) then
          a0 =a0+ez(i)
          x0 =x0+xz(i)*ez(i)
          y0 =y0+yz(i)*ez(i)
        endif
      enddo
      x0 = x0/a0
      y0 = y0/a0
      do i=1,ncellz                !  second momenta
        if(errz(i).gt.10.) cycle
        d1 =(xz(i)-x1)**2+(yz(i)-y1)**2
        d2 =(xz(i)-x2)**2+(yz(i)-y2)**2
        if(d1.lt.dist2cut.or.d2.lt.dist2cut) then
          a=ez(i)/a0
          xx=xx+a*(xz(i)-x0)**2
          yy=yy+a*(yz(i)-y0)**2
          yx=yx+a*(xz(i)-x0)*(yz(i)-y0)
        endif
      enddo
c      print *,'*momenta*  a0,x0,y0,xx,yy,yx=',a0,x0,y0,xx,yy,yx
      return
      end
c
c-----------------------------------------------------------------------
      subroutine killgamma(egc,disc,xm2c)
      implicit none
#include "parameters.inc"
#include "cells.inc"
#include "debug.inc"
#include "gamma.inc"
#include "geom.inc"
#include "misc.inc"
      integer i, ig, ig1, ig1min, ig2, ig2min, ix, iy, ngam0
      real dis, dismin, dist, dist2g, distig, dnear, xm12, xm2g, xm2min, ym12
      real egc, disc, xm2c
      if(ldebug) print *,'*killgam started* egcut,discut,xm2cut =',egc,disc,xm2c
      ngam0=ngam
      ngam=0
      do ig=1,ngam0
        if(egam(ig).ge.egc) then
          ngam=ngam+1
          egam(ngam)=egam(ig)
          xgam(ngam)=xgam(ig)
          ygam(ngam)=ygam(ig)
          chigam(ngam)=chigam(ig)
          gtime(ngam)=gtime(ig)
          ndf(ngam)=ndf(ig)
          lnclust(ngam)=lnclust(ig)
          dnear=1.E30
          do i=1,ncell
            ix=iaddr(i)/48
            iy=iaddr(i)-ix*48
            if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle
            distig=(xgam(ig)-X(i))**2+(ygam(ig)-Y(i))**2
            if(distig.lt.dnear.and.E(i).gt.0.) then
              dnear=distig
              indmc(ngam)=iaddr(ig)
              gtime(ngam)=time(ig)
            endif
          enddo
        endif
      enddo
    1 if(ngam.lt.2) return
      dismin=999999.
      do ig1=1,ngam-1
        do ig2=ig1+1,ngam
          dist=sqrt((xgam(ig1)-xgam(ig2))**2+(ygam(ig1)-ygam(ig2))**2)
          dis=2.*sqrt(egam(ig1)*egam(ig2))*dist/(egam(ig1)+egam(ig2))
          if(dis.lt.dismin) then
            ig1min=ig1
            ig2min=ig2
            dismin=dis
          endif
        enddo
      enddo
      if(ldebug) print *,'*killgamma*  dismin =',dismin
      if(dismin.le.disc) then
        xm12=xgam(ig1min)*egam(ig1min)+xgam(ig2min)*egam(ig2min)
        ym12=ygam(ig1min)*egam(ig1min)+ygam(ig2min)*egam(ig2min)
        if(egam(ig2min).gt.egam(ig1min)) indmc(ig1min)=indmc(ig2min)
        egam(ig1min)=egam(ig1min)+egam(ig2min)
        xgam(ig1min)=xm12/egam(ig1min)
        ygam(ig1min)=ym12/egam(ig1min)
        if(chigam(ig2min).gt.chigam(ig1min)) chigam(ig1min)=chigam(ig2min)
        if(egam(ig2min).gt.egam(ig1min)) then
          indmc(ig1min)=indmc(ig2min)
          gtime(ig1min)=gtime(ig2min)
        endif
        ngam=ngam-1
        do ig=ig2min,ngam
          egam(ig)=egam(ig+1)
          xgam(ig)=xgam(ig+1)
          ygam(ig)=ygam(ig+1)
          chigam(ig)=chigam(ig+1)
          gtime(ig)=gtime(ig+1)
          indmc(ig)=indmc(ig+1)
          ndf(ig)=ndf(ig+1)
          lnclust(ig)=lnclust(ig+1)
        enddo
        if(ldebug) then
          print *,'*killgamma*   1 gamma is killed.  ig1,ig2 =',ig1min,ig2min
          print '(a,200f8.1)','egam =',(egam(ig),ig=1,ngam)
          print '(a,200f8.2)','xgam =',(xgam(ig),ig=1,ngam)
          print '(a,200f8.2)','ygam =',(ygam(ig),ig=1,ngam)
        endif
        goto 1
      endif
    2 if(ngam.lt.2) return
      xm2min=999999.
      dist2g=999999.
      do ig1=1,ngam-1
        do ig2=ig1+1,ngam
          dist=sqrt((xgam(ig1)-xgam(ig2))**2+(ygam(ig1)-ygam(ig2))**2)
          xm2g=sqrt(egam(ig1)*egam(ig2))*dist/Zecal2
          if(xm2g.lt.xm2min) then
            ig1min=ig1
            ig2min=ig2
            xm2min=xm2g
            dist2g=dist
          endif
        enddo
      enddo
      if(ldebug) print *,'*killgamma*  dist2g,xm2min =',dist2g,xm2min
      if(dist2g.le.disc*2..and.xm2min.lt.xm2c) then
        xm12=xgam(ig1min)*egam(ig1min)+xgam(ig2min)*egam(ig2min)
        ym12=ygam(ig1min)*egam(ig1min)+ygam(ig2min)*egam(ig2min)
        if(egam(ig2min).gt.egam(ig1min)) indmc(ig1min)=indmc(ig2min)
        egam(ig1min)=egam(ig1min)+egam(ig2min)
        xgam(ig1min)=xm12/egam(ig1min)
        ygam(ig1min)=ym12/egam(ig1min)
        if(chigam(ig2min).gt.chigam(ig1min)) chigam(ig1min)=chigam(ig2min)
        if(egam(ig2min).gt.egam(ig1min)) then
          indmc(ig1min)=indmc(ig2min)
          gtime(ig1min)=gtime(ig2min)
        endif
        ngam=ngam-1
        do ig=ig2min,ngam
          egam(ig)=egam(ig+1)
          xgam(ig)=xgam(ig+1)
          ygam(ig)=ygam(ig+1)
          chigam(ig)=chigam(ig+1)
          gtime(ig)=gtime(ig+1)
          indmc(ig)=indmc(ig+1)
          ndf(ig)=ndf(ig+1)
          lnclust(ig)=lnclust(ig+1)
        enddo
        if(ldebug) then
          print *,'*killgamma*   1 gamma is killed'
          print '(a,200f8.1)','egam =',(egam(ig),ig=1,ngam)
          print '(a,200f8.2)','xgam =',(xgam(ig),ig=1,ngam)
          print '(a,200f8.2)','ygam =',(ygam(ig),ig=1,ngam)
        endif
        goto 2
      endif
      if(ldebug) print *,'ngam0',' gamma was before killgamma ',ngam,' is after'
      return
      end
c
c-----------------------------------------------------------------------
      subroutine inimapxy
      implicit none
#include "parameters.inc"
#include "xymap.inc"
      integer i, ij, j
      real xi
      ij=0
      do i=1,64
        xi=(i-32.5)*dcell
        do j=1,48
         ij=ij+1
         ymap(ij)=(j-24.5)*dcell
         xmap(ij)=xi
        enddo
      enddo
      return
      end
c
c-----------------------------------------------------------------------
      real function xy(xc,LG)
      implicit none
      integer LG
      real xc
#include "parameters.inc"
      real x, x2
c   The function calculates a shower coordinate (x or y) in ECAL2 shashlyk.
c   XC is the center of gravity of the shower.
C   Made from 40 GeV electron calibration data. ( Wide beam.)  Run 2008.
      x=XC/dcell-nint(XC/dcell)
      x2=x*x
      if(LG.eq.1) then
c        xy=XC-104.7*x*(x2*(x2+.2177)+.0611)*(x2-.25)*(x2-.392)*dcell
        xy=XC-119.*x*(x2*(x2+.149)+.0413)*(x2-.25)*(x2-.382)*dcell ! LG
      else        ! Sha
c        xy=XC-102.6*x*(x2*(x2+.1536)+.0589)*(x2-.25)*(x2-.391)*dcell
        xy=XC-109.2*x*(x2*(x2+.168)+.0441)*(x2-.25)*(x2-.383)*dcell ! Sha
      endif
      RETURN
      END
