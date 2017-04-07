c        Subset of programs for shower separation in ECAL2.  A.Lednev 21.10.2008.
c        Subroutines and functions for FULL shower fit.
c        showerfit, constr, chisq, gradient.
c-----------------------------------------------------------------
      subroutine showerfit(nc,var,varfit,efit,chisqfit,iter,ifl)
      implicit none
#include "parameters.inc"
      integer nc, iter, ifl
      real var(numgcl3), varfit(numgcl3), efit(numc), chisqfit
#include "chisq.inc"
#include "constr.inc"
#include "debug.inc"
#include "flags.inc"
#include "runpar.inc"
#include "zn.inc"
#include "zonesize.inc"
      logical lder
      real der(numgcls),dermem(numgcls),solut(numgcl),delta(numgcl)
      integer i, ic, ir(numgcl), ivar, maxiter, ncmem, nfail, niter
      real chisqmem, chisqold, deltalim, emem(numc), factor, step
      real varmem(numgcl3)
      data maxiter/30/,deltalim/0.03/
      data lder/.true./
c
c      ldebug=.true.
      do i=1,ndebug
        if(numbev.eq.listdebug(i)) ldebug=.true.
      enddo
      if(ldebug) then
        print *,'*** showerfit ***  nxzone,nyzone =',nxzone,nyzone
        print *,'*** showerfit ***  input e ='
        call printzone(ez)
c        print *,'*** showerfit ***  input err ='
c        call printzone(errz)
c        print *,'*** showerfit ***  input x ='
c        call printzone(xz)
c        print *,'*** showerfit ***  input y ='
c        call printzone(yz)
c        print *,'*** showerfit ***  input dcell ='
c        call printzone(dcell)
      endif
      iter=0
      nfail=0
      step=1.
      do ivar=1,nc
        varfit(ivar)=var(ivar)
        varmem(ivar)=var(ivar)
      enddo
      chisqfit=chisq(nc,varfit,efit)
      if(ldebug) then
        print *
        print *,'* showerfit * init varfit=',(varfit(i),i=1,nc)
      endif
      if(ldebug) then
        print *,'*** showerfit ***  input efit ='
        call printzone(efit)
        print *,'*** showerfit ***  input chisqfit =',chisqfit
      endif
      call gradient(nc,varmem,varfit,efit,chisqfit,niter,5,ifl)
      if(ldebug) then
        print *,'gradient chisqfit =',chisqfit,'  *** niter =',niter
        print *,'gradient varfit =',(varfit(i),i=1,nc)
      endif
c      if(chisqfit.lt.chilim) goto 10
      do ivar=1,nc
        varmem(ivar)=varfit(ivar)
      enddo
      chisqmem=chisq(nc,varmem,emem)
    1 if(iter.ge.maxiter) goto 10
      iter=iter+1
      chisqold=chisqmem
      call vzero(der,numgcls)
      call vzero(solut,numgcl)
      chisqfit=chisq(nc,varfit,efit)
      do ic=1,nc
        delta(ic)=constr(ic,nc,varfit,efit,der,lder)
      enddo
      do i=1,nc**2
        dermem(i)=der(i)
      enddo
      if(ldebug) then
        print *
        print *,'*** showerfit *** iter=',iter
        print *,'  efit ='
        call printzone(efit)
c        do i=1,ncell
c          test(i)=((e(i)-efit(i))/err(i))**2
c        enddo
c        print *,'  test ='
c        call printzone(test)
c        print *,'  efit ='
c        print '(15f6.2)',(efit(i),i=1,ncell)
c        icell=ncell
c        do iy=iyzone,1,-1
c          print '(20f7.2)',(efit(i),i=iy,icell,iyzone)
c          icell=icell-1
c        enddo
        print *,'*** chisqfit =',chisqfit
        print *,'varfit=',(varfit(i),i=1,nc)
        print *,'delta =',(delta(ic),ic=1,nc)
      endif
    3 do ic=1,nc
        if(abs(delta(ic)).gt.deltalim) goto 2
      enddo
      goto 10
    2 continue
c
      if(ldebug) then
        print *,'* showerfit *  der ='
        do ic=1,nc
          print '(9f9.4)',(der(i),i=ic*nc-nc+1,ic*nc)
        enddo
      endif
c
      do ic=1,nc
        solut(ic)=delta(ic)
      enddo
      call reqn(nc,der,nc,ir,ifl,1,solut)
      do i=1,nc**2        ! der was destroyed by reqn
        der(i)=dermem(i)
      enddo
      if(ifl.ne.0) then
        if(ldebug) print *,'* showerfit * matrix inversion failed'
        do ivar=1,nc
          varmem(ivar)=varfit(ivar)
        enddo
        ncmem=nc
        call gradient(nc,varmem,varfit,efit,chisqfit,niter,20,ifl)
        if(nc.ne.ncmem) goto 10
        do ivar=1,nc
          varmem(ivar)=varfit(ivar)
        enddo
        chisqmem=chisq(nc,varmem,emem)
        do ic=1,nc
          delta(ic)=constr(ic,nc,varmem,emem,der,lder)
        enddo
        if(ldebug) then
          print *,'gradient chisqfit =',chisqmem,'  *** niter =',niter
          print *,'gradient varmem=',(varmem(i),i=1,nc)
          print *,'gradient delta=',(delta(i),i=1,nc)
        endif
        if(ifl.ne.0) then
          if(ldebug) print *,'gradient ifl =',ifl
          goto 10
        endif
        goto 3
      else
        nfail=0
      endif
      if(ldebug) print *,'* showerfit * solut=',(solut(i),i=1,nc)
c
      factor=1.
      do ivar=3,nc,3
        if(0.9*varmem(ivar).le.solut(ivar)) then
          factor=factor*0.9*varmem(ivar)/solut(ivar)
         if(ldebug)print *,'*** showerfit ***  factor=',factor,' !! ***'
        endif
      enddo
      do ivar=1,nc
        varmem(ivar)=varmem(ivar)-solut(ivar)*factor
      enddo
      chisqmem=chisq(nc,varmem,emem)
c
      if(ldebug) then
        print *,'varmem=',(varmem(i),i=1,nc)
        print *,'* showerfit * chisqmem=',chisqmem
c        print *,'*showerfit* efit ='
c         icell=ncell
c         do iy=iyzone,1,-1
c           print '(20f7.2)',(e(i),i=iy,icell,iyzone)
c          icell=icell-1
c         enddo
      endif
c
      if(chisqmem.le.chisqfit) then
        do ivar=1,nc
          varfit(ivar)=varmem(ivar)
        enddo
        chisqfit=chisq(nc,varfit,efit)
      else                ! change minimization algorithm
        factor=.3
        do ivar=3,nc,3
          if(0.9*varfit(ivar).le.-solut(ivar)*factor) then
            factor=0.9*varfit(ivar)/abs(solut(ivar))
           if(ldebug)print *,'*** showerfit ***  factor=',factor,' !! ***'
          endif
        enddo
        do ivar=1,nc
          varmem(ivar)=varfit(ivar)+solut(ivar)*factor
        enddo
        chisqmem=chisq(nc,varmem,emem)
        if(ldebug) then
          print *,'*** new chiqmem =',chisqmem,'  factor =',factor
          print *,'*** varmem =',(varmem(ivar),ivar=1,nc)
        endif
        if(chisqmem.le.chisqfit) then
          do ivar=1,nc
            varfit(ivar)=varmem(ivar)
          enddo
          goto 1
        else
          factor=.1
          do ivar=3,nc,3
            if(0.9*varfit(ivar).le.solut(ivar)*factor) then
              factor=0.9*varfit(ivar)/solut(ivar)
             if(ldebug)print *,'*** showerfit ***  factor=',factor,' !! ***'
            endif
          enddo
          do ivar=1,nc
            varmem(ivar)=varfit(ivar)-solut(ivar)*factor
          enddo
          chisqmem=chisq(nc,varmem,emem)
          if(ldebug) then
            print *,'*** new chiqmem =',chisqmem,'  factor =',factor
            print *,'*** varmem =',(varmem(ivar),ivar=1,nc)
          endif
          if(chisqmem.le.chisqfit) then
            do ivar=1,nc
              varfit(ivar)=varmem(ivar)
            enddo
            goto 1
          endif
        endif
        do ivar=1,nc
          varmem(ivar)=varfit(ivar)
        enddo
        call gradient(nc,varmem,varfit,efit,chisqfit,niter,20,ifl)
        if(ifl.ne.0) then
          if(ldebug) print *,'*gradient* ifl =',ifl
          goto 10
        endif
        do ivar=1,nc
          varmem(ivar)=varfit(ivar)
        enddo
        chisqmem=chisq(nc,varmem,emem)
        do ic=1,nc
          delta(ic)=constr(ic,nc,varmem,emem,der,lder)
        enddo
        if(ldebug) then
          print *,'gradient chisqfit =',chisqmem,'  *** niter =',niter
          print *,'gradient varmem=',(varmem(i),i=1,nc)
          print *,'gradient delta=',(delta(i),i=1,nc)
        endif
      endif
      if(chisqfit.lt.chisqold-.1) go to 1
      go to 1
   10 continue
c
      chisqfit=chisq(nc,varfit,efit)
      if(ldebug) then
        print *,'* showerfit * iter =',iter,'  chisqfit =',chisqfit
        print *,'* showerfit * e ='
        call printzone(ez)
c        icell=ncell
c        do iy=iyzone,1,-1
c          print '(20f7.2)',(e(i),i=iy,icell,iyzone)
c          icell=icell-1
c        enddo
        print *,'* showerfit * efit ='
        call printzone(efit)
c        icell=ncell
c        do iy=iyzone,1,-1
c          print '(20f7.2)',(efit(i),i=iy,icell,iyzone)
c          icell=icell-1
c        enddo
      endif
      if(ldebug) print *,'*showerfit finished*   chisqfit =',chisqfit
      if(ldebug) print *,'* showerfit * delta =',(delta(i),i=1,nc)
      if(ldebug) print *,'* showerfit * varfit =',(varfit(i),i=1,nc)

      return
      end

c-----------------------------------------------------------------
      real function constr(ic,nc,var,efit,der,lder)
      implicit none
#include "parameters.inc"
      integer ic, nc
      logical lder
      real var(numgcl3), efit(numc), der(numgcls)
c      chisq = Sum (efit(i)-e(i))/err(i)**2
c      constr(ic) = d(chisq)/d(var(ic))
#include "cell.inc"
#include "dcelldx.inc"
#include "dcelldy.inc"
#include "d2celldx2.inc"
#include "d2celldy2.inc"
#include "d2celldxdy.inc"
#include "icalad.inc"
#include "isitlg.inc"
#include "zn.inc"
      real dx(numgcl),dy(numgcl),ee(numgcl)
      integer i, i3, icc, ig, ig1, ig2, k, LG
      real c, cell1, d2cdx1, d2cdxdy1, d2cdy1, dcdx1, dcdy1, de, e1, er2
      data LG/0/
c
      c=0.
      if(lder) then
        do i=ic,nc
          der((ic-1)*nc+i)=0.
        enddo
      endif
      icc=ic-(ic-1)/3*3
      ig1=(ic-1)/3+1
      e1=var(3*ig1)
      do 99 i=1,ncellz
        de=efit(i)-ez(i)
        er2=errz(i)**2
        LG=IsItLG(icalad(xz(i),yz(i)))
        do ig=1,nc/3
          dx(ig)=var((ig-1)*3+1)-xz(i)
          dy(ig)=var((ig-1)*3+2)-yz(i)
cc          call corrxy(x(i),y(i),var(1,ig),var(2,ig),xc(ig),yc(ig))
cc          dx(ig)=xc(ig)-x(i)
cc          dy(ig)=yc(ig)-y(i)
          ee(ig)=var(3*ig)*cell(dx(ig),dy(ig),LG)
        enddo
        cell1=cell(dx(ig1),dy(ig1),LG)
        dcdx1=dcelldx(dx(ig1),dy(ig1),LG)
        dcdy1=dcelldy(dx(ig1),dy(ig1),LG)
c
        if(icc.eq.1) c=c+de*e1*dcdx1/er2
        if(icc.eq.2) c=c+de*e1*dcdy1/er2
        if(icc.eq.3) c=c+de*cell1/er2
c
        if(.not.lder) go to 99
        d2cdx1=d2celldx2(dx(ig1),dy(ig1),LG)
        d2cdy1=d2celldy2(dx(ig1),dy(ig1),LG)
        d2cdxdy1=d2celldxdy(dx(ig1),dy(ig1),LG)
        k=(ic-1)*nc
        if(icc.eq.1) then
          i3=ig1*3-2
          der(k+i3)=der(k+i3)+e1*(de*d2cdx1+e1*dcdx1**2)/er2
          der(k+i3+1)=der(k+i3+1)+e1*(de*d2cdxdy1+e1*dcdx1*dcdy1)/er2
          der(k+i3+2)=der(k+i3+2)+dcdx1*(ee(ig1)+de)/er2
c          do ig2=1,nc/3
c            if(ig2.ne.ig1) der(k+ig1*3-2)=der(k+ig1*3-2)+
c     +        e1*de*d2celldx2(dx(ig2),dy(ig2),LG)/er2
c          enddo
          do ig2=ig1+1,nc/3
            i3=ig2*3-2
            der(k+i3)=der(k+i3)+
     +       e1*var(3*ig2)*dcdx1*dcelldx(dx(ig2),dy(ig2),LG)/er2
            der(k+i3+1)=der(k+i3+1)+
     +       e1*var(3*ig2)*dcdx1*dcelldy(dx(ig2),dy(ig2),LG)/er2
            der(k+i3+2)=der(k+i3+2)+
     +       e1*dcdx1*cell(dx(ig2),dy(ig2),LG)/er2
          enddo
        elseif(icc.eq.2) then
          i3=ig1*3-2
          der(k+i3+1)=der(k+i3+1)+e1*(de*d2cdy1+e1*dcdy1**2)/er2
          der(k+i3+2)=der(k+i3+2)+dcdy1*(ee(ig1)+de)/er2
c          do ig2=1,nc/3
c            if(ig2.ne.ig1) der(k+ig1*3-1)=der(k+ig1*3-1)+
c     +        e1*de*d2celldy2(dx(ig2),dy(ig2),LG)/er2
c          enddo
          do ig2=ig1+1,nc/3
            i3=ig2*3-2
            der(k+i3)=der(k+i3)+
     +       e1*var(3*ig2)*dcdy1*dcelldx(dx(ig2),dy(ig2),LG)/er2
            der(k+i3+1)=der(k+i3+1)+
     +       e1*var(3*ig2)*dcdy1*dcelldy(dx(ig2),dy(ig2),LG)/er2
            der(k+i3+2)=der(k+i3+2)+
     +       e1*dcdy1*cell(dx(ig2),dy(ig2),LG)/er2
          enddo
        elseif(icc.eq.3) then
          i3=ig1*3-2
          der(k+i3+2)=der(k+i3+2)+cell1**2/er2
          do ig2=ig1+1,nc/3
            i3=ig2*3-2
            der(k+i3)=der(k+i3)+
     +        var(3*ig2)*dcelldx(dx(ig2),dy(ig2),LG)*cell1/er2
            der(k+i3+1)=der(k+i3+1)+
     +        var(3*ig2)*dcelldy(dx(ig2),dy(ig2),LG)*cell1/er2
            der(k+i3+2)=der(k+i3+2)+
     +        cell1*cell(dx(ig2),dy(ig2),LG)/er2
          enddo
        endif
   99 continue
      constr=c
      if(lder) then
        do i=ic+1,nc                ! the matrix must be symmetric.
          der((i-1)*nc+ic)=der((ic-1)*nc+i)
        enddo
      endif
      return
      end

c-----------------------------------------------------------------
      real function chisq(nc,var,efit)
      implicit none
#include "parameters.inc"
      integer nc
      real var(numgcl3), efit(numc)
#include "cell.inc"
#include "icalad.inc"
#include "isitlg.inc"
#include "zn.inc"
      integer i, k, LG
      real dx, dy
      data LG/0/
c        print *,'*** chisq ***  input var =',(var(i),i=1,nc)
c         print *,'*** chisq ***  input dcell ='
c        call printzone(dcell)
c         print *,'*** chisq ***  input e ='
c        call printzone(e)
      do i=1,ncellz
        efit(i)=0.
        LG=IsItLG(icalad(xz(i),yz(i)))
        do k=1,nc,3
          dx=var(k)-xz(i)
          dy=var(k+1)-yz(i)
cc          call corrxy(x(i),y(i),var(k),var(k+1),xc(k),yc(k))
cc          dx=xc(k)-x(i)
cc          dy=yc(k)-y(i)
          efit(i)=efit(i)+var(k+2)*cell(dx,dy,LG)
        enddo
c        test(i)=(e(i)-efit(i))/err(i)
      enddo
      chisq=0.
      do i=1,ncellz
        chisq=chisq+((efit(i)-ez(i))/errz(i))**2
      enddo
c        print *,'*** chisq ***  input efit ='
c        call printzone(efit)
c        print *,'*** chisq ***  test ='
c        call printzone(test)
      return
      end

c-----------------------------------------------------------------
      subroutine gradient(nc,varmem,varfit,efit,chisqfit,niter,maxi,ifl)
      implicit none
#include "parameters.inc"
      integer nc, niter, maxi, ifl
      real varmem(numgcl3), varfit(numgcl3), efit(numc), chisqfit
#include "chisq.inc"
#include "constr.inc"
#include "debug.inc"
#include "zn.inc"
      real delta(numgcl3),deltaold(numgcl3),der(numgcls),varini(numgcl3),emem(numc)
      integer i, ic, j, k, ncmem
      real absdelta, absdeltacut, absdeltaold, chisqini, chisqmem, cosi
      real dchi, deltasq, dlim, scal, step, stepold, sum1, sum2
      data absdeltacut/0.05/
c
c      ldebug=.true.
      ifl=1
      ncmem=nc
      niter=0
      do ic=1,nc
        varini(ic)=varmem(ic)
        varfit(ic)=varmem(ic)
      enddo
      chisqfit=chisq(nc,varini,efit)
      chisqini=chisqfit
      do ic=1,nc
        delta(ic)=constr(ic,nc,varfit,efit,der,.true.)
        deltaold(ic)=delta(ic)
      enddo
      deltasq=0.
      do i=1,nc
        deltasq=deltasq+delta(i)**2
      enddo
      absdelta=sqrt(deltasq)
      absdeltaold=absdelta
      sum1=0.
      sum2=0.
      do i=1,nc
        sum1=sum1+delta(i)**2
        do j=1,nc
          sum2=sum2+der((i-1)*nc+j)*delta(i)*delta(j)
        enddo
      enddo
      if(ldebug) print *,'sum1,sum2=',sum1,sum2
      if(sum2.le.sum1) sum2=sum1
      step=sum1/sum2
      if(ldebug) then
        print *
        print *,'*** gradient ***    init chisq =',chisqfit
        print *,'init. varfit =',(varfit(i),i=1,nc)
        print *,'absdelta =',absdelta
        print *,'init. delta =',(delta(i),i=1,nc)
        print *,'step=',step,'    der ='
        print '(12f8.3)',(der(i),i=1,nc**2)
        print *
      endif
    1 niter=niter+1
      if(ldebug) then
        print *
        print *,'*** niter =',niter,'  varmem =',(varmem(i),i=1,nc)
        print *,'step=',step,'  delta =',(delta(i),i=1,nc)
        print '(a,(12f10.3))','der=',(der(i),i=1,nc**2)
      endif
      do i=3,nc                ! To prevent too big step
        dlim=dcell
        if(i-i/3*3.eq.0) dlim=varmem(i)*.9
        if(abs(delta(i))*step.ge.dlim) then
          if(ldebug) print *,'i,delta,varmem',i,delta(i),varmem(i)
          do while(abs(delta(i))*step.ge.dlim)
            step=step/2.
          enddo
          if(ldebug) then
            print *,' *** negative energy prevented ***',' step=',step
          endif
          if(delta(i)*step.lt..001) then
            if(ldebug) then
              print *,'*** Zero energy. One gamma is killed ***'
              print *,'varfit =',(varfit(k),k=1,nc)
            endif
            k=0
            do j=1,nc
              if(j.lt.i-2.or.j.gt.i) then
                k=k+1
                varmem(k)=varfit(j)
              endif
            enddo
            nc=nc-3
            if(ldebug) print *,'varmem =',(varmem(k),k=1,nc)
            do k=1,nc
              varfit(k)=varmem(k)
            enddo
            goto 10
          endif
        endif
      enddo
      do i=1,nc
        varmem(i)=varmem(i)-delta(i)*step
cc        varmem(i)=varmem(i)-delta(i)/absdelta*step
      enddo
      chisqmem=chisq(nc,varmem,emem)
      do ic=1,nc
cc        delta(ic)=constr(ic,nc,varmem,emem,der,.false.)
        delta(ic)=constr(ic,nc,varmem,emem,der,.true.)
      enddo
      deltasq=0.
      do i=1,nc
        deltasq=deltasq+delta(i)**2
      enddo
      absdelta=sqrt(deltasq)
      if(absdelta.lt.absdeltacut) go to 10
      if(absdelta*step.lt.absdeltacut.and.niter.ge.5) go to 10
      if(ldebug) then
        print *
        if(chisqmem.ge.chisqfit) print *,'******'
        print *,'absdelta =',absdelta
        print *,'delta =',(delta(i),i=1,nc)
        print *,'varmem =',(varmem(i),i=1,nc)
        print *,'  chisqmem=',chisqmem,'  chisqfit=',chisqfit
      endif
      stepold=step
      if(chisqmem.le.chisqfit) then
        do i=1,nc
         varfit(i)=varmem(i)
        enddo
        chisqfit=chisq(nc,varfit,efit)
        cosi=0.
        do ic=1,nc
          cosi=cosi+delta(ic)*deltaold(ic)
        enddo
        cosi=cosi/(absdelta*absdeltaold)
c        scal=(1.-cosi)*absdeltaold/absdelta
        scal=(2.*absdeltaold*(1.+cosi)+absdelta)/(absdeltaold+absdelta)
        if(scal.lt..1) scal=.25
        if(scal.gt.4.) scal=4.
c        step=step/scal
        step=step*scal
        do ic=1,nc
          deltaold(ic)=delta(ic)
        enddo
        absdeltaold=absdelta
      else
        do i=1,nc
         varmem(i)=varfit(i)
         delta(i)=deltaold(i)
        enddo
        dchi=chisqmem-chisqfit
cc        step=stepold/(dchi/stepold/absdeltaold+1.)/2.
cc        if(ldebug) print *,'dchi=',dchi,'  stepold,step =',stepold,step
cc        if(step.lt.stepold*.1) step=stepold*.1
        step=step/4.
        if(step.gt..001) then
          goto 1
        else
          goto 10
        endif
      endif
      if(ldebug) print *,'cosi =',cosi,'  scal =',scal,'  step =',step

      sum1=0.
      sum2=0.
      do i=1,nc
        sum1=sum1+delta(i)**2
        do j=1,nc
          sum2=sum2+der((i-1)*nc+j)*delta(i)*delta(j)
        enddo
      enddo
      if(ldebug) print *,'sum1,sum2=',sum1,sum2
      if(sum2.le.sum1) sum2=sum1
      step=sum1/sum2

      if(niter.lt.maxi) go to 1
   10 chisqfit=chisq(nc,varfit,efit)
      if(ldebug) print *,'absdelta =',absdelta,'  step =',step,
     ,  '  chisqfit =',chisqfit,'  niter =',niter
      if(chisqini-chisqfit.gt.0.001) then
       ifl=0
      else
        ifl=1
        nc=ncmem
        do ic=1,nc
          varfit(ic)=varini(ic)
        enddo
        chisqfit=chisqini
        if(ldebug) then
          print *,'***gradient failed***'
          print *,'varfit =',(varfit(i),i=1,nc)
        endif
      endif
      return
      end
