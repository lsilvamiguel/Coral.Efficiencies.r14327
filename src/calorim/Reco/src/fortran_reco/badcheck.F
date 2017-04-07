c        Subset of programs for shower separation in ECAL2.  A.Lednev 3.11.2008.
c        Check for dead and noisy channels.
c        badcheck, readbad, writbad.

      subroutine badcheck
      implicit none
#include "parameters.inc"
      logical hexist
      character*24 tit
#include "badlist.inc"
#include "cells.inc"
#include "chkbad.inc"
#include "debug.inc"
#include "flags.inc"
#include "gclb.inc"
#include "geom.inc"
#include "thresholds.inc"
#include "xymap.inc"
      integer i, iad, ini, ix, iy, ic
      real rstcut, rmncut, rsqcut
      data rstcut/5./,rmncut/0.05/,rsqcut/5./,ini/1/

      if(ini.ne.0) then
        do i=1,numc
          emn(i)=0.
          est(i)=0.
          rmn(i)=0.
          rsq(i)=0.
          rst(i)=0.
        enddo
        if(lbad) call readbad
        ini=0
      endif
      ndead=0
      noisy=0
      do i=1,numc
        ix=(i-1)/48
        iy=i-1-ix*48
        if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) then ! hole
          if(ix-iholel.eq.1.or.iholer-ix.eq.1.or.iy-iholed.eq.1.or.iholeu-iy.eq.1) then
            rmn(i)=0.
            rst(i)=0.
            ncell=ncell+1
            e(ncell)=0.
            iaddr(ncell)=i-1
            x(ncell)=xmap(i)
            y(ncell)=ymap(i)
            err(ncell)=.555
            time(ncell)=0.
            cycle
          endif
        endif
        if(ix.lt.44.and.iy.eq.47) cycle                ! empty row
        if(rst(i).gt.rstcut) then
          if(rmn(i).lt.rmncut) then
            ndead=ndead+1
            idead(ndead)=i
            if(lbad) then
              do ic=1,ncell
                if(iaddr(ic).eq.i-1) exit
                ncell=ncell+1
                e(ncell)=0.
                iaddr(ncell)=i-1
                x(ncell)=xmap(i)
                y(ncell)=ymap(i)
                err(ncell)=9.99
                time(ncell)=0.
                exit
              enddo
            endif
          endif
          if(rsq(i).gt.rsqcut) then
            noisy=noisy+1
            inoisy(noisy)=i
          endif
        endif
      enddo
c      print *,'*badcheck*   ncell =',ncell,'  ndead =',ndead,'  noisy =',noisy
      if(lhist) call hf1(111,float(ndead+noisy),1.)
      do i=1,ndead
        iad=idead(i)-1
        ix=iad/48
        iy=iad-ix*48
        if(ix.ge.iholel.and.ix.le.iholer.and.iy.ge.iholed.and.iy.le.iholeu) cycle ! hole
        if(lhist) call hf2(112,ix+.5,iy+.5,1.)
      enddo
      do i=1,noisy
        iad=inoisy(i)-1
        ix=iad/48
        iy=iad-ix*48
        if(lhist) call hf2(113,ix+.5,iy+.5,1.)
      enddo
      RETURN
      END 
c
      subroutine readbad
      implicit none
#include "parameters.inc"
#include "chkbad.inc"
#include "files.inc"
#include "flags.inc"
#include "gclb.inc"
#include "runpar.inc"
      integer i
C 
      if(.not.lbad) return
      write(badcells(9:14),'(i6.6)') NRUN
      OPEN(UNIT=8,FILE=badcells,STATUS='OLD',ERR=10)
      read(8,*) 
      READ(8,*,err=10) (emn(i),i=1,numc)
      read(8,*) 
      READ(8,*,err=10) (est(i),i=1,numc) 
      read(8,*) 
      READ(8,*,err=10) (rmn(i),i=1,numc) 
      read(8,*) 
      READ(8,*,err=10) (rsq(i),i=1,numc) 
      read(8,*) 
      READ(8,*,err=10) (rst(i),i=1,numc) 
      CLOSE(UNIT=8)
cc      do i=1,numc
cc        if(rst(i).gt.rstcut.and.rmn(i).gt.rmncut.and.rmn(i).lt.8.) then
cc          gcal(i)=gcal(i)/rmn(i)
cc          rst(i)=0.
cc        endif
cc      enddo
      RETURN
   10 print '(a,a)', ' *** readbad ***   File ',badcells,' initialyzed'
      do i=1,numc
        emn(i)=0.
        est(i)=0.
        rmn(i)=0.
        rsq(i)=0.
        rst(i)=0.
      enddo
      RETURN
      END 
C 
      subroutine writbad
      implicit none
#include "parameters.inc"
#include "chkbad.inc"
#include "files.inc"
#include "flags.inc"
#include "runpar.inc"
      integer i
C 
      if(.not.lbad) return
      write(badcells(9:14),'(i6.6)') NRUN
      OPEN(UNIT=8,FILE=badcells,status='UNKNOWN',err=10)
      WRITE(8,'(256(/12f8.2))') (emn(i),i=1,numc)
      WRITE(8,'(256(/12f8.1))') (est(i),i=1,numc)
      WRITE(8,'(256(/12f8.3))') (rmn(i),i=1,numc)
      WRITE(8,'(256(/12f8.3))') (rsq(i),i=1,numc)
      WRITE(8,'(256(/12f8.1))') (rst(i),i=1,numc)
      CLOSE(UNIT=8)
      print '(3a)',' ** writbad **  File ',badcells,' is written'
      RETURN
   10 print '(3a)',' ** writbad ** OPEN failure.  File ',badcells
      STOP
      END 
