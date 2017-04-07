c        Subset of programs for shower separation in ECAL2.  A.Lednev 3.11.2008.
c        Subroutines for FAST shower fit.
c        fastreco, GAMS, c3to5, cellxy, chisq1, chisq2, GAM, mom1, mom2, TWOGAM,
c        ORDER, cell0, CPRINT.
c-----------------------------------------------------------------
      subroutine fastreco(ipcl,leng)
      implicit none
#include "parameters.inc"
      integer ipcl, leng
#include "cells.inc"
#include "cluster.inc"
#include "debug.inc"
#include "flags.inc"
#include "gamma.inc"
      integer i, j, iadc(numc), iwrk(200)
      real adc(numc), wrk(10000)
c      print *,'***fastreco*** ipcl,leng=',ipcl,leng
      do i=1,leng
        j=listcells(ipcl+i-1)
        adc(i)=e(j)
        iadc(i)=iaddr(j)
      enddo
c      call cprint(leng,adc,iadc)
      call gams(leng,adc,iadc,wrk,iwrk)
c      print *,'* fastreco found',ngam,' gamma *'
      return
      end
c-----------------------------------------------------------------------
c        A.Lednev  15.4.1992.        The last correction 02.08.2008.
      SUBROUTINE GAMS(NADC,ADC,IADC,WRK,IWRK)
      implicit none
#include "parameters.inc"
      integer NADC, IADC(numc), IWRK(NADC)
      real ADC(numc), WRK(NADC,0:52)
#include "cell0.inc"
#include "icalad.inc"
#include "isitlg.inc"
#include "xy.inc"
#include "ccuts.inc"
#include "debug.inc"
#include "flags.inc"
#include "gamma.inc"
#include "thresholds.inc"
      INTEGER IGMPK(2,100)
      integer i, ic, icell, idx, ig, in, ipk, IPNPK(100), iter, IX, ixpk, ixy, IXYmax, IXYmin, ixypk, iy, iyc, iypk
      integer LENG, LG, ndof(100), NGAMinput, niter, NPK
      real A, AC, CHI2, CHISQ, CHISQ1, dof, DX, DY, EE, EGM, EPK(100), GM(4,100), pkmin, RATIO, WK
      real Xpeak, XPK(100), XX, XX2, Ypeak, YPK(100), YY, YY2
      DATA CHISQ1/1000./,CHI2/20./,niter/3/,idx/48/
c
      IF(NGAM.GE.199)        RETURN
c      print *,'*GAMS*  NADC =',NADC,'  Input adc,iadc ='
c      print '(10(f7.2,i6))',(adc(i),iadc(i),i=1,nadc)
      call order(NADC,ADC,IADC)
      NGAMinput=NGAM
      pkmin=egcut/4.
cc      call cprint(nadc,adc,iadc)
c
c  Peaks search :
c
      NPK=0
      DO 1 IC=1,NADC
        AC=ADC(IC)
        IF(AC.LT.PKMIN)        goto 1
        IXY=IADC(IC)
        IXYmax=IXY+idx+1
        IXYmin=IXY-idx-1
        call cellxy(ixy,icell,XX,YY)
c      print *
c      print *,'ic,IXY,X,Y =',ic,IXY,XX,YY
c      print '(a,i4,f6.1,3i6)','ic,AC,IAC,IXYmin,IXYmax',ic,AC,IADC(IC),IXYmin,IXYmax
        iyc=ixy-ixy/idx*idx
        in=ic
        do while (in.lt.nadc.and.IADC(in).le.IXYmax)
          in=in+1
          iy=IADC(in)-IADC(in)/idx*idx
c          print '(a,i4,f6.1,3i6)','in,ADC,IADC,iyc,iy=',in,ADC(in),IADC(in),iyc,iy
          IF(iabs(iy-iyc).LE.1.AND.ADC(IN).GE.AC) goto 1
        enddo
        in=ic
        do while (in.gt.1.and.IADC(in).ge.IXYmin)
          in=in-1
          iy=IADC(in)-IADC(in)/idx*idx
c          print '(a,i4,f6.1,3i6)','in,ADC,IADC,iyc,iy=',in,ADC(in),IADC(in),iyc,iy
          IF(iabs(iy-iyc).LE.1.AND.ADC(IN).GT.AC) goto 1
        enddo
        NPK=NPK+1           !  Peak found
        IPNPK(NPK)=IC
c        print '(2(a,i4))','***GAMS*** npk =',npk,'  ipnpk =',ipnpk(npk)
        IF((NPK+2)*NADC.GE.10000) return
        IF(NPK.GT.100) exit
    1 CONTINUE
      IF(NPK.EQ.0) then
        RETURN
      ENDIF
c      print '(a,i4)','*** GAMS *** npk =',npk
      do ipk=1,npk
        ixypk=iadc(ipnpk(ipk))
        ixpk=ixypk/idx
        iypk=ixypk-ixpk*idx
        Xpeak=ixpk
        Ypeak=iypk
        if(lhist) call hf2(301,Xpeak,Ypeak,1.)
c        print *,'ipk,ixypk,Xpeak,Ypeak =',ipk,ixypk,Xpeak,Ypeak
      enddo
      if(nadc.gt.500) then
        print *,'*** GAMS *** nadc,npk =',nadc,npk
        return
      endif
c  Peaks have been found  ***

      IF(NPK.EQ.1) THEN                !   Gamma search for one peak :
        NGAM=NGAM+1
        CHISQ=CHI2
        CALL GAM(NADC,ADC,IADC,CHISQ,dof,GM(3,NGAM),GM(1,NGAM),GM(2,NGAM),
     ,                          GM(3,NGAM+1),GM(1,NGAM+1),GM(2,NGAM+1))
        GM(4,NGAM)=CHISQ
        ndof(ngam)=nint(dof)
        IF(GM(3,NGAM+1).GT.0..AND.NGAM.LE.199) THEN    ! Two gamma found
          NGAM=NGAM+1
          GM(4,NGAM)=CHISQ
          ndof(ngam)=nint(dof)
        ENDIF
c      print'(a,i4)', ' * GAMS *   NGAM =',NGAM

      ELSE                        ! Gamma search for several peaks.
c      print *,' * GAMS *   First step.  One gamma in each peak.'
        RATIO = 1.
        DO I=1,NADC
          WRK(I,NPK+1)=1.
          iwrk(i)=iadc(i)
        ENDDO

        DO iter=1,niter
          DO I=1,NADC
            WRK(I,0)=0.
          ENDDO
c
          DO IPK=1,NPK
            IC=IPNPK(IPK)
            IF(iter.gt.1)  RATIO=WRK(IC,IPK)/WRK(IC,NPK+1)
            EGM=ADC(IC)*RATIO
            ixypk=IADC(IC)
            IXPK=ixypk/idx
            IYPK=ixypk-ixpk*idx
            call cellxy(ixypk,icell,XPK(IPK),YPK(IPK))
            EPK(IPK)=EGM
            XPK(IPK)=EGM*XPK(IPK)
            YPK(IPK)=EGM*YPK(IPK)
            DO IN=IC+1,NADC
              ixy=IADC(IN)
              IX=ixy/idx
              IY=ixy-ix*idx
              IF(ixy-ixypk.GT.idx+1) EXIT
              IF(iabs(IY-IYPK).LE.1) THEN
                IF(iter.gt.1)  RATIO=WRK(IN,IPK)/WRK(IN,NPK+1)
                EGM=ADC(IN)*RATIO
                EPK(IPK)=EPK(IPK)+EGM
                call cellxy(ixy,icell,XX2,YY2)
                XPK(IPK)=XPK(IPK)+EGM*XX2
                YPK(IPK)=YPK(IPK)+EGM*YY2
              ENDIF
            ENDDO
            DO IN=IC-1,1,-1
              ixy=IADC(IN)
              IX=ixy/idx
              IY=ixy-ix*idx
              IF(ixypk-ixy.GT.idx+1) EXIT
              IF(iabs(IY-IYPK).LE.1) THEN
                IF(iter.gt.1)  RATIO=WRK(IN,IPK)/WRK(IN,NPK+1)
                EGM=ADC(IN)*RATIO
                EPK(IPK)=EPK(IPK)+EGM
                call cellxy(ixy,icell,XX2,YY2)
                XPK(IPK)=XPK(IPK)+EGM*XX2
                YPK(IPK)=YPK(IPK)+EGM*YY2
              ENDIF
            ENDDO
            XPK(IPK)=XPK(IPK)/EPK(IPK)
            YPK(IPK)=YPK(IPK)/EPK(IPK)
c            print *,'ipk,XPK(ipk),YPK(ipk) =',ipk,XPK(ipk),YPK(ipk)
            icell=icalad(XPK(IPK)*dcell,YPK(IPK)*dcell)
            LG=IsItLG(icell)
            Xpeak=xy(XPK(IPK)*dcell,LG)/dcell
            Ypeak=xy(YPK(IPK)*dcell,LG)/dcell
c            print *,'LG,Xpeak,Ypeak =',LG,Xpeak,Ypeak
            DO I=1,NADC
              ixy=IADC(I)
              call cellxy(ixy,icell,XX2,YY2)
c              DX=abs(XX2-XPK(IPK))
c              DY=abs(YY2-YPK(IPK))
              DX=abs(XX2-Xpeak)
              DY=abs(YY2-Ypeak)
              A=0.
c            IF(DX.LE.2.5.AND.DY.LE.2.5) A=EPK(IPK)*cell(DX*dcell,DY*dcell,dcell)
              LG=IsItLG(icell)
              IF(DX.LE.2.5.AND.DY.LE.2.5) A=EPK(IPK)*cell0(DX*dcell,DY*dcell,LG)
              WRK(I,IPK)=A
              WRK(I,0)=WRK(I,0)+WRK(I,IPK)
            ENDDO
          ENDDO
          DO I=1,NADC
            WK=WRK(I,0)
            if(WK.le.thrcommon) WK=thrcommon
            WRK(I,NPK+1)=WK
          ENDDO
c          print *,'*GAMS* iter =',iter,'  wrk(i,0),iwrk ='
c          print '(10(f7.2,i6))',(wrk(i,0),iwrk(i),i=1,nadc)
        ENDDO
c      print *,'*GAMS* wrk(i,0),iwrk ='
c      print '(10(f7.2,i6))',(wrk(i,0),iwrk(i),i=1,nadc)
c
c      print *,'*GAMS*  npk,nadc=',npk,nadc
        DO IPK=1,NPK
          LENG=0
          DO I=1,NADC
            IF(WRK(I,0).GT.0.) THEN
              EE=ADC(I)*(WRK(I,IPK))/WRK(I,0)
              IF(EE.GT.thrcommon) THEN
                LENG=LENG+1
                WRK(LENG,NPK+1)=EE
                IWRK(LENG)=IADC(I)
              ENDIF
            ENDIF
          ENDDO
c          print *,'*GAMS*  Peak',ipk
c          call cprint(leng,wrk(1,npk+1),iwrk)
c        print *,'* GAMS *  IPK =',ipk,'  LENG =',LENG,'   WRK,IWRK ='
c        print '(10(f7.2,i6))',(wrk(i,npk+1),iwrk(i),i=1,leng)
          IF(NGAM.GE.199)        GO TO 500
          IF(LENG.EQ.0)        CYCLE
          NGAM=NGAM+1
          CHISQ=CHISQ1
          CALL GAM(LENG,WRK(1,NPK+1),IWRK(1),CHISQ,dof,GM(3,NGAM),GM(1,NGAM),
     ,      GM(2,NGAM),GM(3,NGAM+1),GM(1,NGAM+1),GM(2,NGAM+1))
          GM(4,NGAM)=CHISQ
          ndof(ngam)=nint(dof)
          IGMPK(1,IPK)=NGAM
          IGMPK(2,IPK)=NGAM
          IF(GM(3,NGAM+1).GT.egcut.AND.NGAM.LE.199) THEN
            NGAM=NGAM+1
            GM(4,NGAM)=CHISQ
            ndof(ngam)=nint(dof)
            IGMPK(2,IPK)=NGAM
          ENDIF
c          print '(a,2i4,a,f8.2)', '*GAMS* 1st step. IPK,NGAM =',IPK,NGAM,'  chisq =',chisq
        ENDDO
c
c   This is the second step :   ( 1 or 2 gamma in each peak )
c
c      print *,' * GAMS *   Second step'
        DO I=1,NADC
          WRK(I,0)=0.
        ENDDO
        DO IPK=1,NPK
          if(igmpk(1,ipk).lt.1) cycle
          DO I=1,NADC
            WRK(I,IPK)=0.
            ixy=IADC(I)
            call cellxy(ixy,icell,XX2,YY2)
            LG=IsItLG(icell)
            DO IG=IGMPK(1,IPK),IGMPK(2,IPK)
              DX=XX2-GM(1,IG)
              DY=YY2-GM(2,IG)
c              IIA=nint(GM(3,IG)*cell(DX*dcell,DY*dcell,dcell))
              A=GM(3,IG)*cell0(DX*dcell,DY*dcell,LG)
              WRK(I,IPK)=WRK(I,IPK)+A
              WRK(I,0)=WRK(I,0)+A
            ENDDO
          ENDDO
        ENDDO
cc      do ipk=1,npk
cc        print *,'*GAMS* ipk=',ipk,' wrk =',(wrk(i,ipk),i=1,nadc)
cc      enddo
c
        NGAM=NGAMinput
        DO IPK=1,NPK
          LENG=0
          DO I=1,NADC
            IF(WRK(I,0).GT.0.) THEN
              EE=ADC(I)*WRK(I,IPK)/WRK(I,0)
              IF(EE.GT.thrcommon) THEN
                LENG=LENG+1
                WRK(LENG,NPK+1)=EE
                IWRK(LENG)=IADC(i)
              ENDIF
            ENDIF
          ENDDO
          IF(NGAM.GE.199)exit
          IF(LENG.EQ.0)        cycle
          NGAM=NGAM+1
          CHISQ=CHI2
          CALL GAM(LENG,WRK(1,NPK+1),IWRK(1),CHISQ,dof,GM(3,NGAM),GM(1,NGAM),
     ,          GM(2,NGAM),GM(3,NGAM+1),GM(1,NGAM+1),GM(2,NGAM+1))
          GM(4,NGAM)=CHISQ
          ndof(ngam)=nint(dof)
          IF(GM(3,NGAM+1).GT.0..AND.NGAM.LE.199) THEN
            NGAM=NGAM+1
            GM(4,NGAM)=CHISQ
            ndof(ngam)=nint(dof)
          ENDIF
cc        print '(a,2i4)', ' * GAMS * 2step *  IPK,NGAM =',IPK,NGAM
        ENDDO
      ENDIF
c
c  Conversion of gamma coordinates (cells to cm ).
c
  500 CONTINUE
      DO IG=NGAMinput+1,NGAM
        xgam(ig)=GM(1,ig)*dcell
        ygam(ig)=GM(2,ig)*dcell
        egam(ig)=GM(3,ig)
        chigam(ig)=GM(4,ig)
        ndf(ig)=ndof(ig)
      ENDDO
cc      if(numbev.eq.303) then
cc        print *,'*** GAMS *** ngam =',ngam,'   x,y,e,chisq ='
cc        print '(4f10.2)',(xgam(i),ygam(i),egam(i),chigam(i),i=ngaminput+1,ngam)
cc      endif
      RETURN
      END
c-----------------------------------------------------------------------
c  Convert three-dimensional variable (eps,dx,dy) to six-dimensional variable
c  Input values:  e0,x0,y0 (constants) and eps,dx,dy (variable)
c  Output values: e1,x1,y1,e2,x2,y2 - two gamma energies and coordinates.
      subroutine c3to5(e0,x0,y0,eps,dx,dy,e1,x1,y1,e2,x2,y2)
      implicit none
      real e0, x0, y0, eps, dx, dy, e1, x1, y1, e2, x2, y2
c     eps=(e1-e2)/e0, (e0=e1+e2),  x0*e0=x1*e1+x2*e2,  dx=x1-x2.
      e1=e0*(1+eps)/2.
      e2=e0-e1
      x1=x0+dx*(1.-eps)/2.
      y1=y0+dy*(1.-eps)/2.
      x2=x0-dx*(1.+eps)/2.
      y2=y0-dy*(1.+eps)/2.
      return
      end
c-----------------------------------------------------------------------
c  Calculates the cell coordinates (X,Y) and cell index in calib. file (icell)
c  Input value:  ixy - cell address
c  Output values: icell, X, Y.
      subroutine cellxy(ixy,icell,X,Y)
      implicit none
      integer ixy, icell
      real X, Y
#include "icalad.inc"
#include "parameters.inc"
      integer idx, ix, iy
      data idx/48/
c
      ix = ixy/idx
      iy = ixy-ix*idx
      X = (ix-31.5)
      Y = (iy-23.5)
      icell=icalad(X*dcell,Y*dcell)
      return
      end
c-----------------------------------------------------------------------
c  Chi squared calculation subroutine
c  Input values: nadc, adc(numc),iadc(numc),e1,x1,y1
c  Output value: chisq
      subroutine chisq1(nadc,adc,iadc,e1,x1,y1,chisq)
      implicit none
#include "parameters.inc"
      integer nadc, iadc(numc)
      real adc(numc), e1, x1, y1, chisq
#include "cell0.inc"
#include "errcell2.inc"
#include "isitlg.inc"
      integer I, idx, icell, ixy, LG
      real A, CONST, D, xx2, yy2
      data CONST/1.5/,idx/48/
      CHISQ=0.
      DO I=1,NADC
        ixy=iadc(i)
        call cellxy(ixy,icell,xx2,yy2)
        LG=IsItLG(icell)
        A  = E1*cell0((x1-xx2)*dcell,(y1-yy2)*dcell,LG)
cc        D  = CONST * A * (1.-A/E1)
        D  = errcell2(A) * (1.-A/E1)
        IF(D.LT.0.) D = 0.
        D  = D + .0001
        CHISQ=CHISQ+(ADC(I)-A)**2/D
c        print *,'i,adc(i),a,d,chisq =',i,adc(i),a,d,chisq
      ENDDO
      return
      end
c-----------------------------------------------------------------------
c  Chi squared calculation subroutine
c  Input values: nadc, adc(numc),iadc(numc),e1,x1,y1
c  Output value: chisq
      subroutine chisq2(nadc,adc,iadc,e1,x1,y1,e2,x2,y2,chisq)
      implicit none
#include "parameters.inc"
      integer nadc, iadc(numc)
      real adc(numc), e1, x1, y1, e2, x2, y2, chisq
#include "cell0.inc"
#include "errcell2.inc"
#include "isitlg.inc"
      integer I, icell, idx, ixy, LG
      real A, CONST, D, xx2, yy2
      data CONST/1.5/,idx/48/
      CHISQ=0.
      DO I=1,NADC
        ixy=iadc(i)
        call cellxy(ixy,icell,xx2,yy2)
        LG=IsItLG(icell)
        A  = E1*cell0((x1-xx2)*dcell,(y1-yy2)*dcell,LG)
     +     + E2*cell0((x2-xx2)*dcell,(y2-yy2)*dcell,LG)
cc        D  = CONST * A * (1.-A/E1)
        D  = errcell2(A) * (1.-A/(E1+E2))
        IF(D.LT.0.) D = 0.
        D  = D + .0001
        CHISQ=CHISQ+(ADC(I)-A)**2/D
c        print *,'i,adc(i),a,d,chisq =',i,adc(i),a,d,chisq
      ENDDO
      return
      end
c-----------------------------------------------------------------------
c  Correction of coordinates (X,Y) due to PWO-LG boundary
c  Input values: (x0,y0) - center-of-gravity coordinates
c                 (X,Y) -  coordinates to be corrected
c  Output values: X and Y (corrected)
      subroutine corrxy0(x0,y0,X,Y)        ! This subroutine is for GAMS-4pi.
      implicit none
      real x0, y0, X, Y
      logical lpwo,l0pwo
      real ratio
      data ratio/1.4/
c
      if(abs(x0-24.).gt.2..or.abs(y0-24.).gt.2.) then
        l0pwo=.false.
      else
        l0pwo=.true.
      endif
      if(abs(X-24.).gt.2..or.abs(Y-24.).gt.2.) then
        lpwo=.false.
      else
        lpwo=.true.
      endif
      if(l0pwo.and..not.lpwo) then        ! (x0,y0) is in PWO, (X,Y) is out.
        if(X.gt.26.) X=(X-26.)/ratio+26.
        if(X.lt.22.) X=(X-22.)/ratio+22.
        if(Y.gt.26.) Y=(Y-26.)/ratio+26.
        if(Y.lt.22.) Y=(Y-22.)/ratio+22.
        return
      endif
      if(.not.l0pwo.and.lpwo) then        ! (X,Y) is in PWO, (x0,y0) is out.
        if(x0.gt.26.) X=(X-26.)*ratio+26.
        if(x0.lt.22.) X=(X-22.)*ratio+22.
        if(y0.gt.26.) Y=(Y-26.)*ratio+26.
        if(y0.lt.22.) Y=(Y-22.)*ratio+22.
      endif
      return
      end
c----------------------------------------------------------------------
C     1 dec. 1992.     ONE GAMMA FIT PROGRAM.   Modified on 21 aug. 2000.
      SUBROUTINE GAM(NADC,ADC,IADC,CHISQ,dof,E1,X1,Y1,E2,X2,Y2)
      implicit none
#include "parameters.inc"
      integer NADC, IADC(numc)
      real ADC(numc), CHISQ, dof, E1, X1, Y1, E2, X2, Y2
#include "icalad.inc"
#include "isitlg.inc"
#include "xy.inc"
      logical lhist
#include "fcuts.inc"
#include "debug.inc"
      integer icell, it, LG
      real chi0, chi00, chid, chil, CHImem, chir, chisq0, chiu, dis, dof1, dxy, EE, parx, pary, stepmin, stepx, stepy, x0, XX, y0, YY
      data dxy/.03/,stepmin/.0003/
C
c      print *,'** GAM**  NADC =',nadc,'  ADC,IADC ='
c      print '(10(f6.1,i6))',(adc(i),iadc(i),i=1,nadc)
      E2=0.
      X2=0.
      Y2=0.
      CHImem=CHISQ
      chisq=99999.
      call mom1(nadc,adc,iadc,e1,x1,y1)
      IF(NADC.LE.0)        RETURN
      icell=icalad(x1*dcell,y1*dcell)
      LG=IsItLG(icell)
      x1=xy(x1*dcell,LG)/dcell
      y1=xy(y1*dcell,LG)/dcell
      call chisq1(nadc,adc,iadc,e1,x1,y1,chi0)
      chisq0=chi0
      DOF=NADC-2
      IF(DOF.LT.1.) DOF=1.
      CHISQ=chi0
c      print '(a,f8.0,2f8.3,2f8.1)', ' * gam *  e1,x1,y1,chimem,chisq =',
c     ,  e1,x1,y1,chimem,chisq
c      if(chisq.le.chimem/2.) return
      x0=x1
      y0=y1
      it=0
    1 continue
      it=it+1
      call chisq1(nadc,adc,iadc,e1,x0+dxy,y0,chir)
      call chisq1(nadc,adc,iadc,e1,x0-dxy,y0,chil)
      call chisq1(nadc,adc,iadc,e1,x0,y0+dxy,chiu)
      call chisq1(nadc,adc,iadc,e1,x0,y0-dxy,chid)
      if(chi0.gt.chir.or.chi0.gt.chil) then
        stepx=dxy
        if(chir.gt.chil) stepx=-stepx
      else
        stepx=0.
        parx=chir+chil-2.*chi0
        if(parx.gt.0.) stepx=-dxy*(chir-chil)/(2.*parx)
      endif
      if(chi0.gt.chiu.or.chi0.gt.chid) then
        stepy=dxy
        if(chiu.gt.chid) stepy=-stepy
      else
        stepy=0.
        pary=chiu+chid-2.*chi0
        if(pary.gt.0.) stepy=-dxy*(chiu-chid)/(2.*pary)
      endif
c      print '(a,i4,2f8.4)', ' * gam *  it,stepx,stepy =',it,stepx,stepy
      if(abs(stepx).lt.stepmin.and.abs(stepy).lt.stepmin) go to 2
      call chisq1(nadc,adc,iadc,e1,x0+stepx,y0+stepy,chi00)
c      print '(a,i4,f8.0,2f8.3,2f8.1)', ' * gam *  it,e1,x0,y0,chi0,chi00 =',
c     ,  it,e1,x0,y0,chi0,chi00
      if(chi00.ge.chi0) go to 2
      chi0=chi00
      x0=x0+stepx
      y0=y0+stepy
      go to 1
    2 continue
      if(chi0.lt.chisq0) then
        x1=x0
        y1=y0
        chisq=chi0
      endif
c      print '(a,f8.0,2f8.3,2f8.1)',' * gam *  e1,x1,y1,chimem,chisq =',
c     ,  e1,x1,y1,chimem,chisq
      if(chisq.le.chimem) return
      EE=E1
      XX=X1
      YY=Y1
      dof1=dof
      CALL TWOGAM(NADC,ADC,IADC,CHISQ,dof,EE,XX,YY,E2,X2,Y2)
c      print '(a,2(f8.0,2f8.2),f8.2)',' * after twogam *  2(e,x,y),chisq =',
c     ,  ee,xx*dcell,yy*dcell,e2,x2*dcell,y2*dcell,chisq
      IF(E2.GT.0.) then
        dis=2.*sqrt(EE*E2*((xx-x2)**2+(yy-y2)**2))/(EE+E2)*dcell
        if(lhist)call hf1(17,dis,1.)
c        print *,'dis,discut =',dis,discut
        if(dis.gt.discut) then
          E1=EE
          X1=XX
          Y1=YY
        else
          E2=0.
          dof=dof1
         endif
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
c     The first momenta calculation.
c  Input values: nadc,adc(numc),iadc(numc)
c  Output values: a0, x0, y0 - Sum of energies and center-of-gravity coordinates
      subroutine mom1(nadc,adc,iadc,a0,x0,y0)
      implicit none
#include "parameters.inc"
      integer nadc, iadc(numc)
      real adc(numc), a0, x0, y0
      integer i, icell, idx, ixy
      real a, a0lg, a0pwo, ratio, x0lg, x0pwo, xx2, y0lg, y0pwo, yy2
      data idx/48/,ratio/1.4/
c
c      print *,'*mom1* adc ='
c      print '(10f8.3)',(adc(i),i=1,nadc)
      a0=0.
      x0=0.
      y0=0.
      if(nadc.le.0) return
      a0lg=0.
      x0lg=0.
      y0lg=0.
      a0pwo=0.
      x0pwo=0.
      y0pwo=0.
      do i=1,nadc
        a=adc(i)
        ixy=iadc(i)
        call cellxy(ixy,icell,xx2,yy2)
        if(icell.lt.numc) then
          a0lg=a0lg+a
          x0lg=x0lg+a*xx2
          y0lg=y0lg+a*yy2
        endif
      enddo
      a0=a0lg
      if(a0.le.0.) return
      if(a0lg.gt.0.) then
        x0lg=x0lg/a0lg
        y0lg=y0lg/a0lg
      endif
c      if(a0lg.gt.0..and.a0pwo.gt.0.) then        ! Correction near LG-PWO boundary
c        if(x0lg.lt.22.) then
c          if(a0lg.lt.a0pwo) then
c            x0lg=(x0lg-22.)/ratio+22.
c          else
c            x0pwo=(x0pwo-22.)*ratio+22.
c          endif
c        elseif(x0lg.gt.26.) then
c          if(a0lg.lt.a0pwo) then
c            x0lg=(x0lg-26.)/ratio+26.
c          else
c            x0pwo=(x0pwo-26.)*ratio+26.
c          endif
c        endif
c        if(y0lg.lt.22.) then
c          if(a0lg.lt.a0pwo) then
c            y0lg=(y0lg-22.)/ratio+22.
c          else
c            y0pwo=(y0pwo-22.)*ratio+22.
c          endif
c        elseif(y0lg.gt.26.) then
c          if(a0lg.lt.a0pwo) then
c            y0lg=(y0lg-26.)/ratio+26.
c          else
c            y0pwo=(y0pwo-26.)*ratio+26.
c          endif
c        endif
c      endif
      a0=a0lg+a0pwo
      x0=(a0lg*x0lg+a0pwo*x0pwo)/a0
      y0=(a0lg*y0lg+a0pwo*y0pwo)/a0
      return
      end
c-----------------------------------------------------------------------
c     The first and second momenta calculation.
c  Input values:  nadc,adc(numc),iadc(numc)
c  Output values: a0,x0,y0 - Sum of energies and center-of-gravity coordinates
c                 xx,yy,yx - Second momenta
      subroutine mom2(nadc,adc,iadc,a0,x0,y0,xx,yy,yx)
      implicit none
#include "parameters.inc"
      integer nadc, iadc(numc)
      real adc(numc), a0, x0, y0, xx, yy, yx
      integer i, icell, idx, ixy
      real a, xx2, yy2
      data idx/48/
c
c      print *,'*mom2* adc ='
c      print '(10f8.3)',(adc(i),i=1,nadc)
      a0=0.
      x0=0.
      y0=0.
      xx=0.
      yy=0.
      yx=0.
      if(nadc.le.0) return
      call mom1(nadc,adc,iadc,a0,x0,y0)
      do i=1,nadc                !  second momenta
        a=adc(i)/a0
        ixy=iadc(i)
        call cellxy(ixy,icell,xx2,yy2)
c        call corrxy0(x0,y0,xx2,yy2)
        xx=xx+a*(xx2-x0)**2
        yy=yy+a*(yy2-y0)**2
        yx=yx+a*(xx2-x0)*(yy2-y0)
      enddo
      return
      end
c-----------------------------------------------------------------------
C    15 apr.1992.     TWO GAMMA FIT PROGRAM.   Modified on 21 aug. 2000.
      SUBROUTINE TWOGAM(NADC,ADC,IADC,CHISQ,dof,E1,X1,Y1,E2,X2,Y2)
      implicit none
#include "parameters.inc"
      integer NADC, IADC(numc)
      real ADC(numc), CHISQ, dof, E1, X1, Y1, E2, X2, Y2
#include "cell0.inc"
#include "errcell2.inc"
#include "isitlg.inc"
#include "ccuts.inc"
#include "runpar.inc"
      integer i, icell, idx, ITER, ixy, LG
      real A0, A1, A2, cell1, cell2, CHI2, chisqc, CONST, cosi
      real D, dchi, dchi0, dchida, DD, DELCH, DX0, DX1, DX2, DXC, DXY, DY0, DY1, DY2, DYC
      real E0, E1C, E2C, ex, eps, EPS0, eps1, eps2, epsc, epsmax
      real gr, GRC, GRE, GREC, GRX, GRXC, GRY, GRYC, GX1, GX2, GY1, GY2
      real r, rsq, rsq2, scal, step, STPMIN, u, x0, x1c, x2c, xx, xx2, Y2C, y0, y1c, yx, yy, yy2
      DATA idx/48/,STPMIN/.005/,CONST/0.20/,DELCH/6./,epsmax/.9/
C
      call mom2(nadc,adc,iadc,e0,x0,y0,xx,yy,yx)
      E2=0.
      X2=0.
      Y2=0.
      IF(NADC.LE.2)        RETURN
c     Choosing of the starting point
      dxy = XX-YY
      rsq2= dxy**2+4.*YX**2
      if(rsq2.lt.1.e-20) rsq2=1.e-20
      rsq = sqrt(rsq2)
      DXC =-sqrt((rsq+dxy)*2.)
      DYC = sqrt((rsq-dxy)*2.)
      if(YX.ge.0.) DYC=-DYC
      r=sqrt(DXC**2+DYC**2)
      epsc=0.
      DO I=1,NADC
        ixy=iadc(i)
        call cellxy(ixy,icell,xx2,yy2)
c        call corrxy0(x0,y0,xx2,yy2)
        u  = (xx2-x0)*dxc/r+(yy2-y0)*dyc/r
        epsc=epsc-adc(i)*u*abs(u)
      ENDDO
      epsc=epsc/(e0*rsq)
      if(epsc.gt.epsmax) epsc=epsmax
      if(epsc.lt.-epsmax) epsc=-epsmax
      dxc=dxc/sqrt(1.-epsc**2)
      dyc=dyc/sqrt(1.-epsc**2)
cc      if(numbev.eq.303) then
cc        print *,'*twogam* start E0,X0,Y0,EPSC,DXC,DYC=',E0,X0,Y0,EPSC,DXC,DYC
cc      endif
c     Start of iterations
      step = .1
      COSI= 0.
      CHI2=1.E35
      DO ITER=1,100
        CALL C3TO5(E0,X0,Y0,EPSC,DXC,DYC,E1C,X1C,Y1C,E2C,X2C,Y2C)
c        print '(a,i4,2(f7.0,2f8.3))', ' * iter, e1,x1,y1, e2,x2,y2 =',
c     ,  iter,e1c,x1c+1.,y1c+1.,e2c,x2c+1.,y2c+1.
        eps1 = (1.+epsc)/2.
        eps2 = (1.-epsc)/2.
        call chisq2(nadc,adc,iadc,e1c,x1c,y1c,e2c,x2c,y2c,chisqc)
c        DO I=1,NADC        !  Chi**2 calculation
c          EX  = ADC(I)
c          ixy = iadc(i)
c          call cellxy(ixy,icell,xx2,yy2)
cc          call corrxy0(x0,y0,xx2,yy2)
c          DX1 = X1C - xx2
c          DY1 = Y1C - yy2
c          DX2 = X2C - xx2
c          DY2 = Y2C - yy2
c          LG=IsItLG(icell)
c          A0 = E1C*cell0(DX1*dcell,DY1*dcell,LG) + E2C*cell0(DX2*dcell,DY2*dcell,LG)
cc          D  = CONST * A0 * (1.-A0/E0)
c          D  = errcell2(A0) * (1.-A0/E0)
c          IF(D.LT.0.) D = 0.
c          D  = D + .0001
c          CHISQC = CHISQC+(A0-EX)**2/D
c        ENDDO
cc      if(numbev.eq.303) then
cc        print *,'*twogam* numbev=',numbev,'  iter',iter,'  chisqc =',chisqc
cc        print *,' E1C,E2C,X1C,X2C,Y1C,Y2C=',E1C,E2C,X1C,X2C,Y1C,Y2C
cc      endif
        IF(CHISQC.GE.CHI2) THEN  !  New step if no improvement
c          print '(a,3f12.3,f11.3)', ' * twogam * gr(E,X,Y)c, chisq ='
c     ,    ,grec,grxc,gryc,CHISQC
          if(iter.gt.1) then
            dchi = CHISQC-CHI2
            dchi0= gr*step
            step = .5*step/sqrt(1.+dchi/dchi0)
          endif
          step = .5*step
        ELSE               !  Calculation of gradient
          GREC = 0.
          GRXC = 0.
          GRYC = 0.
          DO I=1,NADC
            EX  = ADC(I)
            ixy = iadc(i)
            call cellxy(ixy,icell,xx2,yy2)
            DX1 = X1C - xx2
            DY1 = Y1C - yy2
            DX2 = X2C - xx2
            DY2 = Y2C - yy2
            LG=IsItLG(icell)
            cell1=cell0(DX1*dcell,DY1*dcell,LG)
            cell2=cell0(DX2*dcell,DY2*dcell,LG)
            A1  = E1C*cell1
            A2  = E2C*cell2
            A0 = A1 + A2
cc            D  = CONST * A0 * (1.-A0/E0)
            D  = errcell2(A0) * (1.-A0/E0)
            IF(D.LT.0.) D = 0.
            D  = D + .0001
            DD = (A0 - EX)/D
c            dchida = DD
            dchida = DD*(2.-DD*CONST*(1.-2.*A0/E0))
            GX1 = (E1C*cell0((X1C+.05-xx2)*dcell,DY1*dcell,LG)-A1)*20.
            GX2 = (E2C*cell0((X2C+.05-xx2)*dcell,DY2*dcell,LG)-A2)*20.
            GY1 = (E1C*cell0(DX1*dcell,(Y1C+.05-yy2)*dcell,LG)-A1)*20.
            GY2 = (E2C*cell0(DX2*dcell,(Y2C+.05-yy2)*dcell,LG)-A2)*20.
            GREC = GREC + dchida*((A1/E1C-A2/E2C)*E0
     -           - (GX1+GX2)*DXC-(GY1+GY2)*DYC)/2.
            GRXC = GRXC + dchida*(GX1*eps2-GX2*eps1)
            GRYC = GRYC + dchida*(GY1*eps2-GY2*eps1)
          ENDDO
          GRC  = SQRT(GREC**2+GRXC**2+GRYC**2)
          IF(GRC.LT.1.E-10) GRC=1.E-10
          if(iter.gt.1) then
            cosi = (GRE*GREC+GRX*GRXC+GRY*GRYC)/(GR*GRC)
            scal = abs(gr/grc-cosi)
            if(scal.lt..1) scal=.1
            step=step/scal
          endif
          CHI2=CHISQC
          EPS0= EPSC
          DX0 = DXC
          DY0 = DYC
          GRE = GREC
          GRX = GRXC
          GRY = GRYC
          GR  = GRC
c          print '(a,3f12.3,f11.3,a)',' * twogam * gr(E,X,Y)c, chisq ='
c     ,    ,grec,grxc,gryc,CHISQC,' *'
        ENDIF
        EPSC= EPS0- step*GRE/gr
        do while (abs(epsc).ge.epsmax+.01)
          step=step/2.
          EPSC= EPS0- step*GRE/gr
        enddo
        DXC = DX0 - step*GRX/gr
        DYC = DY0 - step*GRY/gr
c        print *,' * twogam *  iter',iter,'  step,CHI2 =',step,CHI2
        IF(step*gr.LT.STPMIN) GO TO 101
      ENDDO
  101 CONTINUE
      IF(CHISQ-CHI2.LT.DELCH.or.E0*(1-abs(eps))/2..lt.egcut) RETURN
      DOF=NADC-3
      IF(DOF.LT.1.) DOF=1.
      CHISQ=CHI2
      CALL C3TO5(E0,X0,Y0,EPS0,DX0,DY0,E1,X1,Y1,E2,X2,Y2)
c      print '(a,2(f7.0,2f8.3),f7.2)', ' * twogam * 2(e,x,y),chisq =',
c     , e1,x1,y1,e2,x2,y2,chisq
      RETURN
      END
c-----------------------------------------------------------------------
c        A.Lednev        ! to order the GAMS information  (data,address)
      SUBROUTINE ORDER(NADC,ADC,iaddr)
      implicit none
#include "parameters.inc"
      integer NADC, iaddr(numc)
      real ADC(numc)
      logical ALARM
      integer I, iad, ifail, K
      real DAT
      ALARM=.FALSE.
      GO TO 1
      ENTRY CHORD(NADC,ADC,iaddr,ifail)        ! only to check the right order
      ALARM=.TRUE.
      ifail=0
    1 IF(NADC.LT.2)  RETURN
      DO K=2,NADC
      IF(iaddr(K).LE.iaddr(K-1)) THEN
        IF(ALARM) then
          ifail=1
          return
        ENDIF
        DAT=ADC(K)
        iad=iaddr(K)
        DO I=K-1,0,-1
          IF(I.GE.1) THEN
            IF(iad.LT.iaddr(I)) THEN
              ADC(I+1)=ADC(I)
              iaddr(I+1)=iaddr(I)
            ELSE
              ADC(I)=DAT
              iaddr(I)=iad
              EXIT
            ENDIF
          ELSE
            ADC(1)=DAT
            iaddr(1)=iad
          ENDIF
        ENDDO
      ENDIF
      ENDDO
      RETURN
      END
c-----------------------------------------------------------------------
      real function cell0(dx,dy,LG)
      implicit none
      integer LG
      real dx, dy
#include "cell.inc"
      cell0=cell(dx,dy,LG)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CPRINT(LENG,BUF,IBUF)
      implicit none
      integer LENG, IBUF(*)
      real BUF(*)
      integer I, IDX, IX, IXY, IY, J, J1, K, MAXX, MAXY, MINX, MINY, MINYJ1, NX, NY
      real ALINE(15), AMAX, thr, XLINE(15)
      data thr/.2/
      DATA IDX/48/  ! IDX=64 for (0:4095) GAMS addr.  IDX=256 for byte struct.
C
      IXY=IBUF(1)
      IX=IXY/IDX
      IY=IXY-IX*IDX
      MINX=IX
      MAXX=IX
      MINY=IY
      MAXY=IY
      AMAX=BUF(1)
      IF(LENG.GT.1) THEN
        DO I=2,LENG
          IXY=IBUF(I)
          IX=IXY/IDX
          IY=IXY-IX*IDX
          IF(IX.LT.MINX) MINX=IX
          IF(IX.GT.MAXX) MAXX=IX
          IF(IY.LT.MINY) MINY=IY
          IF(IY.GT.MAXY) MAXY=IY
          IF(BUF(I).GT.AMAX) AMAX=BUF(I)
        ENDDO
      ENDIF
      IF(AMAX.LT.THR) RETURN
      NX=MAXX-MINX+1
      IF(NX.GT.15) NX=15
      DO I=1,NX
        XLINE(I)=MINX+I-1.
      ENDDO
      PRINT '(A5,15F5.0)', '  Y X',(XLINE(I),I=1,NX)
      NY=MAXY-MINY+1
      DO J1=1,NY
        J=NY+1-J1
        DO I=1,NX
          ALINE(I)=0.
        ENDDO
        DO I=1,LENG
          IXY=IBUF(I)
          IX=IXY/IDX
          IY=IXY-IX*IDX
          IF(IY.EQ.MINY+J-1) THEN
            K=IX-MINX+1
            IF(K.LE.NX) ALINE(K)=BUF(I)
          ENDIF
        ENDDO
        MINYJ1=MINY+J-1
        PRINT '(1X,I3,A1,15f5.1)', MINYJ1,'.',(ALINE(I),I=1,NX)
      ENDDO
      RETURN
      END
