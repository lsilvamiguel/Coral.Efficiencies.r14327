      subroutine extrap_sm2(extrapol,miter,mlassc,mlassc3,
     &  lassoc,par_rochi,
     &  x_ref2,par_trk2,erpr_trk2,chi2_tab2,chi2_tab2_ndf,listhit3)
*
********************************************************************
*                 
*----- Extrapolation of tracks found in LAS from RICH1 to SM2
*
********************************************************************

#include "PARAM.h"
#include "PLANES.h"
#include "HITS_OUT.h"
#include "HITS_IN.h"
#include "LIN_ALL.h"
#include "MAGNET.h"

     
      dimension lassoc(2,*),par_rochi(2)
      dimension x_ref2(*)
      dimension par_trk2(4,*)
      dimension erpr_trk2(4,*)
      dimension chi2_tab2(*),chi2_tab2_ndf(*)

      dimension listhit3(ntrack,*)

      dimension coor(mxpl),sig(mxpl)
      dimension xpl23(2*mxpl),cs23(2*mxpl)
      dimension ss23(2*mxpl),hit23(2*mxpl),err23(2*mxpl)
      dimension ans(10)
      dimension distm(mxpl),tol(mxpl),jout(mxpl)
      dimension line(mxpl)
      logical extrapol

*----- count rejects, candidates,failures.....             
        nrej=0
        ncand=0
        nfail=0

*------ Min. distance from SM1 center.............
  
         dist_mg1=1.2

*------ End of first step in extrapolation........
     
          zref1=11.0
          eps=1.e-5

      
      road = par_rochi(1)                   ! Road width
      chilim = par_rochi(2)                 ! Max. chi2/nd

      if (mlassc.le.0) return
*------ loop over reconstructed tracks............

      do 150 jt = 1,mlassc
       if(.not.extrapol)go to 149             ! extrapolation switch off
* 
*---- return arguments (par_trk2,erpr_trk2,chi2_tab2) are preset to -999.......
*

   
       chi2_tab2(miter+jt) = -999.
       chi2_tab2_ndf(miter+jt) = -999.
      
      l_dw = lassoc(2,jt)
      if (l_dw.lt.1.or.l_dw.gt.nlin_fdb) go to 150

*----- extrapolate to ALL extra planes............

      do 10 kpl = 1,ntpl_sup
      dxx = xplan_sup(kpl) - xlinb(l_dw)
      y_pl = ylinb(l_dw) + yplinb(l_dw) * dxx
      z_pl = zlinb(l_dw) + zplinb(l_dw) * dxx
      
        jout(kpl)=0

*----test if inside chamber limits..............

        siy = 0.
        siz = 0.
      if (ityp3(kpl).eq.666) then  
          siy = size_y
          siz = size_z
      elseif (ityp3(kpl).eq.444) then
          siy = size_yg
          siz = size_zg
      elseif (ityp3(kpl).eq.333) then
          siy = size_ysf
          siz = size_zsf 
      endif
 
*      print*,xplan_sup(kpl),xlinb(l_dw)
*      print*,y_pl,ylinb(l_dw),yplinb(l_dw)
*      print*,z_pl,zlinb(l_dw),zplinb(l_dw)
      
      if (abs(y_pl).gt.0.5*siy.or.abs(z_pl).gt.0.5*siz) then
       jout(kpl) = 1
       go to 10
      endif

*------------------------------------------------------
      vy = errlib(1,l_dw) + errlib(3,l_dw) * (dxx**2)
      vz = errlib(2,l_dw) + errlib(4,l_dw) * (dxx**2)
      coor(kpl) = sinn_sup(kpl) * y_pl - cosn_sup(kpl) * z_pl
      sig(kpl) = sqrt(sinn_sup(kpl)**2 * vy + cosn_sup(kpl)**2 * vz)

*     print '("Xplane,coor,sig:",f10.3,2f10.6)',xplan_sup(kpl),
*    1  coor(kpl),sig(kpl)

 10    continue


*--- set max. nr. of standard deviations depending on momentum..........
*             (presently not used) 

        tlf = 1

*------ Fetch hits from line2......................

      nplan_u = 0

      do 155 kpl = 1,ntpl2
*---reject if too close to SM1.....................

      if (xplan2(kpl)-x0_m1.lt.dist_mg1)  go to 155
      nplan_u = nplan_u + 1
      xpl23(nplan_u) = xplan2(kpl)
      cs23(nplan_u)  = cosn2(kpl)
      ss23(nplan_u)  = sinn2(kpl)


      line(nplan_u) = line2(kpl,l_dw)

      if(line2(kpl,l_dw).eq.8888)then
       hit23(nplan_u) = 800.
       err23(nplan_u) = 800.
      else
       hit23(nplan_u) = clusters2(line2(kpl,l_dw))
       err23(nplan_u) = clu_errors2(line2(kpl,l_dw))
      endif

155   continue

*------- Same for additionnal planes...............

*----count planes (for test output)...........................
      n_add = 0
      j_plane = 0

      do 156 kpl = 1,ntpl_sup

*---reject if too far for first step...............
      if (xplan_sup(kpl).gt.zref1)  go to 157
      nplan_u = nplan_u + 1
      xpl23(nplan_u) = xplan_sup(kpl)
      cs23(nplan_u)  = cosn_sup(kpl)
      ss23(nplan_u)  = sinn_sup(kpl)

*----- no search if out of chamber........
      if (jout(kpl).eq.1) then
        hit23(nplan_u) = 900.
        line(nplan_u)=8888
        go to 156
      endif

*--- find hit closest to extrapolation.............
      distm(nplan_u) = 1.0E+10
      do kh = 1,jcshits(kpl)

      yzthts=clusters3(numhts(kpl,kh))

      dist = abs(yzthts - coor(kpl))
      if (dist.lt.distm(nplan_u))  then
         distm(nplan_u) = dist
         hit23(nplan_u) = yzthts
         err23(nplan_u) = clu_errors3(numhts(kpl,kh))
        
         line(nplan_u)=numhts(kpl,kh)

***      tol(nplan_u) = sqrt(err23(nplan_u)**2 + sig(kpl)**2)
*--------use preset road width .............................
         tol(nplan_u) = road
      endif
      enddo

*------do not use hit if out of road....

      if (distm(nplan_u).gt.tlf*tol(nplan_u))then
           hit23(nplan_u) = 800.
           line(nplan_u) = 8888
      endif

      if (distm(nplan_u).le.tlf*tol(nplan_u))  then
        n_add = n_add + 1    
        j_plane = kpl
      endif

156   continue
157   continue
*
*----count and eliminate if out of acceptance....................

      jplok = 0
      do   kcc = 1,ntpl_sup
      if (xplan_sup(kcc).le.zref1) then
       if (jout(kcc).eq.0) jplok = jplok + 1
      endif
      enddo

      if (jplok.lt.2)  then
           nrej = nrej + 1
           go to 149
      endif



***   do kk = 12,nplan_u
***   print '("Nr.-x-cs-ss-h-dist-max:",i4,6f10.5)',
***  1  kk,xpl23(kk),cs23(kk),ss23(kk),hit23(kk),
***  2   distm(kk),tlf*tol(kk)
***   enddo

*-------leave if not more than 1 plane with accepted hit........

           
      if (n_add.le.1) then
*      print '("======> leaving  n_add=",i3)',n_add
       nfail = nfail + 1
       go to 149
      endif

*------Get parameters at end  of fitted range.....

      ncand = ncand + 1
      ans(1) =  xplan_sup(j_plane)

      call stfit1(nplan_u,cs23,ss23,xpl23,hit23,err23,nhit,ans)

*      print '(" Extrap. Fit 1--chi/nd=",e10.4)',ans(10)
*      print '(" Z-X-Y:",3f10.4)',ans(1),ans(2),ans(3)
*      print '(" slopes: ",2f10.6)',ans(4),ans(5)
*      print '(" slp.err:",2f10.6)',sqrt(ans(8)),sqrt(ans(9))

*-----Extrapolate to remaining planes.................

      do 165 kpl = 1,ntpl_sup

      if (xplan_sup(kpl).le.zref1)  go to 165
      dxx = xplan_sup(kpl) - ans(1)
      y_pl = ans(2) + ans(4) * dxx
      z_pl = ans(3) + ans(5) * dxx
      vy = ans(6) + ans(8) * (dxx**2)
      vz = ans(7) + ans(9) * (dxx**2)
      coor(kpl) = sinn_sup(kpl) * y_pl - cosn_sup(kpl) * z_pl
      sig(kpl) = sqrt(sinn_sup(kpl)**2 * vy + cosn_sup(kpl)**2 * vz)

*     print '("Xplane,coor,sig:",f10.3,2f10.6)',xplan_sup(kpl),
*    1  coor(kpl),sig(kpl)

 165   continue

*------ Redo the fit using the same hits in previous planes....

      do 166 kpl = 1,ntpl_sup

      if (xplan_sup(kpl).le.zref1)  go to 166
      nplan_u = nplan_u + 1
      xpl23(nplan_u) = xplan_sup(kpl)
      cs23(nplan_u)  = cosn_sup(kpl)
      ss23(nplan_u)  = sinn_sup(kpl)

*------ no search if out of chamber limits..............

      if (jout(kpl).eq.1)  then
         hit23(nplan_u) = 900.
         line(nplan_u)=8888
         go to 166
      endif

*--- find hit closest to extrapolation.............

      distm(nplan_u) = 1.0E+10
      do kh = 1,jcshits(kpl)

         yzthtss=clusters3(numhts(kpl,kh))

      dist = abs(yzthtss - coor(kpl))
      if (dist.lt.distm(nplan_u))  then
         distm(nplan_u) = dist
         hit23(nplan_u) = yzthtss
         err23(nplan_u) = clu_errors3(numhts(kpl,kh))

         line(nplan_u) = numhts(kpl,kh)

***      tol(nplan_u) = sqrt(err23(nplan_u)**2 + sig(kpl)**2)
*------ use preset road width ...........................
         tol(nplan_u) = road
      endif
      enddo

*------do not use hit if out of road.....................

      if (distm(nplan_u).gt.tlf*tol(nplan_u))then
         hit23(nplan_u) = 800.
         line(nplan_u) = 8888
      endif

      if (distm(nplan_u).le.tlf*tol(nplan_u))  then
        n_add = n_add + 1    
        j_plane = kpl
      endif
166   continue
      ans(1) =  xplan_sup(j_plane)

*     do kk = 12,nplan_u
*     print '("Nr.-x-cs-ss-h-dis-max:",i4,6f10.5)',
*    1  kk,xpl23(kk),cs23(kk),ss23(kk),hit23(kk),
*    2  distm(kk),tlf*tol(kk)
*     enddo

      call stfit1(nplan_u,cs23,ss23,xpl23,hit23,err23,nhit,ans)

*----- apply chi2/nd cut.....try recycling if too large (discard highest contribution)
*
       if (ans(10).gt.chilim) then

***     print '("nhit, chi2/nd initial:",i4,1x,,e10.4)',nhit,ans(10)

        jplmc = 0
        tchi = 0.

        do kk=13,nplan_u
        th = ss23(kk) * (ans(2) +(xpl23(kk)-ans(1))*ans(4))
     1      -cs23(kk) * (ans(3) +(xpl23(kk)-ans(1))*ans(5))

        if (hit23(kk).gt.80.0) then
          chi_pl = 0.
        else
          chi_pl = (hit23(kk)-th)**2/err23(kk)**2
        endif
        if (chi_pl.gt.tchi) then
             tchi = chi_pl
             jplmc = kk
        endif
        enddo

        hit23(jplmc) = 800.
        line(jplmc) = 8888

*     do kk = 12,nplan_u
*     print '("Nr.-x-cs-ss-h-dis-max:",i4,6f10.5)',
*    1  kk,xpl23(kk),cs23(kk),ss23(kk),hit23(kk),
*    2  distm(kk),tlf*tol(kk)
*     enddo

        call stfit1(nplan_u,cs23,ss23,xpl23,hit23,err23,nhit,ans)

**      print '("nhit, chi2/nd new:",i4,1x,e10.4)',nhit,ans(10)

        if (ans(10).gt.chilim) then
          nfail = nfail + 1
          go to 149
        endif
       endif
*---------------------------------------------------------------------------
**     print '(" Extrap. Fit 2--chi/nd=",e10.4)',ans(10)
**     print '(" Z-X-Y:",3f10.4)',ans(1),ans(2),ans(3)
**     print '(" Fin.slopes: ",2f10.6)',ans(4),ans(5)
**     print '(" slp.errors:",2f10.6)',sqrt(ans(8)),sqrt(ans(9))

*------store hits used in this fit in HITS_OUT........
      
*---- now rewritting hits pointers list listhit3 for output
            
            do ka=1,ntpl_sup
              listhit3(miter+jt,ka)=8888
            do ja=1,nplan_u
            if(abs(xpl23(ja)-xplan_sup(ka)).lt.eps)then
              listhit3(miter+jt,ka)=line(ja)

*             print*,jt,ka,ja,listhit3(miter+jt,ka)
            endif
           
            enddo
            enddo


*------update line parameters in LIN_ALL..............

      xlinb(l_dw) = ans(1)
      ylinb(l_dw) = ans(2)
      zlinb(l_dw) = ans(3)
      yplinb(l_dw)= ans(4)
      zplinb(l_dw)= ans(5)
      call ucopy(ans(6),errlib(1,l_dw),5)

*------store last point (in front of SM2) in LIN_ALL......

      xlast(l_dw) = xlinb(l_dw)
      ylast(l_dw) = ylinb(l_dw)
      zlast(l_dw) = zlinb(l_dw)
      ey(l_dw) = sqrt(errlib(1,l_dw))
      ez(l_dw) = sqrt(errlib(2,l_dw))

*     print '(" End point at Z=",f10.4)',xlast(l_dw)
*     print '(" X =",f10.5," +/-",f10.5)',ylast(l_dw),ey(l_dw)
*     print '(" Y =",f10.5," +/-",f10.5)',zlast(l_dw),ez(l_dw)


*-- now second helix

      x_ref2(miter+jt) = xlast(l_dw)              ! xref for second helix

      par_trk2(1,miter+jt) = ylast(l_dw)          ! Y
      par_trk2(2,miter+jt) = zlast(l_dw)          ! Z
      par_trk2(3,miter+jt) = ans(4)               ! slope in Y
      par_trk2(4,miter+jt) = ans(5)               ! slope in Z

      erpr_trk2(1,miter+jt) = ey(l_dw)            ! error in Y
      erpr_trk2(2,miter+jt) = ez(l_dw)            ! error in Z
      erpr_trk2(3,miter+jt) = sqrt(ans(8))        ! error for slope Y
      erpr_trk2(4,miter+jt) = sqrt(ans(9))        ! error for slope Z

      chi2_tab2_ndf(miter+jt) = ans(10)               ! chi2/nd
      chi2_tab2(miter+jt) = ans(10)*(nhit-4)          ! chi2


      mlassc3=mlassc3+1
      go to 150

*----take care of rejects and failures......back to results from tk_fit...

149   continue


      
       listhit3(miter+jt,1)=9999                  ! just flag - extrapolation fault
      do ka=2,ntpl_sup
       listhit3(miter+jt,ka)=8888             
      enddo

150   continue

*       print '(" so far rejects-cand.-failures:",3i5)',nrej,ncand,nfail

*       print *,(jt,jf_sup(miter+jt),jt=1,mlassc)
      end









