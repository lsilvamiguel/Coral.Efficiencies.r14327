      subroutine remdup(miter,mlassc,lassoc,
     1 chi2_tab,chi2_tab_ndf,x_ref,par_trk,erpr_trk,
     2 x_ref2,par_trk2,erpr_trk2,traject)

***********************************************************
*                                                         *
*   routine called from assoc to remove tracks using a    *
*   duplicated segment before or behind SM1               *
*   (selection is based on best chi2 from spline fit)     *
*                                                         *
***********************************************************

#include "PARAM.h"


      dimension lassoc(2,*)
      dimension chi2_tab(*),chi2_tab_ndf(*)
      dimension par_trk(5,*),erpr_trk(5,*)
      dimension par_trk2(4,*),erpr_trk2(4,*) 
      dimension x_ref(*),x_ref2(*)
      dimension traject(mxplan,*)

*
*  local variables..................
*
      dimension lvass(2,ntrack)
      dimension vchi2(ntrack),vchi22(ntrack)
      dimension vv_trk(5,ntrack),ev_trk(5,ntrack)
      dimension vv_trk2(4,ntrack),ev_trk2(4,ntrack) 
      dimension vxref(ntrack),vxref2(ntrack)
      dimension vtrject(mxplan,ntrack)

      logical ldup
*
*--- quit if there is only one track found.....
*



      if (mlassc.le.1)return


      jwr = 0
10    continue
      jcmin = 0
      chi2_min = 1.0e+08


      do 11 k = 1,mlassc
      if (chi2_tab_ndf(miter+k).lt.chi2_min) then
*
*---check first this is an acceptable candidate---
*
      ldup = .false.
      if (jwr.gt.0) then
        do kd = 1,jwr
        if (lassoc(1,k).eq.lvass(1,kd))  ldup = .true.
        if (lassoc(2,k).eq.lvass(2,kd))  ldup = .true.
        enddo
      endif
      if (ldup) go to 11

        jcmin = k
        chi2_min = chi2_tab_ndf(miter+jcmin)
      endif
11    continue

*
*  leave if nothing has been found...
*
      if (jcmin.le.0) go to 20
*
      jwr = jwr + 1
*
*----fill temporary variables.............
*
      do ii = 1,5
      vv_trk(ii,jwr) = par_trk(ii,miter+jcmin)
      ev_trk(ii,jwr) = erpr_trk(ii,miter+jcmin)
      enddo

      do ii = 1,4
      vv_trk2(ii,jwr) = par_trk2(ii,miter+jcmin)
      ev_trk2(ii,jwr) = erpr_trk2(ii,miter+jcmin)
      enddo   

      do ii = 1,mxplan
      vtrject(ii,jwr) = traject(ii,miter+jcmin)
      enddo

                        
      vchi2(jwr) = chi2_tab_ndf(miter+jcmin)
      vchi22(jwr) = chi2_tab(miter+jcmin)
      vxref(jwr) = x_ref(miter+jcmin)
      vxref2(jwr) = x_ref2(miter+jcmin)  
      lvass(1,jwr) = lassoc(1,jcmin)
      lvass(2,jwr) = lassoc(2,jcmin)
      chi2_tab_ndf(miter+jcmin) = 1.0e+09
      go to 10
*
*----- end of process......................
*

20    continue

*
*----- copy output list....................
*

      do kd = 1,jwr

      do ii = 1,5
      par_trk(ii,miter+kd) = vv_trk(ii,kd)
      erpr_trk(ii,miter+kd)= ev_trk(ii,kd)
      enddo

      do ii = 1,4
      par_trk2(ii,miter+kd) = vv_trk2(ii,kd)
      erpr_trk2(ii,miter+kd)= ev_trk2(ii,kd)
      enddo  
       
      do ii = 1,mxplan
      traject(ii,miter+kd) = vtrject(ii,kd)
      enddo

      chi2_tab_ndf(miter+kd) = vchi2(kd)
      chi2_tab(miter+kd) = vchi22(kd)
      x_ref(miter+kd) = vxref(kd)
      x_ref2(miter+kd) = vxref2(kd) 

      lassoc(1,kd) = lvass(1,kd)
      lassoc(2,kd) = lvass(2,kd)

      enddo
*
*----- reset counter.........................
*
      mlassc = jwr
      
      return
      end


