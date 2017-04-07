       subroutine proj_fringe(l1,miter,mlassc,mlassc1,mlassc2,
     & lassoc,prolimit_ff,chi2_tab,chi2_tab_ndf,
     & x_ref,x_ref2,par_trk,par_trk2,traject,erpr_trk,erpr_trk2,
     & listhit1,listhit2,listhit3)
********************************************************************
*                                                                  *
*           fringe field  reconstruction (projection method)       *
*                                                                  *
*         looks for continuation  in SDC02 of tracks seen in SDC01 *
*                                                                  *
********************************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "HITS_IN.h"
#include "HITS_OUT.h"
#include "TUNE.h"
#include "LIN_ALL.h"
#include "INI_STAT.h"
         
      common/localny/jdc02,ifirst

      dimension chi2_tab(*),chi2_tab_ndf(*)
      dimension par_trk(5,*),erpr_trk(5,*)
      dimension par_trk2(4,*),erpr_trk2(4,*)
      dimension prolimit_ff(*)
      dimension traject(mxplan,*),traj(mxplan)          

      dimension xpar(5),cov(5,5)                
      dimension lassoc(2,*)
      dimension x_ref(*) ,x_ref2(*)
      dimension par_end(5)
      dimension vv_behind(4)
      dimension htt_sdc(8),jpl_sdc(8),nhi_sdc2(8)
      dimension err_sdc(8)

      dimension coor(10),dth_min(10),jth_min(10)
 
     
       dimension listhit1(ntrack,*),listhit2(ntrack,*),
     & listhit3(ntrack,*) 
       logical fringe,modset,flaga,fla_sdc2
       logical extrapol,druk,ffsel

       data ifirst/1/

*---------------------------------------
*modification 30Nov
 
       cut_ff_nb=frf1       !
       cut_ff_bd=frf2       ! kasia
       cut_ff_sp=frf3       !

*----------------------------------------------------
*--  search first sdc2 plane (jdc02)  ..........
*--     only sdc2
*

        if(ifirst.eq.1) then
*        print '("cuts":,3f8.4)',cut_ff_nb,cut_ff_bd,cut_ff_sp
        jdc02 = -99

        do kpl=1,ntpl2
        if(ityp2(kpl).eq.222)then 
           jdc02 = kpl
           go to 5
        endif
        enddo

5       continue

        ifirst = 0
        if (jdc02.le.0) then
         print *, "sdc02     not found ??",kpl,ntpl2,jdc02
         stop "proj_fringe"
        endif
        endif
*-------------------------------------------------------------------------
        x_behind = 0.5 * (xplan2(jdc02) + xplan2(jdc02+7))  ! center of sdc2 station
   
         nmin_dc = 5
         modset = .false.

         if(l1.eq.0) then

           li_m=0
           rd_y(1)=rd_yf
           rd_z(1)=rd_zf
           rd_th(1)=rd_thf 
           chi_par = chi_par_ff

           flaga=.true.
           call rec_det1_new(modset,fringe,flaga)      ! det. 1 reconstruction (new2 version)
         endif

         mlassc = 0
         mlassc2 = 0
       
       nlin_frfl = 0

*      print *,"nlin_fda=",nlin_fda

       do 20 kt = 1,nlin_fda

       call tkfit_frfl(kt,nmin_dc,xref,xpar,cov,chi2,ndd,
     1   xref_end,par_end,ifail,traj)

          if (ndd.le.5.or.ifail.ne.0) go to 20

          nlin_frfl = nlin_frfl + 1

           ndf=ndd-5
           chi_eff = chi2/(ndd-5)
           
*----- check  chi/nd limit .............................................

             proba=prob(chi2,ndf)   

            if(ffsel)then
      if(proba.lt.prolimit_ff(1).or.proba.gt.prolimit_ff(2).or.
     &       ifail.ne.0)then
             nlin_frfl = nlin_frfl - 1
             go to 20
      endif 
            else
      if (chi_eff.gt.chi_fringe.or.ifail.ne.0) then
             nlin_frfl = nlin_frfl - 1
             go to 20
      endif
            endif

      if(.not.ffsel.and.(proba.lt.prolimit_ff(1).or.
     &        proba.gt.prolimit_ff(2)))then
             nlin_frfl = nlin_frfl - 1
             go to 20
      endif 

          chi2_tab(miter+nlin_frfl) = chi2
          chi2_tab_ndf(miter+nlin_frfl) = chi_eff

          if (abs(xpar(5)).lt.1.0E-10) xpar(5) = 1.0E-10
          p_fit = 1./xpar(5)
          ep_fit = sqrt(cov(5,5))/xpar(5)**2

*--------------------------------------------------------------
          do kk=1,4
          par_trk(kk,miter+nlin_frfl) = xpar(kk)
          erpr_trk(kk,miter+nlin_frfl)= cov(kk,kk)
          par_trk2(kk,miter+nlin_frfl) = par_end(kk)
          erpr_trk2(kk,miter+nlin_frfl)= cov(kk,kk)     ! errors copied from starting point
          enddo

          par_trk(5,miter+nlin_frfl) = p_fit
          erpr_trk(5,miter+nlin_frfl)= ep_fit
          par_end(5) = p_fit

*------ position along the beam where the fit param. are given
          x_ref(miter+nlin_frfl)=xref
          x_ref2(miter+nlin_frfl)=xref_end
*------fitted trajectory (at plane position)...................
          do krw = 1,mxplan
          traject(krw,miter+nlin_frfl) = traj(krw)
          enddo
*
*----try to improve by extrap. behind M1 (sdc2)......
*

*         print '("par_end:",5f10.5)',(par_end(krw),krw=1,5)
*         print '("xref_end:",f10.5)',xref_end
*         print '("x_behind:",f10.5)',x_behind
          call track_fd(par_end,xref_end,x_behind,vv_behind)

*
       ifail_sdc = 0
       fla_sdc2=.false.
*
*--- quit if not inside  SDC02 ..............................
*
*     print '("x-y,xp,yp behind",4f10.4)', (vv_behind(krw),krw=1,4)        
      if (abs(vv_behind(1)).gt.0.7.or.abs(vv_behind(2)).gt.0.62)go to 18
*
*--- look for nearest hits in Y , Z and th planes of SDC02......................
*
         jok = 0
         do kpl = jdc02,jdc02+7
         yy = vv_behind(1) + vv_behind(3)*(xplan2(kpl)-x_behind)
         zz = vv_behind(2) + vv_behind(4)*(xplan2(kpl)-x_behind)
         coor(kpl) = sinn2(kpl)*yy - cosn2(kpl)*zz
          dth_min(kpl) = 1.0e+10
          jth_min(kpl) = -999
          jth_tt = -999
          do kk=1,jc2hits(kpl)
          cccc = clusters2(numht2(kpl,kk))
          if (abs(coor(kpl)-cccc).lt.dth_min(kpl))  then
            dth_min(kpl) = abs(coor(kpl)-cccc)
            jth_tt = kk
          endif
          enddo
         d_cut = cut_ff_nb
         if (jyzt2(kpl).ne.2) d_cut = cut_ff_bd
         if (dth_min(kpl).lt.d_cut)  then
          jok = jok + 1
          jth_min(kpl) = jth_tt
         endif
***      print '("kpl,dth",i5,e13.4)',kpl,dth_min(kpl)
         enddo
*
*--- require that at least 4 planes are OK....
*
*         print '(" jok_sdc ",i4)', jok
         if (jok.lt.4) then     
           ifail_sdc = -80
           go to 18 
         else
*
*
          do 10 kpl = jdc02,jdc02+7
          jpl_sdc(kpl) = -99
          if (jth_min(kpl).le.0) go to 10
           jhh = jth_min(kpl)
*--- mod. 31 july 2002 - kk
          jpl_sdc(kpl) = kpl
          nhi_sdc2(kpl)=numht2(kpl,jhh)
          htt_sdc(kpl) = clusters2(numht2(kpl,jhh))
          err_sdc(kpl) = clu_errors2(numht2(kpl,jhh))
10        continue    

**        print '("jpl-sdc",6i5)',(jpl_sdc(krw),krw=1,8)
**        print '("htt",6f10.4)',(htt_sdc(krw),krw=1,8)

          endif
*
*-------------------------
           ifail = 0

*          print '("OK till tkfit_frext")'

       call tkfit_frext(kt,jok,htt_sdc,err_sdc,jpl_sdc,xref,xpar,cov,
     1  chi2,ndd,xref_end,par_end,ifail,traj)

          ndf=ndd-5
          chi_eff = chi2/(ndd-5)
*      print '("chi_eff,ifail:",f10.5,i4)',chi_eff,ifail

*-----  ignore if chi2 too large...............................
*
            if(ifail.ne.0)go to 18
            proba=prob(chi2,ndf)   
*      print '("chi2,ndf,prob:",f10.5,i4,f10.5)',chi2,ndf,proba

            if(ffsel)then
      if(proba.lt.prolimit_ff(1).or.proba.gt.prolimit_ff(2).or.
     &       ifail.ne.0)go to 18
            
            else
      if (chi_eff.gt.chi_fringe.or.ifail.ne.0) go to 18
            endif

          if(.not.ffsel.and.(proba.lt.prolimit_ff(1).or.
     &        proba.gt.prolimit_ff(2)))go to 18

*
*------otherwise replace parameters by updated values...........
*      ......... make chi2 negative as a flag...................
*      .... not needed - just listhit2 gives info about (kk, Aug.2002)      
*

          chi2_tab(miter+nlin_frfl) = chi2
          chi2_tab_ndf(miter+nlin_frfl) = chi_eff
          if (abs(xpar(5)).lt.1.0E-10) xpar(5) = 1.0E-10
          p_fit = 1./xpar(5)
          ep_fit = sqrt(cov(5,5))/xpar(5)**2

          do kk=1,4
          par_trk(kk,miter+nlin_frfl) = xpar(kk)
          erpr_trk(kk,miter+nlin_frfl)= cov(kk,kk)
          par_trk2(kk,miter+nlin_frfl) = par_end(kk)
          erpr_trk2(kk,miter+nlin_frfl)= cov(kk,kk)     ! errors copied from starting point
          enddo

          par_trk(5,miter+nlin_frfl) = p_fit
          erpr_trk(5,miter+nlin_frfl)= ep_fit
          par_end(5) = p_fit

*------ position along the beam where the fit param. are given
          x_ref(miter+nlin_frfl)=xref
          x_ref2(miter+nlin_frfl)=xref_end
*
*------fitted trajectory (at plane position)...................
          do krw = 1,mxplan
          traject(krw,miter+nlin_frfl) = traj(krw)
          enddo
       
*------ flag for succesful bridged
           fla_sdc2=.true.
           mlassc2=mlassc2+1
*       
*------------------------
*
  18       continue

           lassoc(1,nlin_frfl) = kt             ! line nb in det 1
           lassoc(2,nlin_frfl) = 2000           ! line nb in det 2 - flag

  20      continue

        mlassc = nlin_frfl

*---------------------------end fringe field...........................

        mlassc1= nlin_frfl

         if(mlassc.eq.0)then
            return
         else

            do l=1,mlassc
            do k=1,ntpl
            listhit1(miter+l,k)=line1(k,l)           
            enddo    

            do k=1,ntpl2
              listhit2(miter+l,k)=  8888                                    ! flag  
            if(fla_sdc2)then                    
             do j=1,8
             if(jpl_sdc(j).eq.k)then
              listhit2(miter+l,k)=  nhi_sdc2(j)
              go to 118
             endif
             enddo
 118         continue
            endif
            enddo

       if(.not.fla_sdc2)listhit2(miter+l,1)=  9999                          ! flag
       if(fla_sdc2)listhit2(miter+l,1)=listhit2(miter+l,1)+9999             ! flag for fringe-field 

            do k=1,ntpl_sup
            listhit3(miter+l,k)=  8888                                      ! flag           
            enddo
            enddo

         endif
        return
        end










