
       subroutine assoc(set,miter,mlassc,lassoc,prolimit,
     &  liref,chi2_tab,chi2_tab_ndf,x_ref,par_trk,erpr_trk,
     &  x_ref2,par_trk2,erpr_trk2,traject)

**********************************************************
*                                                        *
*   routine for connecting lines from first and second   *
*          telescope through magnetic field SM1          *
*                                                        *
**********************************************************

#include "PARAM.h"
#include "TUNE.h"
#include "PLANES.h"
#include "MAGNET.h"
#include "HITS_OUT.h"
#include "LIN_ALL.h"





          dimension lassoc(2,*),prolimit(*),liref(*)
          dimension liner(ntrack)
          dimension chi2_tab(*),chi2_tab_ndf(*)
          dimension xpar(5),cov(5,5)
          dimension par_trk(5,*),erpr_trk(5,*)
          dimension par_trk2(4,*),erpr_trk2(4,*)
          dimension x_ref(*),x_ref2(*)
          dimension par_end(4)
          dimension traject(mxplan,*),traj(mxplan)

          logical set,sele,seles,norsel

*          print *,prolimit(1),prolimit(2)
*          print *,' in assoc..nr lines ',nlin_fda,nlin_fdb

 
             loop=0
             if(set)loop=1
         
          do 10 k=1,nlin_fda
      
*             print*,'in asscoc k ',k,(line1(j,k),j=1,ntpl)
 

C ------ loop on lines in 2d detector................

          do 20 l=1,nlin_fdb 

*             print*,'in asscoc l ',l,(line2(j,l),j=1,ntpl2)


          if (abs(zplina(k)-zplinb(l)).gt.dzplin) go to 20

C--------  check on Z at X= x0_m1 ............  centre of magnet

          zmca = zlina(k) + (x0_m1 - xlina(k))*zplina(k)
          zmcb = zlinb(l) + (x0_m1 - xlinb(l))*zplinb(l)

          if (abs(zmca-zmcb).gt.dzpc) go to 20

          mlassc = mlassc + 1               ! # of reconstructed lines
          seles=.true.

C-kk-Nov.15 -- spline fit (RW & KK)
*-------now do track fit.................................................

          if(loop.eq.1)ifail=1000

 13       continue

          call tkfit_pr(k,l,xref,xpar,cov,chi2,ndd,
     &     xref_end,par_end,ifail,thres,traj)


*---------TEMP!!!!!
           if (ifail.eq.998) then
            mlassc = mlassc -1
            go to 20
           endif
*----------------- 
         
* ---- reject if chi2 greater than limit................................

           if(chi2.lt.0.)chi2 = 100000.

           if(ndd.le.5.and.loop.ne.0)then
             loop=2
           elseif(ndd.le.5.and.loop.eq.0)then
             mlassc = mlassc -1
             go to 20
           endif

           ndf=ndd-5
           chi_eff = chi2/(ndd-5) 
*           print*,mlassc,chi_eff,loop,ifail
           
           
*----- select chi/nd limit according to momentum........................
            ichi = 2
            if (abs(xpar(5)).gt.0.5) ichi = 3
            if (abs(xpar(5)).lt.0.0666) ichi = 1
*---- low momentum : ichi =3, medium: ichi=2 , high: ichi=1
*----------------------------------------------------------------------- 
       
         proba=0.0
         sele=.true.    

            if(norsel)then
            if(ifail.eq.0)proba=prob(chi2,ndf)
      if(proba.lt.prolimit(1).or.proba.gt.prolimit(2).or.ifail.ne.0)
     &      sele=.false.
            else
      if(chi_eff.gt.chi_par2(ichi).or.ifail.ne.0)sele=.false.
            endif 

           if(.not.sele)then       
             if(seles)then
               ndff=ndf
               seles=.false.
             endif
            if(loop.eq.2)then
               loop=1
           
            elseif(loop.eq.1.and.chi_eff.lt.100000)then          
            loop=loop+1
*           print*,'reiteration',loop
           go to 13

            endif

             if(loop.eq.1)then
             do kop=1,ntpl
              if(line1(kop,k).gt.8888)line1(kop,k)=line1(kop,k)-8888
             enddo
             do kop=1,ntpl2
              if(line2(kop,l).gt.8888)line2(kop,l)=line2(kop,l)-8888
             enddo
             endif
 
             mlassc = mlassc -1
             go to 20

           else

             if(loop.eq.1)then                      
             do kop=1,ntpl
              if(line1(kop,k).gt.8888)line1(kop,k)=line1(kop,k)-8888
             enddo
             do kop=1,ntpl2
              if(line2(kop,l).gt.8888)line2(kop,l)=line2(kop,l)-8888
             enddo
*
*------    in case of succes just to update clusters    
*
             elseif(loop.eq.2)then     
                loop=1
                mclu1=0
                mclu2=0

             do kop=1,ntpl
               liner(kop)=99                              ! flag
              if(line1(kop,k).gt.8888)then
               liner(kop)=line1(kop,k)-8888
               line1(kop,k)=8888
              endif
               if(line1(kop,k).ne.8888)mclu1=mclu1+1
             enddo
             do kop=1,ntpl2
               liner(ntpl+kop)=99                         ! flag
              if(line2(kop,l).gt.8888)then
               liner(ntpl+kop)=line2(kop,l)-8888
               line2(kop,l)=8888
              endif
               if(line2(kop,l).ne.8888)mclu2=mclu2+1
             enddo
*
*------    if success but patological one - back
*        
             If((ndf+5)*100/(ndff+5).lt.phit)mclu1=0    ! to protect artificial combination after refit

             if(mclu1.lt.liref(1).or.mclu2.lt.liref(2))then      
              do kop=1,ntpl
              if(liner(kop).ne.99)line1(kop,k)=liner(kop)
              enddo
              do kop=1,ntpl2
              if(liner(ntpl+kop).ne.99)line2(kop,l)=liner(ntpl+kop)
              enddo

               mlassc=mlassc-1
               go to 20
             endif
 

             endif
           endif


*
*--- probability test if chi2 selection used -- suppres by choose prob.<0 and.>1
*
          if(.not.norsel)then 
             proba=prob(chi2,ndf)
          if(proba.lt.prolimit(1).or.proba.gt.prolimit(2))then
             mlassc = mlassc -1
             go to 20
          endif
          endif
           
             chi2_tab_ndf(miter+mlassc)=chi_eff
             chi2_tab(miter+mlassc)=chi2

*-----------------------------------------------------------------------
*-------momentum.....................................................

         p_fit = 1./xpar(5)
         ep_fit = sqrt(cov(5,5))/xpar(5)**2

*        print*, mlassc,k,l,p_fit,ep_fit,chi_eff

      
*------ store for output ............................................
         do krw = 1,4
         par_trk(krw,miter+mlassc) = xpar(krw)
         erpr_trk(krw,miter+mlassc)= sqrt(cov(krw,krw))
         enddo

         par_trk(5,miter+mlassc) = p_fit
         erpr_trk(5,miter+mlassc)= ep_fit

*------ position of the first helix
          x_ref(miter+mlassc)=xref
          
*------ end point....second helix.....................................
         do krw = 1,4
         par_trk2(krw,miter+mlassc) = par_end(krw)
         erpr_trk2(krw,miter+mlassc)= erpr_trk(krw,miter+mlassc)
         enddo

*----- fitted trajectory (at plane positions).........................
          do krw = 1,mxplan
          traject(krw,miter+mlassc) = traj(krw)
          enddo

          x_ref2(miter+mlassc) = xref_end
*------------------------------------------------------------------------

          if(mlassc.gt.ntrack) then
            print *,'mlassc too large ',mlassc,ntrack
            return        ! emergency exit
          endif



          lassoc(1,mlassc) = k              ! line nb in det 1
          lassoc(2,mlassc) = l              ! line nb in det 2



20        continue                          ! end loop over l
10        continue                          ! end loop over k    

*
*-------remove tracks with common segments in one of the telescopes...
*
    
           call remdup(miter,mlassc,lassoc,chi2_tab,chi2_tab_ndf,
     1  x_ref,par_trk,erpr_trk,x_ref2,par_trk2,erpr_trk2,traject)

 
    
*          print '("iter = ",i3)',miter
*          do krw =1,mlassc
*          print '("assoc up-dw-chi2-p",2i4,3f10.5)',
*     1    lassoc(1,krw),lassoc(2,krw),chi2_tab(miter+krw),
*     2    chi2_tab_ndf(miter+krw),par_trk(5,miter+krw)
*          enddo
      
          return
          end





