      subroutine recon_out(iv,nr,itrig,jjflags,ntori,
     &          chi2_tab_ndf,chi2_tab,par_trk,erpr_trk)

******************************************************************
*                                                                *
*             routine printing output from recon                 *
*                  for every event (data)                        *
*                                                                *
******************************************************************
#include "PARAM.h"

      dimension chi2_tab_ndf(*),chi2_tab(*)
      dimension jjflags(*),ntori(*)
      dimension par_trk(5,*),erpr_trk(5,*)

      logical fri,extra

               do k =1,nr
               
               ppp=par_trk(5,k)                        ! reconstructed fitted momenta         
               sigp=erpr_trk(5,k)    
               chi2_ndf=chi2_tab_ndf(k)                ! chi2 of the fit
               chi2_line=chi2_tab(k)
 
               chi2=chi2_tab(k)                        ! chi2 global
               ndf=chi2/chi2_ndf  
               proba=prob(chi2,ndf)

           if(jjflags(k).lt.-100000)then               ! no extrapolation
               extra=.false.
               fri=.false.
               jjfl=-jjflags(k)-100000
               lihit3=0
               lihit1=mod(jjfl,1000)                   
               lihit2=(jjfl-lihit1)/1000
            elseif(jjflags(k).lt.0.and.jjflags(k).gt.-100000)then  ! fringe field reconstruction
               fri=.true.
               extra=.false.
               lihit1=-jjflags(k)
               lihit2=0
               lihit3=0
            else                                        ! extrapolation ok
               fri=.false.
               extra=.true.
               lihit1=mod(jjflags(k),1000)         
               lihit2=mod(jjflags(k),100000)/1000
               lihit3=jjflags(k)/100000
            endif
            
               in=ntori(k)
            bl=  sigp/abs(ppp)   

               if(fri)then
      print'("event, trig., p., chi2/ndf, prob. ,nh.1/2  , it.: ",
     & i10,i5,3f8.2,1x,e7.2,2i3,3x,i3," ff rec. ")',
*     & iv,itrig,ppp,chi2_ndg,chi2_line,proba,lihit1,in
     & iv,itrig,ppp,bl,chi2_ndf,proba,lihit1,lihit2,in

               else
               if(extra)then
      print'("event, trig., p., chi2/ndf, prob., nh.1/2/3, it.: ",
     & i10,i5,3f8.2,1x,e7.2,4i3," ext.SM2")',
*     &  iv,itrig,ppp,chi2_ndg,chi2_line,proba,lihit1,lihit2,in
     &  iv,itrig,ppp,bl,chi2_ndf,proba,lihit1,lihit2,lihit3,in

               elseif(.not.extra)then
      print'("event, trig., p., chi2/ndf, prob., nh.1/2  , it.: ",
     & i10,i5,3f8.2,1x,e7.2,2i3,3x,i3)',
*     &  iv,itrig,ppp,chi2_ndg,chi2_line,proba,lihit1,lihit2,in 
     &  iv,itrig,ppp,bl,chi2_ndf,proba,lihit1,lihit2,in 

               endif
               endif
 
             enddo
             
             return
             end














