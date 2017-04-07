       subroutine projections(modset,bad,l,largo,miter,
     & mlassc,lassoc,mlassc1,mlassc3,prolimit,liref,
     & chi2_tab,chi2_tab_ndf,chi2_tab2,chi2_tab2_ndf,
     & x_ref,x_ref2,par_trk,par_trk2,traject,
     & erpr_trk,erpr_trk2,listhit1,listhit2,listhit3)

********************************************************************
*                                                                  *
*               projections reconstruction method                  *
*                                                                  *
********************************************************************
#include "PARAM.h"
#include "TUNE.h"
#include "INI_STAT.h"
#include "PLANES.h"
#include "VERS.h"

          
      dimension lassoc(2,*),prolimit(*),liref(*)   ! for connection between detectors
      dimension chi2_tab(*),chi2_tab2(*)
      dimension chi2_tab_ndf(*),chi2_tab2_ndf(*)
      dimension par_trk(5,*),erpr_trk(5,*)         ! first helix
      dimension par_trk2(4,*),erpr_trk2(4,*)       ! second helix
      dimension x_ref(*),x_ref2(*)                 ! xref for first and second helix
      dimension traject(mxplan,*)                  ! recstrd. track at each detector

      dimension listhit1(ntrack,*),listhit2(ntrack,*),
     & listhit3(ntrack,*)   

      logical modset,fringe,jfringe,extrapol,druk
      logical set,set1,noas,bad,printout,flaga
      logical version,sub 

       set= .false.       
       set1= .false. 
       printout=.false. 
       bad = .false.
       jfringe = .false. 
       noas=.false.

       li_m=0
       mlassc =0
       mlassc1=0
       mlassc3=0
     
       if(.not.sub)set1=.true.

       if(largo.lt.0)then
        set1=.true.
        largo=-largo
       endif
        
       if(largo.ge.10000)then
          largo=largo-10000
          set=.true.
       endif

       if(largo.eq.2)then                          ! reconstruction only behind SM1
         go to 150
       elseif(largo.eq.3)then                      ! only association
         go to 200
       elseif(largo.eq.1)then

        if(version)then                            ! 2001 setups
             flaga=.false.
         call rec_det1_mgmg(modset,jfringe)              
           if(set1.and..not.modset)then              
         call rec_det1_new(modset,jfringe,flaga)   ! new2 ver. det. 1 reconstruction  
         if(.not.modset)call rec_det1_old          ! old det. 1 reconstruction      
           endif

        elseif(.not.version)then       ! 2002 setups
             flaga=.true.

         call rec_det1_new(modset,jfringe,flaga)   ! new2 ver. det. 1 reconstruction 
         if(.not.modset)call rec_det1_old          ! old det. 1 reconstruction      
 
        endif

       endif

         if(li_m.ge.mxlin.or.modset)then
*           print*,'det 1',li_m,modset
           bad=.true.
         endif

*-----repeat for 2d detector ....(behind M1).................


       if(l.eq.1) nevntb= nevntb+ 1 


 150   continue       
        
       li_m=0
             
           call rec_det2(noas)                     ! det. 2 reconstruction 

          mlassc=0
          mlassc1=0 

         if(li_m.ge.mxlin)then
*           print*,'det 2',li_m
           bad=.true.
         endif
      
*------- associate lines in the 2 detectors.............................

 200   continue

*          noas=.false.

      if(.not.noas)then
      call assoc(set,miter,mlassc,lassoc,prolimit,liref,
     &   chi2_tab,chi2_tab_ndf,
     &   x_ref,par_trk,erpr_trk,
     &   x_ref2,par_trk2,erpr_trk2,traject)  
      else
       mlassc=0
      endif
  
        if(mlassc.gt.ntrack) then                   ! bad events: too many lines
            bad = .true.
            print *,'mlassc > limit',mlassc,ntrack
            return
        endif

      if(nevnta.eq.nevntb)then  
  
       mlassc1=mlassc

       if(mlassc.ge.1)
     & call clean_up(miter,mlassc,lassoc,
     & chi2_tab,chi2_tab_ndf,phit,nos,
     & x_ref,par_trk,erpr_trk,
     & x_ref2,par_trk2,erpr_trk2,
     & ntpl,ntpl2,
     & listhit1,listhit2,0)                         ! no clean_up on - only rewrites


*
*----------Extrapolation from RICH to SM2..............................
*

      if (mlassc.ge.1) then
      call extrap_sm2(extrapol,miter,mlassc,mlassc3,lassoc,par_rochi,
     &   x_ref2,par_trk2,erpr_trk2,chi2_tab2,chi2_tab2_ndf,
     &   listhit3)


         do jt = 1,mlassc
*  
*------- second chi2 is same as first one if no extrap. was done... 
*
        if (chi2_tab2(miter+jt).lt.-998.) then
          chi2_tab2(miter+jt) = chi2_tab(miter+jt)
          chi2_tab2_ndf(miter+jt) = chi2_tab_ndf(miter+jt)
         endif

         if (x_ref2(1).gt.0.and.printout)  then
          print '(" jt ",i2," End point at Z=",f10.4)', jt,
     &    x_ref2(miter+jt)
          print '(" X=",f10.5," +/-",f10.5)', 
     &           par_trk2(1,miter+jt),erpr_trk2(1,miter+jt)
          print '(" Y=",f10.5," +/-",f10.5)', 
     &           par_trk2(2,miter+jt),erpr_trk2(2,miter+jt)
          print '(" SL.Y=",f10.5," +/-",f10.5)', 
     &           par_trk2(3,miter+jt),erpr_trk2(3,miter+jt)
          print '(" SL.Z=",f10.5," +/-",f10.5)',
     &           par_trk2(4,miter+jt),erpr_trk2(4,miter+jt)
      print '(" CHI2 =",e13.5)',chi2_tab2(miter+jt)
      print '(" CHI2/ndf =",e13.5)',chi2_tab2_ndf(miter+jt)
         endif
         enddo
  
      endif

        else
        print *,'number of events before and behind SM1 differs!'
        stop
      endif


        return
        end
















