      subroutine recon(
     & ih1,ih2,ih3,ip1,ip2,ip3,phit1,phit2,phit3,
     & err1,err2,err3,
     &    mlassc,nrhit1,nrhit2,nrhit3,
     &    listhit1,listhit2,listhit3,
     &    chi2_tab,x_ref,par_trk,erpr_trk,
     &    chi2_tab2,x_ref2,par_trk2,erpr_trk2,itrg)

**************************************************
*                                                *
*    recon_2001/2002 reconstruction  module      *
*                                                *
**************************************************
#include "PARAM.h"
#include "TUNE.h"
#include "INI_STAT.h"
#include "PLANES.h"
#include "VERS.h"
 
       common/local_rec/rd_yy,rd_zz1,rd_zz2,coralout,
     & rd_tth1,rd_tth2,chi_par_res,itime,
     & nminpl21,nminpl22,nminpl11,nminpl12,epsi

   
       INTEGER ih1,ih2,ih3,itrg
       INTEGER ip1(*), ip2(*), ip3(*)

       DOUBLE PRECISION phit1(*), phit2(*), phit3(*)
       DOUBLE PRECISION err1(*),  err2(*),  err3(*)
        

      dimension lassoc(2,ntrack),ntori(ntrack)         ! for connection between detectors
      dimension jjflags(ntrack),liref(2)
      dimension prolimit(2),prolimit_ff(2)
      dimension chi2_tab(*),chi2_tab2(*)
      dimension chi2_tab_ndf(ntrack),chi2_tab2_ndf(ntrack)
      dimension par_trk(5,*),erpr_trk(5,*)             ! first helix
      dimension par_trk2(4,*),erpr_trk2(4,*)           ! second helix
      dimension x_ref(*),x_ref2(*)                     ! xref for first and second helix
      dimension traject(mxplan,ntrack)                 ! recstd. track position at various det.  
   
      dimension nrhit1(*), nrhit2(*), nrhit3(*)
      dimension listhit1(ntrack,*),listhit2(ntrack,*),
     & listhit3(ntrack,*)   
                                 
      logical modset,bad,coralout,fringe,extrapol,druk          
      logical pat_hit,fullprint,sub,once,once1,once2,once3,better
      logical version

      data itime/1/  

       if(itime.eq.1)then
       if(mc.ne.6)then
       print *
       print *,'=== start of RECON reconstruction ==='
       print *
       elseif(mc.eq.6)then
       print *
       print *,'=== start of RECON cluster profiles mode ==='
       print *,'           ( no reconstruction! )           '
       print *
       endif

         itime = 0
         maxlas=0
         maxev=0

         do i=1,7
         m_track(i)=0
         enddo
         do i=1,100
         lgolden(i)=0
         enddo

         rd_yy=rd_y(1)
         rd_zz1=rd_z(1)
         rd_zz2=rd_z(2)
         rd_tth1=rd_th(1)
         rd_tth2=rd_th(2)

         chi_par_res=chi_par

         nminpl21=nminpl2(1)
         nminpl22=nminpl2(2)
         nminpl11=nminpl1(1)
         nminpl12=nminpl1(2)
         prolimit(1)= chi_par2_tab(4,1)
         prolimit(2)= chi_par2_tab(4,2)
         prolimit_ff(1)= chi_par2_tab(5,1)
         prolimit_ff(2)= chi_par2_tab(5,2)

         prolim1=prolimit(1)
         prolim2=prolimit(2)

         nevnta=0
         coralout=.false.              ! true - printout
         fullprint=.false.
         epsi= 1.e-3        
         
       endif

         if(fullprint.and..not.druk)druk=.true.

**********************************************************
*                                                        *
*----   now clusters from coral into "hits" in recon     *
*                                                        *
**********************************************************

      if(mc.eq.0.or.mc.eq.2.or.mc.eq.4.or.mc.eq.6)
     &  call gluehits(ih1,ih2,ih3,ip1,ip2,ip3,phit1,phit2,phit3,
     &  err1,err2,err3)


        if(pat_hit)then
          print*,'event skipped in recon'
          return
        endif
   

* ---- if mc switch on => clusters transvered in mc_ini

          nevnta=nevnta+1
          nevntb=nevnta-1 

           itrig=itrg
       if(fullprint)print*,'nev,trigger in RECON',nevnta,itrig
       if(mc.eq.6)go to 1                            ! no reconstruction
      
*************************************************************************************
*                                                                                   *
*                                                                                   *
*                             PATTERN RECOGNITION                                   *
*                                                                                   *
*                                                                                   *
*************************************************************************************

    
********************************************************
*                                                      *
*                    PROJECTIONS                       *
*                  inpar2 = 1 - 9                      *
*                                                      *
********************************************************          
 
                                   
          modset = .false.                            ! skip FF if strategy modified.....

          miter=0
          miters=0
          mode1=2
          mode2=3
          mode3=100                                   ! refit switched off

          if(mode11.ge.1)mode3=mode11-1           
          if(mode.eq.0.or.mode.eq.1.or.mode.eq.2)mode=0

*------------------------------------------------------ multi projections

*----- starting parameters
         rd_y(1)=rd_yy
         rd_z(1)=rd_zz1
         rd_z(2)=rd_zz2  
         rd_th(1)=rd_tth1
         rd_th(2)=rd_tth2      
         chi_par = chi_par_res

         prolimit(1)=prolim1
         prolimit(2)=prolim2

         liref(1)=nminpl1(7)
         liref(2)=nminpl2(7)
         once1=.true.   
         once2=.true.
         once3=.true.
         libad=0
         if(.not.sub)then
            mp1=1  
            mp2=100
         endif
*                 
*
*------ first mp1 iterations for subsetup = micros+SDCs+straws
*------ next mp2 for subsetup micros + gems
*------ next for all detectors used
*
*         
* 
*------------------ start projections
*

         do i_ter=1,inpar2                                    ! projections iterations

         if(i_ter.eq.1)legro=1000
         if(.not.sub)i_tere=i_ter   
              once=.false.

         if(i_ter.le.mp1)then
           nminpl1(1)=nminpl11+1
           nminpl1(2)=nminpl12+1
           nminpl2(1)=nminpl21+1
           nminpl2(2)=nminpl22+1
         elseif(i_ter.gt.mp1.and.i_ter.le.mp1+mp2)then
           nminpl1(1)=nminpl1(3)+1
           nminpl1(2)=nminpl1(4)+1
           nminpl2(1)=nminpl2(3)+1
           nminpl2(2)=nminpl2(4)+1
         elseif(i_ter.gt.mp1+mp2)then
           nminpl1(1)=nminpl1(5)+1
           nminpl1(2)=nminpl1(6)+1
           nminpl2(1)=nminpl2(5)+1
           nminpl2(2)=nminpl2(6)+1
         endif
  
 40         continue

        if(sub)then                                           ! subsetups mode

          if(i_ter.le.mp1)then       
              i_tere=i_ter
              once=.true.
       if(i_tere.eq.1.and.once1)then
              call setup(111)                                ! preparing setup: micros + SDCs + straw 
              miter1=1
              once1=.false.
       elseif(i_ter.eq.mp1-1.and.once3)then            
              prolimit(1)=prolimit(1)/10.
              liref(1)=liref(1)-1
              liref(2)=liref(2)-1
              once3=.false.
       elseif(i_ter.eq.mp1.and..not.once3)then
              prolimit(1)=prolimit(1)/10.
              once3=.true.
            endif

          elseif(i_ter.gt.mp1.and.i_ter.le.mp1+mp2)then   
             i_tere=i_ter-mp1
              once=.true.

            if(i_tere.eq.1.and.once2)then                    ! preparing setup: micros + gems
              call setup(444)  
              miter1=miter+1
              prolimit(1)=prolim1
  
             rd_y(1)=rd_yy
             rd_z(1)=rd_zz1
             rd_z(2)=rd_zz2  
             rd_th(1)=rd_tth1
             rd_th(2)=rd_tth2      
             chi_par = chi_par_res

             liref(1)=nminpl1(8)
             liref(2)=nminpl2(8) 
             once2=.false.
            elseif(i_ter.eq.mp1+mp2-1.and..not.once2)then
              prolimit(1)=prolimit(1)/10.
              liref(1)=liref(1)-1
              liref(2)=liref(2)-1
              once2=.true.
            elseif(i_ter.eq.mp1+mp2.and..not.once1)then
              prolimit(1)=prolimit(1)/10.
              once1=.true.
            endif

          elseif(i_ter.gt.mp1+mp2.and.i_ter.le.inpar2)then  

           i_tere=i_ter-mp1-mp2
             once=.true.

           if(i_tere.eq.1.and.once3)then
            call setup(1)                                      ! preparing setup:  all detectors active
            miter1=miter+1

            prolimit(1)=prolim1
            prolimit(2)=prolim2
     
            rd_y(1)=rd_yy
            rd_z(1)=rd_zz1
            rd_z(2)=rd_zz2  
            rd_th(1)=rd_tth1
            rd_th(2)=rd_tth2      
            chi_par = chi_par_res  

            liref(1)=nminpl1(7) 
            liref(2)=nminpl2(7) 
            once3=.false. 
           endif

          endif        

          if((i_tere.le.3.and.i_ter.le.mp1+mp2).or.
     &       (i_tere.le.2.and.i_ter.gt.mp1+mp2))then            ! smaller for sci-fis 
          chi_par2(1) = chi_par2_tab(i_tere,1)
          chi_par2(2) = chi_par2_tab(i_tere,2)
          chi_par2(3) = chi_par2_tab(i_tere,3)
          elseif(i_ter.le.mp1+mp2)then
          chi_par2(1) = chi_par2_tab(3,1)
          chi_par2(2) = chi_par2_tab(3,2)
          chi_par2(3) = chi_par2_tab(3,3)
          endif
        endif                                                   ! end of subset mode

        if(nminpl1(1).eq.1.or.nminpl2(1).eq.1)once=.false.

           nminpl1(1)=nminpl1(1)-1
           nminpl1(2)=nminpl1(2)-1
           nminpl2(1)=nminpl2(1)-1
           nminpl2(2)=nminpl2(2)-1

          dzplin = dzplin_tab(i_tere)                           ! next set of parameters 
          dzpc = dzpc_tab(i_tere)

          if(.not.sub.and.i_ter.le.3)then
          chi_par2(1) = chi_par2_tab(i_ter,1)
          chi_par2(2) = chi_par2_tab(i_ter,2)
          chi_par2(3) = chi_par2_tab(i_ter,3)
          elseif(.not.sub.and.i_ter.gt.3)then
          chi_par2(1) = chi_par2_tab(3,1)
          chi_par2(2) = chi_par2_tab(3,2)
          chi_par2(3) = chi_par2_tab(3,3)
          endif


*----- now increasing path
           if(i_tere.eq.mode1)then
            rd_y(1)=(rd_y(2)+rd_y(1))/2.0
            rd_z(1)=(rd_z(3)+rd_z(1))/2.0
            rd_z(2)=(rd_z(4)+rd_z(2))/2.0
            rd_th(1)=(rd_th(3)+rd_th(1))/2.0
            rd_th(2)=(rd_th(4)+rd_th(2))/2.0

           elseif(i_tere.gt.mode1.and.i_tere.eq.mode2)then
            rd_y(1)=rd_y(2)
            rd_z(1)=rd_z(3)
            rd_z(2)=rd_z(4) 
            rd_th(1)=rd_th(3)
            rd_th(2)=rd_th(4)
           endif

         if(i_ter.eq.1.and.legro.eq.1000)legro=1
         If(i_tere.eq.1)largo=1
       
   
         if(mode.ne.0)then      

          if(i_tere.lt.mode1.and.i_tere.gt.1)then
           largo=2                                           ! reconstruction behind SM1 only
           if(mlassc.eq.0)largo=3                            ! association only           

          elseif(i_tere.ge.mode1)then

         if(.not.sub.and.i_ter.eq.mode)then
      call remover_hits(100,1,miter,mlassc,lassoc,           ! removing segments left (no setup mode)
     & listhit1,listhit2)
         elseif(sub.and.i_tere.eq.mode)then
      call remover_hits(100,miter1,miter,mlassc,lassoc,      ! removing segments left (setup mode)
     & listhit1,listhit2)
         endif


           largo=1
          endif
         endif

           if(i_tere.gt.mode3)largo = largo + 10000          ! just flag for reiteration

*--- version for 2001 data: 
           if(i_ter.gt.mp1+mp2)largo = -largo                ! flag for rec_det1_new and rec_det1_old    
       
       call projections(modset,bad,legro,largo,miter,
     & mlassc,lassoc,mlassc1,mlassc3,prolimit,liref,
     & chi2_tab,chi2_tab_ndf,chi2_tab2,chi2_tab2_ndf,
     & x_ref,x_ref2,par_trk,par_trk2,traject,
     & erpr_trk,erpr_trk2,listhit1,listhit2,listhit3)   
            legro=0
            if(bad)libad=libad+1
       
      if(fullprint)print *,'proj ',i_ter,': ',mlassc  

             do i=1,mlassc
               ntori(miter+i)=i_ter
             enddo
 
           miter=miter+mlassc
           maxt_extrap=maxt_extrap+mlassc3

           maxt_proj(i_ter)=maxt_proj(i_ter) + mlassc     
           mlassc1_proj(i_ter) = mlassc1_proj(i_ter) + mlassc1
        
       
       if(inpar2.eq.i_ter.and..not.fringe) go to 100
    

       if(mlassc.ne.0.and.i_tere.lt.mode)then
         call remover_hits(1,1,miter,mlassc,lassoc,                 ! used hits removed only behind
     &   listhit1,listhit2)
       elseif(mlassc.ne.0.and.i_tere.ge.mode)then
         call remover_hits(0,1,miter,mlassc,lassoc,                 ! all used hits removed
     &   listhit1,listhit2)
       endif

       if(once.and.better)go to 40                                ! once again with different condition
 
        modset = .false.
        enddo                                                     ! end of loop over projections
*             
*------------------------------------------------------ end of reconstruction 
*


********************************************************
*                                                      *
* fringe field reconstruction module:  fringe = .true. *
*                                                      *
*                    PROJECTIONS                       *
*                  j_te = 1 - 9                        *
*               (end if nothing found)                 *
*                                                      *
********************************************************     

       if(fringe.and..not.modset)then                              ! fringe field reconstruction
       if(sub.and.mp1+mp2.eq.inpar2)call setup(1)

         nminpl1(1) = nminpl11
         nminpl1(2) = nminpl12 

          do j_ter=1,inpar2                                        ! iteration loop for ffr
          i_one=0                                                  ! reconstruction

         call proj_fringe(i_one,miter,mlassc,mlassc1,mlassc2,      ! reconstruction
     & lassoc,prolimit_ff,chi2_tab,chi2_tab_ndf,
     & x_ref,x_ref2,par_trk,par_trk2,traject,erpr_trk,erpr_trk2,
     & listhit1,listhit2,listhit3)

      if(fullprint)print *,'fringe field rec. ',j_ter,': ' ,
     &      mlassc 
         
           maxt_fringe = maxt_fringe + mlassc  
           maxt_fringe1 = maxt_fringe1 + mlassc2 
           mlassc1_fringe = mlassc1_fringe + mlassc1

           do i=1,mlassc
           ntori(miter+i)=j_ter
           if(sub)ntori(miter+i)=j_ter + 100                                 
           enddo 
         
           miter=miter+mlassc

      if((mlassc.eq.0.and.j_ter.gt.1).or.j_ter.eq.inpar2)go to 100 ! jump to ntuple 

         call remover_hits(2,1,miter,mlassc,lassoc,                  ! used hits removed
     &   listhit1,listhit2)

           enddo
 
         endif
*
*------------------------- end of fringe field reconstruction
*
        
  100  continue
           
       if(sub)call setup(0)                                        ! restoration od subset 1

       if(libad.eq.inpar2.and.miter.eq.0.and.mlassc.eq.0)then
            l_bad=l_bad+1

          if(druk)print *,'bad event - reconstruction suspicious',
     &       libad
          if(l_bad.gt.0.and.l_bad.le.100)ib(l_bad)=nevnta
       endif

       mlassc1=mlassc   

c---- here should be last clean_up routine ... not needed?

       mlassc=miter
       maxtr=maxtr+mlassc                        ! counts number of entries for ntuple

       iv_count=iv_count+1
            
       if(itrig.eq.1)then
           maxtrIT=maxtrIT+mlassc
           iv_countIT=iv_countIT+1
       elseif(itrig.eq.2)then
           maxtrMT=maxtrMT+mlassc
           iv_countMT=iv_countMT+1
       elseif(itrig.eq.3)then
           maxtrITMT=maxtrITMT+mlassc
           iv_countITMT=iv_countITMT+1
       elseif(itrig.eq.4)then
           maxtrLT=maxtrLT+mlassc
           iv_countLT=iv_countLT+1
       elseif(itrig.eq.5)then
           maxtrLTIT=maxtrLTIT+mlassc
           iv_countLTIT=iv_countLTIT+1
       elseif(itrig.eq.6)then
           maxtrLTMT=maxtrLTMT+mlassc
           iv_countLTMT=iv_countLTMT+1
       elseif(itrig.eq.7)then
           maxtrLTITMT=maxtrLTITMT+mlassc
           iv_countLTITMT=iv_countLTITMT+1
       elseif(itrig.eq.8)then
           maxtrOT=maxtrOT+mlassc
           iv_countOT=iv_countOT+1
       elseif(itrig.eq.9)then
           maxtrOTIT=maxtrOTIT+mlassc
           iv_countOTIT=iv_countOTIT+1
       elseif(itrig.eq.10)then
           maxtrOTMT=maxtrOTMT+mlassc
           iv_countOTMT=iv_countOTMT+1
       elseif(itrig.eq.256)then
           maxtrINCMT=maxtrINCMT+mlassc
           iv_countINCMT=iv_countINCMT+1
       else
           maxtrAT=maxtrAT+mlassc
           iv_countAT=iv_countAT+1
       endif

       if(mc.eq.1.or.mc.eq.3.or.mc.eq.5)then
             mlassc2=0
             mlassc3=0
             mlassc4=0      

           do ja=1,mlassc
            if(par_trk(5,ja).ge.2)mlassc2=mlassc2+1
            if(par_trk(5,ja).ge.2.and.par_trk(5,ja).le.50)
     &         mlassc3=mlassc3+1
            if(par_trk(5,ja).ge.2.and.par_trk(5,ja).le.5)
     &         mlassc4=mlassc4+1
           enddo
          maxtra(1)=maxtra(1)+mlassc2
          maxtra(2)=maxtra(2)+mlassc3
          maxtra(3)=maxtra(3)+mlassc4

       endif

       call check_hit(mlassc,jjflags,ntori,                   ! counting hits
     & listhit1,listhit2,listhit3)

 1     continue

*
*------------ MONTE CARLO ------------*
*
      if(mc.eq.1.or.mc.eq.3.or.mc.eq.5)then

       call check_hit_mc(mlassc,jjflags,phit,                 ! full MC comparison & quality check
     & ntori,listhit1,listhit2,listhit3)

       call ntuple_mc(nevnta,itrig,mlassc,nos,
     & jjflags,ntori,listhit1,listhit2,x_ref,x_ref2,
     & par_trk, par_trk2,traject,erpr_trk,
     &     chi2_tab_ndf,chi2_tab,mc,druk)                     ! books MC ntuple         

*
*-------------- DATA -----------------*
*       
      elseif(mc.eq.0.or.mc.eq.2.or.mc.eq.4.or.mc.eq.6)then
     
       if(druk.and.mlassc.ne.0.and.mc.eq.0)then
       call recon_out(nevnta,mlassc,itrig,jjflags,
     &   ntori,chi2_tab_ndf,chi2_tab,par_trk,erpr_trk)

       elseif(mc.ne.0)then

       call ntuple(nevnta,itrig,mlassc,jjflags,ntori,
     & listhit1,listhit2,x_ref,x_ref2,par_trk,par_trk2,
     & traject,erpr_trk,chi2_tab_ndf,chi2_tab,mc,druk)        ! books ntuple for data
       endif                                                                                 

      endif

       if(mc.eq.6)return
*
*--------------------- now preparing output to coral
*

       if(maxlas.lt.mlassc)then
          maxlas=mlassc
          maxev=1
          lgolden(maxev)=nevnta
       elseif(maxlas.eq.mlassc)then
          maxev=maxev+1
          lgolden(maxev)=nevnta
       endif

       m_track(6)=maxlas
       m_track(7)=maxev
  
       if(mlassc.eq.1)then
           m_track(1)=m_track(1)+1
       elseif(mlassc.eq.2)then
           m_track(2)=m_track(2)+1
       elseif(mlassc.eq.3)then
           m_track(3)=m_track(3)+1
       elseif(mlassc.eq.4)then
           m_track(4)=m_track(4)+1
       elseif(mlassc.gt.4)then
           m_track(5)=m_track(5)+1
       endif

       do itr=1,mlassc

         if(sub)then
         if(ntori(itr).le.mp1)then
          ntpl_t=ntpl_s1
          ntpl2_t=ntpl2_s1
         elseif(ntori(itr).gt.mp1.and.ntori(itr).le.mp1+mp2)then
          ntpl_t=ntpl_s2
          ntpl2_t=ntpl2_s2
         elseif(ntori(itr).gt.mp1+mp2)then
          ntpl_t=ntpl_tmp
          ntpl2_t=ntpl2_tmp
         endif
         elseif(.not.sub)then
          ntpl_t=ntpl
          ntpl2_t=ntpl2
         endif

         pp=par_trk(5,itr)
         par_trk(5,itr)=1./pp
         erpr_trk(5,itr)= (erpr_trk(5,itr)/pp**2)**2           ! covariance matrix!!

*--------- just for a moment
         chi2_tab(itr) = -chi2_tab(itr)

         ind=0
         do k=1,ntpl_t
         if(listhit1(itr,k).ne.8888)then
            ind=ind+1
            listhit1(itr,ind)=listhit1(itr,k)
         endif
         enddo

         nrhit1(itr)=ind

         ind=0
         do k=1,ntpl2_t
         if(listhit2(itr,k).ne.8888)then
            ind=ind+1
            listhit2(itr,ind)=listhit2(itr,k)
         endif
         enddo

         nrhit2(itr)=ind
             
         ind=0
         do k=1,ntpl_sup
         if(listhit3(itr,k).ne.8888)then
            ind=ind+1
            listhit3(itr,ind)=listhit3(itr,k)
         endif
         enddo

         nrhit3(itr)=ind
       enddo


*------------  test print for output to Coral - for check only
       if(coralout)then

       print *,'Recon output to Coral : ev no., no of tracks',
     +   nevnta,mlassc

       do i=1,mlassc
        print *,'before SM1: number of hits:',nrhit1(i)          
        print *,'  hit list',(listhit1(i,j),j=1,nrhit1(i))
        print *,'after  SM1: number of hits:',nrhit2(i)          
        print *,'  hit list',(listhit2(i,j),j=1,nrhit2(i))
        print *,'after RICH before SM2: number of hits:',nrhit3(i)
        print *,'  hit list',(listhit3(i,j),j=1,nrhit3(i))     
        print *,' for first helix: chi2, x-ref',chi2_tab(i),x_ref(i)       
        print *,'    fitted parameters:',(par_trk(j,i),j=1,5)
        print *,'    errors on trk par:',(erpr_trk(j,i),j=1,5) 
        print *,' for second helix: chi2, x-ref',chi2_tab2(i),
     &                              x_ref2(i)
        print *,' fitted parameters (sec.helix):',
     &                           (par_trk2(j,i),j=1,4)
        print *,' errors on trk par (sec.helix):',
     &                          (erpr_trk2(j,i),j=1,4) 
       enddo
       endif

       
*
*-----------next event .............................
*

*            print*,'next event'

             
            return
            end
         




