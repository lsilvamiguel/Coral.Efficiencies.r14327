      subroutine reconini(imc,in,road,proj,
     1  rlengh,chi,tar,nbpl,region,xmag,nblack) 
**************************************************
*                                                *
*       recon_2001 initialization module         *
*          fulfils common blocks                 *
*                                                *
**************************************************
#include "PARAM.h"
#include "TUNE.h"
#include "INI_STAT.h"
#include "PLANES.h"
#include "TARGET.h"
#include "MAGNET.h"
#include "VERS.h"

       REAL  road(*),proj(*),rlengh(*),chi(*),
     1       tar(*),region(*),xmag(*) 


       integer nblack(*)

       dimension  imc(*), in(*) ,nbpl(*)
       logical printout,fringe,fringe2,extrapol,druk
       logical version,sub,ffsel,norsel,better

         eps=1.0e-5

         if(nblack(1).eq.0)then
         i_option2=0                  !  0 - all planes on, 1 - planes from "black list" off
         else
         i_option2=1 
         endif

          n_black1=nblack(2)
          n_black2=nblack(3)


         printout=.false.
         druk=.false.
         better=.false.                                                         ! better but slow
         version=.false.                                                        ! 2002 setup

      if(in(4).eq.8)then
         in(4)=0
         better=.true.
      endif 
      if(in(4).eq.1.or.in(4).eq.3.or.in(4).eq.5.or.in(4).eq.7)then
          version=.true.                                                        ! 2001 setup
      elseif(in(4).eq.0.or.in(4).eq.2.or.in(4).eq.4.or.in(4).eq.6)then
          version=.false.                                                       ! 2002 setup
      else
       print*,'error in choosing setups - 2002 probability switched on'
       in(4)=0
      endif
      

         inpar2=in(1)                 ! options for reconstruction method (see below)
               
            sub=.false.
         if(inpar2.gt.30000.and.inpar2.lt.40000)then      ! subsetups mode switch on:
            inpar2=inpar2-30000                           ! no print-out: print 3 same as 0
            sub=.true.
         elseif(inpar2.gt.40000.and.inpar2.lt.50000)then  ! subsetups mode switch on:   
            inpar2=inpar2-40000                           ! partially print-out print 4 same as 1
            sub=.true.
            druk = .true.                                 ! reconstruction info written for every event
         elseif(inpar2.gt.50000)then                      ! subsetups mode switch on:
            inpar2=inpar2-50000                           ! full print-out print 5 same as 2
            sub=.true.
            printout = .true.                             ! planes & parameters  print-out
            druk=.true.                                   ! reconstruction info written for every event
         elseif(inpar2.gt.20000.and.inpar2.lt.30000)then  ! subsetups mode switch off: all detectors on
            inpar2 = inpar2 - 20000                       ! full print-out print 2 same as 5
            printout = .true.                             ! planes & parameters  print-out
            druk=.true.                                   ! reconstruction info written for every event
         elseif(inpar2.gt.10000.and.inpar2.le.20000)then  ! subsetups mode switch off: all detectors on
            inpar2 = inpar2 - 10000                       ! partially print-out print 1 same as 4
            druk = .true.                                 ! reconstruction info written for every event
         endif

      print*                    
      print*          
      print*
      print*,'    RRRRR     EEEEE   CCCCCC     OOOOOO     NN   NN    '
      print*,'    R   RR   EE      CCC---CC  OOO----OOO  NNNN---NN   '
      print*,'    RRRRR   EEEEE   CCC        OOO----OOO  NN  N  NN   '
      print*,'    RR--RR   EE      CCC---CC  OOO----OOO  NN---NNNN   '
      print*,'    RR---RR   EEEEE   CCCCCC     OOOOOO     NN   NN    '
      print*                            
      print*
      print*,'             VERSION 10.0 for 2002 data               '  
      print*
      print*,' ********************************************************'
      print*,' *                                                      *'
      print*,' *          RECON COMPASS FIRST SPECTROMETER            *'
      print*,' *               reconstruction program                 *'
      print*,' *                                                      *'
               if(.not.version)
     &print*,' *         2002 version - please check option file      *'
               if(version)
     &print*,' *         2001 version - please check option file      *'
      print*,' *                                                      *'
         if(.not.printout)then
      print*,' *   authors:   Roland Windmolders & Warsaw group       *'
      print*,' *   contact person:  Krzysztof Kurek                   *'
      print*,' *                                                      *'
      print*,' ********************************************************'
      print*
         else
      print*,' *                                                      *'
      print*,' *                      authors:                        *'
      print*,' *                                                      *'
      print*,' *     Roland Windmolders (Mons/Bonn)(main author)      *'
      print*,' *                         and                          *'
      print*,' *                    Warsaw Group:                     *'
      print*,' *    Katarzyna Kowalik (kowalik@fuw.edu.pl - C++ part) *'
      print*,' *    Krzysztof Kurek (kurek@fuw.edu.pl - fortran)      *'
      print*,' *             Ewa Rondio & Adam Mielech                *'
      print*,' *                                                      *'
      print*,' ********************************************************'
       print *            
       print *,'  **************************************************'
       print *,'  *                                                *'
       print *,'  *             please note that:                  *'
       print *,'  *   RECON is using legal SI units: meters!       *'
       print *,'  *   recon frame is the same as comgeant frame    *'
       print *,'  *        and differs from CORAL frame !          *'
       print *,'  *   conversions are done automatically in recon  *'
       print *,'  *   hits convention in recon is:  Y*SIN - Z*COS  *'
       print *,'  *                warning:                        *'
       print *,'  *       straws stations planes can be "glued"    *'
       print *,'  *   to simplify combinatorics - see code below   *'
       print *,'  *                                                *'
       print *,'  **************************************************'  
       print *
    
      
*---------------------  now geometry files: target.pos & magnet.pos
    
         endif

     
      xtarg0=tar(1)   !  
      hltarg=tar(2)   !  target parameters
      ytarg=tar(3)    !
      ztarg=tar(4)    !
      
  
          if(printout)then 
      print *
      print *,'*****************************************************'
      print *,'*                                                   *'
      print *,'*   RECON options and  parameters (from opt_file):  *'
      print *,'*                                                   *'
          if(.not.version)
     &print *,'*                  setup 2002                       *'
          if(version)
     &print *,'*                  setup 2001                       *'
      print *,'*                                                   *'
      print *,'*===================================================*'
      print *

      print '("*----- target & solenoid parameters: ",4f7.3)',      
     &    xtarg0,hltarg,ytarg,ztarg
     

        if(ytarg.ne.ztarg)then
          print*, 'non-cylindrical targets cell shape !'
          print*, 'non-cylindrical targets cell shape !'

       endif
       endif

      x0_m1=xmag(1)                                     
      hlm1=xmag(2)                 ! magnet parameters

        if(printout)
     & print '("*---- SM1 magnet parameters: ",2f7.3)', x0_m1,hlm1
    

         mc=imc(1)
         mmc=mc

         if(mc.gt.100)then
          mp1=(mc-mod(mc,100))/100
          mp22=mod(mc-mp1*100,10)
          mp2=(mc-mp1*100-mp22)/10
          mc=mc-mp1*100-mp2*10     
         elseif(mc.gt.10.and.mc.lt.100)then
          mp1=(mc-mod(mc,10))/10
          mc=mc-mp1*10
          mp1=0
          mp2=0
         endif    

         if(mmc.gt.100.and..not.sub)then
       print*,'WARNING: subset option off, mp1/mp2 parameters ignored'
           mp1=0
           mp2=0
           
         elseif(mmc.lt.100.and.sub)then
       print*,'WARNING: subsetups option on, mp1/mp2 parameters ',
     &        'not specified'
       print*,'         subsetups option off ' 
           mp1=0
           mp2=0
           sub=.false.
         elseif(mmc.lt.100.and..not.sub)then
           mp1=0
           mp2=0
         endif
                                          
         mmc=mc                       ! for recon_end
  
        if(sub)then 
        if(mp2.eq.0)then
           print*,'WARNING:'
           print*,' mp2=0 not allowed in this mode - switched ty mp2=1'
           mp2=1
         elseif(mp1.eq.0)then
           print*,'WARNING:'
           print*,' mp1=0 not allowed in this mode - switched ty mp1=1'
           mp1=1
         endif
         endif

        nn_tel1=imc(2)    
        nn_tel2=imc(3)                ! min # of hits in tel.1 and 2 for acc. track (external hits)  

        if(sub.and.printout)then
        print '("*---- subsetups opt. on, mc/mp1/mp2: ",3i3)',mc,mp1,mp2 
        elseif(.not.sub.and.printout)then
        print '("*---- subsetups opt. off,mc/mp1/mp2: ",2i3)',mc,mp1,mp2
        endif   

        if(mc.eq.6)then
        print*
        print*,'*****************************************'
        print*,'*                                       *'
        print*,'* only hit profiles - NO reconstruction *'
        print*,'*                                       *'
        print*,'*****************************************'
        endif
 
        if((mc.eq.1.or.mc.eq.3.or.mc.eq.5).and.printout)
     & print '("*---- min number of hits in tel.1/2: ",2i3)',
     & nn_tel1,nn_tel2

        
          par_rochi(1) = road(14)     ! roadwidth for extrapolation from Rich to SM2
          rd_thf=road(13)             ! roadwidth in theta for det.1 in fringe field
          rd_yf=road(11)
          rd_zf=road(12)

         frf1=road(15)            
         frf2=road(16)                !  fringe field parameters 
         frf3=road(17) 

         distm2_yz = region(3)   
         distm2_th = region(4)

         sagi_tot  = region(5)    ! for a 1 GeV particle for a distance of 1m and B.dl=0.1 T
*--------------  f = 0.3 * B.dl * L / (8 * P) 

         rd_comp = region(6)      ! max.dist* max. slope = 0.11 * 0.1 (curvature negligeable on 11 cm) 
         rd_com2 = region(7)      ! max.dist* max. slope = 0.15 * 0.1 (curvature negligeable on 11 cm)

         distmx_yz = region(1)
         distmx_th = region(2)

         do i=1,2
          rd_y(i)=road(i)             ! roadwidth in y for det.1 (target shot)
          rd_z(i) = road(i+2)         ! roadwidth in z for det.1 & det.2 (target shot)
          rd_z(i+2) = road(i+4)       ! roadwidth in z for det.2 (target shot)
          rd_th(i)=road(i+6)          ! roadwidth in theta for det.1 & det.2
          rd_th(i+2)=road(i+8)        ! roadwidth in theta for det.1 & det.2

          det1_act(i)= 0.0            ! active region (in m) for telescope 1
          det2_act(i)= 0.0            ! active region (in m) for telescope 2

          nminpl1(i)=nbpl(i)          ! min nr of addit. planes in tel.1 (th and all)(micromegas+SDC1) 
          nminpl1(i+2)=nbpl(i+2)      ! min nr of addit. planes in tel.1 (th and all)(micromegas)
          nminpl1(i+4)=nbpl(i+4)      ! min nr of addit. planes in tel.1 (th and all)(all detectors)
          nminpl1(i+6)=nbpl(i+6)      ! min nr of addit. planes in tel.1 (th and all)(for refit set.1/2)
          nminpl2(i)=nbpl(i+8)        ! min nr of addit. planes in tel.2 (th and all)(SDC1/2 + straws)
          nminpl2(i+2)=nbpl(i+10)     ! min nr of addit. planes in tel.2 (th and all)(gems)
          nminpl2(i+4)=nbpl(i+12)     ! min nr of addit. planes in tel.2 (th and all)(all detectors)
          nminpl2(i+6)=nbpl(i+14)     ! min nr of addit. planes in tel.2 (th and all)(for refit set.1/2)
         enddo

          ten1=0.1                    ! one-tenth parameter for rec_det1_new2
          ten2=0.1                    ! one-tenth parameter for rec_det2
          ten3=0.1                    ! one-tenth parameter for rec_det1_mgmg

       if(printout)then
      print '("*---- active region for tel. 1 (0/0 default)(in m.): ",
     &        2f5.2)',det1_act(1),det1_act(2)              
      print '("*---- active region for tel. 2 (0/0 default)(in m.): ",
     &       2f5.2)', det2_act(1),det2_act(2)
      print '("*---- one-tenth parameter for reconstruction routines: ",
     &       3f5.2)', ten1,ten2,ten3
      print '("*---- distmx_yz, distmx_th, distm2_yz,distm2_th ",
     &       4f5.2)',distmx_yz, distmx_th, distm2_yz,distm2_th
      print '("*---- sagi_tot, rd_comp, rd_com2 ",
     &       3f5.3)',sagi_tot, rd_comp, rd_com2
      if(sub)then
      print '("*---- min no of planes in tel.1 (th & all)(mic.+SDC1): "
     &   ,2i3)',     nminpl1(1),nminpl1(2)
      print '("*---- min no of planes in tel.1 (th & all)(mic.): "
     &   ,2i3)',     nminpl1(3),nminpl1(4)
      print '("*---- min no of planes in tel.1 (th & all)(all det.): "
     &   ,2i3)',     nminpl1(5),nminpl1(6)
      else
      print '("*---- min no of planes in tel.1 (th & all): "
     &   ,2i3)',     nminpl1(1),nminpl1(2)
      endif
      print '("*---- min no of planes in tel.1 (for refit - set1/2): "
     &   ,2i3)',     nminpl1(7),nminpl1(8)

      if(sub)then
      print '("*---- min no of planes in tel.2 (th & all)(SDC23+st.): "
     &   ,2i3)',     nminpl2(1),nminpl2(2)
      print '("*---- min no of planes in tel.2 (th & all)(gems): "
     &   ,2i3)',     nminpl2(3),nminpl2(4)
      print '("*---- min no of planes in tel.2 (th & all)(all det.): "
     &   ,2i3)',     nminpl2(5),nminpl2(6)
      else
      print '("*---- min no of planes in tel.2 (th & all): ",2i3)',
     &               nminpl2(1),nminpl2(2)
      endif
      print '("*---- min no of planes in tel.2 (for refit - set1/2): "
     &   ,2i3)',     nminpl2(7),nminpl2(8)


      print '("*---- roads in Y: ",2f7.3)',rd_y(1),rd_y(2)
      print '("*---- roads in Z: ",4f7.3)',(rd_z(i),i=1,4)
      print '("*---- roads in th: ",4f7.3)',(rd_th(i),i=1,4)
      print '("*---- roads in Y for fringe field: ",f7.3)',rd_yf
      print '("*---- roads in Z for fringe field: ",f7.3)',rd_zf
      print '("*---- roads in th for fringe field: ",f7.3)',rd_thf
      print '("*---- roads for extrapolation to SM2: ",f7.3)',
     &        par_rochi(1)
      print '("*---- fringe field par: ",3f7.3)',frf1,
     &   frf2,frf3     
       endif


        do i=1,8
         dzplin_tab(i)=proj(i)          ! tolerances in proj. - slopes
         dzpc_tab(i)=proj(i+8)          ! tolerances in proj. - positions
        enddo

         plow_mcs=proj(8)               ! plow_mcs parameter
         proj(8)=proj(7)

         if(printout)then
      print '("*---- tolerances in slopes: ",8f8.2)',(proj(j),j=1,8)
      print '("*---- tolerances in positions: ",8f8.2)',
     &             (proj(j),j=9,16)
         endif
  
        dmcs1=rlengh(1)   
        dmcs2=rlengh(2)                 ! frac. of rad.lenght for each plane in mic.,st/SDC,sci-fi,gems
        dmcs3=rlengh(3)   
        dmcs4=rlengh(4)   

         if(printout)then
      print *,'*---- radiative lenghts for used detectors: '
      print '("*----   micromega stations: ",f7.4)',dmcs1
      print '("*----   Saclay drifts chambers & straws: ",f7.4)',dmcs2
      print '("*----   sci-fi stations: ",f7.4)',dmcs3
      print '("*----   gem stations: ",f7.4)',dmcs4
         endif

       
        phit=float(in(2))               ! % of hits (proj)  for "correct " asssoc. and for refit
                        
        min3=in(3)                      ! parameter for error scaling,spline refit and fudge for mcs

         call decod(min3,thres,fudge,vrscal)
 

        ffsel=.false.                                             ! chi2 selection in ff
        if(in(4).eq.0.or.in(4).eq.1.or.in(4).eq.4.or.in(4).eq.5)
     &         ffsel=.true.                                       ! probability selection in ff
      
        chi_fringe = chi(16)                         ! chi2 parameter for fringe field reconstruction 
        chi_par2_tab5(1)=chi(14)                     ! probability limits (lower)
        chi_par2_tab5(2)=chi(15)                     ! probability limits (upper)

*--- decoding inpar2 parameter:

        call decod_inpar2(inpar2,mode11,mode,fringe,extrapol)

        fringe2=fringe                  ! for planes3.F - multiple coulomb switch

        if(printout)then
      print '("*---- number of iterations in projections: ",i3)',inpar2
      print '("*---- mode option: ",i3)',mode
       if(fringe)then 
      print* ,'*-fringe field reconstruction on: (',fringe,fringe2,')'
       else
      print* ,'*-fringe field reconstruction off:(',fringe,fringe2,')'
       endif
        if(mode11.ne.0)then
      print '("*---- refit in spline on: ",i3)',mode11
        else
      print '("*---- refit in spline off: ",i3)',mode11
        endif
        if(extrapol)then
      print* ,'*---- extrapolation on: (',extrapol,')'
        else
      print* ,'*---- extrapolation off: (',extrapol,')'
        endif
      print '("*---- min number of hits for good/ghost track: ",f4.1)',
     &         phit
        endif
       
        chi_par=chi(1)                  ! chi^2/nd for chi_min for segment (normal reconstruction)
        chi_par_ff=chi(16)              ! chi^2/nd for chi_min for segment in fringe field

       do i=1,3
        chi_par2_tab1(i)=chi(i+1)  
        chi_par2_tab2(i)=chi(i+4)  
        chi_par2_tab3(i)=chi(i+7)       !  chi^2/nd for tkfit

        chi_par2_tab4(i)=chi(i+10)      ! probability limits & chi2 for extrapolation to SM2
       enddo 

        par_rochi(2) = chi_par2_tab4(3) ! chi2 parameter for extrapolation from Rich to SM2

                norsel=.false.                                    ! chi2 selection cuts
        if(in(4).eq.0.or.in(4).eq.1.or.in(4).eq.6.or.in(4).eq.7)
     &          norsel=.true.                                     ! probability selection cuts

       if(printout)then
      print '("*--- parameter for re-fit in spline (thres) ",
     & f5.2)',thres
      print '("*--- parameter for mcs (fudge) ",
     & f5.2)',fudge
      print '("*--- parameter for error scaling (vrscal) ",
     & 5f5.2)',vrscal
      print '("*--- chi2 parameter for segment reconstruction: ",
     & f7.3)',chi_par
      print '("*--- chi2 parameter for segment reconstruction in ff: ",
     & f7.3)',chi_par_ff
      print '("*--- chi2 parameter for extrapolation to SM2: ",
     & f7.3)',par_rochi(2)

      if(norsel)then
      print '("*--- selection mode: probability cuts               ")'
      else
      print '("*--- selection mode: chi2 cuts                      ")'
      endif
      print '("*--- probability limits: ",
     & 2f9.6)', (chi_par2_tab4(i),i=1,2)
      print '("*--- chi2 for spline; (diff. mom. range), iter. 1/2: ",
     & 3f7.2)', (chi_par2_tab1(i),i=1,3)
      print '("*---                                     iter. 3/4: ",
     & 3f7.2)', (chi_par2_tab2(i),i=1,3)
      print '("*---                                   iter. 5/6/7: ",
     & 3f7.2)', (chi_par2_tab3(i),i=1,3)

      if(ffsel)then
      print '("*--- selection mode in ff rec.: probability cuts    ")' 
      else
      print '("*--- selection mode: chi2 cuts in ff rec.           ")'
      endif
      print '("*--- probability limits for ff: ",
     & 2f9.6)',(chi_par2_tab5(i),i=1,2)
      print '("*--- chi2 parameter for fringe field reconstruction: ",
     & f5.2)',chi_fringe

       print '("*--- plow_mcs parameter: ",f9.2)',plow_mcs
       print *,'*'
       print *,'*=====================================================*'
       print *
       endif

*----------------- parameters initialization

              iv_count=0
              itr_count=0
              maxtr=0
*-- triggers
              iv_countIT=0
              iv_countMT=0
              iv_countLT=0
              iv_countAT=0
              iv_countITMT=0
              iv_countLTIT=0
              iv_countLTMT=0
              iv_countLTITMT=0
              iv_countOT=0
              iv_countOTMT=0
              iv_countOTIT=0
              iv_countINCMT=0

               maxtrIT=0
               maxtrMT=0
               maxtrLT=0
               maxtrAT=0
               maxtrITMT=0
               maxtrLTIT=0
               maxtrLTMT=0
               maxtrLTITMT=0
               maxtrOT=0
               maxtrOTIT=0
               maxtrOTMT=0
               maxtrINCMT=0


*--- end of trigger mask initialization      
               ntr_licz=0
               ntr_ac_licz =0
               ntr_ac_licz_new =0
           
               maxt_fringe =0
               maxt_fringe1=0
               mlassc1_fringe =0

               maxt_extrap=0

            l_bad=0                                  ! # of bad events (hang reconstruction codes)

            do io=1,100
             nos(io)=0
             if(io.lt.10)maxt_proj(io)=0
             if(io.lt.10)mlassc1_proj(io)=0
             if(io.lt.8)ntra_ac(io)=0
             if(io.lt.4)maxtra(io)=0
            enddo
 
      do i=1,3
        chi_par2_tab(1,i) = chi_par2_tab1(i)
        chi_par2_tab(2,i) = chi_par2_tab2(i)
        chi_par2_tab(3,i) = chi_par2_tab3(i)
        chi_par2_tab(4,i) = chi_par2_tab4(i)
        if(i.lt.3)chi_par2_tab(5,i) = chi_par2_tab5(i)
      enddo

             
         lolita=1

       print *
       print *,'reconstruction method: '
 
       if(inpar2.eq.1)then
         print *,lolita,'. single projections (',inpar2,')'
       elseif(inpar2.eq.2)then
         print *,lolita,'. double projections (',inpar2,')'
       elseif(inpar2.eq.3)then
         print *,lolita,'. triple projections (',inpar2,')'
       elseif(inpar2.eq.4)then
         print *,lolita,'. quadra projections (',inpar2,')'
       elseif(inpar2.eq.5)then
         print *,lolita,'. qvinta projections (',inpar2,')'
       elseif(inpar2.eq.6)then
         print *,lolita,'. sexta projections (',inpar2,')'
       elseif(inpar2.eq.7)then
         print *,lolita,'. septa projections (',inpar2,')'
       elseif(inpar2.eq.8)then
         print *,lolita,'. octa projections (',inpar2,')'
       elseif(inpar2.eq.9)then
         print *,lolita,'. nenta projections (',inpar2,')'
        else
         print *,'incorrect value of inpar2 - how many projections?'
         stop
       endif

          lolita=lolita+1 

         if(fringe)then
            print*,lolita,'. fringe field on (',fringe,')'
         elseif(.not.fringe)then
            print*,lolita,'. fringe field off (',fringe,')'
         endif 

          lolita=lolita+1

          if(mode11.eq.0)then
          print *, lolita,'. re-fit in spline off !!'
          elseif(mode11.eq.1)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,   ' iteration'  
          elseif(mode11.eq.2)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,    ' iteration'  
          elseif(mode11.eq.3)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,    ' iteration'  
          elseif(mode11.eq.4)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,    ' iteration'  
          elseif(mode11.eq.5)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,    ' iteration'  
          elseif(mode11.eq.6)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,    ' iteration'  
          elseif(mode11.eq.7)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &      mode11,    ' iteration'  
          elseif(mode11.eq.8)then
       print *, lolita,'. re-fit in spline on - starting from ',
     &     mode11,     ' iteration'  
          else
       print *,'incorrect value of mode11 - refit in spline off !!'
          mode11=0
          endif

          lolita=lolita+1

          if(extrapol)then
       print *, lolita,'. extrapolation to SM2 on (',extrapol,')'
          else
       print *, lolita,'. extrapolation to SM2 off !! (',extrapol,')'
          endif

          lolita=lolita+1

          if(sub.and.version)then                         ! 2001 setups
       print *, lolita,'. subsetups mode on (',sub,')'
       print *, ' subsetups are defined as follows: '
             if(inpar2.gt.mp1+mp2)then
       print *, ' micromegas + SDC + straws -> first (',mp1,')',
     &           'iterations '
       print *, ' micromegas + gems  -> next (',mp2,') iterations '
       print *, ' all detectors   -> next (',inpar2-mp1-mp2,')',
     &           'iterations '
             elseif(inpar2.eq.mp1+mp2.and.mp2.gt.0)then
       print *, ' micromegas + SDC + straws -> first (',mp1,') ',
     &           'iterations '
       print *, ' micromegas + gems  -> next (',mp2,') iterations '
             elseif(inpar2.eq.mp1)then
       print *, ' micromegas + SDC + straws -> (',mp1,') iterations '
             elseif(inpar2.gt.mp1.and.inpar2.lt.mp1+mp2)then
       print *, ' micromegas + SDC + straws -> first (',mp1,') ',
     &           'iterations '
       print *, ' micromegas + gems  -> (',inpar2-mp1,') iterations '
             elseif(inpar2.lt.mp1)then
       print *, ' micromegas + SDC + straws '
             endif
           endif

        if(sub.and..not.version)then                         ! 2002 setups (no straws up to now)                    
       print *, lolita,'. subsetups mode on (',sub,')'
       print *, ' subsetups are defined as follows: '
             if(inpar2.gt.mp1+mp2)then
       print *, ' micromegas/SDC1 + SDC2/3/straws  -> first (',mp1,')',
     &           'iterations '
       print *, ' micromegas + gems -> next (',mp2,') iterations '
       print *, ' all detectors     -> next (',inpar2-mp1-mp2,')',
     &           'iterations '
             elseif(inpar2.eq.mp1+mp2.and.mp2.gt.0)then
       print *, ' micromegas/SDC1 + SDC2/3/straws  -> first (',mp1,') ',
     &           'iterations '
       print *, ' micromegas + gems -> next (',mp2,') iterations '
             elseif(inpar2.eq.mp1)then
       print *, ' micromegas/SDC1 + SDC2/3/straws  -> (',mp1,')',
     &           ' iterations '
             elseif(inpar2.gt.mp1.and.inpar2.lt.mp1+mp2)then
       print *, ' micromegas/SDC1 + SDC2/3/straws  -> first (',mp1,') ',
     &           'iterations '
       print *, ' micromegas + gems -> (',inpar2-mp1,') iterations '
             elseif(inpar2.lt.mp1)then
       print *, ' micromegas/SDC1 + SDC2/3/straws '
             endif
           endif

       print *
       print *
       print *,'===  RECON initialization done !  === '
       print *
       print *

          if(mc.gt.1)call ntuple_ini(mc)                    ! books ntuples
       return
       end



      subroutine decod_inpar2(inpar2,mode11,mode,fringe,extrapol)
*****************************************************************
*                                                               *
*----  just to decod inpar2 parameter                           *
*                                                               *
*****************************************************************

        logical fringe,extrapol

        if(inpar2.gt.10.and.inpar2.le.38)then
         mode11=0         
         mode=0

         inn3=mod(inpar2,10)
         inn4=(inpar2-inn3)/10
         inpar2=inn3

         if(inn4.eq.1)then
         fringe=.true.
         extrapol=.true.
         elseif(inn4.eq.2)then
         fringe=.false.
         extrapol=.true.
         elseif(inn4.eq.3)then
         fringe=.true.
         extrapol=.false.
         else
           print *, 'WARNING: wrong "ffex" combination; switch to 0'
           fringe=.false.
           extrapol=.false.  
         endif

       
      elseif(inpar2.gt.100.and.inpar2.le.838)then
       mode11=0

       inn= mod(inpar2,100)
       mode=(inpar2-inn)/100

       if(inn.gt.10.and.inn.le.38)then
         inn3=mod(inn,10)
         inn4=(inn-inn3)/10
         inpar2=inn3

         if(inn4.eq.1)then
         fringe=.true.
         extrapol=.true.
         elseif(inn4.eq.2)then
         fringe=.false.
         extrapol=.true.
         elseif(inn4.eq.3)then
         fringe=.true.
         extrapol=.false.
         else
           print *, 'WARNING: wrong "ffex" combination; switch to 0'
           fringe=.false.
           extrapol=.false.    
         endif

       elseif(inn.gt.0.and.inn.le.8)then
          inpar2=inn
          fringe=.false.
          extrapol=.false.
       else
          print *,'errors: check inpar2'
          stop
       endif

      elseif(inpar2.gt.1000.and.inpar2.le.8818)then
       inn= mod(inpar2,1000)
       mode11=(inpar2-inn)/1000
   
       if(inn.ge.100)then
        inn2 = mod(inn,100)
        mode=(inn-inn2)/100
       else
        mode=0
        inn2=inn
       endif

       if(inn2.gt.10.and.inn2.le.38)then
         inn3=mod(inn2,10)
         inn4=(inn2-inn3)/10
         inpar2=inn3

         if(inn4.eq.1)then
         fringe=.true.
         extrapol=.true.
         elseif(inn4.eq.2)then
         fringe=.false.
         extrapol=.true.
         elseif(inn4.eq.3)then
         fringe=.true.
         extrapol=.false.
         else
           print *, 'WARNING: wrong "ffex" combination; switch to 0'
           fringe=.false.
           extrapol=.false.
         endif
     
       elseif(inn2.gt.0.and.inn2.le.8)then
          inpar2=inn2
          fringe=.false.
          extrapol=.false.
       else
          print *,'errors: check inpar2'
          stop
       endif

      elseif(inpar2.gt.0.and.inpar2.le.8)then
        fringe =.false.
        mode=0 
        mode11=0
        extrapol=.false.
      else
          print *,'errors: check inpar2'
          stop
      endif


          if(mode.gt.inpar2)then
             mode =0
             print *, 'WARNING: wrong "mode" combination; switch to 0'
          endif
          if(mode11.gt.inpar2)then
             mode11 =0
      print *, 'WARNING: wrong "mode1" combination; refit in spline off'
          endif


          return
          end

      subroutine decod(min3,thres,fudge,vrscal)
*****************************************************************
*                                                               *
*----  just to decod inpar2 parameter                           *
*                                                               *
*****************************************************************

         dimension vrscal(*)

      if(min3.lt.1111111)then
         print*,' warning - no error scaling  parameter specified;'
         print*,' changed for 1 for all detector'
         print*,' warning - no thres  parameter specified;'
         print*,' changed for thres = 10000 - no refit'
         print*,' warning - no fudge parameter specified;'
         print*,' changed for fudge = 1'
          fudge = 1. 
          thres = 10000.
          do i=1,5
          vrscal(i)=1
          enddo
          return
       endif

     
        min33=mod(min3,1000000)
        vrscal(1)=float((min3-min33)/1000000)
        min333=mod(min33,100000)
        vrscal(2)=float((min33-min333)/100000)
        min3333=mod(min333,10000)
        vrscal(3)=float((min333-min3333)/10000)
        min33333=mod(min3333,1000)
        vrscal(4)=float((min3333-min33333)/1000)
        min333333=mod(min33333,100)
        vrscal(5)=float((min33333-min333333)/100)
        min3333333=mod(min333333,10)
        fudge=float((min333333-min3333333)/10)
       
        thres=min3-vrscal(1)*1000000-vrscal(2)*100000-
     &  vrscal(3)*10000-vrscal(4)*1000-vrscal(5)*100 -fudge*10

*----fine tuning of errors.....RW......Oct.10,2002
*
        do k= 1,5
        ivrs = vrscal(k) + 0.001
        if (ivrs.eq.9) then
         vrscal(k) = 0.25
        elseif (ivrs.eq.8)  then
         vrscal(k) = 0.50
        else
         vrscal(k) = 1.25**(ivrs-1)
        endif
        enddo
****************

        return
        end
    







