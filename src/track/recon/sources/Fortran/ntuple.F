
      subroutine ntuple_ini(mc)
**************************************************
*                                                *
*     ntuples initialization routine             *
*                                                *
**************************************************
#include "NTUPLE.h"

*--- ntuple stuff
*
*     mc parameter defines kind of output
*     mc=2 & 4 - MC info is not used at all, run for data
*     mc=3 & 5 - makes correspodence reconstructed - generated track


      character*8 chtags_mc1(ncolum_mc1)
      character*8 chtags_mc2(ncolum_mc2)
      character*8 chtags_mc3(ncolum_mc3)
      character*8 chtags_mc4(ncolum_mc4)
  
      character*8 chtags1(ncolum1)
      character*8 chtags2(ncolum2)
      character*8 chtags3(ncolum3)
      character*8 chtags4(ncolum4) 

      data nprime/1000/
      print *, 'ntuple_ini starts, mc=',mc


      if(mc.le.1) then
        print *,' error in ntuple parameter!!'
        print *,'  ntuple not booked - just ignored'
        return
      endif

      if(mc.eq.3)then                       ! full MC quality check - standard

      data chtags_mc1/'nr_ev','itr_no', 'nghgem','nghsf','nghsup',
     1        'nghd','nghm','phit','ngh1','ngh2','nrh1','nrh2',
     2    'nrh3','sigp','pr','pg','pg_1pl','z_1pl',
     3    'delp','chi2','proba','f_flag','ex_flag','id_track','na_flag',
     4  'xx','yy','xp','yp','xx_1pl','yy_1pl','xp_1pl','yp_1pl',
     5  'sigx','sigy','sigxp','sigyp'/

      data chtags_mc2/'nr_ev','ng','na','nr','i_muon'/  
 
      elseif(mc.eq.5)then                   ! full MC quality check - full

      data chtags_mc3/'nr_ev','itr_no','nghsup',
     1  'phit','ngh1','ngh2','nrh1','nrh2','nrh3','sigp','pr','pg',
     2  'pg_1pl','z_1pl','delp','chi2','proba','f_flag','ex_flag',
     3  'id_track','na_flag','xx_1pl','yy_1pl','xp_1pl','yp_1pl',
     4  'sigx','sigy','sigxp','sigyp','x_hel',
     5  'hel1','hel2','hel3','hel4',
     6 'y_ref1','z_ref1','y_ref1','z_ref2',
     7 'y_ref3','z_ref3','y_ref4','z_ref4',
     3 'y_ref5','z_ref5','y_ref6','z_ref6',
     4 'y_ref7','z_ref7','y_ref8','z_ref8',
     3 'y_ref9','z_ref9','y_ref10','z_ref10',
     7 'trig_mask','iter',
     4 'scifi1_1','scifi1_2','scifi1_3',
     5 'mic1_1','mic1_2','mic1_3','mic1_4',
     5 'mic2_1','mic2_2','mic2_3','mic2_4',
     6 'scifi2_1','scifi2_2','scifi2_3',
     7 'mic3_1','mic3_2','mic3_3','mic3_4',
     8 'sdc1_1','sdc1_2','sdc1_3','sdc1_4',
     8 'sdc1_5','sdc1_6','sdc1_7','sdc1_8',
     9 'sdc2_1','sdc2_2','sdc2_3','sdc2_4',
     9 'sdc2_5','sdc2_6','sdc2_7','sdc2_8',
     & 'gem1_1','gem1_2','gem1_3','gem1_4',
     1 'sdc3_1','sdc3_2','sdc3_3','sdc3_4',
     1 'sdc3_5','sdc3_6','sdc3_7','sdc3_8',
     2 'gem2_1','gem2_2','gem2_3','gem2_4', 
     3 'st1_1','st1_2','st1_3','st2_1','st2_2','st2_3',
     4 'st3_1','st3_2','st3_3','st4_1','st4_2','st4_3',
     5 'st5_1','st5_2','st5_3','st6_1','st6_2','st6_3',
     6 'st7_1','st7_2','st7_3','st8_1','st8_2','st8_3',
     8 'st9_1','st9_2','st9_3','st10_1','st10_2','st10_3',
     & 'st11_1','st11_2','st11_3','st12_1','st12_2','st12_3',
     7 'gem3_1','gem3_2', 'gem3_3','gem3_4',
     1 'scifi3_1','scifi3_2',
     4 'tscifi1_1','tscifi1_2','tscifi1_3',
     5 'tmic1_1','tmic1_2','tmic1_3','tmic1_4',
     5 'tmic2_1','tmic2_2','tmic2_3','tmic2_4',
     6 'tscifi2_1','tscifi2_2','tscifi2_3',
     7 'tmic3_1','tmic3_2','tmic3_3','tmic3_4',
     8 'tsdc1_1','tsdc1_2','tsdc1_3','tsdc1_4',
     8 'tsdc1_5','tsdc1_6','tsdc1_7','tsdc1_8',
     9 'tsdc2_1','tsdc2_2','tsdc2_3','tsdc2_4',
     9 'tsdc2_5','tsdc2_6','tsdc2_7','tsdc2_8',
     & 'tgem1_1','tgem1_2','tgem1_3','tgem1_4',
     1 'tsdc3_1','tsdc3_2','tsdc3_3','tsdc3_4',
     1 'tsdc3_5','tsdc3_6','tsdc3_7','tsdc3_8',
     2 'tgem2_1','tgem2_2','tgem2_3','tgem2_4', 
     3 'tst1_1','tst1_2','tst1_3','tst2_1','tst2_2','tst2_3',
     4 'tst3_1','tst3_2','tst3_3','tst4_1','tst4_2','tst4_3',
     5 'tst5_1','tst5_2','tst5_3','tst6_1','tst6_2','tst6_3',
     6 'tst7_1','tst7_2','tst7_3','tst8_1','tst8_2','tst8_3',
     8 'tst9_1','tst9_2','tst9_3','tst10_1','tst10_2','tst10_3',
     & 'tst11_1','tst11_2','tst11_3','tst12_1','tst12_2','tst12_3',
     9 'tgem3_1','tgem3_2','tgem3_3','tgem3_4',
     1 'tscifi3_1','tscifi3_2'/
     
       data chtags_mc4/'nr_ev','ng','na','nr','trig_mask','i_muon'/


      elseif(mc.eq.2)then                  !   data - standard

       data chtags1/'nr_ev','ir_no','pr','sigp','chi2','proba',
     1               'nrh1','nrh2','nrh3','f_flag','ex_flag',
     2               'x_hel','hel1','hel2','hel3','hel4',
     3    'y_ref1','z_ref1','y_ref2','z_ref2','trig_mask','iter'/

       data chtags2/'nr_ev','nr','trig_mask'/

      elseif(mc.eq.4.or.mc.eq.6)then       !   data - full
 
       if(mc.eq.4)then
       data chtags3/'nr_ev','ir_no','pr','sigp','chi2','proba',
     1 'nrh1','nrh2','nrh3','f_flag','ex_flag','x_hel',
     2 'hel1','hel2','hel3','hel4',
     3 'y_ref1','z_ref1','y_ref2','z_ref2',
     4 'y_ref3','z_ref3','y_ref4','z_ref4',
     3 'y_ref5','z_ref5','y_ref6','z_ref6',
     4 'y_ref7','z_ref7','y_ref8','z_ref8',
     3 'y_ref9','z_ref9','y_ref10','z_ref10',
     4 'trig_mask','iter', 
     4 'scifi1_1','scifi1_2','scifi1_3',
     5 'mic1_1','mic1_2','mic1_3','mic1_4',
     5 'mic2_1','mic2_2','mic2_3','mic2_4',
     6 'scifi2_1','scifi2_2','scifi2_3',
     7 'mic3_1','mic3_2','mic3_3','mic3_4',
     8 'sdc1_1','sdc1_2','sdc1_3','sdc1_4',
     8 'sdc1_5','sdc1_6','sdc1_7','sdc1_8',
     9 'sdc2_1','sdc2_2','sdc2_3','sdc2_4',
     9 'sdc2_5','sdc2_6','sdc2_7','sdc2_8',
     & 'gem1_1','gem1_2','gem1_3','gem1_4',
     1 'sdc3_1','sdc3_2','sdc3_3','sdc3_4',
     1 'sdc3_5','sdc3_6','sdc3_7','sdc3_8',
     2 'gem2_1','gem2_2','gem2_3','gem2_4', 
     3 'st1_1','st1_2','st1_3','st2_1','st2_2','st2_3',
     4 'st3_1','st3_2','st3_3','st4_1','st4_2','st4_3',
     5 'st5_1','st5_2','st5_3','st6_1','st6_2','st6_3',
     6 'st7_1','st7_2','st7_3','st8_1','st8_2','st8_3',
     8 'st9_1','st9_2','st9_3','st10_1','st10_2','st10_3',  
     & 'st11_1','st11_2','st11_3','st12_1','st12_2','st12_3',
     9 'gem3_1','gem3_2','gem3_3','gem3_4',
     1 'scifi3_1','scifi3_2',
     4 'tscifi1_1','tscifi1_2','tscifi1_3',
     5 'tmic1_1','tmic1_2','tmic1_3','tmic1_4',
     5 'tmic2_1','tmic2_2','tmic2_3','tmic2_4',
     6 'tscifi2_1','tscifi2_2','tscifi2_3',
     7 'tmic3_1','tmic3_2','tmic3_3','tmic3_4',
     8 'tsdc1_1','tsdc1_2','tsdc1_3','tsdc1_4',
     8 'tsdc1_5','tsdc1_6','tsdc1_7','tsdc1_8',
     9 'tsdc2_1','tsdc2_2','tsdc2_3','tsdc2_4',
     9 'tsdc2_5','tsdc2_6','tsdc2_7','tsdc2_8',
     & 'tgem1_1','tgem1_2','tgem1_3','tgem1_4',
     1 'tsdc3_1','tsdc3_2','tsdc3_3','tsdc3_4',
     1 'tsdc3_5','tsdc3_6','tsdc3_7','tsdc3_8',
     2 'tgem2_1','tgem2_2','tgem2_3','tgem2_4', 
     3 'tst1_1','tst1_2','tst1_3','tst2_1','tst2_2','tst2_3',
     4 'tst3_1','tst3_2','tst3_3','tst4_1','tst4_2','tst4_3',
     5 'tst5_1','tst5_2','tst5_3','tst6_1','tst6_2','tst6_3',
     6 'tst7_1','tst7_2','tst7_3','tst8_1','tst8_2','tst8_3',
     8 'tst9_1','tst9_2','tst9_3','tst10_1','tst10_2','tst10_3',   
     & 'tst11_1','tst11_2','tst11_3','tst12_1','tst12_2','tst12_3',
     7 'tgem3_1','tgem3_2','tgem3_3','tgem3_4',
     1 'tscifi3_1','tscifi3_2'/
       endif
     
       data chtags4/'nr_ev','nr','trig_mask',
     4 'scifi1_1','scifi1_2','scifi1_3',
     5 'mic1_1','mic1_2','mic1_3','mic1_4',
     5 'mic2_1','mic2_2','mic2_3','mic2_4',
     6 'scifi2_1','scifi2_2','scifi2_3',
     7 'mic3_1','mic3_2','mic3_3','mic3_4',
     8 'sdc1_1','sdc1_2','sdc1_3','sdc1_4',
     8 'sdc1_5','sdc1_6','sdc1_7','sdc1_8',
     9 'sdc2_1','sdc2_2','sdc2_3','sdc2_4',
     9 'sdc2_5','sdc2_6','sdc2_7','sdc2_8',
     & 'gem1_1','gem1_2','gem1_3','gem1_4',
     1 'sdc3_1','sdc3_2','sdc3_3','sdc3_4',
     1 'sdc3_5','sdc3_6','sdc3_7','sdc3_8',
     2 'gem2_1','gem2_2','gem2_3','gem2_4', 
     3 'st1_1','st1_2','st1_3','st2_1','st2_2','st2_3',
     4 'st3_1','st3_2','st3_3','st4_1','st4_2','st4_3',
     5 'st5_1','st5_2','st5_3','st6_1','st6_2','st6_3',
     6 'st7_1','st7_2','st7_3','st8_1','st8_2','st8_3',
     8 'st9_1','st9_2','st9_3','st10_1','st10_2','st10_3',
     & 'st11_1','st11_2','st11_3','st12_1','st12_2','st12_3',
     7 'gem3_1','gem3_2','gem3_3','gem3_4',
     1 'scifi3_1','scifi3_2'/

      endif


      lrecl=1024

      call hropen(2,'ntuple','ntuple.track','n',lrecl,istat)
      print *,' ntuple_ini...open done'

      if(mc.eq.3)then
      call hbookn(10,'ntuple_track_data',ncolum_mc1,'ntuple',nprime
     1                  ,chtags_mc1)

      call hbookn(11,'ntuple_event_data',ncolum_mc2,'ntuple',nprime
     1                ,chtags_mc2)

      elseif(mc.eq.5)then
      call hbookn(10,'ntuple_track_data',ncolum_mc3,'ntuple',nprime
     1                  ,chtags_mc3)
 
      call hbookn(11,'ntuple_event_data',ncolum_mc4,'ntuple',nprime
     1                ,chtags_mc4)         

      elseif(mc.eq.2)then
      call hbookn(10,'ntuple_track_data',ncolum1,'ntuple',nprime
     1                 ,chtags1)

      call hbookn(11,'ntuple_event_data',ncolum2,'ntuple',nprime
     1                ,chtags2)

      elseif(mc.eq.4)then
      call hbookn(10,'ntuple_track_data',ncolum3,'ntuple',nprime
     1                 ,chtags3)
 
      call hbookn(11,'ntuple_event_hit_profile_data',ncolum4,
     1     'ntuple',nprime,chtags4)    
      print *,' ntuple_ini....hbookn done'

      elseif(mc.eq.6)then
  
      call hbookn(11,'ntuple_event_hit_profile_data',ncolum4,
     1     'ntuple',nprime,chtags4)    
      print *,' ntuple_ini....hbookn done'
      endif
  

*---------------------  end of input for ntuple

      return
      end


       subroutine ntuple(iv,itrig,nr,
     & jjflags,ntori,listhit1,listhit2,
     & x_ref,x_ref2,par_trk,par_trk2,traject,erpr_trk,
     & chi2_tab_ndf,chi2_tab,mc,druk)

************************************************************************
*                                                                      *
*   routine to create ntuple with parameters of reconstructed tracks   *
*                     no MC quality check !!                           *
*                                                                      *
************************************************************************
#include "PARAM.h"
#include "HITS_IN.h"
#include "PLANES.h"
#include "NTUPLE.h"


      dimension xtuple1(ncolum1),xtuple3(ncolum3)
      dimension xtuple2(ncolum2),xtuple4(ncolum4) 
      dimension chi2_tab_ndf(*),chi2_tab(*)
      dimension jjflags(*),ntori(*)
      dimension par_trk(5,*),erpr_trk(5,*)
      dimension x_ref(*),x_ref2(*),par_trk2(4,*)
      dimension listhit1(ntrack,*),listhit2(ntrack,*)

      dimension traject(mxplan,*)

      logical fri,extra,druk,sub,confir,mib

      data x_refer1/4.69/        ! first sdc behind (sdc2)
      data x_refer2/5.1/         ! second sdc behind (sdc3)
      data x_refer3/5.38/        ! straw1 (two first double layers YZ)
      data x_refer4/5.56/        ! straw2 (two last double layers YZ)
      data x_refer5/1.51/        ! first mic (mic1)
      data x_refer6/2.01/        ! second mic (mic2)
      data x_refer7/2.5/         ! third mic (mic3)
      data x_refer8/4.8/         ! first gem (gem1)
      data x_refer9/5.21/        ! second gem (gem2)
      data x_refer10/5.67/       ! third gem (gem3)



      confir=.false.             ! if full input clusters - change for true (very big !!!!)

      nttc=ntpl_cor+ntpl2_cor
      nttc1=2*nttc+38 

       if(nttc1.gt.ncolum3)then
          print*,'ntuple size too small - too many planes used'
          stop
       endif

       if(mc.eq.6)then
         nr=0
         confir=.true.
       endif

           if(nr.eq.0)then
           if(mc.eq.2)then

               xtuple1(1)=iv
               xtuple1(2)=0
               xtuple1(21)=itrig
               xtuple1(22)=0

               do ika=3,20
               xtuple1(ika)=888
               enddo   

               call hfn(10,xtuple1)
     
*-- now booking xtuple2 (event ntuple)
               xtuple2(1)=iv
               xtuple2(3)=itrig 
               xtuple2(2)=0

               call hfn(11,xtuple2)
  
            elseif(mc.eq.4.or.mc.eq.6)then
               if(mc.eq.4)then

               xtuple3(1)=iv
               xtuple3(2)=0               

               do ika=3,nttc1
               xtuple3(ika)=888
               if(ika.eq.37)xtuple3(37) =itrig
               if(ika.eq.38)xtuple3(38) =0
               enddo

               call hfn(10,xtuple3)
               endif
*-- now booking xtuple4 (event/hit profile ntuple)                 

*--------------------- *
             if(.not.confir)then
               xtuple4(1) =iv
               xtuple4(2)=0
               xtuple4(3) =itrig    
               
             ki=0
             do ki1=1,ntpl_cor
             if(mib(ki1))then
             ki=ki+1
          xtuple4(3+ki1)=jchits_mc(ki)                            ! # of clusters
             else
          xtuple4(3+ki1)=-333                                     ! flag - plane off
             endif
             enddo

        
             ki=0
             do ki2=1,ntpl2_cor
             if(mib(ntpl_cor+ki2))then
             ki=ki+1
          xtuple4(3+ntpl_cor+ki2)=jc2hits_mc(ki)                  ! # of clusters
             else
          xtuple4(3+ntpl_cor+ki2)=-333                            ! flag - plane off
             endif
             enddo     

             call hfn(11,xtuple4)    
*---------------------- *

             elseif(confir)then     
               xtuple4(1) =iv
               xtuple4(2)=0
               xtuple4(3) =itrig  
  
                do jclu=1,jmax                                    ! over clusters

               ki=0
               do ki1=1,ntpl_cor
               if(mib(ki1))then
                ki=ki+1
               if(jclu.le.jchits_mc(ki))then
          xtuple4(3+ki1)=clusters1(numhit_mc(ki,jclu))            ! hit profile
               else
          xtuple4(3+ki1)=-555                                     ! flag - less clusters
               endif
               else
          xtuple4(3+ki1)=-333                                     ! flag - plane off
               endif
               enddo
 
           
              if(ki.ne.ntpl_tmp)then
               print*,'error in logic - ntuple1',ki,ntpl_tmp
               stop
              endif



               ki=0
               do ki2=1,ntpl2_cor
               if(mib(ntpl_cor+ki2))then
               ki=ki+1
               if(jclu.le.jc2hits_mc(ki))then
         xtuple4(3+ntpl_cor+ki2)=clusters2(numht2_mc(ki,jclu))   ! hit profile
               else
         xtuple4(3+ntpl_cor+ki2)=-555                            ! flag - less clusters
               endif
               else    
         xtuple4(3+ntpl_cor+ki2)=-333                            ! flag - plane off
               endif
               enddo
        
              if(ki.ne.ntpl2_tmp)then
               print*,'error in logic - ntuple2',ki,ntpl2_tmp
               stop
              endif


             call hfn(11,xtuple4)

            enddo
            endif
            endif

               return
            endif
*-------------------------------- if not patological one

            do k =1,nr
               
               ppp=par_trk(5,k)                 ! reconstructed fitted momenta
               pabs=abs(ppp)         
               sigp=erpr_trk(5,k)               ! momenta error

               chi2_ndg=chi2_tab_ndf(k)         ! chi2 of the fit
               chi2_line=chi2_tab(k)
 
               chi2=chi2_tab(k)                 ! chi2 global
               ndf=chi2/chi2_ndg  

               proba=prob(chi2,ndf)
*               print*, chi2,ndf,proba
        
*-- reference planes:

*--- SDC2      
      y_ref1=(x_refer1 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref1=(x_refer1 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- SDC3  
      y_ref2=(x_refer2 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref2=(x_refer2 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- st1 
      y_ref3=(x_refer3 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref3=(x_refer3 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- st2  
      y_ref4=(x_refer4 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref4=(x_refer4 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)

*--- mic1  
      y_ref5=(x_refer5 - x_ref(k))*par_trk(3,k) + par_trk(1,k)
      z_ref5=(x_refer5 - x_ref(k))*par_trk(4,k) + par_trk(2,k)
*--- mic2  
      y_ref6=(x_refer6 - x_ref(k))*par_trk(3,k) + par_trk(1,k)
      z_ref6=(x_refer6 - x_ref(k))*par_trk(4,k) + par_trk(2,k)
*--- mic3  
      y_ref7=(x_refer7 - x_ref(k))*par_trk(3,k) + par_trk(1,k)
      z_ref7=(x_refer7 - x_ref(k))*par_trk(4,k) + par_trk(2,k) 

*--- gem1  
      y_ref8=(x_refer8 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref8=(x_refer8 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- gem2  
      y_ref9=(x_refer9 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref9=(x_refer9 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)  
*--- gem3  
      y_ref10=(x_refer10 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref10=(x_refer10 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k) 

*--- end of reference planes

  
            if(jjflags(k).lt.-100000.and.jjflags(k).gt.-200000)then               ! no extrapolation
               extra=.false.
               fri=.false.
               jjfl=-jjflags(k)-100000
               lihit3=0
               lihit1=mod(jjfl,1000)                   
               lihit2=(jjfl-lihit1)/1000
            elseif(jjflags(k).lt.-200000)then         ! fringe field with bridging 
               extra=.false.
               fri=.true.
               jjfl=-jjflags(k)-200000
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

             
*------------------- mc=2 
               if(mc.eq.2) then

               xtuple1(1)=iv
               xtuple1(2)=k
               xtuple1(3)=ppp
               xtuple1(4)=sigp
               xtuple1(5)=chi2_ndg
               xtuple1(6)=proba
               xtuple1(7)=lihit1
               xtuple1(8)=lihit2
               xtuple1(9)=lihit3
               if(fri)xtuple1(10)=1                   ! ff
               if(.not.fri)xtuple1(10)=0
               if(extra)xtuple1(11)=1                 ! extrapolation ok
               if(.not.extra)xtuple1(11)=0            ! extrapolation fault

               xtuple1(12)=x_ref(k)
               xtuple1(13)=par_trk(1,k)
               xtuple1(14)=par_trk(2,k)
               xtuple1(15)=par_trk(3,k)
               xtuple1(16)=par_trk(4,k)

               xtuple1(17)=y_ref1
               xtuple1(18)=z_ref1
               xtuple1(19)=y_ref2
               xtuple1(20)=z_ref2
               xtuple1(21)=itrig                      ! trigger mask

               xtuple1(22)=ntori(k)                   ! which iteration 
   
               call hfn(10,xtuple1)
* -- now booking xtuple2 (event ntuple)
 
              if(k.eq.1)then                          ! only ones!

               xtuple2(2) =nr
               xtuple2(1)=iv
               xtuple2(3)=itrig
 
               call hfn(11,xtuple2)     
              endif
  
*------------------- mc=4
               elseif(mc.eq.4) then
 
               xtuple3(1)=iv
               xtuple3(2)=k
               xtuple3(3)=ppp
               xtuple3(4)=sigp
               xtuple3(5)=chi2_ndg
               xtuple3(6)=proba
               xtuple3(7)=lihit1
               xtuple3(8)=lihit2
               xtuple3(9)=lihit3
               if(fri)xtuple3(10)=1                    ! ff
               if(.not.fri)xtuple3(10)=0
               if(extra)xtuple3(11)=1                 ! extrapolation ok
               if(.not.extra)xtuple3(11)=0            ! extrapolation fault
 
               xtuple3(12)=x_ref(k)
               xtuple3(13)=par_trk(1,k)
               xtuple3(14)=par_trk(2,k)
               xtuple3(15)=par_trk(3,k)
               xtuple3(16)=par_trk(4,k)

               xtuple3(17)=y_ref1
               xtuple3(18)=z_ref1
               xtuple3(19)=y_ref2
               xtuple3(20)=z_ref2
               xtuple3(21)=y_ref3
               xtuple3(22)=z_ref3
               xtuple3(23)=y_ref4
               xtuple3(24)=z_ref4
               xtuple3(25)=y_ref5
               xtuple3(26)=z_ref5
               xtuple3(27)=y_ref6
               xtuple3(28)=z_ref6
               xtuple3(29)=y_ref7
               xtuple3(30)=z_ref7
               xtuple3(31)=y_ref8
               xtuple3(32)=z_ref8
               xtuple3(33)=y_ref9
               xtuple3(34)=z_ref9
               xtuple3(35)=y_ref10
               xtuple3(36)=z_ref10
               xtuple3(37)=itrig                      ! trigger mask
               xtuple3(38)=ntori(k)                   ! which iteration 

                  in=ntori(k)  

                  kk=0
                  kk1=0 
                  kkk=0

               do i1 =1,ntpl_cor
           
               if(mib(i1))then
               kkk=kkk+1
                clus=8888
                traje=8888

      if(.not.fri.and.sub.and.in.le.mp1.and.
     &         (ityp1(ntpl_cor+kkk).eq.111.                             ! set 1
     &    or.ityp1(ntpl_cor+kkk).eq.222))then
                  kk=kk+1
            traje= traject(kk,k)
       if(listhit1(k,kk).ne.8888)then
           clus=clusters1(listhit1(k,kk)) 
       else
           clus=-444
       endif
      elseif(.not.fri.and.sub.and.in.gt.mp1.and.
     &         in.le.mp1+mp2.and.                                       ! set 2
     &        ityp1(ntpl_cor+kkk).eq.111)then
                  kk1=kk1+1
            traje= traject(kk1,k)
       if(listhit1(k,kk1).ne.8888)then
           clus=clusters1(listhit1(k,kk1)) 
       else
           clus=-444
       endif

      elseif((sub.and.in.gt.mp1+mp2).or..not.sub.or.
     &       (sub.and.fri))then
            traje= traject(kkk,k)
       if(listhit1(k,kkk).ne.8888)then
           clus=clusters1(listhit1(k,kkk))
       else
           clus=-444 
       endif
      endif

          xtuple3(38+i1+nttc) = traje 
          xtuple3(38+i1)=clus
    
                 else                                                    ! plane off
          xtuple3(38+i1)=-333
          xtuple3(38+i1+nttc)=-333

                 endif                  
               enddo



                 kk=0
                 kk1=0 
                 kkk=0

               do i2=1,ntpl2_cor

               if(mib(ntpl_cor+i1))then
                kkk=kkk+1
                 clus2=8888
                 traje2=8888

      if(.not.fri.and.sub.and.in.le.mp1.and.
     &             (ityp2(ntpl2_cor+i2).eq.222.or.                       ! set 1
     &                ityp2(ntpl2_cor+i2).eq.555))then 
                kk = kk+1

            traje2= traject(kk+ntpl_tmp,k)
       if(listhit2(k,kk).ne.8888)then
           clus2=clusters2(listhit2(k,kk))
       else
           clus2=-444
       endif  
      elseif(.not.fri.and.sub.and.in.gt.mp1.and.in.le.mp1+mp2.and.
     &        ityp2(ntpl2_cor+i2).eq.444)then
                kk1=kk1+1 

             traje2= traject(kk1+ntpl_tmp,k)
       if(listhit2(k,kk1).ne.8888)then
           clus2=clusters2(listhit2(k,kk1))
       else
           clus2=-444
       endif

      elseif((sub.and.in.gt.mp1+mp2).or..not.sub.or.
     &       (sub.and.fri))then

             traje2= traject(kkk+ntpl_tmp,k)
       if(listhit2(k,kkk).ne.8888)then
           clus2=clusters2(listhit2(k,kkk))
       else
           clus2=-444
       endif
      endif
      
          xtuple3(38+ntpl_cor+i2+nttc) = traje2   
          xtuple3(38+ntpl_cor+i2)=clus2
   
                 else                                                    ! plane off
          xtuple3(38+ntpl_cor+i2)=-333
          xtuple3(38+ntpl_cor+i2+nttc)=-333

                 endif                 
      
               enddo 
 
               call hfn(10,xtuple3) 

              if(k.eq.1)then                                             ! only ones!!
          
* -- now booking xtuple4 (event/cluster profile  ntuple)
             if(.not.confir)then
               xtuple4(1) =iv
               xtuple4(2)=0
               xtuple4(3) =itrig    
               
       
             do ki1=1,ntpl_cor
             if(mib(ki1))then
          xtuple4(3+ki1)=jchits_mc(ki1)                           ! # of clusters
             else
          xtuple4(3+ki1)=-333                                     ! flag - plane off
             endif
             enddo


             do ki2=1,ntpl2_cor
             if(mib(ntpl_cor+ki2))then
          xtuple4(3+ntpl_cor+ki2)=jc2hits_mc(ki2)                 ! # of clusters
             else
          xtuple4(3+ntpl_cor+ki2)=-333                            ! flag - plane off
             endif
             enddo     

             elseif(confir)then     
               xtuple4(1) =iv
               xtuple4(2)=0
               xtuple4(3) =itrig  
  
                do jclu=1,jmax                                     ! over clusters

               do ki1=1,ntpl_cor
               if(mib(ki1))then
               if(jclu.le.jchits_mc(ki1))then
          xtuple4(3+ki1)=clusters1(numhit_mc(ki1,jclu))            ! hit profile
               else
          xtuple4(3+ki1)=-555                                      ! flag
               endif
               else
          xtuple4(3+ki1)=-333                                      ! flag - plane off
               endif
               enddo
 

               do ki2=1,ntpl2_cor
               if(mib(ntpl_cor+ki2))then
               if(jclu.le.jc2hits_mc(ki2))then
          xtuple4(3+ntpl_cor+ki2)=clusters2(numht2_mc(ki2,jclu))  ! hit profile
               else
          xtuple4(3+ntpl_cor+ki2)=-555                            ! flag
               endif
               else
          xtuple4(3+ntpl_cor+ki2)=-333                            ! flag - plane off
               endif
               enddo
        
             call hfn(11,xtuple4)

            enddo
            endif
            endif

 
        endif
*------ end of mc=2/4   

       
            if(druk)then
              bl=  sigp/abs(ppp)         

               if(fri)then
      print'("event, trig., p., chi2/ndf, prob. ,nh.1/2  , it.: ",
     & i10,i5,3f8.2,1x,e7.2,2i3,3x,i3," ff rec. ")',
*     & iv,itrig,ppp,chi2_ndg,chi2_line,proba,lihit1,in
     & iv,itrig,ppp,bl,chi2_ndg,proba,lihit1,lihit2,in

               else
               if(extra)then
      print'("event, trig., p., chi2/ndf, prob., nh.1/2/3, it.: ",
     & i10,i5,3f8.2,1x,e7.2,4i3," ext.SM2")',
*     &  iv,itrig,ppp,chi2_ndg,chi2_line,proba,lihit1,lihit2,in
     &  iv,itrig,ppp,bl,chi2_ndg,proba,lihit1,lihit2,lihit3,in

               elseif(.not.extra)then
      print'("event, trig., p., chi2/ndf, prob., nh.1/2  , it.: ",
     & i10,i5,3f8.2,1x,e7.2,2i3,3x,i3)',
*     &  iv,itrig,ppp,chi2_ndg,chi2_line,proba,lihit1,lihit2,in 
     &  iv,itrig,ppp,bl,chi2_ndg,proba,lihit1,lihit2,in 

               endif
               endif
            endif

             enddo
             
 
             return
             end


       subroutine ntuple_mc(iv,itrig,mlassc,no,jjflags,ntori,
     & listhit1,listhit2,x_ref,x_ref2,par_trk,
     & par_trk2,traject,erpr_trk,chi2_tab_ndf,chi2_tab,mc,druk)  

************************************************************************
*                                                                      *
*         routine to create ntuple with parameters of tracks           *
*                     full MC quality check here                       *
*                                                                      *
************************************************************************
#include "PARAM.h"
#include "MC.h"
#include "HITS_IN.h"
#include "PLANES.h"
#include "NTUPLE.h"
   
    
      dimension xtuple1(ncolum_mc1),xtuple2(ncolum_mc2)
      dimension xtuple3(ncolum_mc3),xtuple4(ncolum_mc4)
      
      dimension chi2_tab_ndf(*),l_hitg2(ntrack)
      dimension no(100),li_glue(100)
      dimension jjflags(*),ntori(*)
      dimension par_trk(5,*),erpr_trk(5,*)
      dimension traject(mxplan,*)

      logical fri,sub,extra,muon,pileup,accep,druk,confir

      data x_refer1/4.69/        ! first sdc behind (sdc2)
      data x_refer2/5.1/         ! second sdc behind (sdc3)
      data x_refer3/5.38/        ! straw1 (two first double layers YZ)
      data x_refer4/5.56/        ! straw2 (two last double layers YZ)
      data x_refer5/1.51/        ! first mic (mic1)
      data x_refer6/2.01/        ! second mic (mic2)
      data x_refer7/2.5/         ! third mic (mic3)
      data x_refer8/4.8/         ! first gem (gem1)
      data x_refer9/5.21/        ! second gem (gem2)
      data x_refer10/5.67/       ! third gem (gem3)
  
      confir=.false.             ! if full input clusters - change for true (very big !!!!)

      nttc=ntpl_cor+ntpl2_cor
      nttc1=2*nttc+56

       if(nttc1.gt.ncolum_mc3)then
          print*,'ntuple size too small - too many planes used'
          stop
       endif
  

*--- counters reset
        li_miss=0
        nonac=0
        nonac1=0
        nonac2=0
        nonac3=0
        nonac4=0
        nonac5=0
        nonac6=0
        nonac7=0
        nonac8=0
        nonac9=0
        nonac10=0
        nonac11=0
        nonac12=0
        nonac13=0
        nonac14=0
        nonac15=0
        nonac16=0
        nonac17=0
        nonac18=0
        nonac19=0

        itr_count1=0
        itr_count2=0 
        itr_count3=0
        itr_count4=0
        itr_count5=0 
        itr_count6=0
        itr_count7=0

        itr_count1_new=0
        itr_count2_new=0 
        itr_count3_new=0
        itr_count4_new=0
        itr_count5_new=0 
        itr_count6_new=0
        itr_count7_new=0

        n_acc=0
        n_acc8=0

        do jo=1,100
          li_glue(jo)=0
        enddo

        if(ntr.eq.0.and.jjflag(1,1).ne.-10000)then
           print*, 'error in logic; 0 accepted tracks'
           stop
        endif


*---- loop over generated tracks

        if(ntr.gt.0)then
        do 1 j=1,ntr                           ! over gen.

                    
          
          if(label(j).eq.-1) go to 1           ! no accepted track
          if(label(j).eq.-100) nonac13=nonac13+1
          if(label(j).eq.-200) nonac14=nonac14+1
          if(label(j).eq.-300) nonac15=nonac15+1

          accep=.true.

          if(ped_tab(j).gt.5000.)then
              nonac19=nonac19+1
              accep=.false.
              ped_tab(j)=ped_tab(j)-10000.
          endif
            
           if(accep) n_acc=n_acc+1

           ngl=0
           licz_glue=0


            lhit1g=mod(jj_hit(1,j),1000)             ! generated hits
            lhit2g=mod(jj_hit(1,j),100000)/1000
            lhitmm=jj_hit(1,j)/100000
            lhitd=mod(jj_hit(2,j),1000)
            lhitsf=mod(jj_hit(2,j),100000)/1000
            lhitgem=jj_hit(2,j)/100000
            lhitsup=jj_hit(3,j)

            if(lhit2g.lt.8) n_acc8=n_acc8+1          ! tracks with < 8 hits after but acc.
*        
* --- looking for reconstructed (if exists) track
*---- now reconstructed tracks
*
          do k =1,nr                                 ! over reconstructed

            ngl=jjflag(k,2)                          ! no of gen. assoc. track with rec.

               chi2_ndf=chi2_tab_ndf(k)              ! chi2 fit
               chi2=chi2_tab(k)                      ! chi2 global
               ndf=chi2/chi2_ndf
               proba=prob(chi2,ndf)

            if(ngl.eq.j) then                        ! this one!
               kas=k
*-- reference planes:
*-- reference planes:

*--- SDC2      
      y_ref1=(x_refer1 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref1=(x_refer1 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- SDC2  
      y_ref2=(x_refer2 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref2=(x_refer2 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- st1 
      y_ref3=(x_refer3 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref3=(x_refer3 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- st2  
      y_ref4=(x_refer4 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref4=(x_refer4 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)

*--- mic1  
      y_ref5=(x_refer5 - x_ref(k))*par_trk(3,k) + par_trk(1,k)
      z_ref5=(x_refer5 - x_ref(k))*par_trk(4,k) + par_trk(2,k)
*--- mic2  
      y_ref6=(x_refer6 - x_ref(k))*par_trk(3,k) + par_trk(1,k)
      z_ref6=(x_refer6 - x_ref(k))*par_trk(4,k) + par_trk(2,k)
*--- mic3  
      y_ref7=(x_refer7 - x_ref(k))*par_trk(3,k) + par_trk(1,k)
      z_ref7=(x_refer7 - x_ref(k))*par_trk(4,k) + par_trk(2,k) 

*--- gem1  
      y_ref8=(x_refer8 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref8=(x_refer8 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)
*--- gem2  
      y_ref9=(x_refer9 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref9=(x_refer9 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k)  
*--- gem3  
      y_ref10=(x_refer10 - x_ref2(k))*par_trk2(3,k) + par_trk2(1,k)
      z_ref10=(x_refer10 - x_ref2(k))*par_trk2(4,k) + par_trk2(2,k) 

*--- end of reference planes       


*---- first track is muon from the beam in coral, second - scattered muon
            muon=.false.
            if(ngl.eq.2)muon = .true.                ! muon

            lflag=jjflag(k,1)                        ! flag if good associate

*---- clusters for reconstructed

            if(jjflags(k).lt.-100000)then            ! no extrapolation
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
            else                                     ! extrapolation ok
               fri=.false.
               extra=.true.
               lihit1=mod(jjflags(k),1000)         
               lihit2=mod(jjflags(k),100000)/1000
               lihit3=jjflags(k)/100000
            endif
        

             l_hitg2(k)=800                          ! nb of hits after SM1 for ffr track
             if(fri)l_hitg2(k)=lhit2g

             nphit=jjflag(k,5)                       ! percent of hits


             if((lflag.eq.0.or.lflag.eq.100).and.licz_glue.eq.0) then

               xx = par_trk(1,k)
               yy = par_trk(2,k)
               xp = par_trk(3,k)
               yp = par_trk(4,k)

               ppp=abs(par_trk(5,k))              ! reconstructed fitted momenta

               sigp=erpr_trk(5,k)                 ! momenta error
               sigx=erpr_trk(1,k)
               sigy=erpr_trk(2,k)
               sigxp=erpr_trk(3,k)
               sigyp=erpr_trk(4,k)

             endif
               licz_glue=licz_glue+1
             endif

          enddo                                   ! end of loop over rec. tracks

*
*---------- not reconstructed
*

               if(licz_glue.eq.0)then             ! not reconstructed

               ppp= -999                       
               sigp=-999
               chi2_ndg=-999
               lihit1=-999
               lihit2=-999
               lihit3=-999
               nphit=-999
               extra=.false.
               fri=.false.
               sigp=-999                 ! momenta error
               sigx=-999
               sigy=-999
               sigxp=-999
               sigyp=-999
       
             y_ref1=-999
             z_ref1=-999  
             y_ref2=-999
             z_ref2=-999
    

               elseif(licz_glue.gt.1)then         ! double (multiple) reconstructed

                li_glue(licz_glue-1)=li_glue(licz_glue-1)+1

                li_miss=li_miss+licz_glue-1

               if(lflag.ne.0.and.lflag.ne.100)then         ! no proper number of hits
                print *,'no proper number of hits - error!'
                stop
               endif

               endif                              ! end of non reconstructed

*
*------------- end of not reconstructed --------------------------------------*
*------------- good reconstructed track --------------------------------------*
*

       p_gen =ped_tab(j)                          ! momentum of gen.track
       p_gen2 =ped_tab2(j,1)                      ! momentum of gen.track in first plane
       xp_1pl = ped_tab2(j,2)                     ! angle in first plane
       yp_1pl = ped_tab2(j,3)                     ! angle in first plane
       x_gen2 =xped_tab2(j,1)                     ! x (z) position of first plane of momentum
       x_gen3 =xped_tab2(j,2)                     ! y (x) position of ...
       x_gen4 =xped_tab2(j,3)                     ! z (y) position of ...


       pileup=.false.
       if(p_gen.lt.0)then
          p_gen=-p_gen
          pileup=.true.
       endif

       if(ppp.eq.-999) then
          delp_p=-999
       else
       delp_p=100.*(p_gen-ppp)/p_gen              ! momentum resolution
       endif

*------------------ ntuples mc=3 (medium)
          if(mc.eq.3) then

               xtuple1(1)=iv
               xtuple1(2)=j
               xtuple1(3)=lhitgem
               xtuple1(4)=lhitsf
               xtuple1(5)=lhitsup                  
               xtuple1(6)=lhitd
               xtuple1(7)=lhitmm
               xtuple1(8)=nphit
               xtuple1(9)=lhit1g
               xtuple1(10)=lhit2g
               xtuple1(11)=lihit1
               xtuple1(12)=lihit2
               xtuple1(13)=lihit3
               xtuple1(14)= sigp
               xtuple1(15)=ppp
               xtuple1(16)=p_gen
               xtuple1(17)=p_gen2
               xtuple1(18)=x_gen2
               xtuple1(19)=delp_p
               xtuple1(20)=chi2_ndg
               xtuple1(21)=proba
               if(fri)xtuple1(22)=1                ! fringe field rec. flag
               if(.not.fri)xtuple1(22)=0           ! normal reconstruction
               if(extra)xtuple1(23)=1              ! extrapolation ok
               if(.not.extra)xtuple1(23)=0         ! extrapolation fault
               if(pileup)xtuple1(24)=1             ! pileup track
               if(.not.pileup)xtuple1(24)=-1       ! good generated track
               if(accep)xtuple1(25)=0              ! accepted
               if(.not.accep)xtuple1(25)=-1        ! no in acceptance (target shot)
                                    
               xtuple1(26) = xx
               xtuple1(27) = yy
               xtuple1(28) = xp
               xtuple1(29) = yp
               xtuple1(30) = x_gen4
               xtuple1(31) = x_gen3
               xtuple1(32) = xp_1pl
               xtuple1(33) = yp_1pl
               xtuple1(34) = sigx  
               xtuple1(35) = sigy   
               xtuple1(36) = sigxp  
               xtuple1(37) = sigyp   



               call hfn(10,xtuple1)

*--------------- ntuple - mc=5 - very large 
          elseif(mc.eq.5) then

 
               xtuple3(1)=iv
               xtuple3(2)=j
               xtuple3(3)=lhitsup
               xtuple3(4)=nphit
               xtuple3(5)=lhit1g
               xtuple3(6)=lhit2g
               xtuple3(7)=lihit1
               xtuple3(8)=lihit2
               xtuple3(9)=lihit3
               xtuple3(10)= sigp
               xtuple3(11)=ppp
               xtuple3(12)=p_gen
               xtuple3(13)=p_gen2
               xtuple3(14)=x_gen2
               xtuple3(15)=delp_p
               xtuple3(16)=chi2_ndg
               xtuple3(17)=proba
               if(fri)xtuple3(18)=1                ! fringe field rec. flag
               if(.not.fri)xtuple3(18)=0           ! normal reconstruction
               if(extra)xtuple3(19)=1              ! extrapolation ok
               if(.not.extra)xtuple3(19)=0         ! extrapolation fault  
               if(pileup)xtuple3(20)=1             ! pileup track
               if(.not.pileup)xtuple3(20)=-1       ! good generated track
               if(accep)xtuple3(21)=0              ! accepted
               if(.not.accep)xtuple3(21)=-1        ! no in acceptance (target shot)
 
             
               xtuple3(22) = x_gen4
               xtuple3(23) = x_gen3
               xtuple3(24) = xp_1pl
               xtuple3(25) = yp_1pl
               xtuple3(26) = sigx
               xtuple3(27) = sigy
               xtuple3(28) = sigxp
               xtuple3(29) = sigyp

               xtuple1(30)=x_ref(kas)
               xtuple3(31)=par_trk(1,kas)
               xtuple3(32)=par_trk(2,kas)
               xtuple3(33)=par_trk(3,kas)
               xtuple3(34)=par_trk(4,kas)

               xtuple3(35)=y_ref1
               xtuple3(36)=z_ref1
               xtuple3(37)=y_ref2
               xtuple3(38)=z_ref2
               xtuple3(39)=y_ref3
               xtuple3(40)=z_ref3
               xtuple3(41)=y_ref4
               xtuple3(42)=z_ref4
               xtuple3(43)=y_ref5
               xtuple3(44)=z_ref5
               xtuple3(45)=y_ref6
               xtuple3(46)=z_ref6
               xtuple3(47)=y_ref7
               xtuple3(48)=z_ref7
               xtuple3(49)=y_ref8
               xtuple3(50)=z_ref8
               xtuple3(51)=y_ref9
               xtuple3(52)=z_ref9
               xtuple3(53)=y_ref10
               xtuple3(54)=z_ref10

               xtuple3(55)=itrig                      ! trigger mask
               xtuple3(56)=ntori(kas)                 ! which iteration 

                  in=ntori(kas)  

*  -------------------------- not finished yet !!!!!!

                  kk=0 
               do i1 =1,ntpl_tmp

                do jj=1, jchits_mc(i1)
                if(l_hit1(1,i1,jj).ne.0)then
                 do nj=1,l_hit1(1,i1,jj)
           if(l_hit1(nj+1,i1,jj).eq.j)then
          xtuple3(40+i1)=clusters1(numhit_mc(i1,jj))+10000 
                go to 10
                else
                xtuple3(40+i1)=18888
                endif
                enddo
               else
           xtuple3(40+i1)=28888
                endif
                enddo
 10              continue

       if(sub.and.in.le.mp1+mp2.and.ityp1(ntpl_cor+i1).eq.111)then
                  kk=kk+1
                 if(listhit1(kas,kk).ne.8888)then   
          xtuple3(40+i1)=clusters1(listhit1(kas,kk))
                 else
          xtuple3(40+i1)=8888
                 endif
       elseif((sub.and.in.gt.mp1+mp2).or..not.sub)then
*             print*,kas,i1,kk,listhit1(kas,i1)  
                 if(listhit1(kas,i1).ne.8888)then
          xtuple3(40+i1)=clusters1(listhit1(kas,i1))
                 else
          xtuple3(40+i1)=8888
                 endif                  
       endif
               enddo

                 kk=0
                 kk1=0 
               do i2=1,ntpl2_tmp

               do jj=1, jc2hits_mc(i2)
                if(l_hit2(1,i2,jj).ne.0)then
                 do nj=1,l_hit2(1,i2,jj)
           if(l_hit2(nj+1,i2,jj).eq.j)then
          xtuple3(40+ntpl_tmp+i2)=clusters2(numht2_mc(i2,jj))+10000 
                go to 20
                else
           xtuple3(40+ntpl2_tmp+i2)=18888  
                endif
                enddo
               else
           xtuple3(40+ntpl2_tmp+i2)=28888
                endif
                enddo
 20              continue

      if(sub.and.in.le.mp1.and.(ityp2(ntpl2_cor+i2).eq.222.or.
     &                ityp2(ntpl2_cor+i2).eq.555))then 
                kk = kk+1
               if(listhit2(kas,kk).ne.8888)then
          xtuple3(40+ntpl_tmp+i2)=clusters2(listhit2(kas,kk))
               else
          xtuple3(40+ntpl_tmp+i2)=8888   
               endif
      elseif(sub.and.in.gt.mp1.and.in.le.mp1+mp2.and.
     &        ityp2(ntpl2_cor+i2).eq.444)then
                kk1=kk1+1     
              if(listhit2(kas,kk1).ne.8888)then   
          xtuple3(40+ntpl_tmp+i2)=clusters2(listhit2(kas,kk1))
               else
          xtuple3(40+ntpl_tmp+i2)=8888
               endif
      elseif((sub.and.in.gt.mp1+mp2).or..not.sub)then
               if(listhit2(kas,i2).ne.8888)then
          xtuple3(40+ntpl_tmp+i2)=clusters2(listhit2(kas,i2))
               else
          xtuple3(40+ntpl_tmp+i2)=8888
               endif          
      endif
               enddo 

               liro=li_pla1+li_pla2-ntpl_tmp-ntpl2_tmp


               if(liro.lt.0)then
            liro=0
            print*,'too many planes used - change size of ntuple'
               endif 
     
               do i12=1,liro
             xtuple3(40+ntpl_tmp+ntpl2_tmp+i12)= 9999              ! flag-not used planes
               enddo
 
               call hfn(10,xtuple3)         
           endif

  1     continue                                  ! end of gen. tracks loop

*
*---------- end of generated tracks
*
           endif


           do i=1,50
           no(i)= li_glue(i) + no(i)
           enddo

           no(99)= li_miss + no(99)


* -- now loop over reconstructed tracks to find ghosts

          if(druk)
     &print*,'no of rec/acc/acc(<8 behind) tracks: ',nr,n_acc,n_acc8

              if(nr.gt.0)then
              do j=1,nr                           ! once again over rec. tracks
              
           lflag=jjflag(j,1)
           ngl=jjflag(j,2)

*                print*,j,lflag,ngl

               if(ntr.eq.0)lflag=-1               ! only ghosts possible

               ppp=abs(par_trk(5,j))              ! reconstructed fitted momenta
               sigp=erpr_trk(5,j)                 ! momenta error

               chi2_ndf=chi2_tab_ndf(j)              ! chi2 fit
               chi2=chi2_tab(j)                      ! chi2 global
               ndf=chi2/chi2_ndf
               proba=prob(chi2,ndf)

*-- reference planes:

      y_straw=(x_straw-x_ref2(j))*par_trk2(3,j)+par_trk2(1,j)
      z_straw=(x_straw-x_ref2(j))*par_trk2(4,j)+par_trk2(2,j)

      y_sdc=(x_sdc - x_ref2(j))*par_trk2(3,j) + par_trk2(1,j)
      z_sdc=(x_sdc - x_ref2(j))*par_trk2(4,j) + par_trk2(2,k)

*--- end of reference planes

               sigp=erpr_trk(5,j)                        ! momenta error
               sigx=erpr_trk(1,j)
               sigy=erpr_trk(2,j)
               sigxp=erpr_trk(3,j)
               sigyp=erpr_trk(4,j)

*---- clusters for reconstructed
            if(jjflags(j).lt.-100000)then                ! no extrapolation
               extra=.false.
               fri=.false.
               jjfl=-jjflags(j)-100000
               lihit3=0
               lihit1=mod(jjfl,1000)                   
               lihit2=(jjfl-lihit1)/1000
            elseif(jjflags(j).lt.0.and.jjflags(j).gt.-100000)then  ! fringe field reconstruction
               fri=.true.
               extra=.false.
               lihit1=-jjflags(j)
               lihit2=0
               lihit3=0
            else                                         ! extrapolation ok
               fri=.false.
               extra=.true.
               lihit1=mod(jjflags(j),1000)         
               lihit2=mod(jjflags(j),100000)/1000
               lihit3=jjflags(j)/100000
            endif

           lflag12=jjflag(j,3)

*
*---------------- ghosts ----------------------------------*
*

                   if(lflag.le.-1)then                    ! reconstructed but ghosts

             if(ngl.ne.0)then
                 print *,'error in ghost definition!'
                 stop
             endif

           if(lflag.eq.-1.or.lflag.eq.-2)then
           if(fri)then
             print*,'error in logic 0 - fringe/no fringe rec.'
           endif

            in=ntori(j)
           
            if(extra)then
      if(druk)print '("event,ghost, p_rec, chi2, it: ",
     &i5,f8.2,i3,f8.2,i2,
     &            " ex. to SM2")',iv,ppp,ngl,chi2_ndg,in
            nonac17=nonac17+1

            elseif(.not.extra)then
      if(druk)
     & print '("event,ghost, p_rec, chi2, it: ",i5,f8.2,i3,f8.2)',
     &          iv,ppp,ngl,chi2_ndg,in
            endif
     

         elseif(lflag.eq.-100)then                  !  (fringe field)
          if(.not.fri)then
            print*,'error in logic 1 - fringe/no fringe rec.'
          endif

         if(druk)
     & print '("event,ghost-ffr., p_rec, chi2: ",i5,f8.2,i3,f8.2)',
     &          iv,ppp,ngl,chi2_ndg
          endif


            if(mc.eq.3) then
               xtuple1(1)=iv
               xtuple1(2)=-888
               xtuple1(3)=-888
               xtuple1(4)=-888
               xtuple1(5)=-888                   
               xtuple1(6)=-888
               xtuple1(7)=-888
               xtuple1(8)=-888
               xtuple1(9)=-888
               xtuple1(10)=-888
               xtuple1(11)=lihit1
               xtuple1(12)=lihit2
               xtuple1(13)=lihit3
               xtuple1(14)=sigp
               xtuple1(15)=ppp
               xtuple1(16)=-888                    ! ghost
               xtuple1(17)=-888
               xtuple1(18)=-888
               xtuple1(19)=-888
               xtuple1(20)=-888

               xtuple1(21)=-888
               xtuple1(22)=chi2_ndg
               xtuple1(23)=proba
               if(fri)xtuple1(24)=1
               if(.not.fri)xtuple1(24)=0
               if(extra)xtuple1(25)=1
               if(.not.extra)xtuple1(25)=0
               xtuple1(26)= xx
               xtuple1(27)= yy
               xtuple1(28) = xp
               xtuple1(29) = yp
               xtuple1(30)=-888
               xtuple1(31)=-888 
               xtuple1(32) =-888 
               xtuple1(33) =-888 
               xtuple1(34)= sigx
               xtuple1(35)= sigy
               xtuple1(36) = sigxp
               xtuple1(37) = sigyp

               call hfn(10,xtuple1)

            elseif(mc.eq.5) then

               xtuple3(1)=iv
               xtuple3(2)=-888
               xtuple3(3)=-888
               xtuple3(4)=-888
               xtuple3(5)=-888
               xtuple3(6)=-888
               xtuple3(7)=lihit1
               xtuple3(8)=lihit2
               xtuple3(9)=lihit3
               xtuple3(10)=sigp
               xtuple3(11)=ppp
               xtuple3(12)=-888                    ! ghost
               xtuple3(13)=-888
               xtuple3(14)=-888
               xtuple3(15)=-888
               xtuple3(16)=chi2_ndg
               xtuple3(17)=proba
               if(fri)xtuple3(18)=1                ! fringe field rec. flag
               if(.not.fri)xtuple3(18)=0           ! normal reconstruction
               if(extra)xtuple3(19)=1              ! extrapolation ok
               if(.not.extra)xtuple3(19)=0         ! extrapolation fault
               xtuple3(20)=-888                    ! pileup track
               xtuple3(21)=-888                    ! accepted

               xtuple3(22) =-888 
               xtuple3(23) =-888
               xtuple3(24) =-888 
               xtuple3(25) =-888
               xtuple3(26) = sigx
               xtuple3(27) = sigy
               xtuple3(28) = sigxp
               xtuple3(29) = sigyp

               xtuple1(30)=x_ref(j)
               xtuple3(31)=par_trk(1,j)
               xtuple3(32)=par_trk(2,j)
               xtuple3(33)=par_trk(3,j)
               xtuple3(34)=par_trk(4,j)

               xtuple3(35)=y_ref1
               xtuple3(36)=z_ref1
               xtuple3(37)=y_ref2
               xtuple3(38)=z_ref2
               xtuple3(39)=y_ref3
               xtuple3(40)=z_ref3
               xtuple3(41)=y_ref4
               xtuple3(42)=z_ref4
               xtuple3(43)=y_ref5
               xtuple3(44)=z_ref5
               xtuple3(45)=y_ref6
               xtuple3(46)=z_ref6
               xtuple3(47)=y_ref7
               xtuple3(48)=z_ref7
               xtuple3(49)=y_ref8
               xtuple3(50)=z_ref8
               xtuple3(51)=y_ref9
               xtuple3(52)=z_ref9
               xtuple3(53)=y_ref10
               xtuple3(54)=z_ref10

               xtuple3(55)=itrig                      ! trigger mask
               xtuple3(56)=ntori(j)                   ! which iteration

                  in=ntori(j) 
                  kk=0
               do i1 =1,ntpl_tmp

       if(sub.and.in.le.mp1+mp2.and.ityp1(ntpl_cor+i1).eq.111)then
                  kk=kk+1
                 if(listhit1(j,kk).ne.8888)then
          xtuple3(40+i1)=clusters1(listhit1(j,kk))
                 else
          xtuple3(40+i1)=8888
                 endif
       elseif((sub.and.in.gt.mp1+mp2).or..not.sub)then
*             print*,j,i1,kk,listhit1(j,i1)
                 if(listhit1(j,i1).ne.8888)then
          xtuple3(40+i1)=clusters1(listhit1(j,i1))
                 else
          xtuple3(40+i1)=8888
                 endif
       endif
               enddo

                 kk=0
                 kk1=0
               do i2=1,ntpl2_tmp

      if(sub.and.in.le.mp1.and.(ityp2(ntpl2_cor+i2).eq.222.or.
     &                ityp2(ntpl2_cor+i2).eq.555))then
                kk = kk+1
               if(listhit2(j,kk).ne.8888)then
          xtuple3(40+ntpl_tmp+i2)=clusters2(listhit2(j,kk))
               else
          xtuple3(40+ntpl_tmp+i2)=8888
               endif
      elseif(sub.and.in.gt.mp1.and.in.le.mp1+mp2.and.
     &        ityp2(ntpl2_cor+i2).eq.444)then
                kk1=kk1+1
              if(listhit2(j,kk1).ne.8888)then
          xtuple3(40+ntpl_tmp+i2)=clusters2(listhit2(j,kk1))
               else
          xtuple3(40+ntpl_tmp+i2)=8888
               endif
      elseif((sub.and.in.gt.mp1+mp2).or..not.sub)then
               if(listhit2(j,i2).ne.8888)then
          xtuple3(40+ntpl_tmp+i2)=clusters2(listhit2(j,i2))
               else
          xtuple3(40+ntpl_tmp+i2)=8888
               endif
      endif
               enddo

               liro=li_pla1+li_pla2-ntpl_tmp-ntpl2_tmp


               if(liro.lt.0)then
            liro=0
            print*,'too many planes used - change size of ntuple'
               endif

               do i12=1,liro
             xtuple3(40+ntpl_tmp+ntpl2_tmp+i12)= 9999              ! flag-not used planes
               enddo

               call hfn(10,xtuple3)                
            endif


               if(lflag12.eq.-1)nonac1=nonac1+1
               if(lflag12.eq.-2)nonac2=nonac2+1

               nonac=nonac+1          ! # of non accepted tracks due to hits

               if(ppp.lt.0.5) nonac3=nonac3+1
               if(ppp.lt.1.0) nonac4=nonac4+1
               if(ppp.lt.1.5) nonac5=nonac5+1
               if(ppp.lt.2.0) nonac6=nonac6+1
               if(ppp.lt.2.5) nonac7=nonac7+1
               if(ppp.gt.100.)nonac8=nonac8+1
               if(ppp.gt.80.)nonac9=nonac9+1

               if(lflag.eq.-100)nonac10 = nonac10 + 1

*        
*---------------------- good reconstructed tracks
*

        else                          ! good reconstructed and associated tracks


        pileup=.false.
        p_gen =ped_tab(ngl)
        p_gen2=ped_tab2(ngl,1)                      ! momentum of gen.track in first plane 
        if(p_gen.lt.0)then
           p_gen=-p_gen
           pileup=.true.
        endif

        if(lflag.eq.0)then
          if(fri)then
              print*,'error in logic 2- fringe/no fringe rec.'
          endif

        accep = .true.
        if(label(ngl).eq.10000)accep = .false.

        
         if(extra.and.pileup)then
      if(druk)print '("event, p_rec, p_gen, p_gen_1pl, chi2, it: ",
     &   i5,3f8.2,i4,f8.2,i2,
     &   " ex. to SM2   pileup track")',iv,ppp,p_gen,p_gen2,
     &  ngl,chi2_ndg,in
               nonac16 = nonac16+1
               nonac18 = nonac18+1

         elseif(extra.and..not.pileup)then
      if(druk)print '("event, p_rec, p_gen, p_gen_1pl, chi2, it: ",
     &   i5,3f8.2,i4,f8.2,i2,
     &   " ex. to SM2 ")',iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,in
               nonac16 = nonac16+1
            
         elseif(.not.extra.and.pileup)then
      if(druk)print '("event, p_rec, p_gen, p_gen_1pl, chi2, it: ",
     &   i5,3f8.2,i4,f8.2,i2,
     &   " pileup track ")', iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,in
               nonac18 = nonac18+1

         elseif(.not.extra.and..not.pileup)then
      if(druk)print '("event, p_rec, p_gen, p_gen_1pl, chi2, it: ",
     &   i5,3f8.2,i4,f8.2,i2)',
     &         iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,in
         endif
        endif



*---- fringe field rec. statistic and printing
        if(lflag.eq.100)then

        if(.not.fri)then
         print*,'error in logic 3- fringe/no fringe rec.'
        endif

             nonac12=nonac12 +1                        ! # of tracks rec. in ff (total)

        if(l_hitg2(j).lt.n_tel2.and.pileup) then

      if(druk)print '("event,ffr.,p_rec,p_gen,p_gen_1pl,chi2: ",
     & i4,3f8.2,i3,f6.2," no cont.-ok ",i3," pileup track ")',
     & iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,l_hitg2(j)

        elseif(l_hitg2(j).lt.n_tel2.and..not.pileup) then

      if(druk)print '("event,ffr.,p_rec,p_gen,p_gen_1pl,chi2: ",
     & i4,3f8.2,i3,f6.2," no cont.-ok ",i3)',
     & iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,l_hitg2(j)

        elseif(l_hitg2(j).ge.n_tel2.and.pileup) then

      if(druk)print '("event,ffr.,p_rec,p_gen,p_gen_1pl,chi2: ",
     & i4,3f8.2,i3,f6.2," track cont.!",i3," pileup track ")',
     & iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,l_hitg2(j)

        nonac11=nonac11+1                              ! # of tracks rec. in ff but cont.

        elseif(l_hitg2(j).ge.n_tel2.and..not.pileup) then

      if(druk)print '("event,ffr.,p_rec,p_gen,p_gen_1pl,chi2: ",
     & i4,3f8.2,i3,f6.2," track cont.!",i3)',
     & iv,ppp,p_gen,p_gen2,ngl,chi2_ndg,l_hitg2(j)

        nonac11=nonac11+1                              ! # of tracks rec. in ff but cont
  
        endif
        endif

*------ end of fringe field rec. staff

        itr_count=itr_count + 1
        if(.not.pileup)itr_count4=itr_count4 + 1
        if(accep)itr_count_new=itr_count_new + 1
        if(accep.and..not.pileup)
     &   itr_count4_new=itr_count4_new + 1


        if(p_gen.ge.2)then
           itr_count1=itr_count1+1
          if(.not.pileup)itr_count5=itr_count5 + 1
          if(accep)itr_count1_new=itr_count1_new + 1
          if(.not.pileup.and.accep)
     &       itr_count5_new=itr_count5_new + 1

        endif

        if(p_gen.ge.2.and.p_gen.le.50)then
           itr_count2=itr_count2+1
          if(.not.pileup)itr_count6=itr_count6 + 1
        if(accep)itr_count2_new=itr_count2_new + 1
        if(.not.pileup.and.accep)
     &     itr_count6_new=itr_count6_new + 1

        endif

        if(p_gen.ge.2.and.p_gen.le.5)then
           itr_count3=itr_count3+1
          if(.not.pileup)itr_count7=itr_count7 + 1
          if(accep)itr_count3_new=itr_count3_new + 1
          if(.not.pileup.and.accep)
     &     itr_count7_new=itr_count7_new + 1
        endif

       endif
      enddo                                            ! end of rec.tracks loop
       endif

* -- now filling xtuple2 (event ntuple)

         if(mc.eq.3) then
           xtuple2(1) =iv
           xtuple2(2) =ntr
           xtuple2(3) =n_acc
           xtuple2(4) =nr
           if(muon)xtuple2(5)=1                        ! reconstructed muon
           if(.not.muon)xtuple2(5)=0                   ! no reconstracted muon in event 

           call hfn(11,xtuple2)
         elseif(mc.eq.5) then
           xtuple4(1) =iv
           xtuple4(2) =ntr
           xtuple4(3) =n_acc
           xtuple4(4) =nr
           xtuple4(5) = itrig
           if(muon)xtuple4(6)=1                        ! reconstructed muon
           if(.not.muon)xtuple4(6)=0                   ! no reconstracted muon in event
 
           call hfn(11,xtuple4)        
         endif

             no(100) = no(100)+ nonac
             no(98) = no(98) + nonac1
             no(97) = no(97) + nonac2
             no(96) = no(96) + nonac3
             no(95) = no(95) + nonac4
             no(94) = no(94) + nonac5
             no(93) = no(93) + nonac6
             no(92) = no(92) + nonac7
             no(91) = no(91) + nonac8
             no(90) = no(90) + nonac9
             no(89) = no(89) + nonac10
             no(88) = no(88) + nonac11
             no(87) = no(87) + nonac12
             no(86) = no(86) + nonac13
             no(85) = no(85) + nonac14
             no(84) = no(84) + nonac15
             no(83) = no(83) + nonac16
             no(82) = no(82) + nonac17
             no(81) = no(81) + nonac18
             no(66) = no(66) + nonac19

             no(80) = no(80) + itr_count1
             no(65) = no(65) + itr_count1_new

             no(79) = no(79) + itr_count2
             no(64) = no(64) + itr_count2_new

             no(78) = no(78) + itr_count3
             no(63) = no(63) + itr_count3_new

             no(77) = no(77) + itr_count4
             no(62) = no(62) + itr_count4_new

             no(76) = no(76) + itr_count5
             no(61) = no(61) + itr_count5_new

             no(75) = no(75) + itr_count6
             no(60) = no(60) + itr_count6_new

             no(74) = no(74) + itr_count7
             no(59) = no(59) + itr_count7_new
      return
      end










