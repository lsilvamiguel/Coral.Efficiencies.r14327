         subroutine mcini(ih1,ih2,ih3,ip1,ip2,ip3,
     +   phit1,phit2,phit3,err1,err2,err3,
     + itr1,itr2,itr3,n_tracks,ped_tr,ped_hit,lab_tr)

**************************************************
*                                                *
*       some MC info  created                    *
* modified Sept 17.2000 (kk)                     *
* modified Jan  28.2001 (AM)                     *
*                                                *
**************************************************

#include "PARAM.h"
#include "INI_STAT.h"
#include "PLANES.h"
#include "HITS_IN.h"
#include "MC.h"

 
       INTEGER ih1,ih2,ih3, n_tracks
       INTEGER ip1(*), ip2(*), ip3(*)
       INTEGER itr1(*), itr2(*), itr3(*)
       INTEGER lab_tr(*)

       DOUBLE PRECISION phit1(*), phit2(*), phit3(*)
       DOUBLE PRECISION err1(*), err2(*), err3(*)
       DOUBLE PRECISION ped_tr(*),ped_hit(6,*)

       logical fringe,extrapol,pileup,ety,printout,druk
       logical accep,pat_hit,sub

       printout=.false.

         npa=0
         npa2=0
        if(sub)then
         npa=ntpl_cor
         npa2=ntpl2_cor
        endif

       call gluehits(ih1,ih2,ih3,ip1,ip2,ip3,
     & phit1,phit2,phit3,err1,err2,err3)

       if(pat_hit)return
           
      ntr=n_tracks
      ntr_ef=0
      ntr_ef1=0
      ntr_ef2=0
      ntr_ef3=0
      ntr_ef4=0
      ntr_ef5=0
      ntr_ef6=0
      ntr_ef7=0

      ntr_ef_new=0
      ntr_ef1_new=0
      ntr_ef2_new=0
      ntr_ef3_new=0
      ntr_ef4_new=0
      ntr_ef5_new=0
      ntr_ef6_new=0
      ntr_ef7_new=0

      n_tel1=nn_tel1
      n_tel2=nn_tel2


      if(ntr.eq.0)then
         print*,'no acceptable tracks !!! - patological event'
         return
      endif

       if(ntr.gt.ntrack)then
       print*,'errors in dimension - ntrack too small (mc)'
       print*, ntr,ntrack
       stop
       endif
            
      do i=1,ntr
        
        if(ped_tr(i).lt.0)then
          ped_tab(i)=-ped_tr(i)
          pileup=.true.
        else
          ped_tab(i)=ped_tr(i)
          pileup=.false.
        endif

          ped_tab2(i,1)=
     &   sqrt(ped_hit(1,i)**2+ped_hit(2,i)**2+ped_hit(3,i)**2)

          ped_tab2(i,2)= ped_hit(1,i)/ped_hit(3,i)
          ped_tab2(i,3)= ped_hit(2,i)/ped_hit(3,i)

*           print*,i,(ped_tab2(i,j),j=1,3)

          xped_tab2(i,1)= ped_hit(4,i)/1000.
          xped_tab2(i,2)= ped_hit(5,i)/1000.
          xped_tab2(i,3)= ped_hit(6,i)/1000.
*          print*,i,(xped_tab2(i,j),j=1,3)

        label(i)=0

        nsum1=0
        nsum2=0
        nsumm=0
        nsumy1=0
        nsumy2=0
        nsumz1=0
        nsumz2=0
        nsumd=0
        nsumsfi=0
        nsumg=0
        nsumsup=0

        ksd = 0
        nd2 = 0
        ndz=0
        ndy=0

        if(lab_tr(i).gt.600) then
           label(i)=-1              ! # of hits < n_tel1
           go to 3
        endif

*       print*, 'track',i,ped_tab(i)

        ih=0
        ij=1
        iad=2
        do k=1,ntpl_tmp
         if (ityp1(npa+k).eq.222)ksd=ksd+1
         if(jchits(k).ge.1) then
           do j =1, jchits(k)          ! over clusters in plane k

              ih=ih+1
              if(itr1(iad-1).ne.ih) goto 33
              if(numhit(k,j).gt.ih) then

*   some hits not used, update pointer in itr1 (iad)

                 do iup=ih,numhit(k,j)-1
                    istep=itr1(iad)
                    iad=iad+2+istep
                 enddo

               ih=numhit(k,j)
               endif

               if(itr1(iad-1).ne.ih) go to 33

*              print *,'ih,iad,no of tr.',ih,iad,itr1(iad)
*              print*,'total',k,j,itr1(iad)

               inte=0
               l_hit1(1,k,j)=0
               do jop=1, itr1(iad)     ! over tracks in cluster numhit(k,j)

                if((iad+jop).gt.montecarlo)then
                 print*,'dimension ERROR 1',montecarlo,iad,jop
                 stop
                endif
*               print*,itr1(iad+jop)

         
               if(itr1(iad+jop).ne.200)then
    
*                 print*,inte  
                    ety=.true.
                  do ii=2,inte+1
                  if(itr1(iad+jop).eq.l_hit1(ii,k,j))ety=.false.
                  enddo

                 if(ety)then
                   inte=inte+1

                  if(l_hit1(1,k,j).ge.maxiclus)then
                   nos(67)=nos(67)+1
*                   print*,k,j,inte,itr1(iad+jop),l_hit1(1,k,j)
*                   print*,(itr1(iad+jopp),jopp=1,itr1(iad))
                   l_hit1(1,k,j)=maxiclus
                   go to 11
                  endif

                 l_hit1(1,k,j)=inte
                 l_hit1(inte+1,k,j)=itr1(iad+jop)

*                 print *,inte, l_hit1(inte+1,k,j)                
*                 print *,'track no.',itr1(iad+jop)

                  if(l_hit1(inte+1,k,j).eq.i.and.
     &      abs(clusters1(numhit(k,j))).gt.1.e-5) then
*                print*, 'cluster-track',jop,i,l_hit1(inte+1,k,j)  
                 
                  nsum1=nsum1+1
                if (ityp1(npa+k).eq.111) nsumm=nsumm+1
      if (ityp1(npa+k).eq.111.and.jyzt(npa+k).eq.1.and.
     &         nsumm.le.4)then
                  nsumy1=nsumy1+1
                  clu1=clusters1(numhit(k,j))
                  kpla1=k
      elseif(ityp1(npa+k).eq.111.and.jyzt(npa+k).eq.2.and.
     &         nsumm.le.4)then
                  nsumz1=nsumz1+1
                  clu2=clusters1(numhit(k,j))
                  kpla2=k
      elseif(ityp1(npa+k).eq.111.and.jyzt(npa+k).eq.1.and.
     &         nsumm.gt.4)then
                  nsumy2=nsumy2+1
                  clu3=clusters1(numhit(k,j))
                  kpla3=k
      elseif(ityp1(npa+k).eq.111.and.jyzt(npa+k).eq.2.and.
     &         nsumm.gt.4)then
                  nsumz2=nsumz2+1
                  clu4=clusters1(numhit(k,j))
                  kpla4=k
      endif

                if (ityp1(npa+k).eq.222) nsumd=nsumd+1
*------ just for fringe field acceptance
            if (ityp1(npa+k).eq.222.and.ksd.gt.8) nd2=nd2+1
            if (ityp1(npa+k).eq.222.and.jyzt(npa+k).eq.2.and.
     &          ksd.gt.8)ndz=ndz+1
            if (ityp1(npa+k).eq.222.and.jyzt(npa+k).eq.1.and.
     &          ksd.gt.8)ndy=ndy+1
            if (ityp1(npa+k).eq.333) nsumsfi=nsumsfi+1
                  endif
                  endif
                  endif
 11               continue
               enddo

*               print*,i,(l_hit1(jjj,k,j),jjj=1,l_hit1(1,k,j)+1)

              iad=iad+2+itr1(iad)
           enddo
          endif
        enddo
                iad_before = iad

        if(nsum1.gt.ntpl)label(i)=-100

        if(nsum1.lt.nn_tel1.and.label(i).eq.0)then
           label(i)=-1            
           go to 3
        endif

*--- special treatment these of generated tracks which have clusters in micros

          accep=.true.

          if(nsumy1.ge.1.and.nsumz1.ge.1.and.
     &       nsumy2.ge.1.and.nsumz2.ge.1)
     &   call accep_tar(accep,kpla1,kpla2,kpla3,kpla4,
     &                   clu1,clu2,clu3,clu4)


        ih=0
        ij=2
        iad=2

        
        do k=1,ntpl2_tmp
         if ( jc2hits(k).ge.1) then
           do j =1, jc2hits(k)

              ih=ih+1

              if(numht2(k,j).gt.ih) then

*    some clusters not used, update pointer in itr2 (iad)

                 do iup=ih,numht2(k,j)-1
                    if(iup.ne.itr2(iad-1)) go to 33

*                   print *,'ih.ia,itr...',ih,iad,itr2(iad-1), 
*     +                 itr2(iad), itr2(iad+1)

                    istep=itr2(iad)
                    iad=iad+2+istep
                 enddo
              ih=numht2(k,j)
              endif

              if(itr2(iad-1).ne.ih) goto 33

*              print *,'ih,numht2,itr:',ih,numht2(k,j),itr2(iad-1),
*     +              'iad, no of tr.',  iad,itr2(iad),'  pl,hit',k,j

   
              inte=0
              l_hit2(1,k,j)=0
*              print*,'total',k,j,itr2(iad)
              do jop=1, itr2(iad)                 ! over tracks in cluster numhit(k,j)
*               print*,itr2(iad+jop)

                if((iad+jop).gt.montecarlo)then
                 print*,'dimension ERROR 2',montecarlo,iad,jop
                 stop
                endif

               if(itr2(iad+jop).ne.200)then
            
*                 print*,inte           
                    ety=.true.
                  do ii=2,inte+1
                  if(itr2(iad+jop).eq.l_hit2(ii,k,j))ety=.false.
                  enddo

                 if(ety)then
                   inte=inte+1 
                  if(inte.ge.maxiclus)then
                   nos(68)=nos(68)+1
                   l_hit2(1,k,j)=maxiclus
                   go to 12
                  endif

                 l_hit2(1,k,j)=inte
                 l_hit2(inte+1,k,j)=itr2(iad+jop)

*                 print *,inte, l_hit2(inte+1,k,j)                  



                 if(l_hit2(inte+1,k,j).eq.i.and.
     &       abs(clusters2(numht2(k,j))).gt.1.e-5) then
                    nsum2=nsum2+1
                    if (ityp2(k).eq.444) nsumg=nsumg+1 
                 endif
                 endif
                 endif
 12              continue
              enddo
         
*               print*,(l_hit2(jjj,k,j),jjj=1,l_hit2(1,k,j)+1)

              iad=iad+2+itr2(iad)
            enddo
          endif
         enddo
               
                iad_behind=iad    

         if(nsum2.gt.ntpl2)
     &        label(i)=-200

         if(.not.fringe.and.nsum2.lt.nn_tel2)then
            label(i) = -1
            go to 3
 
         elseif(.not.fringe.and.nsum2.ge.nn_tel2)then

            jj_hit(1,i) = nsum1 + 1000*nsum2 + 100000*nsumm
            jj_hit(2,i) = nsumd + 1000*nsumsfi+100000*nsumg
       
            ntr_ef=ntr_ef +1
            if(accep)ntr_ef_new=ntr_ef_new+1

            if(.not.pileup)ntr_ef4=ntr_ef4+1
         if(.not.pileup.and.accep)ntr_ef4_new=ntr_ef4_new+1

            if(ped_tab(i).gt.2.)then
               ntr_ef1=ntr_ef1+1
            if(accep)ntr_ef1_new=ntr_ef1_new+1
               if(.not.pileup)ntr_ef5=ntr_ef5+1
         if(.not.pileup.and.accep)ntr_ef5_new=ntr_ef5_new+1
            endif
         
            if(ped_tab(i).ge.2.and.ped_tab(i).le.50.)then
               ntr_ef2=ntr_ef2+1
            if(accep)ntr_ef2_new=ntr_ef2_new+1
               if(.not.pileup)ntr_ef6=ntr_ef6+1
         if(.not.pileup.and.accep)ntr_ef6_new=ntr_ef6_new+1
            endif

            if(ped_tab(i).ge.2.and.ped_tab(i).le.5.)then
               ntr_ef3=ntr_ef3+1
            if(accep)ntr_ef3_new=ntr_ef3_new+1
               if(.not.pileup)ntr_ef7=ntr_ef7+1
         if(.not.pileup.and.accep)ntr_ef7_new=ntr_ef7_new+1
            endif

         elseif(fringe)then

         if(nsum2.ge.nn_tel2.or.
     &   (ped_tab(i).gt.0.25.and.ndz.ge.1.and.ndy.ge.1))then        

              jj_hit(1,i) = nsum1 + 1000*nsum2 + 100000*nsumm
              jj_hit(2,i) = nsumd + 1000*nsumsfi+100000*nsumg       
         
              ntr_ef=ntr_ef +1
          if(accep)ntr_ef_new=ntr_ef_new+1
              if(.not.pileup)ntr_ef4=ntr_ef4+1
          if(.not.pileup.and.accep)ntr_ef4_new=ntr_ef4_new+1
              if(ped_tab(i).gt.2)then
                 ntr_ef1=ntr_ef1+1
          if(accep)ntr_ef1_new=ntr_ef1_new+1
                 if(.not.pileup)ntr_ef5=ntr_ef5+1
          if(.not.pileup.and.accep)ntr_ef5_new=ntr_ef5_new+1
              endif
         
              if(ped_tab(i).ge.2.and.ped_tab(i).le.50)then
                 ntr_ef2=ntr_ef2+1
            if(accep)ntr_ef2_new=ntr_ef2_new+1
                 if(.not.pileup)ntr_ef6=ntr_ef6+1
            if(.not.pileup.and.accep)ntr_ef6_new=ntr_ef6_new+1
              endif

              if(ped_tab(i).ge.2.and.ped_tab(i).le.5)then
                 ntr_ef3=ntr_ef3+1
            if(accep)ntr_ef3_new=ntr_ef3_new+1
                 if(.not.pileup)ntr_ef7=ntr_ef7+1
            if(.not.pileup.and.accep)ntr_ef7_new=ntr_ef7_new+1
              endif
           else
              label(i)=-1
           endif

        else
           label(i)=-1
        endif

*--- now extrapolation to SM2 (behind Rich)
        ih=0
        ij=3
        iad=2
        do k=1,ntpl_sup
         if(jcshits(k).ge.1) then
           do j =1, jcshits(k)

              ih=ih+1

              if(numhts(k,j).gt.ih) then
*        some hits not used, update iadess in itr3 array
                 do iupd=ih,numhts(k,j)-1
              if(iupd.ne.itr3(iad-1)) go to 33      

                    istep=itr3(iad)
                    iad=iad+2+istep
                 enddo
                 ih=numhts(k,j)
              endif

              if(itr3(iad-1).ne.ih) goto 33


*              print*,'total',k,j,itr3(iad)
              inte=0
              l_hit3(1,k,j)=0
              do jop=1, itr3(iad)                 ! over tracks in cluster numhit(k,j)

*               print*,itr3(iad+jop)

                if((iad+jop).gt.montecarlo)then
                 print*,'dimension ERROR 3',montecarlo,iad,jop
                 stop
                endif

               if(itr3(iad+jop).ne.200)then
            
*                 print*,inte           
                    ety=.true.
                  do ii=2,inte+1
                  if(itr3(iad+jop).eq.l_hit3(ii,k,j))ety=.false.
                  enddo

                 if(ety)then
                    inte=inte+1
                  if(inte.ge.maxiclus)then
                   nos(69)=nos(69)+1
                   l_hit3(1,k,j)=maxiclus
                   go to 13
                  endif

                 l_hit3(1,k,j)=inte
                 l_hit3(inte+1,k,j)=itr3(iad+jop)
*                   print*,inte, l_hit3(inte+1,k,j)                

                 if(l_hit3(inte+1,k,j).eq.i.and.
     &     abs(clusters3(numhts(k,j))).gt.1.e-5)then
                    nsumsup = nsumsup + 1
                    jj_hit(3,i) = nsumsup
                 endif
                 endif
                 endif
 13              continue
              enddo
           
*               print*,(l_hit3(jjj,k,j),jjj=1,l_hit3(1,k,j)+1)
              iad=iad+2+itr3(iad)


           enddo
         endif
        enddo
              
               iad_rich=iad   
    

        if(nsumsup.gt.ntpl_sup)label(i)=-300
 
 3      continue
         
        if(pileup)ped_tab(i)=-ped_tab(i)
        if(.not.accep)ped_tab(i)=10000. + ped_tab(i)

      enddo                                              ! end loop over tracks

*           print*,'before SM1: ', iad_before
*           print*,'behind SM1: ', iad_behind
*           print*,'behind RICH: ', iad_rich


      ntr_licz=ntr_licz+ntr                              ! # of generated tracks         
      ntr_ac_licz = ntr_ac_licz + ntr_ef                 ! # of accepted tracks
      ntr_ac_licz_new =ntr_ac_licz_new + ntr_ef_new      ! # of accepted tracks (new) 

      ntra_ac(1) = ntra_ac(1) + ntr_ef1       ! # of accepted tracks p>2
      ntra_ac(2) = ntra_ac(2) + ntr_ef2       ! # of accepted tracks 2<p<50
      ntra_ac(3) = ntra_ac(3) + ntr_ef3       ! # of accepted tracks 2<p<5

      ntra_ac(8) = ntra_ac(8) + ntr_ef1_new       ! # of accepted tracks p>2
      ntra_ac(9) = ntra_ac(9) + ntr_ef2_new       ! # of accepted tracks 2<p<50
      ntra_ac(10) = ntra_ac(10) + ntr_ef3_new     ! # of accepted tracks 2<p<5
*---- now the same for pilup_out calculation

      ntra_ac(4) = ntra_ac(4) + ntr_ef4       ! # of accepted tracks full range
      ntra_ac(5) = ntra_ac(5) + ntr_ef5       ! # of accepted tracks p>2
      ntra_ac(6) = ntra_ac(6) + ntr_ef6       ! # of accepted tracks 2<p<50
      ntra_ac(7) = ntra_ac(7) + ntr_ef7       ! # of accepted tracks 2<p<5

      ntra_ac(11) = ntra_ac(11) + ntr_ef4_new       ! # of accepted tracks full range
      ntra_ac(12) = ntra_ac(12) + ntr_ef5_new       ! # of accepted tracks p>2
      ntra_ac(13) = ntra_ac(13) + ntr_ef6_new       ! # of accepted tracks 2<p<50
      ntra_ac(14) = ntra_ac(14) + ntr_ef7_new       ! # of accepted tracks 2<p<5

*------------- internal check out
      if(printout)then
      print *,'before SM1'
      do k=1,ntpl
      do j=1,jchits(k)
      print*,k,j,numhit(k,j),l_hit1(1,k,j),
     &       (l_hit1(ii,k,j),ii=2,l_hit1(1,k,j)+1)
      enddo
      enddo
      print *,'behind SM1'
      do k=1,ntpl2
      do j=1,jc2hits(k)
      print*,k,j,numht2(k,j),l_hit2(1,k,j),
     &       (l_hit2(ii,k,j),ii=2,l_hit2(1,k,j)+1)
      enddo
      enddo
      print *,'behind Rich'
        do k=1,ntpl_sup
      do j=1,jcshits(k)
      print*,k,j,numhts(k,j),l_hit3(1,k,j),
     &       (l_hit3(ii,k,j),ii=2,l_hit3(1,k,j)+1)
      enddo
      enddo
      endif

      return
 33   continue

*   error in input structure detected - program stopped

      print *, 'MC_INI - wrong input in itr',ij
      print *, 'plane',k,' cluster no.',j,' index in itr',iad,' ih=',ih

      if(ij.eq.2) print *,'  itr2(iad-1)=',itr2(iad-1),numht2(k,j),
     + '  itr2(iad)(no of tr.)=', itr2(iad)
      if(ij.eq.1) print *,'  itr1(iad-1)=',itr1(iad-1),numhit(k,j),
     + '  itr1(iad)(no of tr.)=', itr1(iad)
      if(ij.eq.3) print *,'  itr3(iad-1)=',itr3(iad-1),numhts(k,j),
     + '  itr3(iad)(no of tr.)=', itr3(iad)
      print *,' program STOPPED in MC-INI '
      STOP 33

      end 



      subroutine accep_tar(accep,jply1,jplz1,jply2,jplz2,
     &                  y1,z1,y2,z2)
*****************************************************
*                                                   *
* to accept/no accept tracks which are not comming  *
* from the target (the same criteria as for recon.) *
*                                                   *
*****************************************************
#include "PARAM.h"
#include "TARGET.h"
#include "PLANES.h"
#include "TUNE.h"
        
        common/local_accep/ rd_yy,rd_zz

        logical accep,sub

        data fac/3.0/
        data itime/1/

        if(itime.eq.1)then
          rd_yy=rd_y(2)
          rd_zz=rd_z(2)
          itime=0
        endif
        if(sub)npa=ntpl_cor

        accep=.false.              ! non accepted at the beginning
        tar_rad = ytarg

*   target limits.........................
      x_st = xtarg0 - hltarg
      x_nd = xtarg0 + hltarg

c-kk -- target shot
*------- take into account target radius.............................


    
        if (z2.gt.0.0)  then

          zup  = (xplane(jplz1)-x_st)/(xplane(jplz2)-x_st) * z2
     1   + (xplane(jplz2)-xplane(jplz1))/(xplane(jplz2)-x_st) * tar_rad
          zdw  = (xplane(jplz1)-x_nd)/(xplane(jplz2)-x_nd) * z2
     1   - (xplane(jplz2)-xplane(jplz1))/(xplane(jplz2)-x_nd) * tar_rad

        else

          zup = (xplane(jplz1)-x_nd)/(xplane(jplz2)-x_nd) * z2
     1   + (xplane(jplz2)-xplane(jplz1))/(xplane(jplz2)-x_nd) * tar_rad
          zdw = (xplane(jplz1)-x_st)/(xplane(jplz2)-x_st) * z2
     1   - (xplane(jplz2)-xplane(jplz1))/(xplane(jplz2)-x_st) * tar_rad
  
        endif

        zup = zup + fac*rd_zz
        zdw = zdw - fac*rd_zz





        if (z1.gt.zdw.and.z1.lt.zup)then                     ! chance 


*------- take into account target radius.............................

        if (y2.gt.0.0)  then

       yup  = (xplane(npa+jply1)-x_st)/(xplane(npa+jply2)-x_st) * y2
     1 + (xplane(npa+jply2)-xplane(npa+jply1))/(xplane(npa+jply2)-x_st) 
     2 * tar_rad
       ydw  = (xplane(npa+jply1)-x_nd)/(xplane(npa+jply2)-x_nd) * y2
     1 - (xplane(npa+jply2)-xplane(npa+jply1))/(xplane(npa+jply2)-x_nd) 
     2 * tar_rad

        else

       yup = (xplane(npa+jply1)-x_nd)/(xplane(npa+jply2)-x_nd) * y2
     1 + (xplane(npa+jply2)-xplane(npa+jply1))/(xplane(npa+jply2)-x_nd)
     2 * tar_rad
       ydw = (xplane(npa+jply1)-x_st)/(xplane(npa+jply2)-x_st) * y2
     1 - (xplane(npa+jply2)-xplane(npa+jply1))/(xplane(npa+jply2)-x_st) 
     2 * tar_rad
  
        endif
        yup = yup +  fac*rd_yy
        ydw = ydw -  fac*rd_yy

            if (y1.gt.ydw.and.y1.lt.yup)then                 ! track ok!
           
            accep = .true.
            endif
         endif


          return
          end

