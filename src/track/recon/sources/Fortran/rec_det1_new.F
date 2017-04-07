      subroutine rec_det1_new(modset,jfringe,flaga)
*************************************************************
*                                                           *
*                det.1    reconstruction                    *
*                                                           *
*    New version for 2002 configuration (June 13,2002)      *
*    ---Modified for low momentum tracks---(R.W.)----       *
*                                                           *
*************************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "TARGET.h"
#include "HITS_IN.h"
#include "HITS_OUT.h"
#include "TUNE.h"
#include "LIN_ALL.h"

      dimension lref(4),ivref(4)
      dimension indpl(mxpl),ncand(mxpl),lrem(mxpl)
      dimension ivv(mxpl,nmhit)
      logical imodtg,jmodtg,last_chance
      logical modset,jfringe,flaga

         last_chance=.false.
         ncth=ncth1

         if(flaga)nlin_fda=0
       
 1000    continue                               ! last chance
         if(last_chance) ncth=int(ncth1/2)
 
         eps_det =1.e-5

         do jo=1,4
          lref(jo)=0
          ivref(jo)=0
         enddo

         do jo=1,ntpl+ntpl2
          ncand(jo)=0
          lrem(jo)=0

         do io=1,nmhit
          ivv(jo,io)=0
          jflag(jo,io)=0
         enddo
         enddo
    
*   target radius 

      tar_rad = ytarg
    
*   center of telescope
    
      if(abs(det1_act(1)).lt.eps_det.
     &   and.abs(det1_act(2)).lt.eps_det)then
 
      x_mid = 0.5 *(xplane(1) + xplane(ntpl))

*   take 10 % of length of telescope
      x_ten = ten1 *(xplane(ntpl) - xplane(1))

      else
      x_mid = 0.5 *(det1_act(1)+ det1_act(2))

*   take 10 % of length of telescope
      x_ten = ten1 *(det1_act(2) - det1_act(1))
      
      endif
* --- end of center 

*   target limits.........................
      x_st = xtarg0 - hltarg
      x_nd = xtarg0 + hltarg

*    Nr of lines found so far .....................
      nlin_fd = 0

*    flag usable hits.............................

*   z segment.............................

      do 10 kz2 = nzpl,2,-1                       ! loop over Z planes

       jplz2 = jzpl(kz2)                          ! real # pl w det. 1 (ref.1)
      
*  in FF, use only DC planes (type 222)

       if (jfringe.and.ityp1(jplz2).ne.222)  go to 10     
     
*  escape for FF............................
 
       if (jfringe.and.jchits(jplz2).gt.2*ncth1) go to 10 

*   leave if center of telescope is reached
       if (xplane(jplz2).lt.(x_mid - x_ten)) go to 99

        if(jchits(jplz2).eq.0)go to 10
        do 11 jhz = 1,jchits(jplz2)         
*          print*,jchits(jplz2)

       if (jflag(jplz2,jhz).ne.0)  go to 11 
      
         z2 = clusters1(numhit(jplz2,jhz))
       

        do 20 kz1 = 1,kz2-1         
                           
        jplz1 = jzpl(kz1)                         ! real # pl w det. 1 (ref.2)


*  escape for FF............................
 
       if (jfringe.and.jchits(jplz1).gt.2*ncth1) go to 20

*   stop loop at center of telescope
        if (xplane(jplz1).gt.(x_mid + x_ten)) go to 11


*    modify target pointing for Mmegas if high multiplicity............

        imodtg = .false.
      if(jchits(jplz2).gt.ncth.and.jchits(jplz1).gt.ncth) then
*           print *, " modified tg.pointing Z planes",jplz2,jplz1
           imodtg = .true.
        endif
        


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

        zup = zup + rd_z(1)
        zdw = zdw - rd_z(1)

*------- find optimal z1 if modified target pointing...........

        if (imodtg) then
         dtgmin = 9999.
         itgopt = -99
         zpred = z2 * (xplane(jplz1) - xtarg0)/(xplane(jplz2) - xtarg0)

         if(jchits(jplz1).eq.0)go to 20
         do khz = 1,jchits(jplz1)
         if (jflag(jplz1,khz).eq.0) then
         z1 = clusters1(numhit(jplz1,khz))  
         dtg = abs(z1 - zpred)
         if (dtg.lt.dtgmin) then
            dtgmin = dtg
            itgopt = khz
         endif
         endif
         enddo
        endif


        ncnt_z2 = 0                          

         do 21 khz = 1,jchits(jplz1)                ! loop over hits in Z pl 2

         z1 = clusters1(numhit(jplz1,khz))  


         if (jflag(jplz1,khz).ne.0) go to 21        ! repetition flag

         if (z1.gt.zdw.and.z1.lt.zup) then

*   more selective if modified target pointing......

         if (imodtg.and.khz.ne.itgopt) go to 21



*  z proj. is OK .................................................
         indpl(jplz1) = khz               
         indpl(jplz2) =  jhz
         zpr = (z2-z1)/(xplane(jplz2)-xplane(jplz1))
         ncnt_z2 = ncnt_z2 + 1
         lref(3) = jplz1
         lref(4) = jplz2
         ivref(3) = khz
         ivref(4) = jhz

*  now loop over Y planes.........................................

      do 30 ky2 = nypl,2,-1
       jply2 = jypl(ky2)
*----require same typ as jplz2......................

       if (ityp1(jply2).ne.ityp1(jplz2))  go to 30

*----require (Z,Y) not too far away (max. 2 cm)..................

       if (abs(xplane(jply2)-xplane(jplz2)).gt.distmx_yz) go to 30

*   in FF use only DC planes (type 222)

       if (jfringe.and.ityp1(jply2).ne.222)  go to 30
    

*  escape for FF............................
 
       if (jfringe.and.jchits(jply2).gt.2*ncth1) go to 30

      
*  leave when center of telescope is reached

       if (xplane(jply2).lt.(x_mid - x_ten)) go to 21

       if(jchits(jply2).eq.0)go to 30
       do 31 jhy = 1,jchits(jply2)
     

       if (jflag(jply2,jhy).ne.0) go to 31         ! repetition flag
       y2 = clusters1(numhit(jply2,jhy))

*  check compatibility of Y2 and Z2 using nearest theta plane

       xy2z2 = 0.5 * (xplane(jply2) + xplane(jplz2))
       iok2 = 0

      do 330  jth2 = 1,nthpl
        jplth2 = jthpl(jth2)
        if (abs(xplane(jplth2)-xy2z2).gt.distmx_th) go to 330
        if (ityp1(jplth2).ne.ityp1(jplz2)) go to 330
        ttth2 = sinn(jplth2)*y2 - cosn(jplth2)*z2
  
      if(jchits(jplth2).eq.0)go to 330

       do jrr = 1,jchits(jplth2)
        thhit = clusters1(numhit(jplth2,jrr))
        if (abs(thhit-ttth2).lt.rd_comp) iok2 = 1
       enddo

 330  continue
       if (iok2.eq.0) go to 31


*       print *,"jply2,jplz2,iok2",jply2,jplz2,iok2
*---------------------------------------------------------

        do 40 ky1 = 1,ky2-1       
        jply1 = jypl(ky1)
 
*---- require same type as jplz1

        if (ityp1(jply1).ne.ityp1(jplz1)) go to 40

*---- require Z1, Y1 not too far away..............

        if (abs(xplane(jply1)-xplane(jplz1)).gt.distmx_yz) go to 40


*  escape for FF............................
 
       if (jfringe.and.jchits(jply1).gt.2*ncth1) go to 40 

*  stop at center of telescope

        if (xplane(jply1).gt. (x_mid + x_ten))  go to 31
    
*    modify target pointing for Mmegas if high multiplicity............

        jmodtg = .false.
      if(jchits(jply2).gt.ncth.and.jchits(jply1).gt.ncth) then
*           print *, " modified tg.pointing Y planes",jply2,jply1
           jmodtg = .true.
        endif

*------- take into account target radius.............................

        if (y2.gt.0.0)  then

          yup  = (xplane(jply1)-x_st)/(xplane(jply2)-x_st) * y2
     1   + (xplane(jply2)-xplane(jply1))/(xplane(jply2)-x_st) * tar_rad
          ydw  = (xplane(jply1)-x_nd)/(xplane(jply2)-x_nd) * y2
     1   - (xplane(jply2)-xplane(jply1))/(xplane(jply2)-x_nd) * tar_rad

        else

          yup = (xplane(jply1)-x_nd)/(xplane(jply2)-x_nd) * y2
     1   + (xplane(jply2)-xplane(jply1))/(xplane(jply2)-x_nd) * tar_rad
          ydw = (xplane(jply1)-x_st)/(xplane(jply2)-x_st) * y2
     1   - (xplane(jply2)-xplane(jply1))/(xplane(jply2)-x_st) * tar_rad
  
        endif
        yup = yup +  rd_y(1)
        ydw = ydw -  rd_y(1)

*------- find optimal y1 if modified target pointing...........

        if (jmodtg) then
         dtgmin = 9999.
         jtgopt = -99
         ypred = y2 * (xplane(jply1) - xtarg0)/(xplane(jply2) - xtarg0)

         if(jchits(jply1).eq.0)go to 40

         do khy = 1,jchits(jply1)
         y1 = clusters1(numhit(jply1,khy)) 
         if (jflag(jply1,khy).eq.0) then  
         dtg = abs(y1 - ypred)
         if (dtg.lt.dtgmin) then
            dtgmin = dtg
            jtgopt = khy
         endif
         endif
         enddo
        endif

        ncnt_y2 = 0

         do 41 khy = 1,jchits(jply1) 

         y1 = clusters1(numhit(jply1,khy))


         if (jflag(jply1,khy).ne.0) go to 41    ! repetition flag


         if (y1.gt.ydw.and.y1.lt.yup)then

*  y proj. is OK .................................................


*   more selective if modified target pointing......

         if (jmodtg.and.khy.ne.jtgopt) go to 41

* ---Check compatibility between y1 and z1 using nearest theta plane......

         xy1z1 = 0.5 * (xplane(jply1) + xplane(jplz1))
         iok1 = 0

        do 440 jth1 = 1,nthpl     
          jplth1 = jthpl(jth1)
          if (abs(xplane(jplth1)-xy1z1).gt.distmx_th) go to 440
          if (ityp1(jplth1).ne.ityp1(jplz1))  go to 440
          tth1 = sinn(jplth1)*y1 - cosn(jplth1)*z1

         if(jchits(jplth1).eq.0)go to 440
         do jrr = 1,jchits(jplth1)
          thhit = clusters1(numhit(jplth1,jrr))
          if (abs(thhit-tth1).lt.rd_comp) iok1 = 1
         enddo
 440    continue
         if (iok1.eq.0)  go to 41

*          print *, "jply1,jplz1,iok1:",jply1,jplz1,iok1
*-----------------------------------------
 
         indpl(jply2) = jhy      ! index = cluster in  pl 1
         indpl(jply1) = khy      ! index = cluster in  pl 2
        
         ypr = (y2-y1)/(xplane(jply2)-xplane(jply1))
C-kk  --     ypr - slope from Y ref.planes  
       

         ncnt_y2 = ncnt_y2 + 1    !  index of "good" rec.
         lref(1) = jply1         
         lref(2) = jply2
         ivref(1) = khy           ! hit # in 1 and 2
         ivref(2) = jhy 

*  now define roads for theta  planes and remaining Y,Z planes................
         jrem = 0

         yy0 = clusters1(numhit(jply1,indpl(jply1)))
         zz0 = clusters1(numhit(jplz1,indpl(jplz1)))

         do 141 jplth = 1,ntpl             ! loop over remaining pl in tel.1

         if (jplth.eq.lref(1).or.jplth.eq.lref(2).or.jplth.eq.lref(3)
     1    .or.jplth.eq.lref(4))go to 141

*----type must correspond to type of one of the ref. points............

         jrem = jrem + 1
         lrem(jrem) = jplth
         ncand(jrem) = 0

         if (ityp1(jplth).ne.ityp1(jplz2).and.             
     1      ityp1(jplth).ne.ityp1(jplz1)) go to 141

         yyth = yy0 + ypr * (xplane(jplth)-xplane(jply1))
         zzth = zz0 + zpr * (xplane(jplth)-xplane(jplz1))
         ttth = sinn(jplth) * yyth - cosn(jplth) * zzth

         
*----in bending plane, include max. curvature in road.................

         rd_cor = sagi_tot
         if (abs(xplane(jplz2)-x_mid).lt.x_ten.or.
     1    abs(xplane(jplz1)-x_mid).lt.x_ten) rd_cor = 2 * sagi_tot

         yyup = yyth + rd_cor
         yylw = yyth - rd_cor
         zzup = zzth + rd_th(1)
         zzlw = zzth - rd_th(1)
         th1 = sinn(jplth) * yyup - cosn(jplth) * zzup
         th2 = sinn(jplth) * yyup - cosn(jplth) * zzlw
         th3 = sinn(jplth) * yylw - cosn(jplth) * zzup
         th4 = sinn(jplth) * yylw - cosn(jplth) * zzlw

         thup = max(th1,th2,th3,th4)
         thdw = min(th1,th2,th3,th4)
 
*  select theta hits within roads............................................

             if(jchits(jplth).eq.0)then
                 ncand(jrem)=-100                ! empty plane          
                 go to 141
             endif

          do jrr=1,jchits(jplth)                 ! loop over hits in remaining planes

              yzthit = clusters1(numhit(jplth,jrr))
         if (yzthit.gt.thdw.and.yzthit.lt.thup) then

              ncand(jrem) = ncand(jrem) + 1   
*---- test on max.nr.of hits....................
              if (ncand(jrem).ge.nmhit) then
                print *,' Det.1--plane ',jplth,' too many hits seld.'
                go to 141
              endif

              ivv(jrem,ncand(jrem))= jrr     
           
          endif
          enddo                                   ! end of theta hits


*          print *,"jplth,jrem,ncand",jplth,jrem,ncand(jrem)
*------------------------------------------------------------
141       continue                                ! end of theta planes

*  now find best line using reference points ivref....lines to LIN_ALL.h.....
*
          nlin_pr = nlin_fd                      ! # of lines (0)
          
           numdet=101
          if(jfringe)numdet=201
          call find_best(numdet,lref,ivref,lrem,ncand,ivv)
                 
               if(li_m.gt.mxlin.and..not.last_chance)then
                last_chance=.true.
                rd_yy=rd_y(1)
                rd_zz=rd_z(1)
                rd_tth=rd_th(1)
                rd_y(1)=rd_y(1)/10.
                rd_th(1)=rd_th(1)/10.
                rd_z(1)=rd_z(1)/10.
                li_m=0
                go to 1000               
              
               elseif(li_m.gt.mxlin.and.last_chance)then
*          print *,'first telescope-NEW ',li_m,mxlin
                modset = .true.
                go to 99
               endif

*  in case of success ,do not try to reuse any hits...i.e.jump to 11.......
          if (nlin_fd.gt.nlin_pr) go to 11

          endif
41        continue                              ! endloop over Z in 2
40        continue                              ! endloop over Z (backward) (2)
31        continue                              ! endloop over clusters  Z in 1
30        continue                              ! endloop over pl 1 (Z)
         endif
21        continue                              ! endloop over hits Y in 2
20        continue                              ! endloop over Y (backward) (2)

11        continue                              ! endloop over clusters  Y in 1
10        continue                              ! endloop over pl 1 (Y)

99        continue

          nlin_fda = nlin_fd + nlin_fda
*          print*,'new2 ',nlin_fd,nlin_fda 
    

          li_m=0
          if(last_chance)then
                rd_y(1)=rd_yy
                rd_th(1)=rd_tth
                rd_z(1)=rd_zz
          endif
          end
