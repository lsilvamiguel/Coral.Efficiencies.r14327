      subroutine rec_det1_mgmg(modset,jfringe)
*************************************************************
*                                                           *
*                det.1    reconstruction                    *
*    Nov.14, 6.30 pm                                        *
*    modified Aug.31, 2000 by kk                            *              
*------- modified May 7,2001 for reconstruction in MM-MM    *
*        ref. planes are (2 u) and (2 v) at +/- 45 degrees  *
*        association checked using 2 y planes               *
*         (that is all we have this year!!)                 *
*                                                           *
* modified in March 19 2002 by kk                           *
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
      logical modset,jfringe

        last_chance=.false.
        ncth=ncth1
        nlin_fda=0
     

 1000    continue                                    ! last chance
         if(last_chance) ncth=int(ncth1/2)
         
         eps_det =1.e-1

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
    
      if(abs(det1_act(1)).lt.eps_det.and.
     &   abs(det1_act(2)).lt.eps_det)then

      x_mid = 0.5 *(xplane(1) + xplane(ntpl))

*   take (ten3 e.g.10 %) of length of telescope
      x_ten = ten3 *(xplane(ntpl) - xplane(1))

      else
      x_mid = 0.5 *(det1_act(1)+ det1_act(2))

*   take (ten3 e.g.10 %) of length of telescope
      x_ten = ten3 *(det1_act(2) - det1_act(1))
      
      endif
* --- end of center 

*   target limits.........................
      x_st = xtarg0 - hltarg
      x_nd = xtarg0 + hltarg

*    Nr of lines found so far .....................
      nlin_fd = 0

*    flag usable hits.............................

*   z segment.............................


      do 10 kz2 = nthpl,2,-1                       ! loop over Z planes (theta here!!)

       jplz2 = jthpl(kz2)                          ! real # of pl before SM1
   
       if (jfringe.and.ityp1(jplz2).ne.222)  go to 10
       if (jfringe.and.jchits(jplz2).gt.2*ncth1) go to 10

*   leave if center of telescope is reached             
       if (xplane(jplz2).lt.(x_mid - x_ten)) go to 99

        do 11 jhz = 1,jchits(jplz2)         
**----These are now theta hits.........("v").......................
*          print*,jchits(jplz2)

       if (jflag(jplz2,jhz).ne.0)  go to 11 
      
         z2 = clusters1(numhit(jplz2,jhz))

        do 20 kz1 = 1,kz2-1         
                  
        jplz1 = jthpl(kz1)                         ! real # of pl before SM1

** -----check jplz1 has same orientation as jplz2 .......("v")
            
        if (abs(sinn(jplz2)-sinn(jplz1)).gt.eps_det.or.
     1  abs(cosn(jplz2)-cosn(jplz1)).gt.eps_det)  go to 20          
                   
        if (jfringe.and.jchits(jplz1).gt.2*ncth1) go to 20          

*   stop loop at center of telescope
        if (xplane(jplz1).gt.(x_mid + x_ten)) go to 11
        
*   modify target pointing for Mmegas if high multiplicity............

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

        zup = zup + sqrt(rd_y(1)**2 + rd_z(1)**2)
        zdw = zdw - sqrt(rd_y(1)**2 + rd_z(1)**2)

*------- find optimal z1 if modified target pointing...........

        if (imodtg) then
         dtgmin = 9999.
         itgopt = -99
         zpred = z2 * (xplane(jplz1) - xtarg0)/(xplane(jplz2) - xtarg0)

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
**---remember these are actually theta hits...("v" plane)
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


      do 30 ky2 = nthpl,2,-1

       jply2 = jthpl(ky2)

***---Now check that direction is orthogonal to (z1 or z2).........
       if (abs(cosn(jply2)*cosn(jplz2)+sinn(jply2)*sinn(jplz2)).
     1 gt.eps_det)   go to 30
        
       if (jfringe.and.ityp1(jply2).ne.222)  go to 30
       if (jfringe.and.jchits(jply2).gt.2*ncth1) go to 30       
      
*  leave when center of telescope is reached
          
       if (xplane(jply2).lt.(x_mid - x_ten)) go to 21
          
        
       do 31 jhy = 1,jchits(jply2)
     
*  remember these are actually theta hits  ("u" direction)
       if (jflag(jply2,jhy).ne.0) go to 31         ! repetition flag
       y2 = clusters1(numhit(jply2,jhy))


        do 40 ky1 = 1,ky2-1       

        jply1 = jthpl(ky1)

**---check this plane has same orientation as jply2....

       if (abs(sinn(jply2)-sinn(jply1)).gt.eps_det.or.
     1  abs(cosn(jply2)-cosn(jply1)).gt.eps_det) go to 40
    
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
        yup = yup + sqrt(rd_y(1)**2 + rd_z(1)**2)
        ydw = ydw - sqrt(rd_y(1)**2 + rd_z(1)**2)

*------- find optimal y1 if modified target pointing...........

        if (jmodtg) then
         dtgmin = 9999.
         jtgopt = -99
         ypred = y2 * (xplane(jply1) - xtarg0)/(xplane(jply2) - xtarg0)

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

** Remember these are actually theta hits.......("u" direction)

         if (jflag(jply1,khy).ne.0) go to 41    ! repetition flag


         if (y1.gt.ydw.and.y1.lt.yup)then

*  y proj. is OK .................................................


*   more selective if modified target pointing......

         if (jmodtg.and.khy.ne.jtgopt) go to 41

         indpl(jply2) = jhy      ! index equal to hit in pl 1
         indpl(jply1) = khy      ! index equal to hit in pl 2

*-kk  --     ypr - slope from Y ref.planes          
         ypr = (y2-y1)/(xplane(jply2)-xplane(jply1))

  
         ncnt_y2 = ncnt_y2 + 1    !  index of "good" rec.
         lref(1) = jply1         
         lref(2) = jply2
         ivref(1) = khy           ! hit # in 1 and 2
         ivref(2) = jhy 

*  now define roads for theta  planes and remaining Y,Z planes................
         jrem = 0

         yy0 = clusters1(numhit(jply1,indpl(jply1)))
         zz0 = clusters1(numhit(jplz1,indpl(jplz1)))

         do 141 jplth = 1,ntpl                 ! loop over remaining pl in tel.1

         if (jplth.eq.lref(1).or.jplth.eq.lref(2).or.jplth.eq.lref(3)
     1    .or.jplth.eq.lref(4))go to 141


         jrem = jrem + 1
         lrem(jrem) = jplth
         ncand(jrem) = 0
         uuth = yy0 + ypr * (xplane(jplth)-xplane(jply1))
         vvth = zz0 + zpr * (xplane(jplth)-xplane(jplz1))

**Now recalculate true (Y,Z) from (u,v) pair...................
         yyth = sinn(jply1) * uuth + cosn(jply1) * vvth
         zzth = -cosn(jply1)* uuth + sinn(jply1) * vvth


** and finally prediction in additionnal plane........
         ttth = sinn(jplth) * yyth - cosn(jplth) * zzth

         thup = ttth +  rd_th(1)
         thdw = ttth -  rd_th(1)


*  select theta hits within roads............................................
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

141       continue                                ! end of theta planes

*  now find best line using reference points ivref....lines to LIN_ALL.h.....
*
          nlin_pr = nlin_fd                      ! # of lines (0)
                           
          call find_best(101,lref,ivref,lrem,ncand,ivv)

               
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
*          print*,'micros  ',nlin_fd,nlin_fda

          li_m=0
          if(last_chance)then
                rd_y(1)=rd_yy
                rd_th(1)=rd_tth
                rd_z(1)=rd_zz
          endif
          end
  





