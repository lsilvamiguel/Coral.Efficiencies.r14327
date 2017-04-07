
      subroutine rec_det1_old
*************************************************************
*                                                           *
*        det.1    reconstruction (old version)              *
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
      logical imodtg,jmodtg
     
         ncth=ncth1

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
  

*   target limits.........................
      x_st = xtarg0 - hltarg
      x_nd = xtarg0 + hltarg

*    Nr of lines found so far .....................
      nlin_fd = 0

*    flag usable hits.............................


C-kk -- to find line in first aprox. 4 pl are needed: 2 Y and 2 Z
  

*   y segment.............................


      do 10 ky1 = 1,nypl-1                       ! loop over Y planes
       jply1 = jypl(ky1)                         ! real # pl w det. 1 (ref.1)

      if(jchits(jply1).eq.0)go to 10
        do 11 jhy = 1,jchits(jply1)              ! loop over hits in pl 1


       if (jflag(jply1,jhy).ne.0)  go to 11 
      
C-kk -- if hit was used (flaga ne 0) go out
                                      
       y1 = clusters1(numhit(jply1,jhy))
             

        do 20 ky2 = nypl,ky1+1,-1         

C-kk --  backward loop from last to first

                           
        jply2 = jypl(ky2)                         ! real # pl w det. 1 (ref.2)


*    modify target pointing for Mmegas if high multiplicity............

        jmodtg = .false.
      if(jchits(jply2).gt.ncth.and.jchits(jply1).gt.ncth) then
*           print *, " modified tg.pointing Y planes",jply2,jply1
           jmodtg = .true.
        endif

c-kk -- target shot
        y2_e = (xplane(jply2)-x_nd)/(xplane(jply1)-x_nd) * y1
        y2_s = (xplane(jply2)-x_st)/(xplane(jply1)-x_st) * y1

        if (y2_e.lt.y2_s) then
           yup = y2_s
           ydw = y2_e
        else
           yup = y2_e
           ydw = y2_s
        endif

        yup = yup + rd_y(1)
        ydw = ydw - rd_y(1)

*            print*,'target estimation ',yup,ydw

        ncnt_y2 = 0                      
  

*------- find optimal y1 if modified target pointing...........

        if (jmodtg) then
         dtgmin = 9999.
         jtgopt = -99
         ypred = y1 * (xplane(jply2) - xtarg0)/(xplane(jply1) - xtarg0)

         if(jchits(jply2).eq.0)go to 20
         do khy = 1,jchits(jply2)
         y2 = clusters1(numhit(jply2,khy))
         if (jflag(jply2,khy).eq.0) then  
         dtg = abs(y2 - ypred)
         if (dtg.lt.dtgmin) then
            dtgmin = dtg
            jtgopt = khy
         endif
         endif
         enddo
        endif

         do 21 khy = 1,jchits(jply2)       ! loop over hits in Y in pl 2

          y2 = clusters1(numhit(jply2,khy))

         if (jflag(jply2,khy).ne.0) go to 21     !  repetition flag

                
         if (y2.gt.ydw.and.y2.lt.yup) then

*  y proj. is OK .................................................
  
* ------- more restrictive if modified Target pointing.....

         if (jmodtg.and.jply2.ne.jtgopt)  go to 21
 

         indpl(jply1) = jhy      ! index = hit in pl 1
         indpl(jply2) = khy      ! index = hit in pl 2
        
         ypr = (y2-y1)/(xplane(jply2)-xplane(jply1))
C-kk  --     ypr - slope from Y ref.planes  

         ncnt_y2 = ncnt_y2 + 1    
         lref(1) = jply1         
         lref(2) = jply2
         ivref(1) = jhy          
         ivref(2) = khy 

*  now loop over Z planes.........................................

      do 30 kz1 = 1,nzpl-1
       jplz1 = jzpl(kz1)

       if(jchits(jplz1).eq.0)go to 30

       do 31 jhz = 1,jchits(jplz1)
     

       if (jflag(jplz1,jhz).ne.0) go to 31        ! repetition flag
       z1 = clusters1(numhit(jplz1,jhz))


        do 40 kz2 = nzpl,kz1+1,-1
        jplz2 = jzpl(kz2)

*    modify target pointing for Mmegas if high multiplicity............

        imodtg = .false.
      if(jchits(jplz2).gt.ncth.and.jchits(jplz1).gt.ncth) then
*           print *, " modified tg.pointing Z planes",jplz2,jplz1
           imodtg = .true.
        endif

        z2_e = (xplane(jplz2)-x_nd)/(xplane(jplz1)-x_nd) * z1
        z2_s = (xplane(jplz2)-x_st)/(xplane(jplz1)-x_st) * z1
        if (z2_e.lt.z2_s) then
           zup = z2_s
           zdw = z2_e
        else
           zup = z2_e
           zdw = z2_s
        endif
        zup = zup +  rd_z(1)
        zdw = zdw -  rd_z(1)

*------- find optimal z1 if modified target pointing...........

        if (imodtg) then
         dtgmin = 9999.
         itgopt = -99
         zpred = z1 * (xplane(jplz2) - xtarg0)/(xplane(jplz1) - xtarg0)

         if(jchits(jplz2).eq.0)go to 40

         do khz = 1,jchits(jplz2)
         z2 = clusters1(numhit(jplz2,khz))  
         if (jflag(jplz2,khz).eq.0) then 
         dtg = abs(z2 - zpred)
         if (dtg.lt.dtgmin) then
            dtgmin = dtg
            itgopt = khz
         endif
         endif
         enddo
        endif

        ncnt_z2 = 0


         do 41 khz = 1,jchits(jplz2)

             z2 = clusters1(numhit(jplz2,khz))

         if (jflag(jplz2,khz).ne.0) go to 41    ! repetition flag


         if (z2.gt.zdw.and.z2.lt.zup) then

*   more selective if modified target pointing......

         if (imodtg.and.khz.ne.itgopt) go to 41

*  z proj. is OK .................................................
         indpl(jplz1) = jhz               
         indpl(jplz2) =  khz
         zpr = (z2-z1)/(xplane(jplz2)-xplane(jplz1))
         ncnt_z2 = ncnt_z2 + 1
         lref(3) = jplz1
         lref(4) = jplz2
         ivref(3) = jhz
         ivref(4) = khz

C-kk -- now looking for hits in all remainning planes fitted do line between 1 and 2

*  now define roads for theta  planes and remaining Y,Z planes................
         jrem = 0

         yy0 = clusters1(numhit(jply1,indpl(jply1)))
         zz0 = clusters1(numhit(jplz1,indpl(jplz1)))


         do 141 jplth = 1,ntpl                       ! loop over remaining planes
            
                              
         if (jplth.eq.lref(1).or.jplth.eq.lref(2).or.jplth.eq.lref(3)
     1    .or.jplth.eq.lref(4)) go to 141

C-kk -- if any before go to 141

         jrem = jrem + 1
         lrem(jrem) = jplth
         ncand(jrem) = 0
         yyth = yy0 + ypr * (xplane(jplth)-xplane(jply1))
         zzth = zz0 + zpr * (xplane(jplth)-xplane(jplz1))
         ttth = sinn(jplth) * yyth - cosn(jplth) * zzth

         thup = ttth +  rd_th(1)
         thdw = ttth -  rd_th(1)


*  select theta hits within roads............................................
          if(jchits(jplth).eq.0)then
             ncand(jrem)=-100                    ! flag, empty plane
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
          enddo                                   ! end lopp over theta hits

141       continue                                ! end of theta pl loop

*  now find best line using reference points ivref.........................
*
          nlin_pr = nlin_fd                      ! # of lines found (0)


          call find_best(101,lref,ivref,lrem,ncand,ivv)

               if(li_m.gt.mxlin)then
*                 print*,'first telescope-OLD ',li_m,mxlin
                 go to 99
               endif

*
*  in case of success ,do not try to reuse any hits...i.e.jump to 11.......
          if (nlin_fd.gt.nlin_pr) go to 11

          endif
41        continue                              ! endloop over Z  (2)
40        continue                              ! endloop over Z (backward) (2)
31        continue                              ! endloop over pl 1 (Z)
30        continue                              ! end loop over Z planes 1
         endif
21        continue                              ! end loop over hits  Y in 2
20        continue                              ! end loop over Y planes (backward loop)

11        continue                              ! end loop over hits  Y in 1
10        continue                              ! end loop over Y planes 1

 99       continue

          nlin_fda = nlin_fda+nlin_fd
*          print*,'old ',nlin_fd,nlin_fda,li_m 
        
 
          return
          end


