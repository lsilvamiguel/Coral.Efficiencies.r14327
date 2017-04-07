
          subroutine rec_det2(noas)

**********************************************************
*                                                        *
*             det.2   reconstruction                     *
*                                                        *                                            
*    New version for 2002 configuration (June 13,2002)   *
*    ---Modified for low momentum tracks---(R.W.)----    *
*                                                        *                                            
**********************************************************

C-kk --reverse order - Z first then Y - due to SM1 bending

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
      logical noas,last_chance


         last_chance=.false.
         nlin_fdb=0        
         noas= .false.
 
         if(nlin_fda.eq.0)then
         noas= .true.
         return
         endif

 2       continue
         eps_det=1.e-5

 1000    continue                               ! last chance 
         do jo=1,4
          lref(jo)=0
          ivref(jo)=0
         enddo

         do jo=1,ntpl+ntpl2
          ncand(jo)=0
          lrem(jo)=0

         do io=1,nmhit
          ivv(jo,io)=0
          jflg2(jo,io)=0
         enddo
         enddo


*   target radius 

      tar_rad = ytarg

*   center of telescope
    
      if(abs(det2_act(1)).lt.eps_det.
     &   and.abs(det2_act(2)).lt.eps_det)then
      x_mid2 = 0.5 *(xplan2(1) + xplan2(ntpl2))

*   take 10 % of length of telescope
      x_ten2 = ten2 *(xplan2(ntpl2) - xplan2(1))

      else
      x_mid2 = 0.5 *(det2_act(1)+ det2_act(2))

*   take 10 % of length of telescope
      x_ten2 = ten2 *(det2_act(2) - det2_act(1))
      
      endif
* --- end of center 

 
*   target limits.........................
      x_st = xtarg0 - hltarg
      x_nd = xtarg0 + hltarg
*
*    Nr of lines found so far .....................
      nlin_fd = 0
*
*    flag usable hits.............................
   
*

*
*   z segment.............................
*


      do 10 kz1 = 1,nzpl2-1

       jplz1 = jzpl2(kz1)

*  leave if center of telescope is reached
       if (xplan2(jplz1).gt.(x_mid2+x_ten2)) go to 99

       if(jc2hits(jplz1).eq.0)go to 10
       do 11 jhz = 1,jc2hits(jplz1)

       if (jflg2(jplz1,jhz).ne.0)  go to 11

       z1 = clusters2(numht2(jplz1,jhz))


        do 20 kz2 = nzpl2,kz1+1,-1

        jplz2 = jzpl2(kz2)

      if(ityp2(jplz1).eq.555.and.ityp2(jplz2).eq.555.              ! to protect same straw module
     &   and.modul(jplz1).eq.modul(jplz2))go to 20 

*  stop at center of telescope
        if (xplan2(jplz2).lt. (x_mid2 - x_ten2))  go to 11
*------- take into account target radius.............................

        if (z1.gt.0.0)  then

          zup  = (xplan2(jplz2)-x_nd)/(xplan2(jplz1)-x_nd) * z1
     1 + (xplan2(jplz2)-xplan2(jplz1))/(xplan2(jplz1)-x_nd) * tar_rad
          zdw  = (xplan2(jplz2)-x_st)/(xplan2(jplz1)-x_st) * z1
     1 - (xplan2(jplz2)-xplan2(jplz1))/(xplan2(jplz1)-x_st) * tar_rad

        else

          zup = (xplan2(jplz2)-x_st)/(xplan2(jplz1)-x_st) * z1
     1 + (xplan2(jplz2)-xplan2(jplz1))/(xplan2(jplz1)-x_st) * tar_rad
          zdw = (xplan2(jplz2)-x_nd)/(xplan2(jplz1)-x_nd) * z1
     1 - (xplan2(jplz2)-xplan2(jplz1))/(xplan2(jplz1)-x_nd) * tar_rad
  
        endif
        zup = zup +  rd_z(2)
        zdw = zdw -  rd_z(2)
        ncnt_z2 = 0

         if(jc2hits(jplz2).eq.0)go to 20

         do 21 khz = 1,jc2hits(jplz2)

         z2 = clusters2(numht2(jplz2,khz))

         if (jflg2(jplz2,khz).ne.0) go to 21


         if (z2.gt.zdw.and.z2.lt.zup) then

*  z proj. is OK .................................................

         indpl(jplz1) = jhz                
         indpl(jplz2) = khz               
         zpr = (z2-z1)/(xplan2(jplz2)-xplan2(jplz1))
         ncnt_z2 = ncnt_z2 + 1
         lref(1) = jplz1
         lref(2) = jplz2
         ivref(1) = jhz
         ivref(2) = khz 

*  now loop over Y planes.........................................

     
      do 30 ky1 = 1,nypl2-1

       jply1 = jypl2(ky1) 

*---- require same type as jplz1

        if (ityp2(jply1).ne.ityp2(jplz1)) go to 30

*---- require Z1, Y1 not too far away..............
     
        if (abs(xplan2(jply1)-xplan2(jplz1)).gt.distm2_yz) go to 30

*  leave when center of telescope is reached
       if (xplan2(jply1).gt. (x_mid2+x_ten2)) go to 21
     
       if(jc2hits(jply1).eq.0)go to 30
       do 31 jhy = 1,jc2hits(jply1)

       if (jflg2(jply1,jhy).ne.0) go to 31
       y1 = clusters2(numht2(jply1,jhy))

*
*--- Check compatibility between y1 and z1 using nearest theta plane......
*
         xy1z1 = 0.5 * (xplan2(jply1) + xplan2(jplz1))
         iok1 = 0

       do 440 jth1 = 1,nthpl2
         jplth1 = jthpl2(jth1)
         if (abs(xplan2(jplth1)-xy1z1).gt.distm2_th) go to 440
         if (ityp2(jplth1).ne.ityp2(jplz1))  go to 440
         tth1 = sinn2(jplth1)*y1 - cosn2(jplth1)*z1

        if(jc2hits(jplth1).eq.0)go to 440
        do jrr = 1,jc2hits(jplth1)
          thhit = clusters2(numht2(jplth1,jrr))
 
          if (abs(thhit-tth1).lt.rd_com2) iok1 = 1
        enddo
 440     continue

*
*---
*
       if (iok1.eq.0)  go to 31

       do 40 ky2 = nypl2,ky1+1,-1

      jply2 = jypl2(ky2)

      if(ityp2(jply1).eq.555.and.ityp2(jply2).eq.555.              ! to protect same straw module
     &   and.modul(jply1).eq.modul(jply2))go to 40 

*----require same typ as jplz2......................

       if (ityp2(jply2).ne.ityp2(jplz2))  go to 40

*----require (Z,Y) not too far away (max. 2 cm)..................

       if (abs(xplan2(jply2)-xplan2(jplz2)).gt.distm2_yz) go to 40
*  stop at center of telescope 
      if (xplan2(jply2).lt. (x_mid2 - x_ten2)) go to 31

*----- no vertex pointing, since direction has been modified in M1 .....
        ncnt_y2 = 0

         if(jc2hits(jply2).eq.0)go to 40
         do 41 khy = 1,jc2hits(jply2)

         if (jflg2(jply2,khy).ne.0) go to 41
         y2 = clusters2(numht2(jply2,khy))

*
*---  check compatibility of Y2 and Z2 using nearest theta plane
*

       xy2z2 = 0.5 * (xplan2(jply2) + xplan2(jplz2))
       iok2 = 0


      do 330  jth2 = 1,nthpl2
       jplth2 = jthpl2(jth2)
       if (abs(xplan2(jplth2)-xy2z2).gt.distm2_th) go to 330
       if (ityp2(jplth2).ne.ityp2(jplz2)) go to 330
       ttth2 = sinn2(jplth2)*y2 - cosn2(jplth2)*z2

        if(jc2hits(jplth2).eq.0)go to 330
        do jrr = 1,jc2hits(jplth2)
         thhit = clusters2(numht2(jplth2,jrr))
         if (abs(thhit-ttth2).lt.rd_com2) iok2 = 1
        enddo
 330   continue
     
*
*---
*
       if (iok2.eq.0) go to 41

*  y proj. is OK  anyway..........................................
   
         indpl(jply1) = jhy               
         indpl(jply2) =  khy
         ypr = (y2-y1)/(xplan2(jply2)-xplan2(jply1))
         ncnt_y2 = ncnt_y2 + 1
         lref(3) = jply1
         lref(4) = jply2
         ivref(3) = jhy
         ivref(4) = khy


*  now define roads for theta  planes and remaining Y,Z planes................
         jrem = 0

         yy0 = clusters2(numht2(jply1,indpl(jply1)))
         zz0 = clusters2(numht2(jplz1,indpl(jplz1)))

         do 141 jplth = 1,ntpl2

      if (jplth.eq.lref(1).or.jplth.eq.lref(2).or.jplth.eq.lref(3)
     1    .or.jplth.eq.lref(4))go to 141            


*----type must correspond to type of one of the ref. points............
 

         jrem = jrem + 1
         lrem(jrem) = jplth
         ncand(jrem) = 0

         if (ityp2(jplth).ne.ityp2(jplz2).and.             
     1      ityp2(jplth).ne.ityp2(jplz1)) go to 141

         yyth = yy0 + ypr * (xplan2(jplth)-xplan2(jply1))
         zzth = zz0 + zpr * (xplan2(jplth)-xplan2(jplz1))
         ttth = sinn2(jplth) * yyth - cosn2(jplth) * zzth

         
*----in bending plane, include max. curvature in road.................

         rd_cor = sagi_tot
         if (abs(xplan2(jplz2)-x_mid).lt.x_ten.or.
     1    abs(xplan2(jplz1)-x_mid).lt.x_ten) rd_cor = 2 * sagi_tot

         yyup = yyth + rd_cor
         yylw = yyth - rd_cor
         zzup = zzth + rd_th(2)
         zzlw = zzth - rd_th(2)
         th1 = sinn2(jplth) * yyup - cosn2(jplth) * zzup
         th2 = sinn2(jplth) * yyup - cosn2(jplth) * zzlw
         th3 = sinn2(jplth) * yylw - cosn2(jplth) * zzup
         th4 = sinn2(jplth) * yylw - cosn2(jplth) * zzlw
         thup = max(th1,th2,th3,th4)
         thdw = min(th1,th2,th3,th4)

*  select theta hits within roads............................................

              if(jc2hits(jplth).eq.0)then
                ncand(jrem)=-100                  ! flag, empty plane 
                go to 141
              endif  

              do jrr=1,jc2hits(jplth)
                        
              yztht2 = clusters2(numht2(jplth,jrr))

         if (yztht2.gt.thdw.and.yztht2.lt.thup) then

              ncand(jrem) = ncand(jrem) + 1   
              if (ncand(jrem).ge.nmhit) then
                print *,' Det. 2.. plane ',jplth,' too many hits seld.'
                go to 141
              endif
              ivv(jrem,ncand(jrem))= jrr

          endif
          enddo

141       continue
*  now find best line using reference points ivref.........................
*
     

          nlin_pr = nlin_fd

          call find_best(2,lref,ivref,lrem,ncand,ivv)

               if(li_m.gt.mxlin.and..not.last_chance)then
                last_chance=.true.
                rd_zz=rd_z(2)
                rd_tth=rd_th(2)
                rd_th(2)=rd_th(2)/10.
                rd_z(2)=rd_z(2)/10.
                li_m=0
                go to 1000               
              
               elseif(li_m.gt.mxlin.and.last_chance)then
*          print *,'second telescope ',li_m,mxlin

                go to 99
               endif
         

*
*  in case of success ,do not try to reuse any hits...i.e.jump to 11.......
          if (nlin_fd.gt.nlin_pr) go to 11
*
41        continue
40        continue
31        continue
30        continue
          endif
21        continue
20        continue
11        continue
10        continue

99        continue
         
          nlin_fdb = nlin_fd + nlin_fdb

          if(nlin_fdb.eq.0)noas = .true.
*          print*,'det2 ',nlin_fdb

          if(last_chance)then
                rd_th(2)=rd_tth
                rd_z(2)=rd_zz
          endif

          return
          end




