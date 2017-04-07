         subroutine check_hit(mlassc,jjflags,ntori,
     & listhit1,listhit2,listhit3)
                      

**************************************************
*                                                *
*           for checking how many hits           *
*                                                *
**************************************************

#include "PARAM.h"
#include "PLANES.h"

       dimension listhit1(ntrack,*),listhit2(ntrack,*),
     & listhit3(ntrack,*)  

       dimension jjflags(*),ntori(*)
       logical sub,ff
  
        
            ff=.false.
            if(mlassc.eq.0)return

              ntor=ntori(1)
            do li=1, mlassc            ! found tracks loop
             if(ntori(li).lt.100.and.ntori(li).lt.ntor)then
                print*,'ERROR in logic in check-hit'
             endif   
             ntor=ntori(li)   

         if(sub)then
         if(ntori(li).le.mp1)then
          ntpl_t=ntpl_s1
          ntpl2_t=ntpl2_s1
         elseif(ntori(li).gt.mp1.and.ntori(li).le.mp1+mp2)then
          ntpl_t=ntpl_s1
          ntpl2_t=ntpl2_s2
         elseif(ntori(li).gt.mp1+mp2)then
          ntpl_t=ntpl_tmp
          ntpl2_t=ntpl2_tmp
         endif
         elseif(.not.sub)then
          ntpl_t=ntpl
          ntpl2_t=ntpl2
         endif

              if(ntori(li).gt.100)then
                  ntori(li)=ntori(li)-100      ! fringe field flag
                  ff=.true.
              endif

              jjflags(li)=0
              lihit1=0
              lihit2=0
              lihit3=0
                                 
               do k=1,ntpl_t
       if(listhit1(li,k).ne.8888.and.listhit1(li,k).ne.0) 
     &     lihit1=lihit1+1                    
               enddo
  
                  
*----- now check if fringe field/full spectrometer reconstruction

       if(ff.and.listhit2(li,1).lt.9999)then
           print*,'error in logic - listhi2 in ff'
           stop
       endif

       if(listhit2(li,1).eq.9999)then            ! fringe field reconstruction - nobridged
               jjflags(li) =-lihit1
               listhit2(li,1)=8888
               lihit2=0
               lihit3=0
               go to 30  
       elseif(listhit2(li,1).gt.9999)then        ! fringe field reconstruction - bridged 
               listhit2(li,1)=listhit2(li,1)-9999
               do k=1,ntpl2_t 
       if(listhit2(li,k).ne.8888.and.listhit2(li,k).ne.0)
     & lihit2=lihit2+1                   
               enddo  
       jjflags(li) =-(lihit1+ 1000*lihit2+200000)        ! always < -200000
               lihit3=0
               go to 30  

       endif

*------ end of check

  
  
               do k=1,ntpl2_t 
       if(listhit2(li,k).ne.8888.and.listhit2(li,k).ne.0)
     & lihit2=lihit2+1                   
               enddo 


       if(listhit3(li,1).eq.9999)then
              lihit3=0
          listhit3(li,1)=8888
       jjflags(li) =-(lihit1+ 1000*lihit2+100000)        ! always < -100000
              go to 30
       endif
      
              do k=1,ntpl_sup
       if(listhit3(li,k).ne.8888.and.listhit3(li,k).ne.0)
     &  lihit3=lihit3+1
               enddo 

      jjflags(li) =lihit1+ 1000*lihit2+100000*lihit3

  30           continue


*             print*,'debug1 ',li,lihit1              
             enddo
          
             return
             end


         subroutine check_hit_mc(mlassc,jjflags,phit,
     & ntori,listhit1,listhit2,listhit3)
                      
**************************************************
*                                                *
*           for checking how many hits           *
*         and from which generated lines         *
*         belongs to reconstructed line:         *
*         phit is a criteria parameter in %      *
*                                                *
**************************************************

#include "PARAM.h"
#include "MC.h"
#include "HITS_IN.h"
#include "PLANES.h"

        common/local/lipa,itime

       dimension listhit1(ntrack,*),listhit2(ntrack,*),
     & listhit3(ntrack,*)  
 
        dimension lik(ntrack),liks(ntrack)
        dimension jjflags(*),ntori(*)
        dimension ink(mxpl)

        logical fri,extra,sub

        data itime/1/


        nto=0
        if(sub)then
         nto=ntpl_cor
         nto2=ntpl2_cor
        endif

            if(mlassc.eq.0)return            

            if(itime.eq.1)then
              lipa=0
              itime=0
            endif

            eps=1.e-5

             if(ntr.eq.0)then
              jjflag(1,1)= -10000      ! flag
              return
             endif

             if(ntr.gt.ntrack.or.mlassc.gt.ntrack)then
          print*,'errors in dimension - ntrack too small'
             stop
             endif

*--- real MC comparison
            do li=1, mlassc            ! found tracks loop  
         
             if(sub)then
         if(ntori(li).le.mp1)then
          ntpl_t=ntpl_s1
          ntpl2_t=ntpl2_s1

          nt=0
          do k=1,ntpl_tmp
          if(ityp1(nto+k).eq.111)then
             nt=nt+1
             ink(nt)=k
          endif
          enddo
          do k=1,ntpl2_tmp
          if(ityp2(nto2+k).eq.222.or.ityp2(nto2+k).eq.555)then
             nt=nt+1
             ink(ntpl_t+nt)=k
          endif
          enddo

         elseif(ntori(li).gt.mp1.and.ntori(li).le.mp1+mp2)then
          ntpl_t=ntpl_s1
          ntpl2_t=ntpl2_s2
          nt=ntpl_t

          do k=1,ntpl2_tmp
          if(ityp2(nto2+k).eq.444)then
             nt=nt+1
             ink(ntpl_t+nt)=k
          endif
          enddo
 
         elseif(ntori(li).gt.mp1+mp2)then
          ntpl_t=ntpl_tmp
          ntpl2_t=ntpl2_tmp
         endif
         elseif(.not.sub)then
          ntpl_t=ntpl
          ntpl2_t=ntpl2
         endif
        
             do io=1,5
               jjflag(li,io)=0
             enddo
            
             do io=1,ntr
               lik(io)=0
               liks(io)=0
             enddo
                      
            if(jjflags(li).lt.-100000)then               ! no extrapolation
               extra=.false.
               fri=.false.
               jjfl=-jjflags(li)-100000
               lihit3=0
               lihit1=mod(jjfl,1000)                   
               lihit2=(jjfl-lihit1)/1000
               lihit=lihit1+lihit2
            elseif(jjflags(li).lt.0.and.jjflags(li).gt.-100000)then  ! fringe field reconstruction
               fri=.true.
               extra=.false.
               lihit1=-jjflags(li)
               lihit2=0
               lihit3=0
               lihit=lihit1
            else                                          ! extrapolation ok
               fri=.false.
               extra=.true.
               lihit1=mod(jjflags(li),1000)         
               lihit2=mod(jjflags(li),100000)/1000
               lihit3=jjflags(li)/100000
               lihit=lihit1+lihit2+lihit3
            endif
              
                if(lihit1.lt.n_tel1.and.lihit2.lt.n_tel2) then
                    jjflag(li,3)=-1
                elseif(lihit2.lt.n_tel2) then
                    jjflag(li,3)=-2
                elseif(lihit1.lt.n_tel1) then
                    jjflag(li,3)=-1
                endif

*             print*,li,lihit1,lihit2,lihit3

                 
               do k=1,ntpl_t
                 
               k1=k
               if(sub.and.ntori(li).le.mp1+mp2)k1=ink(k)
               if(listhit1(li,k).ne.8888)then

               do jj=1,jchits_mc(k1)
               
               if(l_hit1(1,k1,jj).gt.0)then
               do jop=1,l_hit1(1,k1,jj)

         if(numhit_mc(k1,jj).eq.listhit1(li,k).and.
     &   l_hit1(jop+1,k1,jj).gt.0.and.l_hit1(jop+1,k1,jj).ne.200
     &   .and.label(l_hit1(jop+1,k1,jj)).ne.-1)   
     &  then

          lik(l_hit1(jop+1,k1,jj))=lik(l_hit1(jop+1,k1,jj))+1

*           print*,li,k,l_hit1(jop+1,k1,jj),lik(l_hit1(jop+1,k1,jj))

            endif

             
               enddo
             endif
               enddo

              endif
               enddo

    
*          print*,'before',li,(i,lik(i),i=1,ntr)
   

           if(.not.fri)then

               do k=1,ntpl2_t
               k2=k
               if(sub.and.ntori(li).le.mp1+mp2)k2=ink(ntpl_t+k)
                       
               if(listhit2(li,k).ne.8888)then
 
               do jj=1, jc2hits_mc(k2)

               if(l_hit2(1,k2,jj).gt.0)then
               do jop=1,l_hit2(1,k2,jj)

      if(numht2_mc(k2,jj).eq.listhit2(li,k).and.
     &  l_hit2(jop+1,k2,jj).gt.0.and.l_hit2(jop+1,k2,jj).ne.200
     &  .and.label(l_hit2(jop+1,k2,jj)).ne.-1) 
     &  then

            lik(l_hit2(jop+1,k2,jj))=lik(l_hit2(jop+1,k2,jj))+1
     
*         print*,li,l_hit2(jop+1,k2,jj),lik(l_hit2(jop+1,k2,jj))

             endif

               enddo
             endif
               enddo


               endif
               enddo

         
            if(extra)then

             do k=1,ntpl_sup                       
               if(listhit3(li,k).ne.8888)then 

               do jj=1, jcshits(k)

              if(l_hit3(1,k,jj).gt.0)then
               do jop=1,l_hit3(1,k,jj)

      if(numhts(k,jj).eq.listhit3(li,k).and.
     &  l_hit3(jop+1,k,jj).gt.0.and.l_hit3(jop+1,k,jj).ne.200
     &  . and.label(l_hit3(jop+1,k,jj)).ne.-1) 
     &   then
         
          liks(l_hit3(jop+1,k,jj))=liks(l_hit3(jop+1,k,jj))+1

*          print*,li,l_hit3(jop+1,k,jj),liks(l_hit3(jop+1,k,jj))


             endif
               enddo
             endif
               enddo

               endif
               enddo
             endif
           endif
          
          
*             print*,'after',li,(i,lik(i),i=1,ntr)
*--- now # of hits:
            
             alicz=0.
             hit=0.0
             l=0
             lika=0
             likat=0           
            
            if(ntr.gt.0)then           
            do j=1,ntr

             if(lihit.eq.0)go to 11 

             lika = lik(j)
 
             if(extra)lika=lik(j)+liks(j)

             if(float(lika)/float(lihit).gt.alicz)then                             
               alicz=float(lika)/float(lihit)
               l=j
               likat=lika
             endif

  11        continue
              enddo
            else
               alicz=0.0
            endif

                hit=alicz*100.
                jjflag(li,5)=int(hit)

                if((hit-100.0).gt.eps)then
                   lipa = lipa +1
                   print*,'error 100',lipa,li

       print *,li,l,hit,lihit,likat,lik(l),liks(l),lihit1,lihit2,lihit3
       print *,li,(jj,lik(jj),jj=1,ntr)
       print *,li,(jj,liks(jj),jj=1,ntr)
       print*
                  endif

*-------------- now true generated track
  
                  if(hit.ge.phit)then
                   jjflag(li,1)=0
                  if(fri)jjflag(li,1)=100                  ! fringe field rec. flag
                   jjflag(li,2)=l
                   jjflag(li,4)=liks(l)

                  else
                   jjflag(li,1)=-1
                  if(fri)jjflag(li,1)=-100                 ! fringe field rec. flag
                   jjflag(li,2)= 0 
                   jjflag(li,4)= 0
                  endif
  10           continue

             enddo
          
             return
             end







