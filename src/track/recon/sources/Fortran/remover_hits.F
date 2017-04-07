      subroutine remover_hits(ip,miter1,miter2,mlassc,lassoc,
     & listhit1,listhit2)

********************************************************************
*                                                                  *
*    routine removes hits from hit's table for reconstructed       *
*                tracks found in previous iterations               *
*               ip = 0 - all used hits removed                     *
*               ip = 1 - hits used behind removed                  *
*               ip = 2 - hits used before removed                  *
*               ip = 100 - special mode - removes segments before  *
*                          SM1 found in "mode"   iterations        *
*                                                                  *
********************************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "HITS_OUT.h"
#include "HITS_IN.h"

      
            dimension lassoc(2,*)
            dimension listhit1(ntrack,*),listhit2(ntrack,*)


          if(mlassc.eq.0.and.ip.lt.100)return    ! nothing to be removed

*
*-------------------------------------* first step - mark used clusters
*


          if(ip.eq.100) then                     ! removing hits before SM1 after first mode iterations      
             do l=miter1,miter2
              do k=1,ntpl
 
              if(jchits(k).ne.0)then                   
              do j=1,jchits(k)         
          if(numhit(k,j).eq.listhit1(l,k))numhit(k,j)=8888         !just flag for marking clusters 
              enddo
              endif

              enddo         
            enddo         

          elseif(ip.lt.100)then                  ! normal procedure for marking clusters

             do i =1,mlassc
               im=lassoc(1,i)
               in=lassoc(2,i)

              if(ip.eq.1)go to 1                  ! hits before will not be removed        

              do k=1,ntpl
               if(jchits(k).ne.0)then               
              do j=1,jchits(k) 
        if(numhit(k,j).eq.line1(k,im))numhit(k,j)=8888
              enddo
               endif
              enddo

 1            continue
         
              if(ip.eq.2)go to 2                  ! hits behind will not be removed
 
             if(in.ne.2000)then                   ! no fringe field reconstruction
              do k=1,ntpl2               
               if(jc2hits(k).ne.0)then           
              do j=1,jc2hits(k)
        if(numht2(k,j).eq.line2(k,in))numht2(k,j)=8888
              enddo
               endif 
              enddo
      
             endif
 2           continue

             enddo                               ! end loop over tracks
           endif

*
*----------- used clusters just marked - second step -> removing procedure
*

           if(ip.eq.1)go to 11 
 
            do k=1,ntpl             
              jch=0                
            if(jchits(k).ne.0)then

              do j=1,jchits(k)               
               if(numhit(k,j).ne.8888)then               
                jch=jch+1
                numhit(k,jch) = numhit(k,j)              
               endif              
              enddo

               jchits(k) = jch

             endif
             enddo           

           if(ip.eq.2)return

 11          continue
            do k=1,ntpl2
              jch=0  
            if(jc2hits(k).ne.0)then 
                        
             do j=1,jc2hits(k)
              if(numht2(k,j).ne.8888)then      
                jch=jch+1
                numht2(k,jch) = numht2(k,j)               
              endif
             enddo

              jc2hits(k) = jch

             endif
             enddo
 
             return
             end









