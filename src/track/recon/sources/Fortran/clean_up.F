      subroutine clean_up(miter,mlassc,lassoc,
     & chi2_tab,chi2_tab_ndf,phit,nos,
     & x_ref,par_trk,erpr_trk,
     & x_ref2,par_trk2,erpr_trk2,
     & ntpl,ntpl2,
     & listhit1,listhit2,nonclean) 

**************************************************
*                                                *
*     routine for cleaning up glued tracks       *
*                                                *
**************************************************
#include "PARAM.h"
#include "HITS_OUT.h"


         dimension lassoc(2,*),nos(100)
         dimension chi2_tab(*),chi2_tab_ndf(*)
         dimension x_ref(*),x_ref2(*)
         dimension par_trk2(4,*),erpr_trk2(4,*)
         dimension par_trk(5,*),erpr_trk(5,*)
        
         dimension listhit1(ntrack,*),listhit2(ntrack,*)

         if(nonclean.eq.0)then                  ! just no cleaning - only rewrititng listhits
            do i=1,mlassc
            do k=1,ntpl
            listhit1(miter+i,k)=line1(k,lassoc(1,i))           
            enddo

            do k=1,ntpl2
            listhit2(miter+i,k)=line2(k,lassoc(2,i))           
            enddo

            enddo
            return
         endif
  

         if(mlassc.eq.1) then

            do k=1,ntpl
   
            listhit1(miter+1,k)=line1(k,lassoc(1,1))           
            enddo

            do k=1,ntpl2
            listhit2(miter+1,k)=line2(k,lassoc(2,1))           
            enddo
       
            return
         endif



*---- cleaning if needed

            eps=1.e-5
            
            licznik_g=0
            licznik_b=0  

            phit1=phit
            phit2=phit   


          do i=1,mlassc-1
               mi=lassoc(1,i)
               ni=lassoc(2,i)
          
          if(mi.eq.-900.and.ni.eq.-900) go to 12

          do j=i+1,mlassc
               mj=lassoc(1,j)
               nj=lassoc(2,j)
              
          if(mj.eq.-900.and.nj.eq.-900) go to 11
               
                 lihit11=0
                 lihit12=0
                 lic1=0

               do k=1,ntpl
          if(line1(k,mi).ne.8888)lihit11=lihit11+1
          if(line1(k,mj).ne.8888)lihit12=lihit12+1
          if(line1(k,mi).ne.8888.and.line1(k,mj).ne.8888)then
 
          if(line1(k,mi).eq.line1(k,mj))lic1=lic1+1
          endif
               enddo
                       
             if(lihit11.gt.lihit12)then
                   lihit1=lihit11
             else
                   lihit1=lihit12
             endif
   
                   if(lihit1.ne.0)then     
                    alicz1=float(lic1)*100/float(lihit1)
                   else
                    alicz1=0.0
                   endif

*                 print *,i,j,lihit1,lic1,alicz1

                 lihit21=0
                 lihit22=0
                 lic2=0

    
               do k=1,ntpl2
          if(line2(k,ni).ne.8888)lihit21=lihit21+1
          if(line2(k,nj).ne.8888)lihit22=lihit22+1
          if(line2(k,ni).ne.8888.and.line2(k,nj).ne.8888)then

          if(line2(k,ni).eq.line2(k,nj))lic2=lic2+1
          endif
               enddo
 
             if(lihit21.gt.lihit22)then
                   lihit2=lihit21
             else
                   lihit2=lihit22
             endif

                    if(lihit2.ne.0)then  
                     alicz2=float(lic2)*100/float(lihit2)
                    else
                     alicz2=0.0
                    endif

*                 print *,i,j,lihit2,lic2 ,alicz2
  
                  lic=lic1+lic2
                  lihit=lihit1+lihit2

                  if(lihit.ne.0)then  
                   alicz=float(lic)*100/float(lihit)
                  else
                   alicz=0.0
                  endif

*                 print *,'alicz-phit',i,j, alicz,phit

                 if(alicz.ge.phit)then                  ! really glued tracks
                     licznik_g=licznik_g +1 
     
                  if(chi2_tab_ndf(i).lt.chi2_tab_ndf(j))then                  
                     lassoc(1,j)=-900
                     lassoc(2,j)=-900
                     go to 11
                  elseif(chi2_tab_ndf(i).gt.chi2_tab_ndf(j))then 
                     lassoc(1,i)=-900
                     lassoc(2,i)=-900
                     go to 12
                  else                                   ! same chi2 fit - another criteria needed
                     
                     chi1=chi_proj(mi,1)+chi_proj(ni,2) 
                     chi2=chi_proj(mj,1)+chi_proj(nj,2)

                   if(chi1.lt.chi2)then
                     lassoc(1,j)=-900
                     lassoc(2,j)=-900
                     go to 11
                   else
                     lassoc(1,i)=-900
                     lassoc(2,i)=-900
                     go to 12
                   endif
                  endif
             

C--- kk --- below criteria are reverse (lines chi2 first). It seems work better?
           
                 elseif(alicz1.ge.phit1.and.alicz2.lt.phit2)then       ! first segment glued
                     licznik_b=licznik_b +1
                       go to 5
                 elseif(alicz1.lt.phit1.and.alicz2.ge.phit2)then       ! second segment glued
                     licznik_b=licznik_b +1

 5                  continue

                     chi1=chi_proj(mi,1)+chi_proj(ni,2) 
                     chi2=chi_proj(mj,1)+chi_proj(nj,2)

                   if(chi1.lt.chi2)then
                     lassoc(1,j)=-900
                     lassoc(2,j)=-900
                     go to 11
                   elseif(chi1.gt.chi2)then
                     lassoc(1,i)=-900
                     lassoc(2,i)=-900
                     go to 12
                   else                                              ! total chi2 criteria


                  if(chi2_tab_ndf(i).lt.chi2_tab_ndf(j))then                  
                     lassoc(1,j)=-900
                     lassoc(2,j)=-900
                     go to 11
                  elseif(chi2_tab_ndf(i).gt.chi2_tab_ndf(j))then 
                     lassoc(1,i)=-900
                     lassoc(2,i)=-900
                     go to 12
                  else                                               ! same chi2 fit - go on so far
                     
                  endif
                 endif
                endif

                              
  11         continue
 
             enddo
  12         continue

             enddo

  
            if(licznik_g.eq.0.and.licznik_b.eq.0) then      ! no glued/bad tracks

             do i=1,mlassc           

            do k=1,ntpl
            listhit1(miter+i,k)=line1(k,lassoc(1,i))           
            enddo

            do k=1,ntpl2
            listhit2(miter+i,k)=line2(k,lassoc(2,i))           
            enddo

            enddo

            return
            endif

C----- if yes  - now rewritting new reconstructed lines:

                 ind=0

                 do i=1,mlassc

                  if(lassoc(1,i).ne.-900.and.lassoc(2,i).ne.-900)then
                    ind=ind+1                                        
                   
                    chi2_tab(miter+ind)= chi2_tab(miter+i)          ! chi2
                    chi2_tab_ndf(miter+ind)= chi2_tab_ndf(miter+i)
                    x_ref(miter+ind)=x_ref(miter+i)
                    x_ref2(miter+ind)=x_ref2(miter+i)

                    do j=1,5
                    par_trk(j,miter+ind)= par_trk(j,miter+i)
                    erpr_trk(j,miter+ind)=erpr_trk(j,miter+i)
                    enddo

                    do j=1,4
                    par_trk2(j,miter+ind)= par_trk2(j,miter+i)
                    erpr_trk2(j,miter+ind)=erpr_trk2(j,miter+i)
                    enddo

                    do k=1,ntpl
                       
                    listhit1(miter+ind,k)=line1(k,lassoc(1,i))           
                    enddo

                    do k=1,ntpl2
                    listhit2(miter+ind,k)=line2(k,lassoc(2,i))           
                    enddo
                   

                  endif
                 enddo

                   mlassc = ind

                    nos(73)=nos(73) + licznik_g
                    nos(72)=nos(72) + licznik_b

            return
            end


      subroutine clean_segments(line,ntpl,phit,ans_sav,chimin,
     &                              new_line)

**************************************************
*                                                *
*      routine for cleaning up glued segments    *
*                     before SM1                 *
*                                                *
**************************************************
#include "PARAM.h"
#include "HITS_OUT.h"
#include "LIN_ALL.h"

          dimension ans_sav(10),line(mxpl)

              new_line=0
              procent1=0
              linuma=0

                do lisofar=1,nlin_fda
                 liclus=0
                 licz_clus1=0
                 licz_clus2=0

                do i=1,ntpl
      if(line1(i,lisofar).ne.8888.and.line1(i,lisofar).ne.0)
     &    licz_clus1=licz_clus1+1
      if(line(i).ne.8888.and.line(i).ne.0)licz_clus2=licz_clus2+1
   
      if(line1(i,lisofar).ne.8888.and.line(i).ne.8888.and.
     &    line1(i,lisofar).ne.0.and.line(i).ne.0)then
      if(line1(i,lisofar).eq.line(i))liclus=liclus+1
   
      endif          
                enddo

      licz_clus=licz_clus2 
      if(licz_clus1.gt.licz_clus2)licz_clus=licz_clus1

                procent=float(liclus)/float(licz_clus)*100.
             if(procent.gt.procent1)then
                procent1=procent
                linuma=lisofar
             endif
                enddo
               
              if(linuma.eq.0)return

             if(procent1.gt.phit)then
                if(chi_proj(linuma,1).gt.chimin)then
                 xlina(linuma) = ans_sav(1)            ! x - telescope center
                 ylina(linuma) = ans_sav(2)            ! intercept y
                 zlina(linuma) = ans_sav(3)            ! intercept z
                 yplina(linuma) = ans_sav(4)           ! dy/dx
                 zplina(linuma) = ans_sav(5)           ! dz/dx
                 chi_proj(linuma,1) = chimin           ! chi2 for line  
                do ise=1,ntpl
                   line1(ise,linuma)=line(ise)
                enddo
                 call ucopy(ans_sav(6),errlia(1,linuma),5)
                endif 
               
                new_line=linuma
  
             endif

              return
              end







