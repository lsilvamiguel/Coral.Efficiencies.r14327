
          subroutine find_best(jdett,lref,ivref,lrem,ncand,ivv)

***************************************************************
*                                                             *
*                routine for best line fitting                *
*                                                             *
***************************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "HITS_IN.h"
#include "HITS_OUT.h"
#include "TUNE.h"
#include "LIN_ALL.h"

         

          dimension lref(4),ivref(4)
          dimension lrem(mxpl),ncand(mxpl),ivv(mxpl,*)
      
          dimension ans(12),ans_sav(12)
       
          dimension ipointer(mxpl),line(mxpl)
          dimension jhitsav(mxpl),jr(mxpl),jhit(mxpl)
          dimension hits(mxpl),errors_hits(mxpl)
          logical lscifi,ff

          ff=.false.
          if(jdett.eq.201)then
             jdett=101
             ff=.true.
          endif 

          do i=1, ntpl+ntpl2
           hits(i)=800.                 ! fit flag 
           line(i)=8888
          enddo

* check min number of planes
          lscifi = .false.

          nthet = 0
          n_add = 0
          if(jdett.eq.1.or.jdett.eq.101)then
             jdet=1
          elseif(jdett.eq.2)then
             jdet=2
          endif


          if (jdet.eq.1) then
          nps = ntpl - 4                   ! ok, 4 just used

          do jt =1,nps
          if ((jyzt(lrem(jt)).eq.3.or.jyzt(lrem(jt)).eq.4).
     1            and.ncand(jt).gt.0)  
     1            nthet = nthet + 1
          if (ncand(jt).gt.0) n_add = n_add + 1
          enddo

          elseif(jdet.eq.2)then
          nps = ntpl2 - 4
          do jt =1,nps     
           
          if ((jyzt2(lrem(jt)).eq.3.or.jyzt2(lrem(jt)).eq.4).
     1           and.ncand(jt).gt.0) 
     1            nthet = nthet + 1
          if (ncand(jt).gt.0) n_add = n_add + 1
          enddo
          endif


*          print *,"fi-best:jdet,th,all:",jdet,nthet,n_add
*--------------------------------------------------
* -- bug corrected (March 2000 RW&KK---completed Jan.2001)

          if (jdet.eq.1.and.ityp1(lref(1)).eq.333.
     1     and.ityp1(lref(2)).eq.333.and.ityp1(lref(3)).eq.333.
     2     and.ityp1(lref(4)).eq.333)  lscifi = .true.
 
          if (lscifi)  then
            if (nthet.lt.2) return


          elseif(jdet.eq.1) then
        
           if(nthet.lt.nminpl1(1))return
           if(n_add.lt.nminpl1(2))return
        
          elseif(jdet.eq.2)then

           if(nthet.lt.nminpl2(1))return
           if(n_add.lt.nminpl2(2))return
          endif
*
          chimin = 1.0E+6

          do l = 1,4
          lll = lref(l)
          if (jdet.eq.1) then
          hits(lll) = clusters1(numhit(lll,ivref(l)))
          errors_hits(lll) = clu_errors1(numhit(lll,ivref(l)))
          line(lll)=numhit(lll,ivref(l))

          elseif(jdet.eq.2)then

          hits(lll) = clusters2(numht2(lll,ivref(l)))
          errors_hits(lll) = clu_errors2(numht2(lll,ivref(l)))
          line(lll)=numht2(lll,ivref(l))
          endif
          enddo

          do 55 lp =1,nps    
          lll = lrem(lp)
          nh = ncand(lp)
          ipointer(lp)=0

          if (nh.le.0) then
*  if there are no hits in this plane,create an artificial one (> 80.0) to make 
*                            stfit1 happy...
            
           
            ncand(lp) = 1                         
            ivv(lp,1) = 1
            ipointer(lp)=lll         
         endif

55        continue


          do ll = 1,nps - 1
          jr(ll) = 1
          enddo

          jr(nps) = 0                   


*------ loop over clusters sequences......................

10        continue                          

*---- get next sequence of hits...........................
          call next_sq(nps,lrem,ncand,ivv,jr,jhit)

*----ended if jr(1) = 0...................................

          if (jr(1).eq.0) go to 15

*----take clusters........................................

          do ll=1,nps
          if (jdet.eq.1) then
       
            if(lrem(ll).eq.ipointer(ll))then
             hits(lrem(ll)) = 800.
             line(lrem(ll)) = 8888
            else  
          hits(lrem(ll)) = clusters1(numhit(lrem(ll),jhit(ll)))
          errors_hits(lrem(ll))=clu_errors1(numhit(lrem(ll),jhit(ll)))
          line(lrem(ll))=numhit(lrem(ll),jhit(ll))
            endif

          elseif(jdet.eq.2)then

            if(lrem(ll).eq.ipointer(ll))then
             hits(lrem(ll)) = 800.
             line(lrem(ll)) = 8888
            else  
             hits(lrem(ll)) = clusters2(numht2(lrem(ll),jhit(ll)))
             f=1.
         errors_hits(lrem(ll))=clu_errors2(numht2(lrem(ll),jhit(ll)))
             line(lrem(ll))=numht2(lrem(ll),jhit(ll))
            endif

          endif
          enddo

*-------- now do the fits..............................

          if (jdet.eq.1) then
          ans(1) =xplane(1) 
 

*-------- now routine from phoenix.....................
            

          call stfit5(ntpl,cosn,sinn,fd,xplane,hits,      
     &                errors_hits,nhit,ans)


          elseif(jdet.eq.2)then                                 ! second telescope

          ans(1) =xplan2(1) 

          call stfit5(ntpl2,cosn2,sinn2,fd2,xplan2,hits,
     &                errors_hits,nhit,ans)
 
          endif
            if (ans(10).lt.chimin) then
              chimin = ans(10)
              chichi = ans(10)*(nhit-5)

              do jo=1,nps
                jhitsav(jo)=jhit(jo)
              enddo
              do jo=1,12
                ans_sav(jo)=ans(jo)
              enddo
  

            endif

          go to 10              ! as long as satisfy condition go to 15

15        continue
 

*          print *,"fi-best,jdet,ch/nd,pb:",jdet,ans_sav(10),ans_sav(11)
*---------------------------------------------------------


*---- if fit OK, keep line, kill hits which have been used.........
C-kk -- chimin - chi^2/df

      
           li_m=li_m+1

       if(li_m.gt.mxlin-mxlin1)return

       if(li_m.le.(mxlin-mxlin1).and.chimin.gt.chi_par)return

       
          if (jdet.eq.1)  then
          do k=1,nps
            if(lrem(k).eq.ipointer(k))then
             hits(lrem(k)) = 800.
             line(lrem(k)) = 8888
            else 
         line(lrem(k))=numhit(lrem(k),jhitsav(k)) 
        
         if(line(lrem(k)).eq.0)line(lrem(k))=8888            ! 2001 safety exit
            endif
          enddo

*          print *,' fit ok.det.1...ans(10)=',chimin

           
           nlin=0
           if(jdett.eq.101)then
              nlin=nlin_fda

              proc=phit                               ! % of clusters for clean segments

           call clean_segments(line,ntpl,proc,ans_sav,chimin,
     &                         new_line)
           
           if(new_line.ne.0)go to 77

           endif

          nlin_fd = nlin_fd + 1                       ! lines found

          if(nlin+nlin_fd.eq.mxlin-mxlin1)li_m=mxlin+1

          xlina(nlin+nlin_fd) = ans_sav(1)            ! x - telescope center
          ylina(nlin+nlin_fd) = ans_sav(2)            ! intercept y
          zlina(nlin+nlin_fd) = ans_sav(3)            ! intercept z
          yplina(nlin+nlin_fd) = ans_sav(4)           ! dy/dx
          zplina(nlin+nlin_fd) = ans_sav(5)           ! dz/dx
          chi_proj(nlin+nlin_fd,jdet) = chimin        ! chi2 for line
                
       
          do ise=1,ntpl
          line1(ise,nlin+nlin_fd)=line(ise)
          enddo

          call ucopy(ans_sav(6),errlia(1,nlin+nlin_fd),5)   
C-kk - errlia(5,..) = ans_sav(10) - chi^2/df


 77       continue
          do k=1,nps
          jflag(lrem(k),jhitsav(k)) = 1
*            print *,' flag set..plane,ord;',lrem(k),jhitsav(k)
          enddo

          do kr=1,4
          jflag(lref(kr),ivref(kr)) = 1
*            print *,' flag set..plane,ord;',lref(kr),ivref(kr)
          enddo

*         print*,' in find ',nlin_fd,(line1(is,nlin_fd),is=1,ntpl)

         elseif(jdet.eq.2)then                                  ! second telescope

          do k=1,nps
            if(lrem(k).eq.ipointer(k))then
             hits(lrem(k)) = 800.
             line(lrem(k)) = 8888
            else  
           line(lrem(k))=numht2(lrem(k),jhitsav(k)) 
            endif
          enddo

*          print *,' fit ok.det.2...ans(10)=',chimin 


          nlin_fd = nlin_fd + 1

          if(nlin_fd.eq.mxlin-mxlin1)li_m = mxlin+1 

          xlinb(nlin_fd) = ans_sav(1)
          ylinb(nlin_fd) = ans_sav(2)
          zlinb(nlin_fd) = ans_sav(3)
          yplinb(nlin_fd) = ans_sav(4)
          zplinb(nlin_fd) = ans_sav(5)
          chi_proj(nlin_fd,jdet) = chimin       ! chi2/ndf for line
          errlib(6,nlin_fd)=chichi

          do ise=1,ntpl2
          line2(ise,nlin_fd)=line(ise)
          enddo

*         print*,nlin_fd,(line2(is,nlin_fd),is=1,ntpl2)

          call ucopy(ans_sav(6),errlib(1,nlin_fd),5)


          do k=1,nps
          jflg2(lrem(k),jhitsav(k)) = 1
          enddo
          do kr=1,4
          jflg2(lref(kr),ivref(kr)) = 1
          enddo
          endif
     

          return
          end




          subroutine next_sq(nps,lrem,ncand,ivv,jr,jhit)

*************************************************************
*                                                           *
*   routine which gives set of hits (sequence) to stfit1    *
*                                                           *
*************************************************************
#include "PARAM.h"

          dimension lrem(*),ncand(*),ivv(mxpl,*),jr(*),jhit(*)
       
          do k= nps,1,-1                   ! backward loop
          jr(k) = jr(k) + 1

          if (jr(k).le.ncand(k)) go to 11
          jr(k) = 1
          enddo
          jr(1) = 0
          return
11        do lp = 1,nps
          jhit(lp) = ivv(lp,jr(lp))        

          enddo
          end

