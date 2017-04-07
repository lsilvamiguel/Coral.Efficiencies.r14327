      subroutine setup(ip)

********************************************************************
*                                                                  *
*    routine which defines subsetups if this mode is active        *
*                          (sub=.true.)                            *        
*                                                                  *
*          ip = 111 - SDC+straws  behind, micromegas before (2001) *
*                        + SDC before (2002)                       *
*          ip = 444 - gems behind, micromegas before               *
*          ip = 1   - all detectors active - restoration           *
*          ip = 0   - end - restoration of setup 1                 *
*                                                                  *
********************************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "HITS_IN.h"
#include "VERS.h"

      common/local_rem/jchits_tmp(mxpl),jc2hits_tmp(mxpl),
     &        numhit_tmp(mxpl,nmhit),numht2_tmp(mxpl,nmhit)  

            logical sub,pri_set,pri_set1,pri_set2
            logical version

            pri_set=.false.                   ! internal debugging
            pri_set1=.false.
            pri_set2=.false.

* ------------------------------------------------- *
*
*    if(sub)  vmat and vresol are organized as follows:
*    vmat(i,j) -> for set 1
*    vmat(ntt+i,ntt+j) -> all detectors
*    vmat(2*ntt+i,2*ntt+ntt+j) -> set 2
*    vmat(2*ntt+ntt2+i,2*ntt+ntt2+j) -> set 1 once again
*                                        
* ------------------------------------------------- *

            if(version)then                     ! 2001 setup
                ifac=2
            elseif(.not.version)then            ! 2002 setup
                ifac=3
            endif

             if(ip.eq.0)then
*------- last time in setup routine: copy planes info
     
            ntpl=ntpl_tmp
            ntpl2=ntpl2_tmp 

            nto1=ifac*ntpl_cor
            nto2=3*ntpl2_cor
  
             do k=1,ntpl_s1
              xplane(k)=xplane(nto1+k)
              sinn(k)=sinn(nto1+k)
              cosn(k)=cosn(nto1+k)
              jyzt(k)=jyzt(nto1+k)
              ityp1(k)=ityp1(nto1+k)
              sgpln(k)=sgpln(nto1+k)
              fd(k)=fd(nto1+k)
             enddo

              do k=1,nypl_s1
               jypl(k)=jypl(nto1+k)
              enddo
              do k=1,nzpl_s1
               jzpl(k)=jzpl(nto1+k)
              enddo
              do k=1,nthpl_s1
               jthpl(k)=jthpl(nto1+k)
              enddo

             do k=1,ntpl2_s1
              xplan2(k)=xplan2(nto2+k)
              sinn2(k)=sinn2(nto2+k)
              cosn2(k)=cosn2(nto2+k)
              jyzt2(k)=jyzt2(nto2+k)
              ityp2(k)=ityp2(nto2+k)
              sgpln2(k)=sgpln2(nto2+k)
              fd2(k)=fd2(nto2+k)
              modul(k)=modul(nto2+k)
             enddo


              do k=1,nypl2_s1
               jypl2(k)=jypl2(nto2+k)
              enddo
              do k=1,nzpl2_s1
               jzpl2(k)=jzpl2(nto2+k)
              enddo
              do k=1,nthpl2_s1
               jthpl2(k)=jthpl2(nto2+k)
              enddo

*-------- copy multi_coulomb info
             ntt1=ntpl_s1+ntpl2_s1
             if(ifac.eq.2)then
              ntt2=ntpl_s1+ntpl2_s2
             elseif(ifac.eq.3)then
              ntt2=ntpl_s2+ntpl2_s2
             endif
             ntt=ntpl_tmp+ntpl2_tmp
             ntot=2*ntt+ntt2

             do i=1,ntt1
             do j=1,ntt1
               vmat(i,j) = vmat(ntot+i,ntot+j)
             enddo
             enddo

             return

         elseif(ip.eq.111)then                   ! first subsetup used 

*------- define working parameters for planes
            
             ntpl=ntpl_s1
             nthpl=nthpl_s1
             nzpl=nzpl_s1
             nypl=nypl_s1

             ntpl2=ntpl2_s1
             nthpl2=nthpl2_s1
             nzpl2=nzpl2_s1
             nypl2=nypl2_s1

             nto=ntpl_cor
             nt=0

         do k=1,ntpl_tmp
           if(ityp1(nto+k).eq.111.or.ityp1(nto+k).eq.222)then
              nt=nt+1
              jchits(nt)=jchits_mc(k)
             do j=1,jchits(nt)
               numhit(nt,j) = numhit_mc(k,j)
             enddo
           endif   
         enddo
        
           if(nt.ne.ntpl_s1)
     &  print*,'error in subset 1: 1',nt,ntpl_s1,ntpl,ntpl_cor

             nt=0
             nto=ntpl2_cor

        do k=1,ntpl2_tmp
          if(ityp2(nto+k).eq.555.or.ityp2(nto+k).eq.222)then
            nt=nt+1
            jc2hits(nt)=jc2hits_mc(k)
           do j=1,jc2hits(nt)
             numht2(nt,j) = numht2_mc(k,j)
           enddo
          endif                   
        enddo

           if(nt.ne.ntpl2_s1)
     &   print*,'error in subset 1: 2',nt,ntpl2_s1,ntpl2,ntpl2_cor

       if(pri_set1)then
          do i=1, ntpl
      print*,i,ityp1(i),jyzt(i),xplane(i),sinn(i),cosn(i),sgpln(i)
      print*,'hits ',i,(numhit(i,j),j=1,jchits(i))
      print*,i,ityp1(i),jchits(i)
          enddo

        print*,'Y planes ',(jypl(i),i=1,nypl)
        print*,'Z planes ',(jzpl(i),i=1,nzpl)
        print*,'th planes ',(jthpl(i),i=1,nthpl)

        do i=1, ntpl2
      print*,i,ityp2(i),jyzt2(i),xplan2(i),sinn2(i),cosn2(i),sgpln2(i)
      print*,'hits2 ',i,(numht2(i,j),j=1,jc2hits(i))
      print*,i,ityp2(i),jc2hits(i)
          enddo

        print*,'Y planes ',(jypl2(i),i=1,nypl2)
        print*,'Z planes ',(jzpl2(i),i=1,nzpl2)
        print*,'th planes ',(jthpl2(i),i=1,nthpl2)
       endif

                return

         elseif(ip.eq.444)then                   ! second subsetup used 
       
*------- define working parameters for planes

           if(ifac.eq.3)then                     ! 2002 (NOT 2001!)
             nto=ntpl_cor

             ntpl=ntpl_s2
             nthpl=nthpl_s2
             nzpl=nzpl_s2
             nypl=nypl_s2

             do j=1, nthpl
              jthpl(j)=jthpl(2*nto+j)
             enddo
             do j=1,nzpl
              jzpl(j)=jzpl(2*nto+j)
             enddo 
             do j=1,nypl
              jypl(j)=jypl(2*nto+j)
             enddo

          do k=1,ntpl_s2
              xplane(k)=xplane(2*nto+k)
              sinn(k)=sinn(2*nto+k)
              cosn(k)=cosn(2*nto+k)
              jyzt(k)=jyzt(2*nto+k)
              ityp1(k)=ityp1(2*nto+k)
              sgpln(k)=sgpln(2*nto+k)
              fd(k)=fd(2*nto+k)
             enddo

           endif
*
*------- 
*
             ntpl2=ntpl2_s2
             nthpl2=nthpl2_s2
             nzpl2=nzpl2_s2
             nypl2=nypl2_s2
    
             do j=1, nthpl2
              jthpl2(j)=jthpl2(2*ntpl2_cor+j)
             enddo
             do j=1,nzpl2
              jzpl2(j)=jzpl2(2*ntpl2_cor+j)
             enddo 
             do j=1,nypl2
              jypl2(j)=jypl2(2*ntpl2_cor+j)
             enddo


             nt=0
             nto=ntpl2_cor

      do k=1,ntpl2_tmp
       if((ityp2(k).eq.222.or.ityp2(k).eq.555).and.k.le.ntpl2_s1)then
           jc2hits_tmp(k)=jc2hits(k)
         do j=1,jc2hits_tmp(k)
           numht2_tmp(k,j) = numht2(k,j)
         enddo
       endif

       if(ityp2(nto+k).eq.444)then
*.or.ityp2(nto+k).eq.333)then
           nt=nt+1
           jc2hits(nt)=jc2hits_mc(k)
         do j=1,jc2hits_mc(k)
           numht2(nt,j) = numht2_mc(k,j)
         enddo
       endif     
      enddo             

            if(nt.ne.ntpl2_s2)
     & print*,'error in subset 2: 2',nt,ntpl2_s2,ntpl2,ntpl2_cor   

             do k=1,ntpl2_s2
              xplan2(k)=xplan2(2*ntpl2_cor+k)
              sinn2(k)=sinn2(2*ntpl2_cor+k)
              cosn2(k)=cosn2(2*ntpl2_cor+k)
              jyzt2(k)=jyzt2(2*ntpl2_cor+k)
              ityp2(k)=ityp2(2*ntpl2_cor+k)
              sgpln2(k)=sgpln2(2*ntpl2_cor+k)
              fd2(k)=fd2(2*ntpl2_cor+k)
             enddo


       if(pri_set2)then
          do i=1, ntpl
      print*,i,ityp1(i),jyzt(i),xplane(i),sinn(i),cosn(i),sgpln(i)
      print*,'hits ',i,(numhit(i,j),j=1,jchits(i))
          enddo

        print*,'Y planes ',(jypl(i),i=1,nypl)
        print*,'Z planes ',(jzpl(i),i=1,nzpl)
        print*,'th planes ',(jthpl(i),i=1,nthpl)

        do i=1, ntpl2
      print*,i,ityp2(i),jyzt2(i),xplan2(i),sinn2(i),cosn2(i),sgpln2(i)
      print*,'hits2 ',i,(numht2(i,j),j=1,jc2hits(i))
          enddo

        print*,'Y planes ',(jypl2(i),i=1,nypl2)
        print*,'Z planes ',(jzpl2(i),i=1,nzpl2)
        print*,'th planes ',(jthpl2(i),i=1,nthpl2)
       endif

*-------- copy multi_coulomb info
             ntt2=ntpl_s1+ntpl2_s2
          if(ifac.eq.3)ntt2=ntpl_s2+ntpl2_s2
             ntt=ntpl_tmp+ntpl2_tmp
             ntot=2*ntt
 
             do i=1,ntt2
             do j=1,ntt2
               vmat(i,j) = vmat(ntot+i,ntot+j)
             enddo
             enddo                  

             return
     

         elseif(ip.eq.1)then               ! third subsetup - all detectors active
*------- define working parameters for planes

             ntpl=ntpl_tmp
             nthpl=nthpl_tmp
             nzpl=nzpl_tmp
             nypl=nypl_tmp
 
             ntpl2=ntpl2_tmp
             nthpl2=nthpl2_tmp
             nzpl2=nzpl2_tmp
             nypl2=nypl2_tmp

               do j=1,nthpl
                jthpl(j)=jthpl(ntpl_cor+j)
               enddo
               do j=1,nzpl
                jzpl(j)=jzpl(ntpl_cor+j)
               enddo
               do j=1,nypl
                jypl(j)=jypl(ntpl_cor+j)
               enddo

               do j=1,nthpl2
                jthpl2(j)=jthpl2(ntpl2_cor+j)
               enddo
               do j=1,nzpl2
                jzpl2(j)=jzpl2(ntpl2_cor+j)
               enddo
               do j=1,nypl2
                jypl2(j)=jypl2(ntpl2_cor+j)
               enddo
    
        

*----------- restore remainnig clusters before SM1
 
             nt=0
             nto=ntpl_cor

        do k=1,ntpl_tmp
            
         if((ityp1(k).eq.111.or.ityp1(k).eq.222).and. 
     &       k.le.ntpl_s1)then
            jchits_tmp(k)=jchits(k)
           do j=1,jchits_tmp(k)
             numhit_tmp(k,j) = numhit(k,j)
           enddo
         endif

         if(ityp1(nto+k).eq.333)then
            jchits(k)=jchits_mc(k)
           do j=1,jchits_mc(k)
             numhit(k,j) = numhit_mc(k,j)
           enddo
         elseif(ityp1(nto+k).eq.111.or.ityp1(nto+k).eq.222)then
             nt=nt+1
             jchits(k)=jchits_tmp(nt)
           do j=1,jchits(k)
             numhit(k,j) = numhit_tmp(nt,j)
           enddo      
         endif          
        enddo

            if(nt.ne.ntpl_s1)
     & print*,'error in subset 3: 1',nt,ntpl_s1,ntpl,ntpl_cor

* 
*----------- restore remainnig clusters behind SM1
*
*
             ntg=0
             nt=0
             nto=ntpl2_cor

      do k=1,ntpl2_tmp                 
         if(ityp2(k).eq.444.and.k.le.ntpl2_s2)then
           jc2hits_tmp(k+ntpl2_s1)=jc2hits(k)
          do j=1,jc2hits(k)
            numht2_tmp(k+ntpl2_s1,j) = numht2(k,j)
          enddo
         endif

      if(ityp2(nto+k).eq.333)then
           jc2hits(k)=jc2hits_mc(k)
         do j=1,jc2hits_mc(k)
           numht2(k,j) = numht2_mc(k,j)
         enddo
      elseif(ityp2(nto+k).eq.222.or.ityp2(nto+k).eq.555)then
          nt=nt+1
          jc2hits(k)=jc2hits_tmp(nt)
         do j=1,jc2hits_tmp(nt)
           numht2(k,j) = numht2_tmp(nt,j)
         enddo
      elseif(ityp2(nto+k).eq.444)then
          ntg=ntg+1
          jc2hits(k)=jc2hits_tmp(ntpl2_s1+ntg)
         do j=1,jc2hits_tmp(ntpl2_s1+ntg)
           numht2(k,j) = numht2_tmp(ntpl2_s1+ntg,j)
         enddo
      endif          
      enddo

            if(nt.ne.ntpl2_s1)
     & print*,'error in subset 3: 3',nt,ntpl2_s1,ntpl2,ntpl2_cor

            if(ntg.ne.ntpl2_s2)
     & print*,'error in subset 3: 4',ntg,ntpl2_s2,ntpl2,ntpl2_cor


*------------------- end of restoration procedure

             do k=1,ntpl
              xplane(k)=xplane(ntpl_cor+k)
              sinn(k)=sinn(ntpl_cor+k)
              cosn(k)=cosn(ntpl_cor+k)
              jyzt(k)=jyzt(ntpl_cor+k)
              ityp1(k)=ityp1(ntpl_cor+k)
              sgpln(k)=sgpln(ntpl_cor+k)
              fd(k)=fd(ntpl_cor+k)
             enddo


             do k=1,ntpl2
              xplan2(k)=xplan2(ntpl2_cor+k)
              sinn2(k)=sinn2(ntpl2_cor+k)
              cosn2(k)=cosn2(ntpl2_cor+k)
              jyzt2(k)=jyzt2(ntpl2_cor+k)
              ityp2(k)=ityp2(ntpl2_cor+k)
              sgpln2(k)=sgpln2(ntpl2_cor+k)
              fd2(k)=fd2(ntpl2_cor+k)
              modul(k)=modul(ntpl2_cor+k)
             enddo

*-------- copy multi_coulomb info
             ntt=ntpl_tmp+ntpl2_tmp
             ntot=ntt 
 
             do i=1,ntt
             do j=1,ntt
               vmat(i,j) = vmat(ntot+i,ntot+j)
             enddo
             enddo                  

       if(pri_set)then
          do i=1, ntpl
      print*,i,ityp1(i),jyzt(i),xplane(i),sinn(i),cosn(i),sgpln(i)
      print*,'hits ',i,(numhit(i,j),j=1,jchits(i))
          enddo

        print*,'Y planes ',(jypl(i),i=1,nypl)
        print*,'Z planes ',(jzpl(i),i=1,nzpl)
        print*,'th planes ',(jthpl(i),i=1,nthpl)

        do i=1, ntpl2
      print*,i,ityp2(i),jyzt2(i),xplan2(i),sinn2(i),cosn2(i),sgpln2(i)
      print*,'hits2 ',i,(numht2(i,j),j=1,jc2hits(i))
          enddo

        print*,'Y planes ',(jypl2(i),i=1,nypl2)
        print*,'Z planes ',(jzpl2(i),i=1,nzpl2)
        print*,'th planes ',(jthpl2(i),i=1,nthpl2)
       endif

               return
         else
            print*,'which sub_setup do you want? incorrect ip index'
            stop
         endif
*
*-------------- end of subsetups mode
*

                return
                end




