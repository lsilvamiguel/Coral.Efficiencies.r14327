
      subroutine multi_coulomb
*************************************************************
*                                                           *
*                covariance matrix for MCS                  *
*        mod. RW and KK - March/April 2002                  * 
*                                                           *
*************************************************************

#include "PARAM.h"
#include "PLANES.h"
#include "VERS.h"

       dimension mul(2*mxpl),vmat_temp(2*mxpl,2*mxpl)   
       logical sub,mib,printe
       logical version
*
* ---------------- multicoulomb scattering effect (cont.)
*
*    if(sub)  vmat and vresol are organized as follows:
*    vmat(i,j) -> for set 1
*    vmat(ntt+i,ntt+j) -> all detectors
*    vmat(2*ntt+i,2*ntt+j) -> set 2
*    vmat(2*ntt+ntt2+i,2*ntt+ntt2+j) -> set 1 once again
*
*

        printe=.false.

        nttc=ntpl_cor+ntpl2_cor
 
*---- All elements for the MCS part.......................................
 
      do k = 1,nttc-1
        do i = k+1,nttc
          do j = k+1,nttc
 
         vmat(i,j) = vmat(i,j)
     1    + (tetdet(k)**2) * (xxx(i)-xxx(k))
     2    *(xxx(j)-xxx(k))
     2    * (ccc(i)*ccc(j) + sss(i)*sss(j))
        enddo
       enddo
      enddo                    

        ntt=ntpl_tmp+ntpl2_tmp 
        mm=0

*           print*,ntt,ntpl_tmp,ntpl2_tmp,nttc
       do i=1,nttc
        if(mib(i))then
            mm=mm+1
            mul(mm)=i
        endif
       enddo

        if(mm.ne.ntt)then
           print*,'error in mcs ',mm,ntt
           stop
        endif

*
* -------- normal 
*
       if(.not.sub)then

        do i=1,ntt
        do j=1,ntt
          vmat_temp(i,j)=vmat(mul(i),mul(j))
        enddo
        enddo


*
* ------- subsetups mode
*
       elseif(sub)then

       ntt1=ntpl_s1+ntpl2_s1
       ntt2=ntpl_s2+ntpl2_s2
       ntt3=2*ntt+ntt2
 
        m1=0
        m2=0                                  
        mx1=0
        mx2=0 

      do i=1,nttc
       if(mib(i))then

         if(i.le.ntpl_cor)then
            mx1=mx1+1

         if(version)then                              ! 2001 setup
         if(ityp1(ntpl_cor+mx1).eq.111)then
         m1=m1+1
         m2=m2+1 
         mul(ntt+m1)=i
         mul(ntt+ntt1+m2)=i 
         endif

         elseif(.not.version)then                     ! 2002 setup
         if(ityp1(ntpl_cor+mx1).eq.111.or.
     &      ityp1(ntpl_cor+mx1).eq.222)then
         m1=m1+1
         mul(ntt+m1)=i
         if(ityp1(ntpl_cor+mx1).eq.111)then
         m2=m2+1
         mul(ntt+ntt1+m2)=i 
         endif
         endif
         endif

         elseif(i.gt.ntpl_cor)then
            mx2=mx2+1

       if(ityp2(ntpl2_cor+mx2).eq.222.
     & or.ityp2(ntpl2_cor+mx2).eq.555)then
         m1=m1+1
         mul(ntt+m1)=i
       endif
     
       if(ityp2(ntpl2_cor+mx2).eq.444)then
*     & or.ityp2(ntpl2_cor+mx2).eq.333)then
         m2=m2+1
         mul(ntt+ntt1+m2)=i
       endif
        endif
        endif
       enddo

      if(m1.ne.ntt1)then
        print*,' error in mcs-sub m1 ',m1,ntt1
        stop
      endif
      if(m2.ne.ntt2) then
        print*,' error in mcs-sub m2 ',m2,ntt2
        stop
      endif

       do i=1,ntt1
       do j=1,ntt1
         vmat_temp(i,j)=vmat(mul(ntt+i),mul(ntt+j))
         vmat_temp(ntt3+i,ntt3+j)=vmat_temp(i,j)
       enddo
       enddo                  

       do i=1,ntt
       do j=1,ntt
         vmat_temp(ntt+i,ntt+j)=vmat(mul(i),mul(j))
       enddo
       enddo          
 
       do i=1,ntt2
       do j=1,ntt2
         vmat_temp(2*ntt+i,2*ntt+j)=
     &   vmat(mul(ntt+ntt1+i),mul(ntt+ntt1+j))
       enddo
       enddo      
     
* ----------- rewritting
       do i=1,2*ntt+ntt2+ntt1
       do j=1,2*ntt+ntt2+ntt1 
        vmat(i,j)=vmat_temp(i,j)
       enddo
       enddo
             
       endif
* ------------------ end normal/sub-setup    

       if(printe)then
       print*,'setup 1'
       do i=1,ntpl_s1+ntpl2_s1
       print*,i,(vmat(i,j),j=1,ntpl_s1+ntpl2_s1)
       enddo
       print*,'setup 2'
       do i=1,ntpl_s2+ntpl2_s2
       print*,i,(vmat(2*ntt+i,2*ntt+j),j=1,ntpl_s2+ntpl2_s2)
       enddo
       print*,'all'
       do i=1,ntpl_tmp+ntpl2_tmp
       print*,i,(vmat(ntt+i,ntt+j),j=1,ntpl_tmp+ntpl2_tmp)
       enddo
       print*,'setup 2 again'
       do i=1,ntpl_s1+ntpl2_s1
       print*,i,(vmat(2*ntt+ntt2+i,2*ntt+ntt2+j),j=1,ntpl_s1+ntpl2_s1)
       enddo
       print*,' too long printout - debugimg mode MCS only - stop'
       stop
       endif

       return
       end







