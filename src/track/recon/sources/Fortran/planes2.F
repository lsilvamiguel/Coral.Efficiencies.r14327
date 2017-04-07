      subroutine planes2(list2_black,npla2,
     &zpla2,angel2,sinan2,cosan2,spresol2,type2,pleff2,uid2)                  
**************************************************
*                                                *
*     transports planes description behind SM1   * 
*             from coral to recon                *
*                                                *
**************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "VERS.h"
#include "MAGNET.h"

       INTEGER npla2,type2(*),uid2(*),list2_black(*)

       REAL zpla2(*),angel2(*),spresol2(*),pleff2(*),
     &      sinan2(*),cosan2(*)
       dimension ind_fd2(mxpl),mad(6)
       double precision  xx(3),h_field(3)
       logical printout,sub,fringe2
       logical version,str,lay
       logical e2_black,mib

         eps=1.0e-2
         eps_det=1.0e-5
*--------- second telescope (behind SM1)


*------------- if subsetups mode: ----------------*
*
* the organization of planes vectors as follows:
*       1 - ntpl2_s1 -> for subsetup 1
*       npla2 - npla2+ntpl2 -> for all
*       2*npla2 - 2*npla2+ntpl2_s2 -> for subsetup 2
*       3*npla2 - once again subsetup 1 (for restoring)
*
*--------------------------------------------------* 
         inda=0
         nypl2= 0
         nzpl2= 0
         nthpl2= 0
         ntpl2=0
         numlay=3                                  ! one layer
         npa=0

         if(sub)then

*---- for subsetups
         npa=npla2
         nypl2= npa
         nzpl2= npa
         nthpl2=npa
         ntpl2=npa
         inda=npa  
   
         if(mxpl.lt.4*npa)then
           print*,'too small dimensions in arrays'
           stop
          endif

         endif

         ntpl2_s1=0
         nthpl2_s1=0
         nzpl2_s1 =0
         nypl2_s1 =0
            
         ntpl2_s2=0
         nthpl2_s2=0
         nzpl2_s2 =0
         nypl2_s2 =0

          if(printout)then     
           print *
           print *,'   first spectrometer - behind SM1'
           print *
          endif

               
            if(i_option2.eq.1.and.n_black2.ne.0)then
          print*
          print*,'====>> WARNING: planes listed below are OFF '
          print*,'                 from reconstruction !  '
          print*
         
          print*,(list2_black(i),i=1,n_black2)
          print*
             endif

          
         ifalse=0
         str=.false.
         lay=.true.
         kstr=0
         nm=1

         ntpl2_cor=npla2  

*----------- loop over planes        
         do i=1,npla2
          
          e2_black=.true.
          mib(ntpl_cor+i)=.true.
 
          if(i_option2.eq.0)go to 4

          do ibra=1,n_black2
           if(uid2(i).eq.list2_black(ibra))then
           e2_black=.false.
           mib(ntpl_cor+i)=.false.
           mib_lab(ntpl_cor+i)=uid2(i)
           ifalse=ifalse+1
           endif
          enddo
 4        continue

*---- multiple coulomb scattering
 
          if(type2(i).eq.555)dmcs=dmcs2
          if(type2(i).eq.333)dmcs=dmcs3
          if(type2(i).eq.222)dmcs=dmcs2
          if(type2(i).eq.444)dmcs=dmcs4 

            
           tetdet(ntpl_cor+i)= 0.0136 * sqrt(dmcs)*fudge  
           xxx(ntpl_cor+i) = zpla2(i)/1000. 
           difi=xxx(ntpl_cor+i) - xxx(ntpl_cor+i-1) 
           if(type2(i).eq.444.and.abs(difi).lt.esp_det)then
           xxx(ntpl_cor+i-1)=zpla2(i-1)/1000.- 0.1*spresol2(i-1)/1000
           xxx(ntpl_cor+i)=zpla2(i)/1000.+ 0.1*spresol2(i)/1000 
           endif
  
            ccc(ntpl_cor+i)=-sinan2(i)      ! coral hit convention -> recon hit convention
            sss(ntpl_cor+i)=cosan2(i)
          
*---- end of multiple c.s.                            

          if(e2_black)then  
          if(zpla2(i)/1000.0.lt.x0_m1)go to 102

          angi2=angel2(i)

         c2=-sinan2(i)       ! coral hit convention -> recon hit convention
         s2=cosan2(i)
         ity_p = 1000
       
     
         if((abs(c2).lt.eps).or.(s2.gt.(1.-eps).and.s2.lt.1.).or.
     &      (s2.lt.(-1.+eps).and.s2.gt.-1.)) ity_p = 1            ! Y plane
        
         if((abs(s2).lt.eps).or.(c2.gt.(1.-eps).and.c2.lt.1.).or. 
     &      (c2.lt.(-1.+eps).and.c2.gt.-1.)) ity_p = 2            ! Z plane
       
           iuid2(i)=uid2(i)

             inda=inda+1
             i_s2(inda)=i  
             xplan2(inda)=zpla2(i)/1000
             sgpln2(inda)=spresol2(i)/1000

             ityp2(inda)=type2(i)
             eff_tel2(inda)=pleff2(i)
             sinn2(inda)=s2
             cosn2(inda)=c2


*--- straws special
         if(str)then                                    ! reset
           kstr=0
           str=.false.
           lay=.true.
           nm=nm+1
         endif
         if(type2(i).eq.555)then
             kstr=kstr +1
             modul(inda)=nm
         endif     
         if(kstr.eq.3)then                              ! straw layer
            lay=.false.
            nm=nm+1
         endif
         if(kstr.eq.5)lay=.true.
         if(kstr.eq.6)str=.true.                        ! double layer
            
*
*---
*
            kgem=0
         if(ity_p.eq.1)then
           jyzt2(inda) = 1                              ! Y plane
           nypl2 = nypl2 + 1
           jypl2(nypl2) = inda-npa
           kgem=1
C--- to protect the same planes position to make spline happy
           if(type2(i).eq.444)
     &     xplan2(inda)=xplan2(inda)- 0.1*sgpln2(inda)

           if(type2(i).eq.555.and.lay)then
             if(kstr.eq.1.or.kstr.eq.5)signu=-1.0
             if(kstr.eq.2.or.kstr.eq.6)signu=1.0
           xplan2(inda)=xplan2(inda)- signu*0.1*sgpln2(inda)
           endif

         elseif(ity_p.eq.2)then
           jyzt2(inda) = 2                              ! Z plane
           nzpl2 = nzpl2 + 1
           jzpl2(nzpl2) = inda-npa
           kgem=2
C--- to protect the same planes position to make spline happy
           if(type2(i).eq.444)
     &     xplan2(inda)=xplan2(inda)+ 0.1*sgpln2(inda) 

           if(type2(i).eq.555.and.lay)then
             if(kstr.eq.1.or.kstr.eq.5)signu=-1.0
             if(kstr.eq.2.or.kstr.eq.6)signu=1.0
           xplan2(inda)=xplan2(inda)- signu*0.1*sgpln2(inda)
           endif

         elseif(c2.gt.0.and.c2.lt.1)then
           jyzt2(inda) = 4                              ! -th plane
           nthpl2 = nthpl2 + 1
           jthpl2(nthpl2) = inda-npa
           kgem=4
C--- to protect the same planes position to make spline happy
           if(type2(i).eq.444)
     &     xplan2(inda)=xplan2(inda)- 0.1*sgpln2(inda)

           if(type2(i).eq.555.and.lay)then
             if(kstr.eq.1.or.kstr.eq.5)signu=-1.0
             if(kstr.eq.2.or.kstr.eq.6)signu=1.0
           xplan2(inda)=xplan2(inda)- signu*0.1*sgpln2(inda)
           endif

         elseif(c2.lt.0.and.c2.gt.-1)then
           jyzt2(inda) = 3                              ! th plane
           nthpl2 = nthpl2 + 1
           jthpl2(nthpl2) = inda-npa
           kgem=3
C--- to protect the same planes position to make spline happy
           if(type2(i).eq.444)
     &     xplan2(inda)=xplan2(inda)+ 0.1*sgpln2(inda)

           if(type2(i).eq.555.and.lay)then
             if(kstr.eq.1.or.kstr.eq.5)signu=-1.0
             if(kstr.eq.2.or.kstr.eq.6)signu=1.0
           xplan2(inda)=xplan2(inda)- signu*0.1*sgpln2(inda)
           endif
         endif

 
C--- now coulomb energy losses
           ind_fd2(inda-npa)=0

           indapri=inda
           if(sub)indapri=indapri-npa

        if(ityp2(inda).eq.555)then 
*--------------- subsetup 1 definition
        if(sub)then
      
             ntpl2_s1=ntpl2_s1+1
             ind_fd2(inda-npa)=1

             xplan2(ntpl2_s1)=xplan2(inda)
             sinn2(ntpl2_s1)=sinn2(inda)
             cosn2(ntpl2_s1)=cosn2(inda)
             jyzt2(ntpl2_s1)=jyzt2(inda)
             ityp2(ntpl2_s1)=ityp2(inda)
             i_s2(ntpl2_s1)=i_s2(inda)
             sgpln2(ntpl2_s1) = sgpln2(inda)
             modul(ntpl2_s1)=modul(inda)
 
             xplan2(3*npa+ntpl2_s1)=xplan2(inda)
             sinn2(3*npa+ntpl2_s1)=sinn2(inda)
             cosn2(3*npa+ntpl2_s1)=cosn2(inda)
             jyzt2(3*npa+ntpl2_s1)=jyzt2(inda)
             ityp2(3*npa+ntpl2_s1)=ityp2(inda) 
             i_s2(3*npa+ntpl2_s1)=i_s2(inda) 
             sgpln2(3*npa+ntpl2_s1) = sgpln2(inda)
             modul(3*npa+ntpl2_s1)=modul(inda)

             if(kgem.eq.3.or.kgem.eq.4)then
                nthpl2_s1=nthpl2_s1+1
                jthpl2(nthpl2_s1)= ntpl2_s1
                jthpl2(3*npa+nthpl2_s1)= ntpl2_s1

             elseif(kgem.eq.2)then
                nzpl2_s1 = nzpl2_s1 + 1
                jzpl2(nzpl2_s1)= ntpl2_s1
                jzpl2(3*npa+nzpl2_s1)= ntpl2_s1

             elseif(kgem.eq.1)then
                nypl2_s1 = nypl2_s1 + 1
                jypl2(nypl2_s1)= ntpl2_s1
                jypl2(3*npa+nypl2_s1)= ntpl2_s1

             endif
       
        
      if(printout)then
      print '("plane no. ",2i3," - ",4f9.5,i2,f9.6,"  straw ",
     & i4,"    sub_set1")',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),iuid2(i)
      endif
         endif
*----------- end of set1

  
        if(printout.and..not.sub)then
      print '("plane no. ",2i3," - ",4f9.5,i2,f9.6,"  straw ",i4)',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),iuid2(i)    

        endif


        elseif(ityp2(inda).eq.444)then
*------------ subsetup 2 definition
        if(sub)then
             ntpl2_s2=ntpl2_s2+1
             ind_fd2(inda-npa)=2

             xplan2(2*npa+ntpl2_s2)=xplan2(inda)
             sinn2(2*npa+ntpl2_s2)=sinn2(inda)
             cosn2(2*npa+ntpl2_s2)=cosn2(inda)
             jyzt2(2*npa+ntpl2_s2)=jyzt2(inda)
             ityp2(2*npa+ntpl2_s2)=ityp2(inda)
             sgpln2(2*npa+ntpl2_s2) = sgpln2(inda)
             i_s2(2*npa+ntpl2_s2) = i_s2(inda)

             if(kgem.eq.3.or.kgem.eq.4)then
                nthpl2_s2=nthpl2_s2+1
                jthpl2(2*npa+nthpl2_s2)= ntpl2_s2
             elseif(kgem.eq.2)then
                nzpl2_s2 = nzpl2_s2 + 1
                jzpl2(2*npa+nzpl2_s2)= ntpl2_s2
             elseif(kgem.eq.1)then
                nypl2_s2 = nypl2_s2 + 1
                jypl2(2*npa+nypl2_s2)= ntpl2_s2
             endif

       
      
        if(printout.and.ityp2(inda).eq.444)then
      print '("plane no. ",2i3," - ",4f9.5,i2,f9.6,"  gem ",
     & i4,"  sub_set2")',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),
     & iuid2(i)
        elseif(printout.and.ityp2(inda).eq.333)then
      print '("plane no. ",2i3," - ",4f9.5,i2,f9.6," sci-fi ",
     & i4,"  sub_set2")',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),
     & iuid2(i)                                     
        endif 
         endif
*----------- end of set2


          if(printout.and..not.sub)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  gem ",i4)',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),
     & iuid2(i)
          endif

        elseif(ityp2(inda).eq.333)then

          if(printout)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  sci-fi ",i4)',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),
     & iuid2(i)   
           endif

        elseif(ityp2(inda).eq.222)then 
       
*--------------- subsetup 1 definition
         if(sub)then

             ntpl2_s1=ntpl2_s1+1
             ind_fd2(inda-npa)=1

             xplan2(ntpl2_s1)=xplan2(inda)
             sinn2(ntpl2_s1)=sinn2(inda)
             cosn2(ntpl2_s1)=cosn2(inda)
             jyzt2(ntpl2_s1)=jyzt2(inda)
             ityp2(ntpl2_s1)=ityp2(inda)
             sgpln2(ntpl2_s1) = sgpln2(inda)
             i_s2(ntpl2_s1)=i_s2(inda)  

             xplan2(3*npa+ntpl2_s1)=xplan2(inda)
             sinn2(3*npa+ntpl2_s1)=sinn2(inda)
             cosn2(3*npa+ntpl2_s1)=cosn2(inda)
             jyzt2(3*npa+ntpl2_s1)=jyzt2(inda)
             ityp2(3*npa+ntpl2_s1)=ityp2(inda)    
             sgpln2(3*npa+ntpl2_s1) = sgpln2(inda) 
             i_s2(3*npa+ntpl2_s1)=i_s2(inda)  

             if(kgem.eq.3.or.kgem.eq.4)then
                nthpl2_s1=nthpl2_s1+1
                jthpl2(nthpl2_s1)= ntpl2_s1
                jthpl2(3*npa+nthpl2_s1)= ntpl2_s1

             elseif(kgem.eq.2)then
                nzpl2_s1 = nzpl2_s1 + 1
                jzpl2(nzpl2_s1)= ntpl2_s1
                jzpl2(3*npa+nzpl2_s1)= ntpl2_s1

             elseif(kgem.eq.1)then
                nypl2_s1 = nypl2_s1 + 1
                jypl2(nypl2_s1)= ntpl2_s1
                jypl2(3*npa+nypl2_s1)= ntpl2_s1

             endif


        if(printout)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  Saclay DC ",
     & i4,"   sub_set1")',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),
     & iuid2(i)
        endif

         endif
*------------------ end od set1
  

          if(printout.and..not.sub)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  Saclay DC ",i4)',
     & i,indapri,eff_tel2(inda),xplan2(inda),
     & sinn2(inda),cosn2(inda),jyzt2(inda),sgpln2(inda),
     & iuid2(i)
  
          endif
          
        endif

 
       
 102     continue

          endif 

         enddo
*----------------------- end of loop over coral planes

         ntpl2=inda

         if(sub)then
           ntpl2=inda-npa
           nypl2=nypl2-npa
           nzpl2=nzpl2-npa
           nthpl2=nthpl2-npa

           nypl2_tmp=nypl2
           nzpl2_tmp=nzpl2
           nthpl2_tmp=nthpl2
         endif

         ntpl2_tmp=ntpl2

      print *
      print '(" no of planes behind SM1 in coral: ",i5)',
     & npla2
      print '(" no of planes behind SM1 off: " ,i5)',ifalse                                               

      print '(" effective no of planes behind SM1 in recon :",
     & 2i4," Y; ",i4," Z; ",i4," th ")',
     & ntpl2,nypl2,nzpl2,nthpl2

      if(sub)then
      print *,'nb of Y planes in set.1: ', nypl2_s1
      print *,'nb of Y planes in set.2: ', nypl2_s2
      print *,'nb of Y planes : ', nypl2_tmp,nypl2
      print *,'nb of Z planes in set.1: ', nzpl2_s1
      print *,'nb of Z planes in set.2: ', nzpl2_s2
      print *,'nb of Z planes : ', nzpl2_tmp,nzpl2
      print *,'nb of th planes in set.1: ', nthpl2_s1
      print *,'nb of th planes in set.2: ', nthpl2_s2
      print *,'nb of th planes : ', nthpl2_tmp,nthpl2
      endif
      print *
      print *,'Y planes: ',(jypl2(npa+i),i=1,nypl2)
      print *,'Z planes: ',(jzpl2(npa+i),i=1,nzpl2)
      print *,'th planes: ',(jthpl2(npa+i),i=1,nthpl2)

  
      print *
      print *,'real coral nb of planes for gems,sci-fi ',
     &        'SDCs and straws: '
      print *,(i,i=1,ntpl2)
      print *,(i_s2(npa+i),i=1,ntpl2)
      print *,(iuid2(i_s2(npa+i)),i=1,ntpl2)
      print *

      print *,'straws modules: '
         modu=0
         k=0
      do i=1,ntpl2
      if(ityp2(i+npa).eq.555)then
      if(k.eq.0)modu=modul(i+npa)
         if(modu.eq.modul(i+npa))then
              k=k+1
              mad(k)=i
              mad(3+k)=iuid2(i_s2(i+npa))
         else
              if(k.ne.3)then
               print*, 'error in straw module logic: ',k
               stop
              endif
              print*,modul(npa+i-1),(mad(l),mad(3+l),l=1,3)
              modu=modul(i+npa)
              k=1
              mad(k)=i
              mad(3+k)=iuid2(i_s2(i+npa))
         endif     
      endif        
      enddo
             print*,modu,(mad(l),mad(3+l),l=1,3)
      print *,'nb. of active straws modules: ',modu

*
*--------------- mag.field double integral........
*                mod. RW. June 2002
*
      Bdl = 0.
      Bdldl = 0.

         fd2(1)=0.0
      if(sub)then
         fd2(3*npa+1)=0.0
         fd2(2*npa+1)=0.0
         fd2(npa+1)=0.0
      endif

      xx(3) = xplan2(1)
      ifds1=1
      ifds2=1

      do  50 k= 2,ntpl2
 49    xx(3) = xx(3) + 0.005
       xx(1) = 0.010 * (xx(3) - xplan2(1))
       xx(2) = 0.005 * (xx(3) - xplan2(1))

      call mag_field(xx,h_field)
 
       Bdl = bdl + 0.3 * h_field(2) * 0.005
       Bdldl = Bdldl + Bdl * 0.005
       if (xx(3).lt.xplan2(k)) go to 49

*      print *,' z,H123:',xx(3),h_field
*      print *,'k,xplan2,hy,bdl,bdl2:',k,xplan2(k),h_field(2),
*     1 bdl,bdldl

       if(sub)then

       if(ind_fd2(k).eq.1)then
        ifds1=ifds1+1
        fd2(ifds1)=bdldl
        fd2(3*npa+ifds1)=bdldl
       elseif(ind_fd2(k).eq.2)then
        ifds2=ifds2+1
        fd2(2*npa+ifds2)=bdldl
       endif

       fd2(npa+k) = bdldl
      elseif(.not.sub)then
       if(ind_fd2(k).ne.0)then
         print*,'error in index for fd2 -> planes2 '
         stop
       endif
   
       fd2(k) = bdldl
      
      endif
 50    continue 

*-- just crosscheck

       if(sub)then
       if(ifds1.ne.ntpl2_s1)then
       print*,'error in ifds1 -> planes2 ',ifds1,ntpl2_s1
         stop
       endif
       if((ifds2-1).ne.ntpl2_s2)then
       print*,'error in ifds2 -> planes2 ',ifds2,ntpl2_s2
         stop
       endif
       endif
*
* ---------------- multicoulomb scattering effect (cont.)    
*

            call  multi_coulomb 

       return
       end


















