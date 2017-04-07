      subroutine planes1(list1_black,npla,
     & zpla,angel,sinan,cosan,spresol,type,pleff,uid)   
**************************************************
*                                                *
*     transports planes description before SM1   * 
*             from coral to recon                *
*                                                *
**************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "TARGET.h"
#include "VERS.h"


       INTEGER npla,type(*),uid(*),list1_black(*)
       REAL zpla(mxpl),angel(mxpl),spresol(mxpl), pleff(*),
     &      sinan(*),cosan(*)
       dimension ind_fd(mxpl)
       double precision xx(3),h_field(3)
       logical printout,version
       logical e1_black,mib,sub 

         eps=1.0e-3

*--------    now detectors and geometry description 
*--------    rewriting variables: formal parameters to commons PLANES.h

         if(n_black1.eq.1.and.list1_black(1).eq.100000)go to 30
             if(i_option2.eq.1.and.n_black1.ne.0)then
          print*
          print*,'====>> WARNING: planes listed below are OFF '
          print*,'                 from reconstruction !  '
          print*
         
          print*,(list1_black(i),i=1,n_black1)
          print*
             endif
        
 30       continue

*--------- first telescope (before SM1,behind PT)

         ntpl=0
         nypl = 0
         nzpl = 0
         nthpl = 0
         npa=0

        if(sub)then
         npa=npla
         ntpl=npa
         nypl =npa
         nzpl =npa
         nthpl = npa
       
          if(mxpl.lt.4*npa)then
           print*,'too small dimensions in arrays'
           stop
          endif

         endif

         ntpl_s1=0
         nypl_s1 = 0
         nzpl_s1 = 0
         nthpl_s1 = 0

         ntpl_s2=0
         nypl_s2 = 0
         nzpl_s2 = 0
         nthpl_s2 = 0

*---------  control printing
          if(printout)then

        print *,'**************************************'
        print *,'*                                    *'
        print *,'*    Planes description in RECON:    *'
        print *,'*                                    *'
        print *,'**************************************'
        print *
        print *,'   first spectrometer - before SM1  '
        print *
          endif

         ifalse = 0                           
         ntpl_cor=npla 


*-------initialize vmat and vresol to zero..........
         do i=1,2*mxpl
         do j=1,2*mxpl
         vmat(i,j) = 0.
         enddo
         enddo
*----------------------------------------------------
         do i=1,npla
*
*---- "black list" of detectors 
*
 
*---- multiple coulomb scattering

          if(type(i).eq.111)dmcs=dmcs1
          if(type(i).eq.333)dmcs=dmcs3
          if(type(i).eq.222)dmcs=dmcs2    

          tetdet(i) =  0.0136  * sqrt(dmcs)*fudge 
         

           xxx(i) = zpla(i)/1000.     
           ccc(i)=-sinan(i)      ! coral hit convention -> recon hit convention
           sss(i)=cosan(i)

           mib(i)=.true.

*---- end of multiple c.s.
 
         e1_black=.true.

         if(i_option2.eq.0)go to 3 
         
          do ibra=1,n_black1
           if(uid(i).eq.list1_black(ibra))then
           e1_black=.false.
           mib(i)=.false.
           mib_lab(i)= uid(i)
           ifalse=ifalse+1
           endif  
          enddo

 3       continue



*--- planes on/off
          if(e1_black)then

             i_tp=1000
             ntpl=ntpl+1
             i_s1(ntpl)=i        
          xplane(ntpl)=zpla(i)/1000                ! position of plane
      

            if(xplane(ntpl).le.xtarg0)then
              print *
              print *,' planes no: ',i,' IN or BEFORE TARGET !!!!!'
              stop
            endif

          sgpln(ntpl)=spresol(i)/1000              ! resolution
 

          ityp1(ntpl)=type(i)                      ! type of detector (111-micros,222-SDC,333-scifi)
          eff_tel1(ntpl)=pleff(i)                  ! plane efficiency 

          angi=angel(i)

            c=-sinan(i)      ! coral hit convention -> recon hit convention
            s=cosan(i)

      if((s.lt.(-1.+eps).and.s.ge.-1.).or.(abs(c).lt.eps).or.
     &   (s.gt.(1.-eps).and.s.le.1.))i_tp=1
 
      if((c.gt.(1.-eps).and.c.le.1.).or.(abs(s).lt.eps).or.
     &   (c.lt.(-1.+eps).and.c.ge.-1.))i_tp=2
       
          sinn(ntpl)=s
          cosn(ntpl)=c

          iuid(ntpl)=uid(i)
           kmic=0

        if(i_tp.eq.1)then
           jyzt(ntpl) = 1                      ! Y plane
           nypl = nypl + 1
           jypl(nypl) = ntpl-npa
           kmic=1
        elseif(i_tp.eq.2)then
           jyzt(ntpl) = 2                      ! Z plane
           nzpl = nzpl + 1
           jzpl(nzpl) = ntpl-npa
           kmic=2
        elseif(c.gt.0.and.c.lt.1.and.i_tp.ne.1.and.i_tp.ne.2)then
           jyzt(ntpl) = 4                      ! -th plane
           nthpl = nthpl + 1
           jthpl(nthpl) = ntpl-npa
           kmic=4
        elseif(c.lt.0.and.c.gt.-1.and.i_tp.ne.1.and.i_tp.ne.2)then
           jyzt(ntpl) = 3                      ! th plane
           nthpl = nthpl + 1
           jthpl(nthpl) = ntpl-npa
           kmic=3
        endif      
    
           ind_fd(ntpl-npa)=0

           ntplpri=ntpl
       if(sub)ntplpri=ntpl-npa

       if(ityp1(ntpl).eq.111.or.ityp1(ntpl).eq.222)then

*----------- subsetup definition:
           if(sub)then
             ntpl_s1=ntpl_s1+1
             ind_fd(ntpl-npa)=1

           if(ityp1(ntpl).eq.111)then              ! second setup (micromegas only)
             ntpl_s2=ntpl_s2+1
             ind_fd(ntpl-npa)=2

             xplane(2*npa+ntpl_s2)=xplane(ntpl)
             sinn(2*npa+ntpl_s2)=sinn(ntpl)
             cosn(2*npa+ntpl_s2)=cosn(ntpl)
             jyzt(2*npa+ntpl_s2)=jyzt(ntpl)
             ityp1(2*npa+ntpl_s2)=ityp1(ntpl)
             iuid(2*npa+ntpl_s2)=iuid(ntpl) 
             sgpln(2*npa+ntpl_s2)=sgpln(ntpl)
           endif

             xplane(ntpl_s1)=xplane(ntpl)
             sinn(ntpl_s1)=sinn(ntpl)
             cosn(ntpl_s1)=cosn(ntpl)
             jyzt(ntpl_s1)=jyzt(ntpl)
             ityp1(ntpl_s1)=ityp1(ntpl)
             iuid(ntpl_s1)=iuid(ntpl)
             sgpln(ntpl_s1)=sgpln(ntpl)
 
             xplane(3*npa+ntpl_s1)=xplane(ntpl)
             sinn(3*npa+ntpl_s1)=sinn(ntpl)
             cosn(3*npa+ntpl_s1)=cosn(ntpl)
             jyzt(3*npa+ntpl_s1)=jyzt(ntpl)
             ityp1(3*npa+ntpl_s1)=ityp1(ntpl)
             iuid(3*npa+ntpl_s1)=iuid(ntpl) 
             sgpln(3*npa+ntpl_s1)=sgpln(ntpl)

             if(kmic.eq.3.or.kmic.eq.4)then
                nthpl_s1=nthpl_s1+1
                jthpl(nthpl_s1)=ntpl_s1
                jthpl(3*npa+nthpl_s1)=ntpl_s1
              if(ityp1(ntpl).eq.111)then
                nthpl_s2=nthpl_s2+1
                jthpl(2*npa+nthpl_s2)=ntpl_s2
              endif

             elseif(kmic.eq.2)then
                nzpl_s1 = nzpl_s1 + 1
                jzpl(nzpl_s1)=ntpl_s1
                jzpl(3*npa+nzpl_s1)=ntpl_s1
              if(ityp1(ntpl).eq.111)then
                nzpl_s2 = nzpl_s2 + 1
                jzpl(2*npa+nzpl_s2)=ntpl_s2
              endif

             elseif(kmic.eq.1)then
                nypl_s1 = nypl_s1 + 1
                jypl(nypl_s1)=ntpl_s1
                jypl(3*npa+nypl_s1)=ntpl_s1
              if(ityp1(ntpl).eq.111)then
                nypl_s2 = nypl_s2 + 1
                jypl(2*npa+nypl_s2)=ntpl_s2
              endif

             endif

        if(printout.and.ityp1(ntpl).eq.111)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  micro-a ",
     & i4," sub_set1/2")',
     & i,ntplpri,eff_tel1(ntpl),xplane(ntpl),
     & sinn(ntpl),cosn(ntpl),jyzt(ntpl),sgpln(ntpl),iuid(ntpl)
        endif

        if(printout.and.ityp1(ntpl).eq.222)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  Saclay DC",
     & i4," sub_set1")',
     & i,ntplpri,eff_tel1(ntpl),xplane(ntpl),
     & sinn(ntpl),cosn(ntpl),jyzt(ntpl),sgpln(ntpl),iuid(ntpl)
        endif
           endif
*--------- end of sub_set1

        if(printout.and..not.sub.and.ityp1(ntpl).eq.111)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  micro-a ",i4)',
     & i,ntplpri,eff_tel1(ntpl),xplane(ntpl),
     & sinn(ntpl),cosn(ntpl),jyzt(ntpl),sgpln(ntpl),iuid(ntpl)
      
        elseif(printout.and..not.sub.and.ityp1(ntpl).eq.222)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  Saclay DC ",i4)',
     & i,ntplpri,eff_tel1(ntpl),xplane(ntpl),
     & sinn(ntpl),cosn(ntpl),jyzt(ntpl),sgpln(ntpl),iuid(ntpl)
        endif

           elseif(ityp1(ntpl).eq.333)then

        if(printout)then
      print '("plane no. ",2i3," - ",4f8.4,i2,f9.6,"  Sci-fi ",i4)',
     & i,ntplpri,eff_tel1(ntpl),xplane(ntpl),
     & sinn(ntpl),cosn(ntpl),jyzt(ntpl),sgpln(ntpl),iuid(ntpl)
        endif

         endif


        endif
 
         enddo

         if(sub)then
           ntpl=ntpl-npa
           nypl=nypl-npa
           nzpl=nzpl-npa
           nthpl=nthpl-npa
           nypl_tmp=nypl
           nzpl_tmp=nzpl
           nthpl_tmp=nthpl

         endif

         ntpl_tmp=ntpl

*
*--------------- mag.field double integral........
*              (mod.RW. June 2002)
*

      Bdl = 0.
      Bdldl = 0.

          fd(1)=0.0
      if(sub)then
         fd(3*npa+1)=0.0
         fd(2*npa+1)=0.0
         fd(npa+1)=0.0
      endif

      xx(3) = xplane(1)
      ifds1=1
      ifds2=1

      do  50 k= 2,ntpl 

 49    xx(3) = xx(3) + 0.005
       xx(1) = 0.010 * (xx(3) - xplane(1))
       xx(2) = 0.005 * (xx(3) - xplane(1))

      call mag_field(xx,h_field)

       Bdl = bdl + 0.3 * h_field(2) * 0.005
       Bdldl = Bdldl + Bdl * 0.005
       if (xx(3).lt.xplane(k)) go to 49

*      print *,' z,H123:',xx(3),h_field
*      print *,'k,xplane,hy,bdl,bdl2:',k,xplane(k),h_field(2),
*     1 bdl,bdldl


      if(sub)then

       if(ind_fd(k).eq.1.or.ind_fd(k).eq.2)then
        ifds1=ifds1+1
        fd(ifds1)=bdldl
        fd(3*npa+ifds1)=bdldl
       endif
       if(ind_fd(k).eq.2)then
        ifds2=ifds2+1
        fd(2*npa+ifds2)=bdldl
       endif

       fd(npa+k) = bdldl
      elseif(.not.sub)then
       if(ind_fd(k).ne.0)then
         print*,'error in index for fd -> planes1 '
         stop
       endif
   
       fd(k) = bdldl
      
      endif

50    continue 

*-- just crosscheck

       if(sub)then
       if((ifds1-1).ne.ntpl_s1)then
       print*,'error in ifds1 -> planes1 ',ifds1,ntpl_s1
         stop
       endif
       if((ifds2-1).ne.ntpl_s2)then
       print*,'error in ifds2 -> planes1 ',ifds2,ntpl_s2
         stop
       endif
       endif

      print *
      print '(" no of planes before SM1 in coral: ",i5)',
     & npla
      print '(" no of planes behind SM1 off: " ,i5)',ifalse
      print '(" effective no of planes before SM1 in recon :",
     & 2i4," Y; ",i4," Z; ",i4," th ")',
     & ntpl,nypl,nzpl,nthpl
 
      if(sub)then
      print *,'nb of Y planes in set.1: ', nypl_s1
      print *,'nb of Y planes in set.2: ', nypl_s2
      print *,'nb of Y planes : ', nypl_tmp,nypl
      print *,'nb of Z planes in set.1: ', nzpl_s1
      print *,'nb of Z planes in set.2: ', nzpl_s2
      print *,'nb of Z planes : ', nzpl_tmp,nzpl
      print *,'nb of th planes in set.1: ', nthpl_s1
      print *,'nb of th planes in set.2: ', nthpl_s2
      print *,'nb of th planes : ', nthpl_tmp,nthpl
      endif
      print *
      print *,'Y planes: ',(jypl(npa+i),i=1,nypl)
      print *,'Z planes: ',(jzpl(npa+i),i=1,nzpl)
      print *,'th planes: ',(jthpl(npa+i),i=1,nthpl)
 
      print *,'real coral nb of planes for sci-fi ',
     &        'and micromegas : '
      print *,(i,i=1,ntpl) 
      print *,(i_s1(npa+i),i=1,ntpl)   
      print *,(iuid(npa+i),i=1,ntpl)
      print *                

     

         return
         end















