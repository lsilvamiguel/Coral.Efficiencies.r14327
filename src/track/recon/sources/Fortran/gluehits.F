      subroutine gluehits(ih1,ih2,ih3,ip1,ip2,ip3,
     &   phit1,phit2,phit3,err1,err2,err3) 

*********************************************************
*                                                       *
*   routine translates clusters from coral to "hits"    *
*                        in recon                       *
*        ( authors: Adam Mielech,Krzysztof Kurek &      *
*                        Ewa Rondio)                    *
*                                                       *
*             rewritten by KK - March 2002              *
*             rewritten by KK - July 2002               *
*                                                       *
*********************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "HITS_IN.h"
#include "VERS.h"


       INTEGER ih1,ih2,ih3
       INTEGER ip1(*), ip2(*), ip3(*)

       DOUBLE PRECISION phit1(*),phit2(*),phit3(*)
       DOUBLE PRECISION err1(*), err2(*), err3(*)

       logical version
       logical pat_hit,fhit,sub
       logical print_coral1,print_coral2
       logical lock,lock1,check1,check2
       logical print_hit1,print_hit2,print_hit3

   
        check1=.false.
        check2=.false.
  
        print_coral1=.false.
        print_coral2=.false.

        print_hit1=.false.
        print_hit2=.false.
        print_hit3=.false.

        if(print_coral1)then
          do i =1,ih1
        print *,'fortran: before', i,ip1(i),phit1(i),err1(i)
          enddo
        endif

        if(print_coral2)then
           do i =1,ih2
        print *,'fortran: after ', i,ip2(i),phit2(i),err2(i)
           enddo
        endif

        jmax=0
        pat_hit=.false.

*--- before SM1
           inhi=1              ! WARNING: assumption that hits are ordered in Z pos.
           fhit=.false.   
           nto=0
           nto2=0

           if(sub)then
             nto=ntpl_cor
             nto2=ntpl2_cor
           endif
*
*---- for debugging ------
*
       
         if(check1)then
           licz_clu1=0
           licz_ih1=0

           do i=1,ih1
           do k=1,ntpl_cor
           if(ip1(i).eq.mib_lab(k))licz_ih1=licz_ih1+1
           enddo
           enddo
           licz_ih1=ih1-licz_ih1
         endif
         if(check2)then
           licz_clu2=0
           licz_ih2=0

           do i=1,ih2
           do k=1,ntpl2_cor
           if(ip2(i).eq.mib_lab(ntpl_cor+k))licz_ih2=licz_ih2+1
           enddo
           enddo
           licz_ih2=ih2-licz_ih2
         endif


        do i=1,ntpl_tmp         ! loop over planes before SM1
         jchits(i)=0
         jchits_mc(i)=0

         if(i.lt.ntpl_tmp)idi=i_s1(nto+i)+1-i_s1(nto+i+1)
         inhi_bis=inhi

         do j=inhi,ih1
    
         if(i.eq.ntpl_tmp)fhit=.false.                     ! special case - last plane

         if(fhit.and.iuid(nto+i).ne.ip1(j).and.idi.eq.0)then
          go to 50
         elseif(fhit.and.iuid(nto+i).ne.ip1(j).and.idi.lt.0)then
          fhit=.false.
          go to 50
         endif

         if(.not.fhit.and.iuid(nto+i).ne.ip1(j))then                       
           inhi=inhi+1
           if(inhi.eq.ih1+1)then
            inhi=inhi_bis                                   ! to protect empty plane (reset counter)
            goto 50
           endif
           go to 51
         endif 

              inhi=inhi+1
              jchits(i)=jchits(i)+1                         ! nb. of clusters in the plane
              jchits_mc(i)=jchits_mc(i)+1

              if(jchits(i).gt.nmhit)then                    ! plane skip
                   jchits(i)=0
                   jchits_mc(i)=0
                   go to 50
              endif
              if(j.gt.nmghit)then
                print*,'dimension ERROR in cluster 1 '
                if(ityp1(nto+i).eq.111)print*,' in micros'
                if(ityp1(nto+i).eq.333)print*,' in sci-fis'
                if(ityp1(nto+i).eq.222)print*,' in sdc1 '
                pat_hit=.true.
                return
              endif

              clusters1(j)=phit1(j)/1000.
              clu_errors1(j)=err1(j)/1000.
              numhit(i,jchits(i))=j 
              numhit_mc(i,jchits(i))=j
 
              fhit=.true.

       if(print_hit1)
     &  print*,' HITS1', j,i,i_s1(nto+i),iuid(nto+i),ip1(j),
     &  clusters1(j),ityp1(nto+i),clu_errors1(j),numhit(i,jchits(i))

 51           continue
           enddo

 50        continue
  
      if(jmax.lt.jchits(i))jmax=jchits(i)

*      if(print_hit1) print*,i,ityp1(nto+i),jchits(i)
*        print*,i,ityp1(nto+i),jchits(i)

        if(check1)licz_clu1=licz_clu1+jchits(i)

           enddo

        if(check1.and.licz_clu1.ne.licz_ih1)then
            print*, 'clusters1 lost',licz_clu1,licz_ih1,ih1
            stop
        endif
*
*----------------- now second telescope
*

             inhi2=1
             fhit=.false.    
             lock=.false.
             inhi2_res=0
             inhi2_gem=0
             inhi2_bis=0
             igem=0
             igem1=0
             igem2=0

             licz_clu2=0                     ! tool for debugging - only if all planes active

       do i=1,ntpl2_tmp
          jc2hits(i)=0 
          jc2hits_mc(i)=0
 
          i_tr=i_s2(nto2+i)                   ! gives the plane number as ordered in coral
          if(i.lt.ntpl2_tmp)idi=i_s2(nto2+i)+1-i_s2(nto2+i+1)
          inhi2_bis=inhi2

          
*--- gems treatment - sort problem-July 5, 2002 kk


      if(ityp2(nto2+i).eq.444)then                                  ! protect gems    

        if(igem.eq.0.or.igem.eq.4.or.igem.eq.8)then
           igem1=iuid2(i_s2(nto2+i))
           igem2=iuid2(i_s2(nto2+i+1))
           inhi2_gem=inhi2
        if(igem.ne.0.and.inhi2_gem.lt.inhi2_res)inhi2_gem=inhi2_res
           lock=.true.
        endif

        if(igem.eq.2.or.igem.eq.6.or.igem.eq.10)then
        if(inhi2_gem.lt.inhi2_res)inhi2_gem=inhi2_res
           igem1=iuid2(i_s2(nto2+i))
           igem2=iuid2(i_s2(nto2+i+1))
        endif

           igem=igem+1
           inhi2=inhi2_gem
           lock1=.true.
      elseif(lock.and.ityp2(nto2+i).ne.444)then                     ! first plane behind gem
           if(inhi2_res.gt.inhi2)inhi2=inhi2_res                  
           lock=.false.
           lock1=.true.
           igem=0
           next=1
    
      endif
*
*---
*
     
*         print*,i,ityp2(nto2+i),inhi2,inhi2_gem,inhi2_res,lock,lock1

         do j=inhi2,ih2 
         
      if(lock.and.ip2(j).ne.igem2.and.ip2(j).ne.igem1)then          ! gems doublet protection
         if(lock1)go to 300
         inhi2_res=j
         go to 200
      endif 

      if(lock.and.lock1)then
         lock1=.false.
         inhi2_gem=j
      endif

      if(i.eq.ntpl2_tmp)fhit=.false.                              ! special case - last plane


         if(lock.and.iuid2(i_tr).ne.ip2(j))go to 300              ! protects gems 

         if(fhit.and.iuid2(i_tr).ne.ip2(j).and.idi.eq.0)then
          go to 200
         elseif(fhit.and.iuid2(i_tr).ne.ip2(j).and.idi.lt.0)then
          fhit=.false.
          go to 200
         endif


        if(iuid2(i_tr).ne.ip2(j).and..not.fhit)then
           inhi2=inhi2+1
           if(inhi2.eq.ih2+1)then
               inhi2=inhi2_bis                              ! to protect empty plane (reset counter)
               go to 200
           endif  
           goto 300
        endif 
            
              
              inhi2=inhi2+1

              jc2hits(i)=jc2hits(i)+1                       ! nb. of clusters in the plane
              jc2hits_mc(i)=jc2hits_mc(i)+1
 
              if(jc2hits(i).gt.nmhit)then                   ! skip plane
                 jc2hits(i)=0
                 jc2hits_mc(i)=0
                 go to 200
              endif
 
              if(j.gt.nmghit)then
                print*,'dimension ERROR in cluster 2  '
              if(ityp2(nto2+i).eq.555)print*,' in straws-ng'
              if(ityp2(nto2+i).eq.222)print*,' in SDC'
              if(ityp2(nto2+i).eq.444)print*,' in gems'
              if(ityp2(nto2+i).eq.333)print*,' in sci-fis'
                pat_hit=.true.
                return
              endif
 
              clusters2(j)=phit2(j)/1000.
              clu_errors2(j)=err2(j)/1000.


              numht2(i,jc2hits(i))=j
              numht2_mc(i,jc2hits(i))=j
              fhit=.true.
        
 
        if(print_hit2)then
        if(ityp2(nto2+i).eq.555)then
       print*,' straws-ng.',j,i,i_s2(nto2+i),iuid2(i_tr),ityp2(nto2+i),  
     & inhi2,jc2hits(i),clusters2(j),clu_errors2(j)
        elseif(ityp2(nto2+i).eq.444)then
          print*,' gems', j,i,i_s2(nto2+i),iuid2(i_tr),ityp2(nto2+i),
     &    inhi2,jc2hits(i),clusters2(j),clu_errors2(j)
        elseif(ityp2(nto2+i).eq.222)then
          print*,' SDCS', j,i,i_s2(nto2+i),iuid2(i_tr),ityp2(nto2+i),
     &    inhi2,jc2hits(i),clusters2(j),clu_errors2(j)
        elseif(ityp2(nto2+i).eq.333)then
          print*,' sci-fis', j,i,i_s2(nto2+i),iuid2(i_tr),ityp2(nto2+i),
     &   inhi2, jc2hits(i),clusters2(j),clu_errors2(j)
        endif
        endif

 
 300      continue
         enddo
 200     continue
           
      if(jmax.lt.jc2hits(i))jmax=jc2hits(i)

*      if(print_hit2) print*,i,ityp2(nto2+i),jc2hits(i)
*       print*,i,ityp2(nto2+i),jc2hits(i)

        if(check2)licz_clu2=licz_clu2+jc2hits(i)


       enddo 
   
        if(check2.and.licz_clu2.ne.licz_ih2)then
         print*, 'clusters2 lost',licz_clu2,licz_ih2,ih2
         stop
        endif


*
*   filling additional planes behind Rich and before SM2
               

             inhi3=1                     ! WARNING: assumption that hits are ordered in Z pos.        
             fhit=.false.    
             lock=.false.
             inhi3_res=0
             inhi3_gem=0
             inhi3_bis=0
             igem=0
             igem1=0
             igem2=0

       do i=1,ntpl_sup
           jcshits(i)=0
           inhi3_bis=inhi3

      if(ityp3(i).eq.444)then                                  ! protect gems    

        if(igem.eq.0.or.igem.eq.4.or.igem.eq.8)then
           igem1=iuid3(i)
           igem2=iuid3(i+1)
           inhi3_gem=inhi
        if(igem.ne.0.and.inhi3_gem.lt.inhi3_res)inhi3_gem=inhi3_res
           lock=.true.
        endif

        if(igem.eq.2.or.igem.eq.6.or.igem.eq.10)then
        if(inhi3_gem.lt.inhi3_res)inhi3_gem=inhi3_res
           igem1=iuid3(i)
           igem2=iuid3(i+1)
        endif

           igem=igem+1
           inhi3=inhi3_gem
           lock1=.true.
      elseif(lock.and.ityp3(i).ne.444)then                     ! first plane behind gem
           if(inhi3_res.ne.0)inhi3=inhi3_res                  
           lock=.false.
           lock1=.true.
           igem=0
    
      endif
*  
       do j=inhi3,ih3

      if(lock.and.ip3(j).ne.igem2.and.ip3(j).ne.igem1)then          ! gems doublet protection
         if(lock1)go to 350
         inhi3_res=j
         go to 360
      endif 

      if(lock.and.lock1)then
         lock1=.false.
         inhi3_gem=j
      endif

          if(i.eq.ntpl_sup)fhit=.false.                             ! special case - last plane
          if(lock.and.iuid3(i).ne.ip3(j))go to 350                  ! protects gems 

          if(fhit.and.iuid3(i).ne.ip3(j))go to 360

          if(iuid3(i).ne.ip3(j).and..not.fhit)then
           inhi3=inhi3+1
           if(inhi3.eq.ih3+1)then
               inhi3=inhi3_bis                              ! to protect empty plane (reset counter)
               go to 360
           endif  
           goto 350
          endif 

              inhi3=inhi3+1 
              jcshits(i)=jcshits(i)+1               ! nb. of clusters in the plane

              if(jcshits(i).gt.nmhit)then           ! skip plane
               jcshits(i)=0
               go to 360
              endif

              if(j.gt.nmghit)then
                print*,'dimension ERROR in cluster 3 '
                pat_hit=.true.
                return
              endif   

              clusters3(j)=phit3(j)/1000.
              clu_errors3(j)=err3(j)/1000.
              numhts(i,jcshits(i))=j  
              fhit=.true.            

       if(print_hit3)
     & print*,' HITS3',j,i,clusters3(j),ityp3(i),clu_errors3(j)   

 350           continue
           enddo
 360     continue

       if(print_hit3)print*,i,ityp3(i),jcshits(i)

       enddo

      return
      end         
       




















