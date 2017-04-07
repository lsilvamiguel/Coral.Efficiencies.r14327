      subroutine planes3(npla3,zpla3,angel3,
     & sinan3,cosan3,spresol3,type3,pleff3,uid3,geom)   
**************************************************
*                                                *
*     transports planes description behind SM1   * 
*     and before RICH (for extrapolation)        *
*             from coral to recon                *
*     and calculate multiple coulomb scattering  *
*                  effect                        *
*                                                *
**************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "VERS.h"


       INTEGER npla3,type3(*),uid3(*)
       REAL zpla3(*),angel3(*),spresol3(*),pleff3(*),
     &      sinan3(*),cosan3(*),geom(3,*)

       logical printout,version,sub,fringe2

         eps=1.0e-3

* -------------  planes after Rich and befor SM2 ---------------
*                 used for track extrapolation 
*--------- chamber sizes needed for extrapolation to SM2

         size_y = geom(1,1)                   ! proportional chamber Y size
         size_z = geom(1,2)                   ! proportional chamber Z size
        size_ysf =geom(2,1)                   ! sci-fi           
        size_zsf =geom(2,2)                   ! sci-fi
        size_yg = geom(3,1)                   ! gem
        size_zg = geom(3,2)                   ! gem

         if(printout)then          
        print *
        print *,'   first spectrometer - between RICH & SM2'
        print *
        print *,' planes used for extrapolation '
        print '(" omega cham. B - sizes: Y & Z: ",2f6.2)',size_y,size_z
        print '(" sci-fi - sizes: Y & Z: ",2f6.2)',size_ysf,size_zsf
        print '(" gems - sizes: Y & Z: ",2f6.2)',size_yg,size_zg
        print *     
         endif


         i=0

         do ii=1,npla3

          i_tp=1000
          i=i+1
          iuid3(i)=uid3(ii)
                  
          xplan_sup(i)=zpla3(ii)/1000             ! position of plane

          sgpln_sup(i)=spresol3(ii)/1000          ! resolution
          ityp3(i)=type3(ii)                      ! type of detector 
          eff_tel3(i)=pleff3(ii)                  ! plane efficiency 

          angi=angel3(ii)



            c=-sinan3(ii)      ! coral hit convention -> recon hit convention
            s=cosan3(ii)

      if((s.lt.(-1.+eps).and.s.ge.-1.).or.(abs(c).lt.eps).or.
     &   (s.gt.(1.-eps).and.s.le.1.))i_tp=1
 
      if((c.gt.(1.-eps).and.c.le.1.).or.(abs(s).lt.eps).or.
     &   (c.lt.(-1.+eps).and.c.ge.-1.))i_tp=2
      

          sinn_sup(i)=s
          cosn_sup(i)=c


          if(i_tp.eq.1)then 
           jyzt3(i) = 1                      ! Y plane
           nypl3 = nypl3 + 1
           jypl3(nypl3) = i
          if(type3(i).eq.444)
     &     xplan_sup(i)=xplan_sup(i)- 0.1*sgpln_sup(i)

          elseif(i_tp.eq.2)then
           jyzt3(i) = 2                      ! Z plane
           nzpl3 = nzpl3 + 1
           jzpl3(nzpl3) = i
          if(type3(i).eq.444)
     &     xplan_sup(i)=xplan_sup(i)+ 0.1*sgpln_sup(i)  

          elseif(c.gt.0.and.c.lt.1)then
           jyzt3(i) = 4                      ! -th plane
           nthpl3 = nthpl3 + 1
           jthpl3(nthpl3) = i
          if(type3(i).eq.444)
     &     xplan_sup(i)=xplan_sup(i)- 0.1*sgpln_sup(i)


          elseif(c.lt.0.and.c.gt.-1)then
           jyzt3(i) = 3                      ! th plane
           nthpl3 = nthpl3 + 1
           jthpl3(nthpl3) = i
        if(type3(i).eq.444)
     &     xplan_sup(i)=xplan_sup(i)+ 0.1*sgpln_sup(i)  

          endif

        if(ityp3(i).eq.666.and.printout)then

        
        print '(" plane no. ",i3," - ",
     &  4f8.4,i2,f9.6," omega chamber ",i4)',
     &  i,eff_tel3(i),xplan_sup(i),
     &  sinn_sup(i),cosn_sup(i),jyzt3(i),sgpln_sup(i),
     &  iuid3(i)
  

        elseif(ityp3(i).eq.444.and.printout)then
    
        print '(" plane no. ",i3," - ",
     &  f8.4,f9.5,2f8.4,i2,f9.6," gem ",i4)',
     &  i,eff_tel3(i),xplan_sup(i),
     &  sinn_sup(i),cosn_sup(i),jyzt3(i),sgpln_sup(i),
     &  iuid3(i)
 

        elseif(ityp3(i).eq.333.and.printout)then
    
        print '(" plane no. ",i3," - ",
     &  4f8.4,i2,f9.6," sci-fi ",i4)',
     &  i,eff_tel3(i),xplan_sup(i),
     &  sinn_sup(i),cosn_sup(i),jyzt3(i),sgpln_sup(i),
     &  iuid3(i)
        endif

      enddo

      ntpl_sup=i

        if(printout)then
      print *,' total number of planes in Rich-SM2 region ',npla3
      print *,' out of them ',i,' will be used to search hits in Recon'
 
        endif

* ---------------- end of additional planes ------------------

      return
      end




