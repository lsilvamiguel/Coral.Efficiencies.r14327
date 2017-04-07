      subroutine tkfit_frext(k_up,nhtt,htt,errr,jpl,xref,xpar,prod,
     1 chi2,np_eff, xref_end,par_end,ifail,traj)
*************************************************************
*                                                           *
*                EXTENDED  VERSION ("frext") of             *
*                 spline fit for projections                *
*    using  upstream detectors and 4 to 8 hits in DC02      *
*                   RW, Aug.30,2002                         *
*           ================================                *
*                                                           *
*   -- uses hits from DC2 which are close enough to track-  *
*       (min. 4 extra hits required)                        *
*                                                           *
*************************************************************
#include "PARAM.h"
#include "PLANES.h"
#include "HITS_OUT.h"
#include "HITS_IN.h"

          parameter (mpl2=mxpl)
          parameter (nps=100)
          parameter (nps2=50)
          parameter (nps3=10)

          dimension htt(8),errr(8),jpl(8)                  ! input arguments
          dimension xpar(5),err(5)
          dimension par_end(4)
          dimension traj(*),vresol(mxpl)
      
          DOUBLE PRECISION xx(3),h_field(3)         ! do not change !!!

          dimension pys(nps),pzs(nps),xnty1(nps),
     1     xnty2(nps),xntz1(nps),xntz2(nps)
          dimension dfy1(nps),dfy2(nps),dfz1(nps),dfz2(nps)
          dimension gy(mpl2,mpl2),a(mpl2,5),at(5,mpl2),cc(mpl2,1)
          dimension prov1(mpl2,5),prod(5,5),work(nps)
          dimension prov2(mpl2,1),prov3(5,1),xxp(5,1)
          dimension xxxx(nps),yyy(nps),zzz(nps)
          dimension dry1(nps),dry2(nps),drz1(nps),drz2(nps)
          dimension xadd(nps),yadd(nps),zadd(nps),dya1(nps),dya2(nps)
          dimension dza1(nps),dza2(nps),xnay1(nps),xnay2(nps),
     1     xnaz1(nps),xnaz2(nps)
          dimension wkspace(nps)  
          dimension kref(4) 
          dimension jphys(nps2),jphy2(nps3),jreal(nps)

*----for chi2 calculation with MCS correlated errors.......................
          dimension vvv(1,mpl2),vvt(mpl2,1),gyv(1,mpl2)
          dimension chi_2(1,1)
   
          external erscal
*-----------------------------------------------------------------
           do i=1,ntpl+ntpl2
           if(i.le.ntpl)then
             vresol(i) = sgpln(i)**2
           else
             vresol(i) = sgpln2(i-ntpl)**2
           endif
           enddo
*--------- above kk June 25 2002


           ifail = 0

          lspline=0              ! internal debugger for spline; 0 off 1 on
      
*----MCS applied only for P < plow_mcs......................................
      
  
*-----Set trajectory positions ("traj") to -8888.  by default...............

            do kp = 1,mxplan
            traj(kp) = -8888.
            enddo

            kref(1)=1                  ! nr. of points in upper part of LAS
            kref(2)=40     
            kref(4) = nps
            kref(3) = kref(4)-nhtt-1   ! account for   hits in SDC2 planes


           do jo=1,nps
            if (jo.le.nps2) jphys(jo)=0
            if (jo.le.nps3) jphy2(jo)=0
            jreal(jo)=0
           enddo
  
**
*---- find nr of used planes in  upper part..............
*
          nplu1 = 0

          do kpl = 1,ntpl

         if(line1(kpl,k_up).eq.8888)then
          tth=800.
         else
          tth = clusters1(line1(kpl,k_up))
         endif

          if (tth.ne.0.and.tth.lt.700.) then
            nplu1 = nplu1 +1
            jphys(nplu1) = kpl
      
          if(lspline.eq.1)then
            print '("np,kpl,x,tth",2i5,2f10.5)',nplu1,jphys(nplu1),
     1    xplane(kpl),tth
          endif

          endif
          enddo

*-------Plane numbers for extension....stored in jpl and already ordered....
           isup = 0
           do ks=1,8    
           if (jpl(ks).gt.0) then
            isup = isup + 1
            jphy2(isup)  = jpl(ks)
           endif 
           enddo
          nplu2 = nhtt




*
*---- x position of ALL  planes....................
*
          xadd(kref(1)) = xplane(1) - 0.001         ! modified RW Febr. 8
          xref=xadd(kref(1))

          xadd(kref(2)) = xplane(ntpl) + 0.01
          neff_up = kref(2) - kref(1) - nplu1
          dx1 = (xadd(kref(2)) - xadd(kref(1))) / neff_up

          xadd(kref(3)) = xplan2(1) - 0.01
          xadd(kref(4)) = xplan2(8) + 0.01
          xref_end = xadd(kref(4))

          dx_end = xadd(kref(4)) - xadd(kref(1))    
          dx2 = (xadd(kref(3)) - xadd(kref(2)))/(kref(3)-kref(2))

            
*
*----upstream part.....................
          do k = 1,nplu1
          xadd(k+1) = xplane(jphys(k))
          jreal(k+1) = jphys(k)

          if(lspline.eq.1)then
          print '("1=k,jph,x",2i5,f10.5)',k+1,
     1      jreal(k+1),xadd(k+1)
          endif

          enddo
*----now nhtt planes behind SM1..............
          do k= 1,nhtt
           xadd(k+kref(3)) = xplan2(jphy2(k))
           jreal(k + kref(3)) = jphy2(k)
          enddo
         
*
*----- now add  intermediate planes...............................

          do k = 1,neff_up  - 1        
          kout = kref(1) + nplu1 + k
          xadd(kout) = xadd(kref(1))+ k*dx1
          jreal(kout) = 0
          enddo

          do k= 1,kref(3)-kref(2) -1
          kout = kref(2) + k
          xadd(kout) = xadd(kref(2)) + k* dx2
          jreal(kout) = 0
          enddo
         
          if(lspline.eq.1)then
          do k=1,kref(2)
          print '("k,jph,x",2i5,f10.5)',k,jreal(k),xadd(k)
          enddo
          endif

*----- order these planes with increasing x.........
*
11        jperm = 0
          do k=1,kref(4)-1
          if (xadd(k+1).lt.xadd(k)) then
            xtt= xadd(k+1)
            jtt= jreal(k+1)
            xadd(k+1) = xadd(k)
            jreal(k+1)= jreal(k)
            xadd(k) = xtt
            jreal(k)= jtt
            jperm = 1
          endif
          enddo
          if (jperm.eq.1) go to 11
*
*---- select y planes for smoothing of slope
*
          nply = 0
          do  k=1,kref(4)
          if (jreal(k).gt.0) then     
            if (xadd(k).le.xplane(ntpl).and.jyzt(jreal(k)).eq.1) then
            nply = nply + 1
            xxxx(nply) = xadd(k)
            yyy(nply) = clusters1(line1(jreal(k),k_up))

           elseif (xadd(k).ge.xplan2(1).and.jyzt2(jreal(k)).eq.1) then
             nply = nply + 1
             xxxx(nply) = xadd(k)
             yyy(nply) = htt(jreal(k)) 
           endif 

           if(lspline.eq.1)then
            print '("nply-k-jp-x-y",3i5,2f10.5)',nply,k,
     1     jreal(k),xxxx(nply),yyy(nply)
           endif

          endif 
          enddo

          
*-----------------------------------------------
          call spnfit(nply,xxxx,yyy,dry1,dry2,xnty1,xnty2,iret)
           if (iret.lt.0) then
            ifail = -10
            return
           endif

          if(lspline.eq.1)then
          print *, 'nply',nply
          do k=1,nply
          print '(" k,xxx,yyy,dy:",i4,3f10.5)',k,xxxx(k),yyy(k)
     1     ,dry1(k)
          enddo
          endif
          

          do kv=1,kref(4)
          call spnpol(nply,xxxx,yyy,dry1,dry2,xadd(kv),yadd(kv))

          if(lspline.eq.1)then
**       if (kv.le.2.or.kv.eq.10.or.kv.ge.58) then
           print '("kv,xadd,yadd",i4,2f10.4)',kv,xadd(kv),yadd(kv)
**       endif
          endif

          enddo
*
*---- select z planes for smoothing of slope
*
          nplz = 0
          do 21 k=1,kref(4)
          if (jreal(k).le.0) go to 21
           if (xadd(k).le.xplane(ntpl).and.jyzt(jreal(k)).eq.2) then
            nplz = nplz + 1
            xxxx(nplz) = xadd(k)
            zzz(nplz) = clusters1(line1(jreal(k),k_up))

           elseif (xadd(k).ge.xplan2(1).and.jyzt2(jreal(k)).eq.2) then
            nplz = nplz + 1
            xxxx(nplz) = xadd(k)
            zzz(nplz) = htt(jreal(k)) 
           endif
           if(lspline.eq.1)then
            print '("nplz-k-jp-x-z",3i5,2f10.5)',nplz,k,
     1     jreal(k),xxxx(nplz),zzz(nplz)
           endif
 21       continue
*
          
          call spnfit(nplz,xxxx,zzz,drz1,drz2,xntz1,xntz2,iret)
            if (iret.lt.0) then
             ifail = -11
             return
            endif

           if(lspline.eq.1)then
          print *, 'nplz',nplz
          do k=1,nplz
          print '(" k,xxx,zzz,dz:",i4,3f10.5)',k,xxxx(k),zzz(k)
     1      ,drz1(k)
          enddo
           endif
*

          do kv=1,kref(4) 
          call spnpol(nplz,xxxx,zzz,drz1,drz2,xadd(kv),zadd(kv))

           if(lspline.eq.1)then
         if (kv.le.2.or.kv.eq.10.or.kv.ge.58) then
           print '("kv,xadd,zadd",i4,2f10.4)',kv,xadd(kv),zadd(kv)
         endif
           endif

          enddo

*--------------Slopes at intermediate points...................................
          call spnfit(kref(4),xadd,yadd,dya1,dya2,xnay1,xnay2,iret)
           if (iret.lt.0) then
             ifail = -12
             return
           endif

          call spnfit(kref(4),xadd,zadd,dza1,dza2,xnaz1,xnaz2,iret)
           if (iret.lt.0) then
             ifail = -13
             return
           endif


           if(lspline.eq.1)then
          do kv=1,kref(4) 
**       if (kv.le.2.or.kv.eq.10.or.kv.ge.58) then
        print '("kv,x-y-dya1-z-dza1",i4,5f10.4)',kv,xadd(kv),
     1    yadd(kv),dya1(kv),zadd(kv),dza1(kv)
**       endif
         enddo
           endif

*---- iteration starts from here................
*
          iterat = 1
*---- there are nhtt hits behind SM1...................................
          nn12 = nplu1 + nhtt

*-------data and errors.............................................
          do k1=1,nn12
          do k2=1,nn12
          gy(k1,k2) = 0.
          enddo
          enddo
          jv = 0
          do kp = 1,kref(4)
          if (jreal(kp).ne.0) then

           jv = jv + 1
           if (jv.le.nplu1)  then

           if (line1(jreal(kp),k_up).ge.8888) then
           cc(jv,1) = 800.
           gy(jv,jv) = 1./sgpln(jreal(kp))**2           
           else
           cc(jv,1) = clusters1(line1(jreal(kp),k_up))
           gy(jv,jv) = 1./(clu_errors1(line1(jreal(kp),k_up))
     1       *erscal(1,jreal(kp)))**2
           endif
             
          elseif (kp.gt.kref(3)) then

*---TEMP....
**         print '("kp,hit,err",i4,2f10.5)',kp,htt(jreal(kp)),
**   1      errr(jreal(kp))
*---end TEMP
           cc(jv,1) = htt(jreal(kp))
           gy(jv,jv) = 1./(errr(jreal(kp))*erscal(2,jreal(kp)))**2
           endif
          
          endif
          enddo
80        continue

*--------------Now spline and double integral ..........................
          do kv = 1,kref(4)
          xx(3) = xadd(kv)                              
          xx(1) = yadd(kv)                 
          xx(2) = zadd(kv)                 
          call mag_field(xx,h_field)
          fac = 0.3*sqrt(1.+dya1(kv)**2 + dza1(kv)**2)

          if(lspline.eq.1)then
         if (kv.le.2.or.kv.eq.10.or.kv.ge.58) then
           print '("kv,xx,hf,fac",i4,7f10.4)',kv,xx,h_field,fac
         endif
          endif

          pys(kv) = fac*(h_field(3)*dza1(kv) + h_field(1) *
     1     dya1(kv) * dza1(kv) - h_field(2) *
     2     (1. + dya1(kv)**2))
          pzs(kv) = fac*(-h_field(3) * dya1(kv) - h_field(2) *
     1     dya1(kv) * dza1(kv) + h_field(1) *
     2     (1. + dza1(kv)**2))

          if(lspline.eq.1)then
         if (kv.le.2.or.kv.eq.10.or.kv.ge.58) then
           print '("kv,pys,pzs",i4,2f10.4)',kv,pys(kv),pzs(kv)
         endif
          endif

          enddo

          call spnfit(kref(4),xadd,pys,dfy1,dfy2,xnty1,xnty2,iret)
           if (iret.lt.0) then
             ifail = -14
             return
           endif

          call spnfit(kref(4),xadd,pzs,dfz1,dfz2,xntz1,xntz2,iret)
           if (iret.lt.0) then
             ifail = -15
             return
           endif

*----------------------------------------------------------------------

          if(lspline.eq.1)then
          do kk=1,kref(4) 
          print  '(i4,"x = ",f5.3,"int2 y,z:",2f10.5)',
     1     kk,xadd(kk),xnty2(kk),xntz2(kk)
          enddo
          endif

*
*------ now do least square fit-------------------------------------
*
*----- here we use the hits in all planes  and     
*       their errors sgpln(2). Values are retrieved using
*       pointer jreal(plane nr).   
*
*---- fill matrices.................................................
*

*-------coefficients ...............................................
*
         do k1=1,nn12
         do k2=1,5
         a(k1,k2) = 0.
         enddo
         enddo

         do k1 = 1,nplu1     
         a(k1,1) = sinn(jphys(k1))
         a(k1,2) = -cosn(jphys(k1))
         enddo

         do k1 = 1,nhtt
         a(k1+nplu1,1) = sinn2(jphy2(k1))
         a(k1+nplu1,2) = -cosn2(jphy2(k1))
         enddo
*

*-------More coefficients ...............................................

         kk = 0

         do k1=1,kref(4)         
         if (jreal(k1).ne.0)  then 
          kk = kk + 1
          work(kk) = xadd(k1)
          a(kk,5)= a(kk,1)*xnty2(k1)
     1        +a(kk,2)*xntz2(k1)

         if(lspline.eq.1)then
        print '("k1,kk,nty-z,a5",2i5,3e13.4)',k1,kk,
     1    xnty2(k1),xntz2(k1),a(kk,5)
         endif

         endif
         enddo

*
         do jk=1,nn12
         work(jk) = work(jk) - xref
         enddo
*------------------------------------------------------
         do jk = 1,nn12           
         a(jk,3) = a(jk,1) * work(jk)
         a(jk,4) = a(jk,2) * work(jk)
         enddo

          if(lspline.eq.1)then
         do jk = 1,nn12
         print '("jk,work,a",i3,6f10.5)',jk,work(jk),
     1   (a(jk,jr),jr=1,5)
         enddo

         do jr = 1,nn12
         print '(" jr,gy",i3,e10.4)',jr,gy(jr,jr)
         enddo
          endif

*-------some multiplications.........................................

         call matmu2(nn12,nn12,5,gy,a,prov1,mpl2,mpl2,mpl2)
         call trans2(nn12,5,a,at,mpl2,5)
         call matmu2(5,nn12,5,at,prov1,prod,5,mpl2,5)

          if(lspline.eq.1)then
         do jr = 1,5
         print '("prod",5e10.4)',(prod(jr,jjr),jjr=1,5)
         enddo
          endif

         ifail = 0
****------- name changed to wkspace to avoid conflict <===============
         call rinv(5,prod,5,wkspace,ifail)
          if (ifail.ne.0)then
              print *,' rinv--ifail:',ifail

      
          if (lspline.eq.1) then
           do jk = 1,nn12
           print '("jk,wkspace,a,cc",i3,7f10.5)',jk,work(jk),
     1     (a(jk,jr),jr=1,5),cc(jk,1)
           enddo
          endif

           return                          
          endif

                         
         call matmu2(nn12,nn12,1,gy,cc,prov2,mpl2,mpl2,mpl2)
         call matmu2(5,nn12,1,at,prov2,prov3,5,mpl2,5)
         call matmu2(5,5,1,prod,prov3,xxp,5,5,5)
          do k=1,5
          xpar(k) = xxp(k,1)
          enddo

*-------results and errors...........................................
          do kk=1,5
           if (prod(kk,kk).le.0.) then 
              print *,' cov. diag.negative',k
              prod(kk,kk) = abs(prod(kk,kk))
           endif
          err(kk) = sqrt(prod(kk,kk))

          if(lspline.eq.1)then
          print '(" Param.",i3,f8.5,"+/-",f8.5)',kk,xpar(kk),err(kk)
          endif

          enddo
*------- make chi2 big if momentum undefined.........................
         if(iterat.eq.2.and.err(5)/xpar(5).gt.0.5) then
           chi_2(1,1) = chi_2(1,1) + 10000.
         endif

*-------chi2.........................................................
         
         v_mean = 0.
         np_eff = 0
         chi2t = 0.

         do  889 kp=1,nn12

         vvv(1,kp) = 0.

         if (a(kp,1).eq.0.0.and.a(kp,2).eq.0.0) go to 889

         vcal = a(kp,1)*xpar(1) + a(kp,2)*xpar(2) + a(kp,3) * xpar(3) +
     1     a(kp,4) * xpar(4) + a(kp,5)*xpar(5) 

*---- kk- error - corrected August 1, 2002
         
         kkp=kp
         if (kp.le.nplu1.and.jphys(kp).ne.0)then
           kkp = jphys(kp)
         elseif(kp.gt.nplu1.and.jphy2(kp).ne.0)then
           kkp = jphy2(kp-nplu1) + ntpl_tmp
         endif
         traj(kkp) = vcal

         vvv(1,kp) = cc(kp,1) - vcal
         np_eff = np_eff + 1
         chi2p = (vvv(1,kp)**2) * gy(kp,kp)
         chi2t = chi2t + chi2p
         v_mean = v_mean + vvv(1,kp)**2

         if(lspline.eq.1)then
         tterr = 1./sqrt(gy(kp,kp))
        print '("kp,a1-2-5,cc,err,vcal,chi2p,chi2T",i4,6f10.5,
     1    2e13.4)',kp,
     2   a(kp,1),a(kp,2),a(kp,5),cc(kp,1),tterr,vcal,chi2p,chi2t
         endif
 889     continue

         call trans2(1,nn12,vvv,vvt,1,mpl2)
         call matmu2(1,nn12,nn12,vvv,gy,gyv,1,mpl2,1)
         call matmu2(1,nn12,1,gyv,vvt,chi_2,1,mpl2,1)
         ppp = 1./xpar(5)

         if(lspline.eq.1)then
          print '("it-CHI2 =",i2,f10.2," for ",i3," planes,p="f6.2)',
     1    iterat,chi_2(1,1),nn12,ppp
         endif

*
*------transfer to chi2 for test in calling routine (s/r assoc)............
*
        chi2 = chi_2(1,1)

*--- get updated dya1,dza1 from fitted trajectory......................
*
          do kv = 1,kref(4)
          dya1(kv) = xpar(3) + xpar(5) * xnty1(kv)
          dza1(kv) = xpar(4) + xpar(5) * xntz1(kv)
          yadd(kv) = xpar(1) + (xadd(kv)-xref)*xpar(3)
     1       + xpar(5)*xnty2(kv)
          zadd(kv) = xpar(2) + (xadd(kv)-xref)*xpar(4)
     1        + xpar(5)*xntz2(kv)
          enddo
           par_end(3) = dya1(kref(4))
           par_end(4) = dza1(kref(4))
*------ iterate with different gy's to account for MCS------------------
        if (iterat.eq.1) then
          iterat = 2
          pb = xpar(5)
          if (abs(pb).lt.1.0e-10) pb = 1.0e-10

          
*------  no iteration if p> plow_mcs ...................................
          if(plow_mcs.lt.1.0e-10)plow_mcs=1.0e-10
          if (abs(pb).lt.1./plow_mcs) go to 700

*----- RW - modif. April 19,2002
          e_pi = sqrt(1./pb**2 + 0.139**2)
          betab = e_pi * pb
*----------- select relevant part of vmat........................
888       continue
*         print *,' jpb =',jpb
          do ii= 1,nn12
           do jj=1,nn12

*--- kk -error- corrected August 1 2002
*--- remodified Aug.30...RW...after redefining nn12 + modif. in proj_fringe
* 
           if (ii.le.nplu1) then
            iii = jphys(ii)
            irs = 1
            irsa = iii
           else
            iii = jphy2(ii-nplu1) + ntpl
            irs = 2
            irsa = iii-ntpl
           endif
           if (jj.le.nplu1)then
            jjj = jphys(jj)
           else
            jjj = jphy2(jj-nplu1) + ntpl
           endif

            gy(ii,jj) = vmat(iii,jjj) * (pb * betab)**2
*------- and add contribution from vresol to diagonal terms
           if(ii.eq.jj) gy(ii,jj) = gy(ii,jj)  
     1      + vresol(iii)*(erscal(irs,irsa)**2)

*           print*,'gy ',ii,jj,gy(ii,jj),iii,jjj,vmat(iii,jjj)

          enddo
          enddo

          
          call rinv(nn12,gy,mpl2,wkspace,ifail)
           if (ifail.ne.0)  then
                print *,' rinv-MCS....ifail 2',ifail
                return
           endif

          go to 80
        endif

*------ end point ........................................................

700     continue
        par_end(1) = xpar(1) + dx_end * xpar(3) + pb * xnty2(kref(4))
        par_end(2) = xpar(2) + dx_end * xpar(4) + pb * xntz2(kref(4))


     
          return
          end


















