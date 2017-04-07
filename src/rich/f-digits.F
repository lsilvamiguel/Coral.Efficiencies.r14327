c-----------------------------------------------------------------------
      SUBROUTINE digits( nhit, ichit, xhit, yhit, ehit, 
c1900     +                   npatt, icto, ixto, iyto, phto )
     +                   npatt, icto, ixto, iyto, phto, ihcco )
c--------------------------------------------------------------
c
      IMPLICIT NONE
c
c------------------------------
#include "finclude/pads.inc"
#include "finclude/clu.inc"
#include "finclude/elpulse.inc"
#include "finclude/enois.inc"
c------------------------------
c
      INTEGER  nhit
      INTEGER  ichit(1)
      REAL  xhit(1), yhit(1), ehit(1)
      INTEGER  npatt
      INTEGER  icto(1), ixto(1), iyto(1)
      REAL  phto(1)
c0104      INTEGER  ihcco(1)                    !   1900
      INTEGER*8  ihcco(1)                    !   0104
c
      INTEGER  khit, jk, kp
      INTEGER  ndph, nmip, nfeed, nphevall
      INTEGER  nphotonsall, ndphall, nmipall, nfeedall
      INTEGER  icath
      REAL  x, y, en
c
      INTEGER  ictw, ixtw, iytw
      INTEGER  kpadin, jpadin, npadin
ccc      INTEGER  knoise
c
      REAL  hcathx, hcathy
ccc      REAL  ddd
      INTEGER  NPHEV
      PARAMETER  (NPHEV = 3000)
      INTEGER  ictin(NPHEV), ixtin(NPHEV), iytin(NPHEV)
c
ccc      INTEGER  kdophgen
c
      INTEGER  kcall
c
      DATA  kcall / 0 /
c
c
      if(kcall.eq.0) then
      kcall = 1
c
      nphotonsall = 0
      ndphall = 0
      nmipall = 0
      nfeedall = 0
c
      hcathx = padx * npadx / 2.
      hcathy = pady * npady / 2.
c
      endif      !   endif kcall
c
c
      do   kp=1,NELMAX
        ihitcc(kp) = 0
      enddo
c
      ndph = 0
      nfeed = 0
      nmip = 0
c
      kpadin = 0
      npadin = 0
c
c
      do   khit=1,nhit
c
c---- process hit (= photon) :
c     ------------------------
      x  = xhit(khit)
      y  = yhit(khit)
      en = ehit(khit)
      icath = ichit(khit)
cpp         print *,khit,ichit(khit),xhit(khit),yhit(khit)
c
c
c---- save incident cerenkov photon hit pads :
c     ----------------------------------------
c---- (no double counting)
      ictw = icath
      ixtw = (x + hcathx) / padx + 1
      iytw = (y + hcathy) / pady + 1
      jpadin = 0
      do   kp=1,kpadin
      if(ictw.eq.ictin(kp) .and. 
     +   ixtw.eq.ixtin(kp) .and. iytw.eq.iytin(kp)) jpadin = 1
         enddo      !   kpadin
      if(jpadin.eq.0) then
      kpadin = kpadin + 1
      ictin(kpadin) = icath
      ixtin(kpadin) = ixtw
      iytin(kpadin) = iytw
      endif      !   jpadin
      npadin = kpadin
c
c
c---- process enoise only :
c--------------------------
      if(knoise.eq.-1) then
      x = x + 10000.
      y = y + 10000.
c@@-----------------
      endif
c
cpp         print *,x,y,1,en,icath,ndph,nmip,nfeed
c
c1900      call phtrack(x,y,1,en,icath,ndph,nmip,nfeed)
      call phtrack(khit,x,y,1,en,icath,ndph,nmip,nfeed)
c------------------------------------------------------
c
cpp      print *,ndph,nmip,nfeed
c
      nphevall = ndph + nfeed + nmip
c
         enddo      !   khit
c
      ndphall = ndphall + ndph
      nmipall = nmipall + nmip
      nfeedall = nfeedall + nfeed
      nphotonsall = nphotonsall + nhit
c
c>>   kdophgen = 0,1
c-------------------
      if(kdophgen.eq.0) then
c
      call elnoise                      !   check rgau = 0.
c-----------------
      call threshold( nphevall, npatt )
c--------------------------------------
c
      else       !   kdophgen
c
      npatt = npadin
      do   jk=1,npatt
      ict(jk) = ictin(jk)
      ixt(jk) = ixtin(jk)
      iyt(jk) = iytin(jk)
      pht(jk) = 50.
      ihcc(jk) = jk                      !   1900
         enddo      !   npatt
c
      endif      !   kdophgen
c
      do   jk=1,npatt
      icto(jk) = ict(jk)
      ixto(jk)  = ixt(jk)
      iyto(jk)  = iyt(jk)
      phto(jk) = pht(jk)
      ihcco(jk) = ihcc(jk)               !   1900
         enddo      !   npatt
c
c
c---- histo pad occupancy :
c     ---------------------
ccc      call hfill(1039,float(npatt),ddd,1.)
c
c
 4000 continue
c
      END
c
c
c--------------------------------------------------------------------
      SUBROUTINE initall(npadxin,padxin,npadyin,padyin,ncathin,
     +                   kdophgenin,knoisein,
ccc     +                   padnamein,
     +                   csiqein, alfain, phsloin, ridfactin,
     +                   dofeedbackin, doreflin, doQEtestin,
ccc     +                   nwiresin, nwirepadin, pitchin,
     +                   pedmin, sigpin, sigmin, sigsin, threshin,
     +                   dosigcompin )
c--------------------------------------------------------------------
c
      IMPLICIT NONE
c
c------------------------------
#include "finclude/pads.inc"
#include "finclude/enois.inc"
#include "finclude/elpulse.inc"
#include "finclude/baspar.inc"
c------------------------------
c
      INTEGER  npadxin, npadyin, ncathin
      REAL  padxin, padyin
c
      INTEGER  kdophgenin, knoisein, dosigcompin
c
      REAL  csiqein(25), alfain(16), phsloin(16), ridfactin(16)
      INTEGER  dofeedbackin, doreflin, doQEtestin
ccc      INTEGER  nwiresin, nwirepadin
ccc      REAL  pitchin
      REAL  pedmin, sigpin, sigmin, sigsin, threshin
ccc      CHARACTER*13  padnamein
ccc      REAL  stepin
ccc      INTEGER  npin
c
ccc      REAL  padim
      INTEGER  np, nstep, ntp, i, j, k
ccc      CHARACTER*13  padname
c
c
      character*128 filename
      character*256 expfilename
      integer       status
      integer       expnm
      external      expnm
c
c
c------ CsI QE for Eph 5.5-7.9 eV (csiqe3.dat file from FRESIM)
c
ccc      DATA  csiqe /0.00025,0.00035,0.0005,0.005,0.0075,0.0125
ccc     +  ,0.01775,0.0235,0.04,0.064,0.094,0.116,0.1585,0.192,0.207
ccc     +  ,0.23,0.2525,0.27,0.2845,0.3035,0.32,0.34,0.345,0.36,0.37/
c
c------ feedback rate in ADC^-1 for 16 cathodes
c
ccc      DATA  alfa /16*0.03/
c
c------ mean Pulse Heihgt
c
ccc      DATA  phslo /16*30./
c
c----- overall QE reduction factor for 16 cathodes
c
ccc      DATA  ridfact /16*1./
c
c-------------------------------------
      kdophgen = kdophgenin
      knoise = knoisein
c
c@@-----------------------------
      dofeedback = dofeedbackin     ! fdbck ph. generation (1=ON)
      dorefl   = doreflin           ! fdbck ph. scat. on CsI surface (1=ON)
      doQEtest = doQEtestin         ! QE test for Cherenkov photons (1=ON)
ccc      dofeedback = 1           ! fdbck ph. generation (1=ON)
ccc      dorefl = 1               ! fdbck ph. scat. on CsI surface (1=ON)
ccc      doQEtest = 0             ! QE test for Cherenkov photons (1=ON)
c@@-----------------------------
c
      do   k=1,25
      csiqe(k) = csiqein(k)
         enddo
      do   k=1,16
      alfa(k) = alfain(k)
      phslo(k) = phsloin(k)
      ridfact(k) = ridfactin(k)
         enddo
c
      padx = padxin
      npadx = npadxin               ! no. of pads along x
      cathx2 = padx * npadx  / 2.   ! cathode half dimension in x
      pady = padyin
      npady = npadyin               ! no. of pads along y
      cathy2 = pady * npady  / 2.   ! cathode half dimension in y   
      ncath = ncathin               ! number of cathodes
ccc      padim = 8        ! pad dimension (mm)
ccc      padx = padim
ccc      pady = padim
ccc      npadx = 72       ! no. of pads along x
ccc      cathx2 = 290.    ! cathode half dimension in x
ccc      npady = 72       ! no. of pads along y
ccc      cathy2 = 290.    ! cathode half dimension in y   
ccc      pitch = 4        ! pitch between wires
ccc      ncath = 16       ! number of cathodes
c
c
ccc      nwires = nwiresin
ccc      nwirepad = nwirepadin
ccc      pitch =  pitchin
      nwires = 144
      nwirepad = 2
      pitch = 4                     ! pitch between wires
c
c
      pedm = pedmin
      sigp = sigpin
      sigm = sigmin
      sigs = sigsin
c@@----------------
ccc      pedm = 64.0   
ccc      sigp = 21.0
ccc      sigm = 1.15
ccc      sigs = 0.30
cmo---------------
c
      thresh =  threshin
      dosigcomp = dosigcompin  ! thresh*sigma compensation (1=ON)
c@@---------------------
ccc      thresh = 3.
c@@------------------
cmo      thresh = 3.            !    <113>  pad/ev enoise
cmo      thresh = 2.5           !    <490>  pad/ev enoise
cmo      thresh = 2.3           !    <890>  pad/ev enoise
cmo      thresh = 2.0           !   <    >  pad/ev enoise
cmo      thresh = 1.8           !   <2980>  pad/ev enoise
cmo------------------
c
c -------------------------------------
c   Read pad table for signal induction
c -------------------------------------
c
ccc      padname = padnamein
ccc      padname = 'padtabl1.8250'
c------------------------------
ccc      character*128 filename
ccc      character*256 expfilename
ccc      integer       status
ccc      integer       expnm
ccc      external      expnm
c
      filename = '$CORAL/src/pkopt/padtabl1.8250'
      expfilename = ''
      status =  expnm( filename, expfilename )
c
ccc      write( 6, '(A)' ) filename  
      write( 6, '(A)' ) expfilename  
c
      open( unit=20, file=expfilename, status='OLD', 
     .      form='FORMATTED', iostat=status, err=777 )
c
ccc      write( 6, '(I10)' ) status
c
      step = 0.25    ! each pad is divided on padx/step steps along X
      np = 5
c
ccc      open(20,file=padname,status='old',err=777)
c
      nstep = nint(padx/step)
      ntp = np**2
c
      do i=1,nwirepad
        do j=1,nstep
          read(20,*)
          do k = 1,ntp
            read(20,*) tableph(i,j,k)
cpp      print *,tableph(i,j,k)
          enddo
        enddo
      enddo
      close(20)
c
      RETURN
c
 777  print *,'file  ',expfilename,'   NOT read!'
      print *,'-------------------------------'
      print *,'ERROR IN RICH MC DIGITIZATION'
ccc      print *,'CORAL  stops'
ccc      stop 777                                //   020902
      do i=1,nwirepad
        do j=1,nstep
          do k = 1,ntp
            tableph(i,j,k) = 0.
          enddo
        enddo
      enddo
      RETURN
c
      END
c
c
c-----------------------------------------------------------------------
c1900      SUBROUTINE phtrack(xp,yp,itype,eph,icath,ndph,nmip,nfeed)
      SUBROUTINE phtrack(khitcc,xp,yp,itype,eph,icath,ndph,nmip,nfeed)
c---------------------------------------------------------------------
c
c --------------------------------------------
c  EVENT generator : MIP, CHERENKOV + FEEDBACK 
c --------------------------------------------
c
c IN:
c     xp - x coordinate of particle hit on cathode plane (mm)
c
c     yp - y coordinate of particle hit on cathode plane (mm)   
c
c     itype - 1 - Cherenkov 141
c           - 2 - MIP
c
c     eph   -  Cherenkov photon energy ( eV ) if itype=1 or
c
c     nsec=int(eph) -  number of secondaries, generated by charged particle
c                      in pad-anode gap if itype=2
c
c     icath - cathode number (1-16)
c
c OUT: 
c       ndph - updated number of detected Cherenkov photons
c
c       nmip - updated number of detected charged particles
c
c       nfeed - updated number of detected feedback photons
c
c *****************************************************************
c
c                  SOME USER GUIDE
c
c This programme is based on RD26 code FRESIM
c Let  np = ndph + nmip + nfeed be a general index of entries.
c Then after call to phtrack for each entry we will have:
c
c typho(np) = 0 ---> MIP
c typho(np) = 1 ---> Cherenkov photon
c typho(np) = 2 ---> feedback from MIP 
c typho(np) = 3 ---> feedback from Cherenkov
c 
c nelcoo(np,1) - hit pad number along x ( 1:NPADX really but 0:NPADX+1 here)
c nelcoo(np,2) - hit pad number along y ( 1:NPADX really but 0:NPADX+1 here)
c nelcoo(np,3) - hit cathode number
c
c pulse(np,25) - pulse heights for hit pad and also for 24 neighbour pads
c
c Basic parameters for pulse generation are :
c
c csiqe(25) - CsI quantum efficiency for 25 energy bins ( Eph=5.5-7.9 eV )
c phslo(16) - chamber gains
c alfa(16)  - feedback rates in ADC^1
c ridfact(16) - QE overall reduction factor for each cathode
c 
c Physical processes flags are: (see INITALL)
c
c dofeedback   feedback photons generation (1=ON, 0=OFF)
c dorefl       feedback photons scattering on CsI surface (1=ON, 0=OFF)
c doQEtest     QE test for Cherenkov photons (1=ON, 0=OFF)
c
c*********************************************************************
c
	IMPLICIT NONE

#include "finclude/pads.inc"
#include "finclude/baspar.inc"
#include "finclude/elpulse.inc"
#include "finclude/dstart.inc"

c----------------
	REAL xp,yp,eph
	INTEGER itype,icath,ndph,nmip,nfeed	
        INTEGER khitcc
c----------------

	INTEGER nsec,jen,index,indey,i,j,indew,n3,n4,k,l,i0,ntp
	INTEGER np,kk,nloop,istart,ifeed,i1,indexw,idummy
	INTEGER ken,nfp,inx,iny,i2,i3

ccc	REAL step,enrg,rj,qe,xqe,random,dmin,d,slope,phgen
	REAL enrg,rj,qe,xqe,random,dmin,d,slope,phgen
	REAL xt,yt,xt1,yt1,phfun,enfp,psi,pi,phi,dcz1,dcx1,dcy1
	REAL xfd,yfd,xxf,yyf,zzf,ro2,fodist,dfeed,trgas,the0,rif0
	REAL the1,ri1,ri2,ri3,ricsi,fp,rif1,rif2,cp,rdm,trg

        parameter (pi=3.141592654)

        ri1=2.95
        ri2=2.75
        ri3=2.5

	nsec = 1

ccc	call hf2(70,xp,yp,1.)

       if(itype.eq.1) then                ! Cherenkov photon

 	 if(doQEtest.eq.1) then           ! do QE test

c*********** find proper energy bin (jen) 

	 enrg=5.5
	 rj=0.
	 jen=1

 20      if(eph.ge.(enrg+rj/10.) .and.eph.lt.(enrg+(rj+1.)/10.) ) then
          goto 25
         else

          if(jen.lt.24) then
           jen=jen+1
           rj=rj+1.
           goto 20
          endif

c!!! 25      endif
         endif
 25      continue


         enrg=enrg+rj/10.

c********** CsI QE test
 
          qe=csiqe(jen)+(csiqe(jen+1)-csiqe(jen))*(eph-enrg)/0.1
          xqe = qe*ridfact(icath)
          call ranlux(random,1)
          if(random.ge.xqe)  go to 777           ! go to end
ccc          call hf2(7,eph,xqe,1.)
	 endif                                   ! end of do QE test
 
	else                                     ! end of Cherenkov specific

	 nsec=int(eph)                           ! for MIP

	endif

ccc        call hf1(600,float(nsec),1.)


c********** active area test

	if(xp.lt.-cathx2.or.xp.gt.cathx2) go to 777  ! go to end	  
	if(yp.lt.-cathy2.or.yp.gt.cathy2) go to 777  ! go to end

c********** type of signal

	if(itype.eq.1) then
	 ndph=ndph+1
	 typho(ndph+nfeed+nmip)=1
	else
	 nmip=nmip+1
	 typho(ndph+nfeed+nmip)=0
	endif

c********* number of hit pad along x and y

        index = 0

       do i = 0,npadx+1
         if (cpadx(i).le.xp.and.
     +  (cpadx(i)+padx).gt.xp) then
        index = i
        nelcoo(ndph+nfeed+nmip,1)=index
        goto 50
        endif
       end do
       
50     indey = 0

       do i = 0,npady+1
         if (cpady(i).le.yp.and.
     + (cpady(i)+pady).gt.yp) then
       indey = i
       nelcoo(ndph+nfeed+nmip,2)=indey
       goto 60
       endif
       end do
         
60     continue

ccc        call hf2(15,float(index),float(indey),1.)

c******* cathode number

	nelcoo(ndph+nfeed+nmip,3)=icath
c
	ihitcc(ndph+nfeed+nmip) = khitcc         !   1900


c -----------------------------------------------------
c  get coordinates for wire hit and position inside pad
c  wire index (n3)
c -----------------------------------------------------
c
c pdv   y -> x
      dmin = abs(cwire(0)-xp)
      indew = 0
      do i = 1,nwires+1
         d = abs(cwire(i)-xp)
         if (d.le.dmin) then
              dmin = d
              indew = i
         endif
      end do
	n3=indew
c     
c ----------------------------------------------------------------------
c  index (n4) inside pad (the pad is subdivided into 32 step of 0.25 mm)
c ----------------------------------------------------------------------
c
c pdv   x -> y
      n4 = int(abs((cpady(nelcoo(ndph+nfeed+nmip,2))-
     +yp))/step) + 1
c
       if(n4.eq.33) then
        n4=1
        nelcoo(ndph+nfeed+nmip,2)=nelcoo(ndph+nfeed+nmip,2) + 1
       endif
c
c
c    detect the wire index inside the pad
c
c pdv   y -> x
      k = int(nwires/npadx)
      indexw = mod(n3,k)
      if (indexw.eq.0) then
       indexw =1
      else
       indexw =2
      endif
c
c ----------------------------------------------------------------
c
c   total PH generation
c

        slope=phslo(icath)
        phgen=0


        do l=1,nsec                      ! for Cherenkov nsec=1
        i0 = 0

        do while (i0.eq.0)
          call ranlux(xt,1)
          call ranlux(yt,1)
c
          xt1 = xt*2000
          yt1 = phfun(xt1,slope)

         if(yt.lt.yt1) then
          phgen = phgen + xt1
          i0 = 1
ccc          call hf1(50,xt1,1.)
         endif
        enddo
      
	enddo

c
c read the ph table for this case
c
      ntp = np**2
c
      do   kk=1,25
      pulse(ndph+nfeed+nmip,kk) = phgen*tableph(indexw,n4,kk)
         enddo

	 
c-------   loops on feedback -------------------------------------

      if(dofeedback.eq.0) go to 777

       nloop=nfeed
       istart=nfeed
       ifeed=0
c pdv   x <-> y
       posx(nfeed+1)=cwire(indew)
       posy(nfeed+1)=yp
       phf(nfeed+1)=phgen

 90    do i1=istart,nloop
c ------- feedback photon --------------------------------------------

        fp=phf(i1+1)*alfa(icath)
        call poissn(fp,nfp,idummy)
c---------------------------------

        do j=1,nfp

         call ranlux(random,1)

        if(random.le.0.57) then
         enfp= 7.5
         ken=21
         ricsi=ri2
         goto 93
        endif

        if(random.gt.0.57 .and. random .le.0.7) then
         enfp= 6.4
         ken=10
         ricsi=ri1
         goto 93
        endif

        if(random.gt.0.7) then       
         enfp= 7.9
         ken=25
         ricsi=ri3
        endif

 93     continue

ccc        call hf1(141,enfp,1.)

        call ranlux(random,1)
        psi=random*2.*pi

        if(psi.lt.1.570796327.or.psi.gt.4.71238898) go to 555
        call ranlux(random,1)

        phi=random*2.*pi
        dcz1=cos(psi)
        dcx1=sin(psi)*cos(phi)
        dcy1=sin(psi)*sin(phi)
        
        xfd=posx(i1+1) -(2.*dcx1/dcz1)
        yfd=posy(i1+1) -(2.*dcy1/dcz1)
         
c********** active area test

	if(xfd.lt.-cathx2.or.xfd.gt.cathx2) go to 555  ! go to ..	  
	if(yfd.lt.-cathy2.or.yfd.gt.cathy2) go to 555  ! go to ..
c     
        xxf=(xfd)**2
        yyf=(yfd)**2
        zzf=4.
        ro2=zzf/(dcz1)**2
        fodist=sqrt(xxf+yyf)
        dfeed=sqrt(ro2)

	trg=trgas(enfp,dfeed/10.)	
        call ranlux(random,1)

ccc	call hf2(10,dfeed/10.,trg,1.)

ccc	call hf1(20,xfd-posx(i1+1),1.)

        if((100.*random).gt.trg) go to 555

ccc	call hf1(30,xfd-posx(i1+1),1.)


c       reflection at the CsI surface

        if(psi.lt.3.14159265) then
         the0= 3.14159265-psi
        else
         the0= psi-3.14159265
        endif

        if(the0.eq.0..or.the0.eq.3.14159265) then
         rif0=1.0
         goto 94
        endif
        
        the1=asin((sin(the0))/ricsi)
        rif1=(tan(the0-the1)/tan(the0+the1))**2
        rif2=(sin(the1-the0)/sin(the0+the1))**2
 
        call ranlux(cp,1)
        rif0=cp*rif1 + (1-cp)*rif2
        if(dorefl.eq.0) rif0=-1.0

 94     continue 
        call ranlux(rdm,1)
        if(rdm.lt.rif0) go to 555   

ccc	call hf1(40,xfd-posx(i1+1),1.)

    
  

c *****      CsI QE

        call ranlux(random,1)
ccc	call hf2(18,enfp,csiqe(ken),1.)	
        if(random.gt.csiqe(ken)*ridfact(icath)) go to 555
        nfeed=nfeed+1
        ifeed=ifeed+1
        posx(nfeed+1)=xfd

c********** type of signal

	if(itype.eq.1) then
	 typho(ndph+nfeed+nmip)=3
	else
	 typho(ndph+nfeed+nmip)=2
	endif


c********* find pad numbers

        inx=0
       do i2 = 0,npadx+1
        if (cpadx(i2).lt.xfd.and.(cpadx(i2)+padx).ge.xfd) then
         inx = i2
         nelcoo(ndph+nfeed+nmip,1)=inx
         goto 97
        endif
       enddo
       
  97   INY = 0
       do i3 = 0,npady+1
        if (cpady(i3).lt.yfd.and.(cpady(i3)+pady).ge.yfd) then
         iny = i3
         nelcoo(ndph+nfeed+nmip,2)=iny
         goto 98
        endif
       end do
         
  98   continue

c******* cathode number

ccc        call hf2(16,float(inx),float(iny),1.)

	nelcoo(ndph+nfeed+nmip,3)=icath
c
ccc	ihitcc(ndph+nfeed+nmip) = khitcc         !   1900


c ----------------------------------------------------------------
c     get coordinates for wire hit and position inside pad
c
c      wire index
c-----------------------------------------------------------------
c pdv  y -> x
      dmin = abs(cwire(0)-xfd)
      indew = 0
      do i = 1,nwires+1
         d = abs(cwire(i)-xfd)
         if (d.le.dmin) then
              dmin = d
              indew = i
         endif
      enddo
      n3 = indew
      posx(nfeed+1)=cwire(indew)


c ----------------------------------------------------------------------
c      index inside pad (the pad is subdivided into 32 steps of 0.25 mm)
c ----------------------------------------------------------------------
c pdv x -> y
      n4 = int(abs((cpady(nelcoo(ndph+nfeed+nmip,2))-
     + yfd))/step) + 1

        if(n4.eq.33) then
         n4=1
         nelcoo(ndph+nfeed+nmip,2)=nelcoo(ndph+nfeed+nmip,2) + 1
        endif
c
c   total PH generation
c

       i0 = 0
       slope = phslo(icath)

      do while (i0.eq.0)
         call ranlux(xt,1)
         call ranlux(yt,1)
c pdv 
         xt1=xt*2000
         yt1 = phfun(xt1,slope)
         if(yt.lt.yt1) then
           phgen = xt1
           phf(nfeed)=phgen
           i0 = 1
ccc         call hf1(51,xt1,1.)
         endif
      enddo
c
c    detect the wire index inside the pad
c
      k = int(nwires/npady)
      indexw = mod(n3,k)
      if (indexw.eq.0) then
       indexw = 1
      else
       indexw = 2
      endif

c
c    read the ph table for this case
c

      ntp = np**2
      do kk = 1,25
      pulse(ndph+nfeed+nmip,kk) = 
     +    phgen*tableph(indexw,n4,kk)
      end do

 
555    continue

       enddo
       enddo

       if(ifeed.ne.nfeed-nloop) print*,'FUCK, ifeed!!'

       if(ifeed.gt.0) then
        istart =nfeed-ifeed+1
        nloop=nfeed
        ifeed=0
        goto 90
       endif

c ----  end of feedback -------------------------- 



777    continue

	RETURN
	END
c
c-----------------------------------------------------------------------
      function phfun(x,y)
c------------------------
c
c------------------------
c PH exponential function
c------------------------
c
      phfun= exp(-x/y)
c
      END
c
c-----------------------------------------------------------------------
      FUNCTION trgas(x,y)
c------------------------
c
c---------------------------------------------------------
c transmission of the gas
c up to 100 percent, with ISO+METH(5-95)
c x = energy in eV
c y = distance in cm
c from curves in Lu-McDonald paper for RICH workshop proc.
c---------------------------------------------------------
c
c  This is the last for CH4
c
         if(x.lt.7.25) am=100000000.
         if(x.ge.7.25.and.x.lt.7.375) am=644.34
         if(x.ge.7.375.and.x.lt.7.5) am=211.4
         if(x.ge.7.5.and.x.lt.8.) am=195.
c
         trgas=100*exp(-y/am)
c
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE geomet
c----------------------
c
c Calculate some chamber geometry values
c (distances from chamber center) for anode wires and pads:
c
c cwire(1:nwires) really, but (0:nwires+1) here,
c this is made for taking into account probable photon convertion
c on grounding (around pad plane).
c 
c cpadx(0:npadx+1)
c
c cpady(0:npady+1)
c
	IMPLICIT NONE
c      
#include "finclude/pads.inc"
c
	REAL y0,x0
	INTEGER i,j
c
c
c
c     WIRES
c-----------
c
c pdv  y -> x
      nwires = npadx*padx/pitch 
c
      x0 =  (npadx / 2) * padx  - pitch / 2
      do i = 0,nwires+1
         cwire(i) = x0 - (i-1)*pitch           ! "centers"
      end do
c
c     PADS
c---------
c
      x0 = ( npadx / 2 ) * padx
      do i = 0, npadx+1
         cpadx(i) = - x0 + (i-1)*padx          ! "left" sides 
      end do
c
      y0 = ( npady / 2 )  * pady
      do j = 0, npady+1
         cpady(j) = - y0 + (j-1)*pady          ! "left" sides
      end do
c
c
c pdv y -> x
      nwirepad = nint(1.0*nwires/npadx)
c
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE pedsigma
c------------------------
c
c ----------------------------------
c     pedestal and sigma calculation
c ----------------------------------
c
c calculates thresholds for each pad for each cathode.
c basic values are: (see INITALL)
c-----------------------------------------------------
c
c (please look at formulas before asking a question)
c pedm - mean pedestal
c sigp - sigma of pedestals
c sigm - mean value of noise
c sigs - mean sigma of noise
c thresh - number of sigmas for thresholding

	IMPLICIT NONE

#include "finclude/pads.inc"
#include "finclude/enois.inc"


	REAL rgau
	INTEGER i,j,k
 
c
c     calculate pedestals
c
 	            do k=1,ncath	
       do i=1,npadx
         do j=1,npady
           ped(k,i,j)=0.0
           call norran(rgau)
ccc	   call hf1(60,rgau,1.)        
           ped(k,i,j) = rgau*sigp + pedm
c!!!            call hf1(220,float(ped(k,i,j)),1.)
ccc            call hf1(220,ped(k,i,j),1.)
         enddo
       enddo
	             enddo
c
c     calculate sigmas + threshold
c
		   do k=1,ncath
       do i=1,npadx
         do j=1,npady
            call norran(rgau)
            sig(k,i,j) = rgau*sigs+sigm  
            thr(k,i,j) = thresh*sig(k,i,j) + ped(k,i,j)
c!!!           call hf1(230,float(sig(k,i,j)),1.)
c!!!           call hf1(210,float(thr(k,i,j)),1.)
ccc           call hf1(230,sig(k,i,j),1.)
ccc           call hf1(210,thr(k,i,j),1.)
           enddo
       enddo
	           enddo
	END
c
c-------------------------------------------------------------------------
      SUBROUTINE elnoise
c-----------------------
c
c For each event calculates noise for each pad
c---------------------------------------------
c
c
	IMPLICIT NONE
c
c----------------------------
#include "finclude/pads.inc"
#include "finclude/enois.inc"
c----------------------------
c
	REAL rgau
	INTEGER i,j,k
c      
	             do k=1,ncath
      do i=1,npadx
         do j=1,npady
         enoise(k,i,j) = 0.000
         call norran(rgau)
         if(knoise.gt.0) rgau = 0.            !!!   18/8/00
cpp         print *,knoise,rgau
ccc         rgau = 0.
c@@-------------------------------
         enoise(k,i,j) = rgau*sig(k,i,j) + ped(k,i,j)
ccc         call hf1(200,enoise(k,i,j),1.)
         enddo
      enddo
		     enddo
c
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE threshold(nphevall,npatt)
c-----------------------------------------
c
c ----------------------------------------------------
c PH on pads and some outputs for NORMAL PAD threshold
c ----------------------------------------------------
c
c Does thresholding for each event
c --------------------------------
c
c  Output is stored in 
c  npatt - number of pads with pulse greater then treshold
c  For each such pad:
c  ict(1:npatt) - number of cathode
c  ixt(1:npatt) - number of pad along x
c  iyt(1:npatt) - number of pad along y
c  pht(1:npatt) - pulse height
c
c and also in padpht(ic,ix,iy) - pulse height for cathode ic,
c x pad number ix, y pad number iy. Zeros for pads without signal
c (below threshold)
c 
c
      IMPLICIT NONE
c
c------------------------------
#include "finclude/pads.inc"
#include "finclude/enois.inc"
#include "finclude/elpulse.inc"
#include "finclude/dstart.inc"
#include "finclude/clu.inc"
c------------------------------
c
      INTEGER n,npho,ix,iy,ii,jj,itab,ipulhei,ipulh,ith,ic
      INTEGER nphevall,npatt,np,ipp
c
      REAL pulhei,totpu
      REAL paph(16,72,72)
ccc      REAL lhitcc(16,72,72)           !   1900
c0104      INTEGER ihitw                   !   1900
c0104      INTEGER lhitcc(16,72,72)        !   091001
c0104      INTEGER lhitccw                 !   091001
      INTEGER*8 ihitw                   !   1900
      INTEGER*8 lhitcc(16,72,72)        !   091001
      INTEGER*8 lhitccw                 !   091001
c
      np = 5
c
      call vzero(padpht,16*72*72)
      call vzero(paph,16*72*72)
ccc      call vzero(lhitcc,16*72*72)         !   1900
      do   ic=1,16
        do   ix=1,72
          do   iy=1,72
            lhitcc(ic,ix,iy) = 0
          enddo
        enddo
      enddo
c
      n = (np-1)/2
c     
c calculates total pulse height for event (totpu)
c------------------------------------------------
c
      do   npho=1,nphevall
       totpu=0
       ic = nelcoo(npho,3)
c
         do   ii=1,np
           ix = nelcoo(npho,1)-n+ii-1
           if(ix.gt.0.and.ix.le.npadx) then
            do   jj=1,np
              iy = nelcoo(npho,2)-n+jj-1
              if(iy.gt.0.and.iy.le.npady) then
               itab = (ii-1)*np+jj
               paph(ic,ix,iy) = paph(ic,ix,iy) + pulse(npho,itab)
c
                if(npho.gt.1) then
                 pulhei = 0
                 pulhei = pulse(npho,itab) + enoise(ic,ix,iy)
                 ipulhei = pulhei
                 ipulh = nint(float(ipulhei)-ped(ic,ix,iy))
                 ith = int(thresh*sig(ic,ix,iy)+1.)
                 if(ipulh.gt.ith) then               
                  pulhei = float(ipulh)
                  totpu = totpu + pulse(npho,itab)
                 endif
                endif
c
              endif
            enddo
          endif
         enddo
c
       if(totpu.gt.0) then
c!!!        call hf1(500,float(totpu),1.)
ccc        call hf1(500,totpu,1.)
       endif
c      
      end do
c
c
c---- added 1900
      do   npho=1,nphevall
       ic = nelcoo(npho,3)
c
         do   ii=1,np
           ix = nelcoo(npho,1)-n+ii-1
           if(ix.gt.0.and.ix.le.npadx) then
            do   jj=1,np
              iy = nelcoo(npho,2)-n+jj-1
              if(iy.gt.0.and.iy.le.npady) then
ccc--- reference to only THREE hits is kept !!!
c0104- reference to only TWO or THREE hits is kept !!!
              lhitccw = lhitcc(ic,ix,iy)
              ihitw = ihitcc(npho)
ccc              if( ihitw .lt. 0 ) print *, '1 ', ihitcc(npho) 
              if(lhitccw.gt.0) then
                ihitw = ihitcc(npho) * MAXPAD         !   091001
ccc                if( ihitw .lt. 0 ) print *, '2 ', ihitcc(npho) 
                if( ihitw .lt. 0 ) ihitw = 0.         !   0104
              endif
              if(lhitccw.gt.MAXPAD) then
                ihitw = ihitcc(npho) * MAXPAD*MAXPAD  !   091001
ccc                if( ihitw .lt. 0 ) print *, '3 ', ihitcc(npho) 
                if( ihitw .lt. 0 ) ihitw = 0.         !   0104
              endif
              if(lhitccw.gt.MAXPAD*MAXPAD) then       !   091001
                ihitw = 0.
              endif
              lhitcc(ic,ix,iy) = lhitccw + ihitw
              endif
            enddo
          endif
         enddo
c      
      enddo
c---- 1900
c
c
c put threshold on pad PH + check photon overlap in fiducial
c-----------------------------------------------------------
c      
       npatt = 0
c
	do   ic=1,ncath	      
         do   ix =1,npadx
           do   iy = 1,npady
c
              paph(ic,ix,iy) = paph(ic,ix,iy) + enoise(ic,ix,iy)
c                             
ccc              call hf1(201,enoise(ic,ix,iy)-thr(ic,ix,iy),1.)
c
              if (paph(ic,ix,iy).gt.thr(ic,ix,iy)) then
                if( dosigcomp.gt.0) then
                  paph(ic,ix,iy) = paph(ic,ix,iy) - ped(ic,ix,iy)
                else
                  paph(ic,ix,iy) = paph(ic,ix,iy) - thr(ic,ix,iy)
                endif
                ipp = nint(paph(ic,ix,iy))
                paph(ic,ix,iy) = float(ipp)
                if(paph(ic,ix,iy).gt.0.) then
                  padpht(ic,ix,iy) = paph(ic,ix,iy)
ccc                  call hf2(17,float(ix),float(iy),1.)
                  npatt = npatt+1
	        ict(npatt) = ic
                ixt(npatt) = ix
                iyt(npatt) = iy
                pht(npatt) = padpht(ic,ix,iy)
                ihcc(npatt) = lhitcc(ic,ix,iy)           !   1900
                endif
              else
               paph(ic,ix,iy) = 0.0
              endif
c
           enddo
          enddo
	 enddo
c
c
      END
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c--------1---------2---------3---------4---------5---------6---------7--
