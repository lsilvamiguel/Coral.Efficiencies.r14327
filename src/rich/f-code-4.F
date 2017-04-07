
c-----------------------------------------------------------------------
      SUBROUTINE rhbook(kof)
c---------------------------
c
c----- date  4/5/99.
c----- version 1/a.
c
      IMPLICIT NONE
c
      INTEGER  kof, khh, k
c
      REAL  h
ccc      INTEGER  NPAWC
ccc      PARAMETER (NPAWC = 1000000)
ccc      COMMON / PAWC / h(NPAWC)
c
ccc      call hlimit(NPAWC)
c-----------------------
c
c
c---- books histograms.
c     -----------------
      khh = 0
c
c---- from CsRCEventClus::getProcClus :
c     ---------------------------------
      call hbook1(kof+ 501,'clu size',50,0.,50.,0.)
      call hbook2(kof+ 511,'clu (max pad)',9,-4.5,4.5,9,-4.5,4.5,0.)
      call hbook2(kof+ 521,'clu*PH (max pad)',9,-4.5,4.5,9,-4.5,4.5,0.)
c     -----------------------------------------------------------------
c
c---- from CsRCRichRec4.cc :
c     ----------------------
      call hbook1(kof+1521,'pads / ev.',100,0.,1000.,0.)
      call hbook1(kof+1525,'clus / ev.',100,0.,1000.,0.)
      call hbook1(kof+1528,'clu / pad',50,0.5,1.,0.)
      call hbook1(kof+1531,'particles / ev.',50,0.,50.,0.)
      call hbook1(kof+1535,'saved clusters / ev.',100,0.,4000.,0.)
c     ------------------------------------------------------------
c
c---- from CsRCEventRings::getEveRings :
c     ----------------------------------
      call hbook1(kof+1545,'ph.s peak / par.',100,0.,100.,0.)
      call hbook1(kof+1547,'ph.s peak / par.-fil',100,0.,100.,0.)
      call hbook1(kof+1548,'ph.s rec. / par.-fil/3sg/w',100,0.,100.,0.)
c     -----------------------------------------------------------------
c
c---- from CsRCPartRing::getRingBf :
c     ------------------------------
      call hbook2(kof+1551,'peak window',64,0.9984,1.,20,0.,20.,0.)
c     -------------------------------------------------------------
c
c---- from CsRCEventRings::getEveRings :
c     ----------------------------------
      call hbook1(kof+1560,'th-threc single ph.',100,-10.,10.,0.)
      call hbook1(kof+1561,'theta single ph.',280,0.,70.,0.)
c     ------------------------------------------------------
c
c---- from CsRCPartRing::histCntRatio :
c     ---------------------------------
      call hbook1(kof+1581,'PHpeakmx',100,0.,100.,0.)        !   phtw = 1.
      call hbook1(kof+1582,'PHpeakav',100,0.,100.,0.)        !   phtw = 1.
      call hbook1(kof+1583,'PHpeakav/mx',100,0.,1.,0.)
      call hbook1(kof+1585,'PHpeak tot',100,0.,100.,0.)
c     -------------------------------------------------

c---- from RCEventAnalysis::partIdent :
c     ---------------------------------
      call hbook1(kof+1568,'qsq vs nph',100,0.,10.,0.)
      call hbook2(kof+1569,'qsq vs nph',100,0.,100.,
     +            24,6.,30.,0.)
      call hbprox(kof+1569,0.)
      call hbslix(kof+1569,12,0.)
c
      call hbook1(kof+1651,'pi mass',100,0.0,1.0,0.)
      call hbook1(kof+1652,'K  mass',100,0.0,1.0,0.)
      call hbook1(kof+1653,'p  mass',100,0.0,1.0,0.)
      call hbook2(kof+1655,'mass vs mass',4,0.,4.,4,0.,4.,0.)
      call hbook2(kof+1656,'p vs dp/p pi',100,-0.25,0.25,
     +            30,0.,30.,0.)
      call hbslix(kof+1656,6,0.)
      call hbook2(kof+1657,'p vs dp/p  K',100,-0.25,0.25,
     +            30,8.,38.,0.)
      call hbslix(kof+1657,6,0.)
      call hbook2(kof+1658,'p vs dp/p  p',100,-0.25,0.25,
     +            30,15.,45.,0.)
      call hbslix(kof+1658,6,0.)
c
      call hbook1(kof+1661,'chi min',100,0.,10.,0.)
      call hbook1(kof+1662,'chi max',100,0.,10.,0.)
      call hbook2(kof+1663,'dchi',100,0.,2.5,
     +            100,0.,10.,0.)
      call hbprox(kof+1663,0.)
      call hbook1(kof+1666,'prob max',100,0.,1.,0.)
      call hbook1(kof+1667,'prob less',100,0.,1.,0.)
      call hbook1(kof+1668,'dprob',100,0.,0.1,0.)
c
      do   k=1,12
      khh = kof +1670 + k
      call hbook2(khh,'mass vs mass',4,0.,4.,4,0.,4.,0.)
         enddo
c     --------------------------------------------------
c
c---- from  CsRCEventAnalysis::MCMonitor :
c     ------------------------------------
      call hbook1(kof+3001,'photon',100,-10.,10.,0.)
      call hbook1(kof+3002,'ring',100,-2.5,2.5,0.)
cxx      call hbook1(kof+3006,'photon spectrum',60,0.,60.,0.)
      call hbook2(kof+3007,'betacal vs thphrec-thphcal',
     +            100,-0.00025,0.00025,80,0.9984,1.,0.)
      call hbslix(kof+3007,8,0.)
      call hbprox(kof+3007,0.)
      call hbook2(kof+3009,'mom vs thphrec-thphcal',100,-5.,5.,
     +            120,0.,60.,0.)
      call hbslix(kof+3009,12,0.)
      call hbook2(kof+3010,'betacal vs thphrec-thphcal',100,-5.,5.,
     +            80,0.9984,1.,0.)
      call hbslix(kof+3010,8,0.)
      call hbprox(kof+3010,0.)
      call hbook2(kof+3011,'ncaphi vs betarec',50,0.,50.,
     +            80,0.9984,1.,0.)
      call hbprox(kof+3011,0.)
      call hbproy(kof+3011,0.)
      call hbslix(kof+3011,8,0.)
      call hbook2(kof+3012,'n pho vs thphrec-thphcal',100,-5.,5.,
     +            40,0.,40.,0.)
      call hbslix(kof+3012,8,0.)
      call hbook2(kof+3013,'betacal vs thphrec-thphcal',100,-5.,5.,
     +            80,0.9984,1.,0.)
      call hbslix(kof+3013,8,0.)
      call hbook2(kof+3014,'betarec vs betacal',80,0.9984,1.,
     +            80,0.9984,1.,0.)
      call hbook1(kof+3015,'nphorec vs beta',64,0.9984,1.,0.)
c
      call hbook1(kof+3016,'nparus vs beta',64,0.9984,1.,0.)
      call hbook1(kof+3017,'nparec vs beta',64,0.9984,1.,0.)
c
      call hbook1(kof+3019,'ring-r-rec vs beta',64,0.9984,1.,0.)
      call hbook1(kof+3020,'ring-r-cal vs beta',64,0.9984,1.,0.)
c
      call hbook1(kof+3021,'average effic',64,0.9984,1.,0.)
c
      call hbook1(kof+3022,'nparus vs mom',30,0.,60.,0.)
      call hbook1(kof+3023,'nparec vs mom',30,0.,60.,0.)
cxx      call hbook1(kof+3024,'effic vs mom',30,0.,60.,0.)
      call hbook1(kof+3025,'nparus vs n phot',30,0.,60.,0.)
      call hbook1(kof+3026,'nparec vs n phot',30,0.,60.,0.)
cxx      call hbook1(kof+3027,'effic vs n phot',30,0.,60.,0.)
      call hbook2(kof+3028,'beta vs rec cut',60,0.,6.,
     +            64,0.9984,1.,0.)
      call hbprox(kof+3028,0.)
      call hbslix(kof+3028,8,0.)
      call hbook2(kof+3029,'p-type vs mom, rejected',60,0.,60.,
     +            15,1.,15.,0.)
      call hbprox(kof+3029,0.)
cxx      call hbook1(kof+3030,'mu mom',100,0.,100.,0.)
      call hbook1(kof+3031,'thphrec-thphcal',100,-2.5,2.5,0.)
c     ---------------------------------------------------------
c
c---- from  CsRCEventAnalysis::dataMonitor :
c     --------------------------------------
      call hbook1(kof+3036,'nparus vs beta',64,0.9984,1.,0.)
      call hbook1(kof+3037,'nparec vs beta',64,0.9984,1.,0.)
c
      call hbook1(kof+3041,'average effic',64,0.9984,1.,0.)
c
      call hbook1(kof+3042,'nparus vs mom',30,0.,60.,0.)
      call hbook1(kof+3043,'nparec vs mom',30,0.,60.,0.)
      call hbook1(kof+3045,'nparus vs n phot',30,0.,60.,0.)
      call hbook1(kof+3046,'nparec vs n phot',30,0.,60.,0.)
      call hbook1(kof+3050,'rec. cluster mult.',20,0.,20.,0.)
      call hbook2(kof+3053,'betacal vs thphrec-thphcaw',100,-5.,5.,
     +            80,0.9984,1.,0.)
      call hbslix(kof+3053,8,0.)
c     -------------------------------------------------------------
c
c---- from  CsRCEventAnalysis::MCMonitor ( cont ) :
c     ---------------------------------------------
      call hbook2(kof+3501,'phiph vs theph-cal',100,-10.,10.,
     +            100,0.,360.,0.)
c     -------------------------------------------------------
c
c---- from CsRCSelMirrors::CsRCSelMirrors :
c     -------------------------------------
      call hbook2(kof+3505,'ypamr vs xpamr',200,-2000.,2000.,
     +            200,-2000.,2000.,0.)
c     -------------------------------------------------------
c
c---- from CsRCPartPhots::getPartPhots :
c     ----------------------------------
      call hbook2(kof+3506,'ydey vs xdey',200,-1600.,1600.,
     +            200,-1600.,1600.,0.)
c     -----------------------------------------------------
c
c---- from  CsRCEventAnalysis::dataMonitor ( cont ) :
c     -----------------------------------------------
      call hbook2(kof+3519,'theta pa-de vs theta pa',80,0.,40.,
     +            90,0.,30.,0.)
      call hbprox(kof+3519,0.)
      call hbproy(kof+3519,0.)
      call hbook2(kof+3520,'theta pa-de vs phi pa',80,0.,40.,
     +            90,0.,360.,0.)
      call hbproy(kof+3520,0.)
c     ---------------------------------------------------------
c
c---- from CsRCEventRings::singlePhoton :
c     -----------------------------------
      call hbook2(kof+3521,'phi vs theph-rec',100,-10.,10.,
     +            60,0.,360.,0.)
      call hbslix(kof+3521,24,0.)
      call hbproy(kof+3521,0.)
      call hbook2(kof+3522,'phi vs theph-rec',100,-10.,10.,
     +            60,0.,360.,0.)
      call hbslix(kof+3522,24,0.)
      call hbproy(kof+3522,0.)
      call hbook2(kof+3523,'phi vs theph-rec',100,-10.,10.,
     +            60,0.,360.,0.)
      call hbslix(kof+3523,24,0.)
      call hbproy(kof+3523,0.)
      call hbook2(kof+3524,'phi vs theph-rec',100,-10.,10.,
     +            60,0.,360.,0.)
      call hbslix(kof+3524,24,0.)
      call hbproy(kof+3524,0.)
c
      call hbook1(kof+3525,'theph-rec',100,-10.,10.,0.)
c
      call hbook2(kof+3526,'phi vs theph-rec',100,-10.,10.,
     +            60,0.,360.,0.)
      call hbslix(kof+3526,24,0.)
      call hbproy(kof+3526,0.)
      call hbook2(kof+3527,'theta vs theph-rec',100,-10.,10.,
     +            60,0.,60.,0.)
      call hbslix(kof+3527,12,0.)
      call hbproy(kof+3527,0.)
      call hbook2(kof+3528,'thpamr vs theph-rec',100,-10.,10.,
     +            72,12.,16.,0.)
      call hbslix(kof+3528,4,0.)
      call hbproy(kof+3528,0.)
      call hbook1(kof+3529,'thpade',80,0.,40.,0.)
c     --------------------------------------------------------
c
c---- from  CsRCEventAnalysis::dataMonitor ( cont ) :
c     -----------------------------------------------
      call hbook2(kof+3601,'ypa vs xpa',100,-500.,500.,
     +            100,-500.,500.,0.)
      call hbprox(kof+3601,0.)
      call hbproy(kof+3601,0.)
      call hbandx(kof+3601,-60.,60.,0.)
      call hbook2(kof+3602,'ampa vs alpa',100,-0.1,0.1,100,-0.1,0.1,0.)
      call hbprox(kof+3602,0.)
      call hbproy(kof+3602,0.)
      call hbandx(kof+3602,-0.012,0.012,0.)
      call hbook2(kof+3603,'alpa vs xpa',100,-500.,500.,100,-0.1,0.1,0.)
      call hbprox(kof+3603,0.)
      call hbproy(kof+3603,0.)
      call hbook2(kof+3604,'ampa vs ypa',100,-500.,500.,100,-0.1,0.1,0.)
      call hbprox(kof+3604,0.)
      call hbproy(kof+3604,0.)
c     ------------------------------------------------------------------
c
c---- from CsRCEventParts::exitWindow :
c     ---------------------------------
      call hbook2(kof+3607,'dtan-b vs dtan-a',100,-0.025,0.025,
     +            100,-0.025,0.025,0.)
      call hbprox(kof+3607,0.)
      call hbproy(kof+3607,0.)
      call hbook2(kof+3608,'exit window',100,-10.,10.,100,-10.,10.,0.)
      call hbprox(kof+3608,0.)
      call hbproy(kof+3608,0.)
      call hbook2(kof+3609,'exit win vs mom',100,-10.,10.,100,0.,60.,0.)
      call hbslix(kof+3609,12,0.)
c     ------------------------------------------------------------------
c
c---- from  CsRCEventAnalysis::dataMonitor ( cont ) :
c     -----------------------------------------------
      call hbook1(kof+3610,'particle mom. acc.',120,0.,120.,0.)
      call hbook1(kof+3611,'particle mom. all ',120,0.,120.,0.)
c     ---------------------------------------------------------
c
c---- from  CsRCEventAnalysis::MCMonitor ( cont ) :
c     ---------------------------------------------
      call hbook1(kof+3615,'particle type all ',30,0.,30.,0.)
c     -------------------------------------------------------
c
c---- from CsRCEventParts::partAnalysis :
c     -----------------------------------
      call hbook1(kof+3617,'particle type rej.',30,0.,30.,0.)
c     -------------------------------------------------------

c---- from  CsRCEventAnalysis::MCMonitor ( cont ) :
c     ---------------------------------------------
      call hbook1(kof+3623,'e+- mom',120,0.,120.,0.)
      call hbook1(kof+3624,'mu+- mom',120,0.,120.,0.)
      call hbook1(kof+3625,'pi+- mom',120,0.,120.,0.)
      call hbook1(kof+3626,'K+- mom',120,0.,120.,0.)
      call hbook1(kof+3627,'ppbar mom',120,0.,120.,0.)
c     ------------------------------------------------
c
c---- from RCEventAnalysis::partIdent ( cont ) :
c     ------------------------------------------
      call hbook2(kof+3701,'c0-c0fit',100,-2.5,2.5,100,-2.5,2.5,0.)
      call hbprox(kof+3701,0.)
      call hbproy(kof+3701,0.)
      call hbook2(kof+3702,'thphrec-thphfit',100,-1.0,1.0,
     +            48,0.,48.,0.)
      call hbprox(kof+3702,0.)
      call hbslix(kof+3702,4,0.)
      call hbproy(kof+3702,0.)
      call hbook1(kof+3705,'chi/nu',100,0.,10.,0.)
      call hbook1(kof+3706,'nit',20,0.,20.,0.)
c     ----------------------------------------
c
c--- from CsRCEventDisplay::doEveDisplay :
c    -------------------------------------
chh            call hbook2(kof+4000...
chh            call hbook2(kof+4200...
chh            call hbook1(kof+4500...
chh            call hbook1(kof+5000...
c    -------------------------------------
c
c
      RETURN
c
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE cchbook1(kh,nx,x1,x2,ww)
c----------------------------------------
c
      call hbook1(kh,' ',nx,x1,x2,ww)
c------------------------------------
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE cchbook2(kh,nx,x1,x2,ny,y1,y2,ww)
c-------------------------------------------------
c
      call hbook2(kh,' ',nx,x1,x2,ny,y1,y2,ww)
c---------------------------------------------
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE histoput
c------------------------
c
c$$      call hrput(0,'/scra/paolo/richrec-cpp.hist','n')
      call hrput(0,'richrec-cpp.hist','n')
c-----------------------------------------------------
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE open_disp
c-------------------------
c
c$$      call hropen(2,'evedisp','/scra/paolo/evedisp-cpp.hist',
      call hropen(2,'evedisp','evedisp-cpp.hist',
c------------------------------------------------------------
     +              'n',1024,ista)
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE cchropen(kk,cc,lcc,file,lff,sta,lss,lrec,ista)
c--------------------------------------------------------------
c
      character*10 cc
      character*100 file
      character*10 sta
      character*10 ccu
      character*100 fileu
      character*10 stau
c
      lcc = lcc - 1
      ccu(1:lcc) = cc
      lff = lff - 1
      fileu(1:lff) = file
      lss = lss - 1
      stau(1:lss) = sta
         print *,kk,lrec
         print *,lcc,'  ',ccu
         print *,lff,'  ',fileu
         print *,lss,'  ',stau
      call hropen(kk,ccu,fileu,stau,lrec,ista)
c---------------------------------------------
      stop 123
c
      RETURN
      END
c
c-----------------------------------------------------------------------
      SUBROUTINE hcdir_disp
c--------------------------
c
      call hcdir('//evedisp',' ')
c--------------------------------
c
      RETURN
      END
