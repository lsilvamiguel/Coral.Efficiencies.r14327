
/*!
   \file    CsRCCluster.cc
   \------------------------
   \brief   CsRCCluster class author.
   \implementation  Paolo Schiavon
   \version 0.01
   \date    October 2000, rev. August 2005
*/

  #include <ostream>
  #include <cstdio>
  #include <iostream>
  #include <fstream>
  #include <cmath>
  #include <cstdlib>

  #include <iomanip>

  #include <unistd.h>
  #include <fcntl.h>

  #include "CsHbookProto.h"

  #include "CsErrLog.h"

//--------------------------
  #include "CsRCCluster.h"

  #include "CsRCExeKeys.h"
  #include "CsRCRecConst.h"
  #include "CsRCHistos.h"
  #include "CsRCDetectors.h"
  #include "CsRCCathode.h"
//--------------------------

  using namespace std;


//==========================================================================
  CsRCCluster::CsRCCluster() {
//----------------------------
    lPads_.clear();
  }

//==========================================================================
  CsRCCluster::CsRCCluster( int kClu, list<CsRCPad*> lPads, int iCat ) {
//----------------------------------------------------------------------
    kClu_ = kClu;
    lPads_ = lPads;
    ic_ = iCat;
    flag_ = true;

    CsRCCathode* cat = CsRCDetectors::Instance()->ptrToCat( ic_ );
    isPMT_ = cat->isPMT();
    isAPV_ = cat->isAPV();
  }

//==========================================================================
  CsRCCluster::CsRCCluster( int kClu, list<CsRCPad*> lPads, int iCat,
//-------------------------------------------------------------------
			    double xClu, double yClu, double PHClu ) {
    kClu_ = kClu;
    lPads_ = lPads;
    ic_ = iCat;
    xc_ = xClu;
    yc_ = yClu;
    PH_ = PHClu;
    flag_ = true;

    CsRCCathode* cat = CsRCDetectors::Instance()->ptrToCat( ic_ );
    isPMT_ = cat->isPMT();
    isAPV_ = cat->isAPV();
  }

//==========================================================================
  CsRCCluster::CsRCCluster( const CsRCCluster &clu ) {
//----------------------------------------------------
    //cout << "RICHONE : CsRCCluster CopyConstructor" << endl;
    kClu_ = clu.kClu_;
    lPads_ = clu.lPads_;
    ic_ = clu.ic_;
    xc_ = clu.xc_;
    yc_ = clu.yc_;
    PH_ = clu.PH_;
    flag_ = clu.flag_;

    isPMT_ = clu.isPMT_;
    isAPV_ = clu.isAPV_;
  }

//==========================================================================
  void CsRCCluster::print() const {
//---------------------------------
    cout << "   Cluster  " << kClu_ << " : " << "  x-pos " << xc_
         << ",  y-pos " << yc_ << ",  cath " << ic_ << ",  PH " << PH_
	 << ",  pads " << lPads_.size()
	 << ",  isMPT " << isPMT_ << ", isAPV " << isAPV_
         << ",  flag " << flag_ << endl;
  }

//==========================================================================
  CsRCCluster::~CsRCCluster() {
//-----------------------------
    lPads_.clear();
  }


//==========================================================================
  bool CsRCCluster::getBackWgtX( double &backWgt ) {
//--------------------------------------------------


//---- Paolo  -  May 2004

    CsRCExeKeys *key = CsRCExeKeys::Instance();

    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    static bool flag = true;
    static float xlimm[4], xlimx[4];
    static float ylimm[4], ylimx[4];
    static float ylim0[4];
    static int j1[4], j2[4];
    static int nsliy = 72;
    static float stepx, stepy;

    const int ncha = 74;
    const int nsli = 238;
    static float hcont[nsli][ncha];

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCluster::getBackWgt" );

//--- histo.s 'back-para.hist' are produced by :
//    /acid17/data06/schiavon/paw/rich0404/back-para.kumac
//    ----------------------------------------------------
      //char* fname = "back-para.hist";
      char* fname = strdup("   ");
      string sflo("   ");
      CsOpt* opt = CsOpt::Instance();
      bool boo = opt->CsOpt::getOpt( "RICHONE", "BackgrParam05", sflo );
      if( boo ) {
	fname = const_cast<char*>(sflo.c_str());
      } else {
	std::cout << " RICHONE, CsRCCluster::getBackWgt : ";
	std::cout << "NO background(05) file available " << std::endl;
	flag = false;
        return  false;
      }
      //int lun = 1;
      //int lrec = 1024;
      //int istat;
      //hropen( lun, " ", fname, " ", lrec, istat );
//@@---------------------------------------------
      //std::cout << istat << "  " << fname << std::endl;
      //if( istat != 0 ) {
      //if( fh < 0 ) {
      //std::cout << "ERROR opening file  " << fname << std::endl;
      //flag = false;
      //return  false;
      //}
      std::cout << " RICHONE, CsRCCluster::getBackWgt : background(05) file  "
      //        << fname << "  in use" << std::endl;
		<< fname << "  requested" << std::endl;

      //hrin_( 0 , 111, 0 );
//    -------------------

      j1[0] = 15;
      //j2[0] = 71;
      j2[0] = 73;
      j1[1] = 15;
      //j2[1] = 71;
      j2[1] = 73;
      j1[2] =  8;
      j2[2] = 67;
      j1[3] =  8;
      j2[3] = 67;

      //stepx =  8.;
      stepx =  16.;
      xlimm[0] = -1230.;
      xlimx[0] = -  50.;
      xlimm[1] =    50.;
      xlimx[1] =  1230.;
      xlimm[2] = -1230.;
      xlimx[2] = -  50.;
      xlimm[3] =    50.;
      xlimx[3] =  1230.;

      ylim0[0] =     0.;
      ylim0[1] =     0.;
      ylim0[2] = -1600.;
      ylim0[3] = -1600.;
      stepy = 20.;
      for( int k=0; k<4; k++ ) {
	ylimm[k] = ylim0[k] + (j1[k]-1)*stepy;
	ylimx[k] = ylim0[k] +  j2[k]   *stepy;
      }

      std::ifstream f( fname, ios::in );
      if( !f ) {
	//std::cout << "NOT f!" << std::endl;
	std::cout << " RICHONE, CsRCCluster::getBackWgt : ";
	std::cout << "Background(05) file NOT available " << std::endl;
	std::string mess = " RICHONE, CsRCCluster::getBackWgt : ";
	std::string err = "Background(05) file NOT available";
	mess.append( err );
	CsErrLog::Instance()->mes( elFatal, mess );
	flag = false;
        return  false;
      } else {
	for( int ks=0; ks<nsli; ks++ ) {
	  for( int kc=0; kc<ncha; kc++ ) f >> hcont[ks][kc];
  	  //std::cout << setprecision( 8 );
	  //for( int kc=0; kc<ncha; kc++ )
	  //  std::cout << kc << "  " << hcont[ks][kc] << " ";
	  //std::cout << std::endl << "  " << ks << std::endl;
          //std::cout << hcont[ks][kc] << " ";
	  //if( (kc+1)%10 == 0 ) std::cout << std::endl;
	}
        //std::cout << std::endl;
        //std::cout << "  " << ks+1 << std::endl;
      }
      //exit( 0 );
//--- added   100201   cppcheck
      //delete  fname;

    }
    backWgt = - 1.;
    if ( !flag ) return  false;

    float xcw = xc_;
    float ycw = yc_;
    int kqua = 0;
    if( xcw < 0.  &&  ycw > 0. ) kqua = 0;
    if( xcw > 0.  &&  ycw > 0. ) kqua = 1;
    if( xcw < 0.  &&  ycw < 0. ) kqua = 2;
    if( xcw > 0.  &&  ycw < 0. ) kqua = 3;

    int nhisto = 0;
    int kyv = 0;
    for( int ky=j1[kqua]; ky<=j2[kqua]; ky++ ) {
      float lim0 = ylim0[kqua] + (ky-1)*stepy;
      float lim1 = lim0 + stepy;
      if( ycw > lim0  &&  ycw <= lim1 ) {
        kyv = ky;
        break;
      }
    }
    if( kyv == 0 ) {
      if( ycw < (ylim0[kqua] + (j1[kqua]-1)*stepy) ) kyv = j1[kqua];
      if( ycw > (ylim0[kqua] +  j2[kqua]   *stepy) ) kyv = j2[kqua];
    }
    kyv -= j1[kqua];
    nhisto = kyv + kqua * 59;
    if( kqua == 3 ) nhisto++;

    float cco = 0.;
    if( nhisto >= 0  &&  nhisto < 238 ) {
      if( xcw < xlimm[kqua] )  xcw = xlimm[kqua];
      if( xcw > xlimx[kqua] )  xcw = xlimx[kqua];
      int kcha = int( (xcw - xlimm[kqua])/stepx );
      if( kcha >= 0  &&  kcha < ncha ) cco = hcont[nhisto][kcha];
      if( cco < 0. ) cco = 0.;
    }
    else {
      std::cout << nhisto << "  NO histo found!  " << ycw << std::endl;
      return  false;
    }
    backWgt = cco;
    //std::cout << setprecision( 8 );
    //std::cout << nhisto << "  " << xcw << "  "
    //      << ycw  << "  " << cco << "  " << backWgt << std::endl;
    //if( backWgt == 0. ) std::cout << nhisto << "  " << xcw << "  "
    //  << ycw  << "  " << cco << "  " << backWgt << std::endl;

    //if( cco > 0. ) {
      // moved to CsRCLikeXXX0Y::likeBackgr   -   040601
      //xh = xcw;
      //yh = ycw;
      //wh = cco;
      //if( hist.hRC1530 ) hist.hRC1530->Fill( xh, yh, wh );
//    ---------------------------------------------------
    //}

    return  true;
  }


//==========================================================================
  bool CsRCCluster::getBackWgt( double &backWgt ) {
//-------------------------------------------------


//- Paolo  -  January 2006

    CsRCExeKeys *key = CsRCExeKeys::Instance();
    CsRCRecConst *cons = CsRCRecConst::Instance();
    static int kOffHi = cons->kOffHi();
    CsRCHistos& hist = CsRCHistos::Ref();
    double xh, yh, wh;

    static bool flag = true;
    static float xlimm[4], xlimx[4];
    static float ylimm[4], ylimx[4];
    static float ylim0[4];
    static int j1[4], j2[4];
    static int nsliy = 72;
    static float stepx, stepy;

    static int format = 0;
    static const int nqua = 4;
    static const int nsli0 = nqua * 60 - 2;
    static int nsli = 0;
    static const int ncha = 74;
    static int nslito = 0;
//@@--------------------------
    static float hcont[nsli0][ncha];

    static std::vector<CsHist2D*> hRCBkgMap;

    static bool firstCall = true;
    if( firstCall ) {
      firstCall = false;
      key->acknoMethod( "CsRCluster::getBackWgt" );

//--- the numeric formatted file 'back-para-...-200x.vector74'
//    FORMAT 0 (old)
//    is produced by :
//    /acid17/data06/schiavon/paw/rich0404/back-para.kumac
//    and related
//    ------------------------------------------------------
//--- the numeric formatted file 'back-para-...-200x.new-vector74'
//    FORMAT 1 (new)
//    is produced by :
//    /acid17/data06/schiavon/back-para-05/paw/
//    back-para-new-vector74.kumac + sum-back24-new.F -> .go
//    ------------------------------------------------------

      char* fname = strdup("   ");
      string sflo("   ");
      int iflo = 0;
      bool boo = false;
      CsOpt* opt = CsOpt::Instance();
      boo = opt->CsOpt::getOpt( "RICHONE", "BackgrParam05Fmt", iflo );
      if( boo ) { format = iflo; }
      boo = opt->CsOpt::getOpt( "RICHONE", "BackgrParam05", sflo );
      if( boo ) {
	fname = const_cast<char*>(sflo.c_str());
      } else {
	std::cout << " RICHONE, CsRCCluster::getBackWgt : ";
	std::cout << "NO background(05) file available " << std::endl;
	flag = false;
        return  false;
      }
      std::cout << " RICHONE, CsRCCluster::getBackWgt : background(05) file  "
		<< fname << "  requested" << std::endl;

      if( format == 0 ) {
//    -------------------
	nsli = 238;
	nslito = nsli;
	j1[0] = 15;
	j2[0] = 73;
	j1[1] = 15;
	j2[1] = 73;
	j1[2] =  8;
	j2[2] = 67;
	j1[3] =  8;
	j2[3] = 67;

	stepx =  16.;
	xlimm[0] = -1230.;
	xlimx[0] = -  50.;
	xlimm[1] =    50.;
	xlimx[1] =  1230.;
	xlimm[2] = -1230.;
	xlimx[2] = -  50.;
	xlimm[3] =    50.;
	xlimx[3] =  1230.;
	stepy = 20.;
	ylim0[0] =     0.;
	ylim0[1] =     0.;
	ylim0[2] = -1600.;
	ylim0[3] = -1600.;
	for( int k=0; k<4; k++ ) {
	  ylimm[k] = ylim0[k] + (j1[k]-1)*stepy;
	  ylimx[k] = ylim0[k] +  j2[k]   *stepy;
	}
      }

      else if( format == 1 ) {
//    ------------------------
	nsli = 37;
	nslito = nqua * nsli;

	stepx =  16.;
	xlimm[0] = -1216.;
	xlimx[0] = -  32.;
	xlimm[1] =    32.;
	xlimx[1] =  1216.;
	xlimm[2] = -1216.;
	xlimx[2] = -  32.;
	xlimm[3] =    32.;
	xlimx[3] =  1216.;
	stepy = 32.;
	ylimm[0] =   272.;
	ylimx[0] =  1456.;
	ylimm[1] =   272.;
	ylimx[1] =  1456.;
	ylimm[2] = -1456.;
	ylimx[2] = - 272.;
	ylimm[3] = -1456.;
	ylimx[3] = - 272.;
      }

//--- 100617 --- //
      else if( format == 2 ) {
//    ------------------------
//    ...
      }

      else { return  false; }

      std::ifstream f( fname, ios::in );
      if( !f ) {
	std::cout << " RICHONE, CsRCCluster::getBackWgt : ";
	std::cout << "Background(05) file NOT available " << std::endl;
	std::string mess = " RICHONE, CsRCCluster::getBackWgt : ";
	std::string err = "Background(05) file NOT available";
	mess.append( err );
	CsErrLog::Instance()->mes( elFatal, mess );
	flag = false;
        return  false;
      } else {
	for( int ks=0; ks<nslito; ks++ ) {
	  for( int kc=0; kc<ncha; kc++ ) f >> hcont[ks][kc];

  	  //std::cout << setprecision( 8 );
	  //for( int kc=0; kc<ncha; kc++ )
	  //  std::cout << kc << "  " << hcont[ks][kc] << " ";
	  //std::cout << std::endl << "  " << ks << std::endl;
	}
      }

      bool printMap = false;
      //bool printMap = true;
      if( printMap ) {
	std::cout << nslito << "  " << ncha << std::endl;
	double hMax = 0.;
	for( int ks=0; ks<nslito; ks++ ) {
	  for( int kc=0; kc<ncha; kc++ ) 
	    if( hcont[ks][kc] > hMax ) hMax = hcont[ks][kc];
	}
	std::cout << setprecision( 8 );
	std::cout << hMax << std::endl;
	//getchar();
	int ent = 0;
	int kq = 0;
	int kl = 0;
	for( int ks=0; ks<nslito; ks++ ) {
	  if( ks%nsli == 0 ) {
	    kq++;
	    std::cout << "quad  " << kq << std::endl;
	    kl = 1;
	    //getchar();
	  }
	  std::cout << "slice  " << kl << std::endl;
	  kl++;
	  std::cout << setprecision( 8 );
	  for( int kc=0; kc<ncha; kc++ ) {
	    std::cout << kc << "  " << hcont[ks][kc]/hMax << " ";
	    if( (kc+1)%7 == 0 ) std::cout << std::endl;
	    ent++;
	  }
	  std::cout << std::endl;
	}
	std::cout << "entries " << ent << std::endl;
	std::cout << std::endl;
	exit(0);
      }

      CsHistograms::SetCurrentPath("/RICH");

      std::string hTitle = "back map";
      int kName = 0;
      int nsli1 = 0;
      int nslih = nslito / 4;
      for( int kH=0; kH<4; kH++ ) {
	kName = kOffHi + 3300 + kH;
	std::stringstream hNBkgMap;
	hNBkgMap << kName;
	hRCBkgMap.push_back( new CsHist2D( hNBkgMap.str(), hTitle,
		    	     ncha, xlimm[kH], xlimx[kH],
			     nsli, ylimm[kH], ylimx[kH] ) );

	int nsli2 = nsli1 + nslih - 1;
	if( format == 0  &&  kH > 1 ) nsli2++;
	float dxx = (xlimx[kH]-xlimm[kH])/float(ncha);
	float dyy = (ylimx[kH]-ylimm[kH])/float(nslih);
	//std::cout << nsli1 << "  " << nsli2 << std::endl;
	int kt = 0;
	for( int ks=nsli1; ks<nsli2; ks++ ) {
	  yh = ylimm[kH] + kt * dyy + dyy/2.;
	  for( int kc=0; kc<ncha; kc++ ) {
	    xh = xlimm[kH] + kc * dxx + dxx/2.;
	    //std::cout << setprecision( 8 );
	    //std::cout << kc << "  " << hcont[ks][kc] << "  ";
	    wh = hcont[ks][kc];
	    if( hRCBkgMap[kH] ) hRCBkgMap[kH]->Fill( xh, yh, wh );
//hh                            ---------------------------------
	  }
	  //std::cout << std::endl << "  " << ks << std::endl;
	  kt++;
	}
	nsli1 = nsli2 + 1;
      }

      CsHistograms::SetCurrentPath("/");
      //exit( 0 );
//--- added   100201   cppcheck
      //delete  fname;

    }
    backWgt = - 1.;
    if ( !flag ) return  false;

    float xcw = xc_;
    float ycw = yc_;
    int kqua = 0;
    if( xcw < 0.  &&  ycw > 0. ) kqua = 0;
    if( xcw > 0.  &&  ycw > 0. ) kqua = 1;
    if( xcw < 0.  &&  ycw < 0. ) kqua = 2;
    if( xcw > 0.  &&  ycw < 0. ) kqua = 3;

    int nhisto = 0;
    int kyv = 0;

    if( format == 0 ) {
//  -----------------
      for( int ky=j1[kqua]; ky<=j2[kqua]; ky++ ) {
	float lim0 = ylim0[kqua] + (ky-1)*stepy;
	float lim1 = lim0 + stepy;
	if( ycw > lim0  &&  ycw <= lim1 ) {
	  kyv = ky;
	  break;
	}
      }
      if( kyv == 0 ) {
	if( ycw < (ylim0[kqua] + (j1[kqua]-1)*stepy) ) kyv = j1[kqua];
	if( ycw > (ylim0[kqua] +  j2[kqua]   *stepy) ) kyv = j2[kqua];
        //std::cout << kyv << "  NO histo found!  " << ycw << std::endl;
      }
      kyv -= j1[kqua];
      nhisto = kyv + kqua * 59;
      if( kqua == 3 ) nhisto++;
    }

    else if( format == 1 ) {
//  ----------------------
      for( int ky=1; ky<=nsli; ky++ ) {
	float lim0 = ylimm[kqua] + (ky-1)*stepy;
	float lim1 = lim0 + stepy;
        if( ycw > lim0  &&  ycw <= lim1 ) {
          kyv = ky;
          break;
        }
      }
      if( kyv == 0 ) {
        if( ycw < ylimm[kqua] ) kyv = 1;
        if( ycw > ylimx[kqua] ) kyv = nsli;
        //std::cout << kyv << "  NO histo found!  " << ycw << std::endl;
      }
      kyv--;
      nhisto = kyv + kqua * nsli;
    }

    else { return  false; }

    float cco = 0.;
    if( nhisto >= 0  &&  nhisto < nslito ) {
      if( xcw < xlimm[kqua] )  xcw = xlimm[kqua];
      if( xcw > xlimx[kqua] )  xcw = xlimx[kqua];
      int kcha = int( (xcw - xlimm[kqua])/stepx );
      if( kcha >= 0  &&  kcha < ncha ) cco = hcont[nhisto][kcha];
      if( cco < 0. ) cco = 0.;
    }
    else {
      std::cout << " RICHONE, CsRCCluster::getBackWgt() : ";
      std::cout << "NO histo found!  " << nhisto << "  " << ycw << std::endl;
      return  false;
    }
    backWgt = cco;

    //std::cout << setprecision( 8 );
    //std::cout << nhisto << "  " << xcw << "  "
    //      << ycw  << "  " << cco << "  " << backWgt << std::endl;
    //if( backWgt == 0. ) std::cout << nhisto << "  " << xcw << "  "
    //  << ycw  << "  " << cco << "  " << backWgt << std::endl;

    return  true;
  }

