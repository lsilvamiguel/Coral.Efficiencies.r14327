// $Id: CsRTRelation.cc,v 1.14 2010/01/28 12:51:25 tnagel Exp $

/*!
   \file    CsRTRelation.cc
   \brief   RT relation for drift-like detectors
   \author  Hugo Pereira
   \version $Revision: 1.14 $
   \date    $Date: 2010/01/28 12:51:25 $
*/

#include "CsRTRelation.h"
#include "CsOpt.h"

using namespace std;

//___BASE_______________________________________________________________________
CsRTRelation::CsRTRelation( CsDetector& det ) :
  dt_(-1),
  d_( &det ),
  hasGrid_(false),
  hasRegularGrid_( false ),
  hasRTPars_( true ),
  _FIT1_( true ),
  _FIT2_( false )
{
  rtGrid_.clear();
  rtGReg_.clear();
  rtPars_.push_back(0);
  rtPars_.push_back(det.getVel()/det.getTiSli());
}

//_______________________________________________________________________________
CsRTRelation& CsRTRelation::operator=( const CsRTRelation& rt )
{
  if( this != &rt ) *this = rt;
	return( *this );
}

//_______________________________________________________________________________
bool CsRTRelation::readFromOpt( )
{
  //=== clean old registered values ===
  rtGrid_.clear();  hasGrid_ = false;
  rtGReg_.clear();  hasRegularGrid_ = false;
  rtPars_.clear();  hasRTPars_ = false;

  list<string> member;          // rtGrid string member. must have at least three items
  list<string>::iterator Im;    

  CsOpt* opt = CsOpt::Instance();

  //=== first try reading RT Relation Parametrisation ===
  //=== FIT1 parametrisation, that is: polynom fit
  if( opt->getOpt( "RT", d_->GetTBName() , member ) ||       // Kept for backward compatibility
      opt->getOpt( "RTFit1", d_->GetTBName() , member ) ) {  // correct key is RTFit1 ( bool _FIT2_ = true )
    if( member.size()<2 ) CsErrLog::mes(elError,"Wrong format for RT relation member.");
    else {
      double rtPar;
      unsigned int i = 0;
      
      _FIT1_ = true;    // Set fit function type
      rtPars_.clear();
      for( Im = member.begin(); Im != member.end(); Im++, i++ ) {
        istringstream( *Im ) >> rtPar;
        rtPars_.push_back( rtPar );
      }
      hasRTPars_ = true;
      cout << "CsRTRelation::readFromOpt - got _FIT1_ for " << d_->GetTBName() <<": ";
      for( unsigned int ip = 0; ip < rtPars_.size(); ip++ ) cout <<" "<< rtPars_[ip] ;
      cout << endl;
      return true;
    }
  } 

  //=== FIT2 parametrisation, that is: Hyperbolic tangeant fit
  if( opt->getOpt( "RTFit2", d_->GetTBName() , member ) ) {   // correct key is RTFit2 ( bool _FIT2_ = true )
    if( member.size()<2 ) CsErrLog::mes(elError,"wrong format for RT relation member.");
    else {
      double rtPar;
      unsigned int i = 0;
      
      _FIT2_ = true;    // Set fit function type
      rtPars_.clear();
      for( Im = member.begin(); Im != member.end(); Im++, i++ ) {
        istringstream( *Im ) >> rtPar;
        rtPars_.push_back( rtPar );
      }
      hasRTPars_ = true;
      cout << "CsRTRelation::readFromOpt - got _FIT2_ for " << d_->GetTBName() <<": ";
      for( unsigned int ip = 0; ip < rtPars_.size(); ip++ ) cout <<" "<< rtPars_[ip] ;
      cout << endl;
      return true;
    }
  } 
  
  //=== second try to read RT Grid ===

  bool grid_in_opt = false;

  while( opt->getOptRec( "RTGrid", d_->GetTBName(), member ) ) {
    double t;
    double r;
    double res;
    unsigned int i=0;
    grid_in_opt = true;
    for( Im = member.begin(); Im != member.end(); Im++, i++ )
    switch ( i ) {
    case 0: istringstream( *Im ) >> t;   break;
    case 1: istringstream( *Im ) >> r;   break;
    case 2: istringstream( *Im ) >> res; break;
    default: break;
    }

    RTGridPoint gp;
    gp.t = t;
    gp.r = r;
    gp.res = res;
    rtGrid_.push_back( gp );
  }

  if (!grid_in_opt) return false;
  
  if( rtGrid_.size() < 2 ) {
    CsErrLog::msg(elError,__FILE__,__LINE__,
		  "%s: Not enough grid point for RT relation.",
		  d_->GetTBName().c_str() );
    return false;
  } else hasGrid_ = true;

  _sortGrid( rtGrid_ );
  if( opt->getOpt( "RT", "useRegularGrids" ) ) _makeRegular();

  cout << "CsRTRelation::readFromOpt. " << rtGrid_.size()
    << ( (hasRegularGrid_) ? " regular":" " )
    << " points read for " << d_->GetTBName().c_str() << endl;

  return true;

}

//_______________________________________________________________________________
double CsRTRelation::getTfromR( const double r, bool &error )
{
  if( !( hasGrid_ || hasRegularGrid_ ) ) {
    CsErrLog::mes( elError, "No Grid set. Abort.");
    error = true;
    return -1;
  }

  vector< RTGridPoint > rt = ( (hasRegularGrid_) ? rtGReg_ : rtGrid_ );

  unsigned int ig = 0;
  double r0 = -1, r1 = -1;
  double t0 = -1, t1 = -1;
  while( ig < rt.size() && rt[ig].r < r ){
    r0 = rt[ig].r;
    t0 = rt[ig].t;
    ig++;
  }

  if( ig == 0 ) return rt.front().t;
  else if( ig == rt.size() ) return rt.back().t;
  //...else
  r1 = rt[ig].r;
  t1 = rt[ig].t;
  return (r0 == r1 ) ? t1:( t0 + (t1-t0)*(r-r0)/(r1-r0) );

}

//_______________________________________________________________________________
double CsRTRelation::getRfromT( const double t, bool &error )
{
  error = false;
  double r = 0;
  
  //===========================
  //=== RESULT FROM THE FIT ===
  if( hasRTPars_ ) {


    if( _FIT1_ ) {

      //=== polynomial version ===
      double tC = t-rtPars_[0];  // offset corrected time
      double v = 0;
      r = 0;
      if (tC>0) {
	for( unsigned int i = 1; i < rtPars_.size(); i++ ) {
	  v += rtPars_[i]*i*pow( tC, double(i-1));
	  r += rtPars_[i]*pow( tC, double(i));
	}
      }
      if( v < 0 ) error = true;

    } else if( _FIT2_ ) {

      //=== hyperbolic tg version ===
      double tC = t - rtPars_[0];  // offset corrected time
      double pol = 0;
      for( unsigned int i = 1; i < rtPars_.size(); i++ ) {
        pol += rtPars_[i]*pow(tC,double(i*2-1));
      }
      double cell = d_->getWirP()*0.5;
      r = cell * tanh( pol/cell );
      
    } 

  //================================
  //=== RESULT FROM REGULAR GRID ===
  } else if( hasRegularGrid_ ) {
    double id = ( t-rtGReg_.front().t )/dt_;
    register unsigned int i = (unsigned int) id;
    if( id < 0 ) { r = 0; error = true; }
    else if( id >= (double) rtGrid_.size()-1 ) { r = rtGReg_.back().r;  error = true; }
  	else r = rtGReg_[i].r + ( rtGReg_[i+1].r-rtGReg_[i].r )*( t - rtGReg_[i].t )/dt_ ;

  //============================
  //=== RESULT FROM RAW GRID ===
  } else if( hasGrid_ ) {
    unsigned int ig = 0;
	  double r0 = -1, r1 = -1;
	  double t0 = -1, t1 = -1;
  	while( ig < rtGrid_.size() && rtGrid_[ig].t < t ){
	    r0 = rtGrid_[ig].r;
	  	t0 = rtGrid_[ig].t;
	  	ig++;
	  }
	
	  if( ig == 0 ){ r = rtGrid_.front().r; error = true; }
	  else if( ig == rtGrid_.size() ){ r = rtGrid_.back().r; error=true; }
	  else {
	  	r1 = rtGrid_[ig].r;
	  	t1 = rtGrid_[ig].t;
	  	r = (t0 == t1 ) ? r1:( r0 + (t-t0)*(r1-r0)/(t1-t0) );
	  }

  //=======================================
  //=== NOTHING KNOWN ABOUT RT RELATION ===
  } else {
    CsErrLog::mes(elError,"CsRTRelation::getRfromT. Nothing known about rt relation.");
    r = 0;
    error = true;
  }

  //=== some checks ===
  if( r < 0 ) {
    r =0;
    error = true;
  } else if( r > 0.5 * d_->getWirP() ) {
    r = 0.5 * d_->getWirP();
    error = true;
  }

  return r;
}

//_______________________________________________________________________________
double CsRTRelation::getRfromT( const double t, const int wire, bool &error )
{
  error = false;
  double r = 0;
  
  //===========================
  //=== RESULT FROM THE FIT ===
  if( hasRTPars_ ) {


    if( _FIT1_ ) {

      //=== polynomial version ===
      double tC = t-rtPars_[0];  // offset corrected time
      double v = 0;
      r = 0;
      if (tC>0) {
	for( unsigned int i = 1; i < rtPars_.size(); i++ ) {
	  v += rtPars_[i]*i*pow( tC, double(i-1));
	  r += rtPars_[i]*pow( tC, double(i));
	}
      }
      if( v < 0 ) error = true;

    } else if( _FIT2_ ) {

      //=== hyperbolic tg version ===
      double tC = t - rtPars_[0];  // offset corrected time
      double pol = 0;
      for( unsigned int i = 1; i < rtPars_.size(); i++ ) {
        pol += rtPars_[i]*pow(tC,double(i*2-1));
      }
      double p1 = d_->Pitch(double(wire-1));
      double p2 = d_->Pitch(double(wire));
      double pmax = ( p1 >= p2 ? p1 : p2 );
      double cell = pmax*0.5;
      cout << "cell = " << cell << endl;
      r = cell * tanh( pol/cell );
      
    } 

  //================================
  //=== RESULT FROM REGULAR GRID ===
  } else if( hasRegularGrid_ ) {
    double id = ( t-rtGReg_.front().t )/dt_;
    register unsigned int i = (unsigned int) id;
    if( id < 0 ) { r = 0; error = true; }
    else if( id >= (double) rtGrid_.size()-1 ) { r = rtGReg_.back().r;  error = true; }
  	else r = rtGReg_[i].r + ( rtGReg_[i+1].r-rtGReg_[i].r )*( t - rtGReg_[i].t )/dt_ ;

  //============================
  //=== RESULT FROM RAW GRID ===
  } else if( hasGrid_ ) {
    unsigned int ig = 0;
	  double r0 = -1, r1 = -1;
	  double t0 = -1, t1 = -1;
  	while( ig < rtGrid_.size() && rtGrid_[ig].t < t ){
	    r0 = rtGrid_[ig].r;
	  	t0 = rtGrid_[ig].t;
	  	ig++;
	  }
	
	  if( ig == 0 ){ r = rtGrid_.front().r; error = true; }
	  else if( ig == rtGrid_.size() ){ r = rtGrid_.back().r; error=true; }
	  else {
	  	r1 = rtGrid_[ig].r;
	  	t1 = rtGrid_[ig].t;
	  	r = (t0 == t1 ) ? r1:( r0 + (t-t0)*(r1-r0)/(t1-t0) );
	  }

  //=======================================
  //=== NOTHING KNOWN ABOUT RT RELATION ===
  } else {
    CsErrLog::mes(elError,"CsRTRelation::getRfromT. Nothing known about rt relation.");
    r = 0;
    error = true;
  }

  //=== some checks ===
  double p1 = d_->Pitch(double(wire-1));
  double p2 = d_->Pitch(double(wire));
  double pmax = ( p1 >= p2 ? p1 : p2 );
  if( r < 0 ) {
    r =0;
    error = true;
  } else if( r > 0.5 * pmax ) {
    r = 0.5 * pmax;
    error = true;
  }

  return r;
}

//_______________________________________________________________________________
double CsRTRelation::getDRDT( const double t, bool &error )
{
  error = false;
  double v = 0;
  
  //===========================
  //=== RESULT FROM THE FIT ===
  if( hasRTPars_ ) {


    if( _FIT1_ ) {

      //=== polynomial version ===
      double tC = t-rtPars_[0];  // offset corrected time
      v = 0;
      for( unsigned int i = 1; i < rtPars_.size(); i++ )
      v += rtPars_[i]*i*pow( tC, double(i-1));

    } else if( _FIT2_ ) {

      //=== hyperbolic tg version ===
      double tC = t - rtPars_[0];  // offset corrected time
      double pol = 0;
      double dpol = 0;
      for( unsigned int i = 1; i < rtPars_.size(); i++ ) {
        pol  += rtPars_[i]*pow(tC,double(i*2-1));
        dpol += rtPars_[i]*(2*i-1)*pow(tC, double(2*(i-1)) );
      }
      double cell = d_->getWirP()*0.5;
      v = dpol/pow( cosh( pol/cell ) , 2.0 );
      
    } 

  //================================
  //=== RESULT FROM REGULAR GRID ===
  } else if( hasRegularGrid_ ) {
    double id = ( t-rtGReg_.front().t )/dt_;
    register unsigned int i = (unsigned int) id;
    if( id < 0 ) { v = 0; error = true; }
    else if( id >= (double) rtGrid_.size()-1 ) { v = 0 ;  error = true; }
  	else v = ( rtGReg_[i+1].r-rtGReg_[i].r )/dt_ ;

  //============================
  //=== RESULT FROM RAW GRID ===
  } else if( hasGrid_ ) {
    unsigned int ig = 0;
	  double r0 = -1, r1 = -1;
	  double t0 = -1, t1 = -1;
  	while( ig < rtGrid_.size() && rtGrid_[ig].t < t ){
	    r0 = rtGrid_[ig].r;
	  	t0 = rtGrid_[ig].t;
	  	ig++;
	  }
	
	  if( ig == 0 ){ v = 0; error = true; }
	  else if( ig == rtGrid_.size() ){ v = 0; error=true; }
	  else {
	  	r1 = rtGrid_[ig].r;
	  	t1 = rtGrid_[ig].t;
	  	v = (t0 == t1 ) ? 0:( (r1-r0)/(t1-t0) );
	  }

  //=======================================
  //=== NOTHING KNOWN ABOUT RT RELATION ===
  } else {
    CsErrLog::mes(elError,"CsRTRelation::getdRdT. Nothing known about rt relation.");
    v = 0;
    error = true;
  }
  
  return v;
}

//_______________________________________________________________________________
double CsRTRelation::getResAtT( const double t, bool& error )
{

  error = false;
  // check grid or regular grid
  if( !( hasGrid_ || hasRegularGrid_ ) ) {
    CsErrLog::mes( elError, "CsRTRelation::getResAtT. No Grid set. Aborted.");
    error = true;
    return -1;
  }

  // result from regular grid
  if( hasRegularGrid_ ) {
    double id = ( t-rtGReg_.front().t )/dt_;
    register unsigned int i = (unsigned int) id;
    if( id < 0 ) return 0;
    if( id >= (double) rtGReg_.size()-1 ) return rtGReg_.back().res;
  	return rtGReg_[i].res + ( rtGReg_[i+1].res-rtGReg_[i].res )*( t - rtGReg_[i].t )/dt_ ;
  }

  // result from raw grid
  unsigned int ig = 0;
	double res0 = -1, res1 = -1;
	double t0 = -1, t1 = -1;
	while( ig < rtGrid_.size() && rtGrid_[ig].t < t ){
	  res0 = rtGrid_[ig].res;
		t0 = rtGrid_[ig].t;
		ig++;
	}
	
	if( ig == 0 ) return rtGrid_.front().res;
	else if( ig == rtGrid_.size() ) return rtGrid_.back().res;
	else {
		res1 = rtGrid_[ig].res;
		t1 = rtGrid_[ig].t;
		return (t0 == t1 ) ? res1:( res0 + (res1-res0)*(t-t0)/(t1-t0) );
	}
}

//_______________________________________________________________________________
void CsRTRelation::dump( FILE* out )
{

	int ig = 0;
	
	if( out == NULL ) return;

  if( hasRTPars_ ) {

    if( _FIT1_ ) fprintf( out, "RTFit1 %8s ", d_->GetTBName().c_str() );
    if( _FIT2_ ) fprintf( out, "RTFit2 %8s ", d_->GetTBName().c_str() );
    for( unsigned int ip = 0; ip < rtPars_.size(); ip++ ) fprintf( out, "%10.5e ", rtPars_[ip] );
    fprintf( out, "\n" );

  } else if( hasGrid_ || hasRegularGrid_ ) {

    vector< RTGridPoint > rt =( ( hasRegularGrid_ ) ? rtGReg_ : rtGrid_ );
	  for( unsigned int ig = 0; ig < rt.size(); ig++ )
    fprintf( out, "RTGrid %8s %10.5f %10.5f %10.5f \n",
		  d_->GetTBName().c_str(),
      rt[ig].t,
      rt[ig].r,
      rt[ig].res );
    fprintf(out, "       %8s %10s %10s %10s \n","DetName","t(ns)","r(mm)","res(mm)");

	} else

  CsErrLog::mes(elError,"CsRTRelation::dump. Nothing known about rt relation.");

  return;
}

//___PRIVATE_METHODS____________________________________________________________
void CsRTRelation::_makeRegular()
{
  if( !( hasGrid_ || hasRTPars_ ) ){
    CsErrLog::mes( elError,"No grid nor parameters set. Abort");
    return;
  }

  unsigned int n = rtGrid_.size();
  if( n < 2 ){
    CsErrLog::mes( elError,"At least 2 Grid points required.");
    return;
  }

  rtGReg_.clear();
  dt_ = (rtGrid_.back().t - rtGrid_.front().t)/((double)n-1);
  for( unsigned int i=0; i< n; i++ ) {
    bool error;
    RTGridPoint pt;
    pt.t   = rtGrid_.front().t + i*dt_;
    pt.r   = getRfromT( pt.t, error );
    pt.res = getResAtT( pt.t, error );

    rtGReg_.push_back( pt );
  }
  hasRegularGrid_ = true;

  return;
}

//_______________________________________________________________________________
struct CsRTRelation::sortRTGridPoints_:
  public binary_function<RTGridPoint, RTGridPoint, bool> {
	bool operator() ( RTGridPoint p1, RTGridPoint p2 ) {
		if( p1.t < p2.t ) return true;
		else return false;
	}
};

//_______________________________________________________________________________
void CsRTRelation::_sortGrid( vector< RTGridPoint > &rt )
{

  list< RTGridPoint > l; l.clear();

  // put rtGrid_ into a list
  for( unsigned int i = 0; i < rt.size(); i++ ) l.push_back( rt[i] );
  l.sort( sortRTGridPoints_() );

  // put the list back into rt
  rt.clear();
  list< RTGridPoint >::iterator Il;
  for( Il=l.begin(); Il!=l.end(); Il++ ) rt.push_back( *Il );

  return;

}

//_______________________________________________________________________________
istream& operator>>(istream& in, CsRTRelation& rt)
{

  // clean old registered values
  rt.rtGrid_.clear();  rt.hasGrid_ = false;
  rt.rtGReg_.clear();  rt.hasRegularGrid_ = false;
  rt.rtPars_.clear();  rt.hasRTPars_ = false;

  CsOpt* opt = CsOpt::Instance();

  string RTtype;
  in >> RTtype;
  if( ! in.good() ) {
    CsErrLog::msg(elError, __FILE__, __LINE__, 
		    "%s: empty stream. << DB.",
		    rt.d_->GetTBName().c_str() );
    return in;
  }
  
  //=== first try reading RT Relation Parametrisation ===
  if( RTtype == "RT" ||       // Kept for backward compatibility <- assuming _FIT1_
      RTtype == "RTFit1" ||   // RTFit1 (bool _FIT1_=true) is polynom fit
      RTtype == "RTFit2" ) {  // RTFit2 (bool _FIT2_=true) is hyperbolic tangeant
    
    //=== Read RTRelation parameters
    rt.rtPars_.clear();
    while( in.good() && !in.eof() ) {
      double rtp;
      in >> rtp;
      if( in.good() ) rt.rtPars_.push_back( rtp );
    }
    
    //=== Check number of parameters
    if (rt.rtPars_.size() < 2) {
      rt.rtPars_.clear();
      rt.hasRTPars_ = false;
      if( in.fail() ) in.clear();  
      CsErrLog::msg(elError, __FILE__, __LINE__, 
		    "%s: Not enough RTParameters read << DB.",
		    rt.d_->GetTBName().c_str() );
      return in;
    }

    //=== Check istream status
    // When in state is ios::fail, this is probably due to the fact that other DB info follows.
    // The state is then reset to ios::good to avoid exceptions ...
    // This is dirty but I don't what else to do.
    if( in.fail() ) in.clear();  
    rt.hasRTPars_ = true;
    CsErrLog::msg(elInfo,__FILE__,__LINE__,"%s: parameters read.",
		  rt.d_->GetTBName().c_str(),rt.rtPars_.size());
    
    //===========================================
    if( RTtype == "RTFit1" || RTtype == "RT"  ) {
      rt._FIT1_ = true;
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
		    "%s: << DB: Parameterisation is _FIT1_.",
		    rt.d_->GetTBName().c_str());
    }
    
    //========================
    if( RTtype == "RTFit2" ) {
      rt._FIT2_ = true;
      CsErrLog::msg(elInfo,__FILE__,__LINE__,
		    "%s: << DB: Parameterisation is _FIT2_.",
		    rt.d_->GetTBName().c_str());
    } 

    return in;
  }

  //=== second try to read RT Grid ===
  else if (RTtype == "RTGrid"){

    rt.rtGrid_.clear();
    while( in.good() && !in.eof() ) {
      RTGridPoint rtg;
      in >> rtg;
      if( in.good() ) rt.rtGrid_.push_back( rtg );
    }

    if (rt.rtGrid_.size() < 2 ) {
      rt.rtGrid_.clear();
      rt.hasGrid_ = false;
      if( in.fail() ) in.clear();  
      CsErrLog::msg(elError,__FILE__,__LINE__, 
		    "%s: Not enough RT grid points << DB.",
		    rt.d_->GetTBName().c_str() );
      return in;
    } else rt.hasGrid_ = true;

    rt._sortGrid( rt.rtGrid_ );
    if (opt->getOpt( "RT", "useRegularGrids" ) ) rt._makeRegular();

    CsErrLog::msg(elInfo, __FILE__, __LINE__,"%s: %d%s points read.",
		  rt.d_->GetTBName().c_str(),rt.rtGrid_.size(),
		  rt.hasRegularGrid_ ? " regular":"");

    // When in state is ios::fail, this is probably due to the fact that other DB info follows.
    // The state is then reseted to ios::good to avoid exceptions
    // This is dirty but I don't what else to do.
    if( in.fail() ) in.clear();     

    return in;
  }
  
  //=== third: RTType is not understood ===
  else {
    CsErrLog::Instance()->msg(elError, __FILE__, __LINE__, 
      "CsRTRelation << from DB ==> %s Wrong RT relation type %s.",
      rt.d_->GetTBName().c_str(), RTtype.c_str() );
    return in;
  }
}

