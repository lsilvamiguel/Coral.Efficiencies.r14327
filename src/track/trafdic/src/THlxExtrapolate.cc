// $Id: THlxExtrapolate.cc 14069 2015-09-17 20:44:46Z lsilva $

/*!
  Interface to different helix extrapolation methods
*/

#include <float.h>
#include <stdio.h>
#include "TOpt.h"
#include "TSetup.h"
#include "TDisplay.h"
#include "THlx.h"
#include "CsGeom.h"

using namespace std;

bool THlx::Extrapolate(THlx& Hout, bool mmap) const
{

  if( &Hout == this){
    cout<<"THlx::Extrapolate ==> Output helix must differ from 'this' helix"<<endl;
    assert(false);
  }

  if(fabs(Hout.Hpar[0] - this->Hpar[0]) < FLT_EPSILON) { // trivial case: extrapolate to the same X
    Hout = *this;
    Hout.path     = fabs(Hout.Hpar[0] - this->Hpar[0]);
    Hout.radLenFr = 0;
    Hout.eLoss = 0;
    return(true);
  }

  double X_final = Hout.Hpar[0];

  if(TOpt::ReMode[20] > 0 && mmap && this->with_mom()) { // "use of material map" is ON and track with momentum 

    const TSetup& setup = TSetup::Ref();
    bool print(false);  if(TOpt::Print[6] > 0) print = true; // debug printout flag
    bool ok(false);
    THlx   Hfrom;
    float  Len, RadLen, Step;
    bool   DIREC   = ( X_final >= this->Hpar[0] ? true : false );

    CsMaterialMap *mp = CsGeom::Instance()->getCsMaterialMap();
    if(mp == NULL){
      cout<<"THlx::Extrapolate ==> Extrapolation with material map had been requested"<<endl
	  <<" but material map is not available!"<<endl;
      assert(false);
    }

    // If start point and and point are out ouf the map
    // and there is no material maps on the way,
    // do ordinary extrapolation and return.
    if(! setup.InMaterialMap(this->Hpar[0]) &&
       ! setup.InMaterialMap(X_final)       &&
       ! setup.IsMaterialBetween(this->Hpar[0], X_final)) {
      return(this->Extrap(X_final, Hout));
    }

    if(print) cout<<"Requested extrapolation from "<<Hpar[0]<<" \t to "<<Hout.Hpar[0]<<endl;

    int i=0;
    Hout = *this;

    double tot_len(0), tot_rad_len(0), tot_eloss(0);
    do{ // iterative extrapolation through material map
 
      // get media properties and recommended step
      //Hout.Print("Input RadLen");
      //cout<<"DIREC = "<<DIREC<<endl;
      mp->getRadLength(Hout,DIREC,RadLen,Step);

      // draw point (if requested) where material map was called
      if(TDisplay::Ptr() != NULL) { // TDisplay object exists
	if(TDisplay::Ref().Mmaptraj()){
	  double r[] = {Hout(0),Hout(1),Hout(2)};
	  TDisplay::Ref().point(r,110);
	}
      }

      if(print) cout<<"   iter. # "<<i<<"  recommended step is "<<Step<<endl;
      // prepare next step
      Hfrom = Hout;
      if(DIREC) Hout.Hpar[0] + Step < X_final ? Hout.Hpar[0] += Step : Hout.Hpar[0] = X_final;
      else      Hout.Hpar[0] + Step > X_final ? Hout.Hpar[0] += Step : Hout.Hpar[0] = X_final;

      // do extrapolation
      ok = Hfrom.Extrap(Hout.Hpar[0], Hout);
      if(!ok) return (false); // propagation error.

      // passed distance
      Len = Hfrom.Dist(Hout); // straight line distance
      Len = Hout.Path();      // trajectory path
      if(print) cout<<"          "<<"  Passed  "<<Len<<" cm.   X final = "<<Hout.Hpar[0]
		    <<"  rad. len. = "<<RadLen<<"  cm."<< endl;
      tot_len     += Len;
      tot_rad_len += Len/RadLen;


      // Add multiple scattering contribution to cov. matrix
      Hout.AddNoise( Len, RadLen);

      // Add energy losses 
      if(TOpt::ReMode[20] > 1) {
	float eloss = mp->getdE(Hfrom,Len); tot_eloss += eloss;
	// Propagate the sigma from straggling in energy loss to the momentum
	// uncertainty in covariance matrix
	if (TOpt::ELossStraggling)
	  Hout.Hcov[14] = Hout.Hpar[5]*Hout.Hpar[5]*Hout.Hpar[5]*Hout.Hpar[5]
	    * ( Hfrom.Hcov[14]/Hfrom.Hpar[5]/Hfrom.Hpar[5]/Hfrom.Hpar[5]/Hfrom.Hpar[5] + mp->getdEStraggling(Hfrom,Len) );
	if (DIREC) Hout.Hpar[5] /= ( 1 - eloss / Hout.Mom() );
	else       Hout.Hpar[5] /= ( 1 + eloss / Hout.Mom() );
      }

    next:
      if( ++i > 5000 ){ 
	cout<<"THlx::Extrapolate ==> The number of iteration exceeded 5000."<<endl;
	cout<<"Last step was "<<Step<<" at X = "<<Hout.Hpar[0]<<"  X final is "<<X_final<<endl;
	break;
      }

    } while ( fabs(Hout.Hpar[0] - X_final) > FLT_EPSILON ); // end of iteration loop

    // save total length passed and rad. len. fraction 
    Hout.path =     tot_len;
    Hout.radLenFr = tot_rad_len;
    Hout.eLoss =    tot_eloss;

    return true;

  } else { // use of material map is OFF or track without momentum.

    return(this->Extrap(X_final, Hout));  
    
  }
  return(false);
}
