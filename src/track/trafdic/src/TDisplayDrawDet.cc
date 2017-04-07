// $Id: TDisplayDrawDet.cc 13508 2012-07-10 17:21:50Z ybedfer $

#include "CsErrLog.h"
#include "CsDetector.h"
#include "CsPixelGEMDetector.h"
#include "CsPixelMumegaDetector.h"
#include "TDisplay.h"
#include "TAlgo.h"
#include "TSetup.h"
#include "TOpt.h"
#include "TConstants.h"
#include "higz.h"

using namespace std;

/*!
  \brief Draw detector planes (TPlane objects) on a projection plane in the event display.
  \arg \a ifl: if !=0, no actual drawing, wait instead for a mouse click locating a detector plane and return its index (in "TSetup::vecPlane").

  - GP P-pixel planes are drawn in all of 4 proj.: X, Y and 2 extras, defined by "TraF" options.
  - Upon option, MP M-pixel planes are drawn in one extra proj., viz. along their v axis.
  - A unique label for the 2 layers of a double layer.
  - Get the angle of the current event display proj. from "TSetup::vecProj".
*/

int TDisplay::DrawDet(int ifl)
{ 

  float x[2], y[2], xx[2], yy[2], xc, yc, dx, dy, dz; 
  const double eps=0.0001;
  int n,istat,nt;
  float px,py;

  const TSetup& setup = TSetup::Ref();

  if(ifl != 0){
    IRQLC(1,20,istat,nt,px,py);
    if(istat == 0) return(-1);
  } 

  //        ********** ANGLE of CURRENT DISPLAYED PROJECTION **********
  int iAngle = (int)floor(setup.vProj()[Proj]/10.+.5);
  double Angle = iAngle*M_PI/180;

  float xbar = 0.002*(Xvf[1]-Xvf[0]), ybar = 0.002*(Yvf[1]-Yvf[0]);
  int zone, ip_picked; float xprv; float dist, dist_min;
  for (zone = 0, xprv = -1e6, ip_picked = -1, dist_min = FLT_MAX;
       zone<(int)setup.vIplFirst().size(); zone++) {

    //      *************** LOOP ON ZONES ***************

    //                 ***** PIXELGEM/MMs *****
    // Pixels could be projected onto any projection. In TraFDic's proj. search,
    // they are projected on 4 different proj's: 0 and 90 degrees, and 2 others,
    // which +/-45 deg. defaults can be overriden by option on a per zone basis.
    double aPixProjs[] = {0,90,45,-45};
    if (TOpt::iCut[29]&1<<zone) { // Zones where special pixelGEM/MM projection
      aPixProjs[2] = TOpt::dCut[82]; aPixProjs[3] = TOpt::dCut[83];
    }
    int drawPPixDet, i03, vProj, iproj; for (i03=drawPPixDet = 0; i03<=3; i03++)
      if (abs(floor(aPixProjs[i03]+.5)-iAngle)<=1) {
	drawPPixDet = 1; break;
      }
    for (iproj = 0, vProj = -1; iproj<(int)setup.vProj().size(); iproj++) {
      int jAngle = (int)floor(setup.vProj()[iproj]/10.+.5);
      if (abs(abs(iAngle-jAngle)-90)<=1) {
	vProj = iproj;
	break; // Breaking here prevents the -90 proj., which contains only BMS, from being retained.
      }
    }

    for (int ipl = 0; ipl<(int)setup.vPlane().size(); ipl++) {

      // *********************************************************************
      // *************** LOOP ON ALL DETECTORS ZONE AFTER ZONE ***************
      // *********************************************************************

      const TPlane &p = setup.vPlane(ipl);
      if (ipl<setup.vIplFirst()[zone]) continue;
      if (setup.vIplLast()[zone]<ipl) break;
      if ((TOpt::Graph[5]&0x1) &&// (Optionally) Do not draw switched OFF planes
	  p.IFlag==0) continue;
      const TDetect &d = setup.vDetect(p.IDetRef);
      int isPixDet = 0; if (d.IType==29) isPixDet = 1; // GP P-pixel
      else if ((TOpt::Graph[5]&0x2) && d.IType==32)
	isPixDet = 2;                    // (upon option) MP M-pixel
      bool isRPD = d.Name.find("RP")==0;
      if (p.IProj==Proj ||
	  (isPixDet==1 && drawPPixDet) || (isPixDet==2 && p.IProj==vProj) ||
	  isRPD) {

	// *************** CURRENT TPLANE in PROJ. "Proj" ***************

	static double orig, pitch; static int nWires;
	if (isPixDet) {                                     // ***** PixelGEM/MM
	  CsDetector *csDet = d.PtrDet();
	  double a = atan2(d.Sa,d.Ca)-Angle, ca = cos(a), sa = sin(a);
	  double uOrig = csDet->getWirD()/10, uPitch = csDet->getWirP()/10;
	  double u0 = uOrig-uPitch/2, u1 = uOrig+uPitch*(csDet->getNWir()+.5);
	  // PixelGEM/MM desciption is simplified (anyway, the rest of the
	  // present detector drawing is also quite simplistic).
	  double vOrig, vPitch; int vNWir; if (d.IType==29) {
	    CsPixelGEMDetector *pg = dynamic_cast<CsPixelGEMDetector*>(csDet);
	    if (!pg) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TDetect \"%s\" w/ type=29 is not a CsPixelGEMDetector",d.Name.c_str());
	    vOrig = pg->getWirDV()/10; vPitch = pg->getWirPV()/10;
	    vNWir = pg->getNWirV();
	  }
	  else {
	    CsPixelMumegaDetector *pm =
	      dynamic_cast<CsPixelMumegaDetector*>(csDet);
	    if (!pm) CsErrLog::msg(elFatal,__FILE__,__LINE__,
 "TDetect \"%s\" w/ type=32 is not a CsPixelMumegaDetector",d.Name.c_str());
	    vOrig = pm->getWirDV()/10; vPitch = pm->getWirPV()/10;
	    vNWir = pm->getNWirV();
	  }
	  double v0 = vOrig-vPitch/2, v1 = vOrig+vPitch*(vNWir+.5);
	  double y0p = u0*ca-v0*sa, y0m = u0*ca-v1*sa;
	  double y1p = u1*ca-v1*sa, y1m = u1*ca-v0*sa;
	  if (fabs(y1p-y0p)>fabs(y1m-y0m)) { y[0] = y0p; y[1] = y1p; }
	  else                             { y[0] = y0m; y[1] = y1m; }
	  x[0]=x[1] = d.X(0);
	}
	else if (isRPD) {                              // ***** RPD
	  x[0] = d.X(0)-d.Siz(2)/2; x[1] = d.X(0)+d.Siz(2)/2;
	  y[0]=y[1] = d.Pitch*d.Nwires/2/M_PI;
	}
	else if (d.IType==30) // Dummy MP plane: skip
	  continue;
	else {                                         // ***** ANY OTHER TPlane
	  orig = d.Uorig; pitch = d.Pitch; nWires = d.Nwires;
	  y[0] = orig-pitch/2.; y[1] = orig+d.Range+pitch/2.;
	  x[0]=x[1] = d.X(0);
	}
	if (ifl==0) {
	  //                 ********** DRAWING MODE **********
	  double pp = fabs(pitch);
	  if (DrOpt[7]>0 &&          // ***** DETAILED DETECTOR DRAWING *****
	      d.IType!=29 && d.IType!=32 /* No details for pixels yet */ &&
	      !isRPD /* No details for RPD yet */) {
	    ISPLCI(1); if (p.IFlag==0) ISPLCI(105); // plane is off
	    for (int iwir = 0; iwir<nWires; iwir++) { // Loop over wires
	      if (x[0]<Xvf[0] || x[0]>Xvf[1]) continue;
	      float ywir(0);
	      if (d.W2Pos.empty()) {
		ywir = orig + pitch*iwir;
	      } else {                 // Variable pitch det.
		ywir = d.W2Pos[iwir];
	      }
	      if (ywir<Yvf[0] || ywir>Yvf[1]) continue;
	      if (d.Kind==1) {                 // ***** DRIFT DETECTOR... *****
		if (d.Name[0]=='S' && d.Name[1]=='T') {         // *****...STRAW
		  IGARC(x[0], ywir, float(pp/2.), float(pp/2.), 0., 0.);
		}
		else {                                           // *****...ELSE
		  float x1=x[0]-pp/2.;  float x2=x[0]+pp/2.;
		  float y1=ywir-pp/2.;  float y2=ywir+pp/2.;
		  IGBOX(x1,x2,y1,y2);
		}
	      }
	      else {                        // ***** STRIP/WIRE DETECTORS *****
		if(d.W2Pos.empty()) { // not a variable pitch detectors
		  float x1=x[0]-pp/10.;  float x2=x[0]+pp/10.;
		  float y1=ywir-pp/4. ;  float y2=ywir+pp/4.;
		  IGBOX(x1,x2,y1,y2);
		} else {               // variable pitch det.
		  float x1=x[0]-d.Siz(0)/2.;  float x2=x[0]+d.Siz(0)/2.;
		  float y1=ywir-pp/4.; float y2=ywir+pp/4.;
		  IGBOX(x1,x2,y1,y2);
		}
	      }
	    } // End of loop over wires
	  }
	  else {                         // ***** SIMPLIFIED DETECTORS DRAWING
	    int icol = p.IFlag==0 ? 105 : 1;
	    ISPLCI(icol); ISLN(1); ISLWSC(3);
	    IPLm(2,x,y); // Draw plane
	    if (isRPD) {
	      y[0] *= -1; y[1] = y[0]; IPLm(2,x,y); y[0] *= -1; y[1] = y[0];
	    }
	    // Have a unique label for the 2 layers of a double layers...
	    // ...in order to avaid overwritting of upstream layer name by
	    // downstream one.
	    string *dlString = NULL;
	    if (d.Name[0]=='D' ||
		d.Name[0]=='M' && d.Name[1]=='A') {
	      if (d.Name[5]=='2') {
		char dlName[] = "DC01X12"; strncpy(dlName,d.Name.c_str(),5);
		dlName[5] = '1'; dlName[6] = '2';
		dlString = new string(dlName);
	      }
	    }
	    else if (d.Name[0]=='S' ||
		     d.Name[0]=='M' && d.Name[1]=='B') {
	      if (d.Name[6]=='d') {
		char dlName[] = "ST03X1uda"; strncpy(dlName,d.Name.c_str(),6);
		dlName[6] = 'u'; dlName[7] = 'd'; dlName[8] = d.Name[7];
		dlString = new string(dlName);
	      }
	    }
	    else if (d.Name[0]=='H' && (d.Name[7]=='u' || d.Name[7]=='d')) {
	      if (d.Name[7]=='d') {
		char dlName[] = "HI04X1_d"; strncpy(dlName,d.Name.c_str(),6);
		dlName[6] = 'u'; dlName[7] = 'd';
		dlString = new string(dlName);
	      }
	    }
	    else dlString = new string(d.Name);
	    if (dlString)
	      DrawString
		(" . "+(*dlString),max(x[0],x[1]),max(y[0],y[1]),icol,1.E15);
	    delete dlString;
	    ISLWSC(1);
	    // Draw "bars"
	    xx[0]=x[0]; xx[1]=x[0]-xbar; yy[0]=y[0]; yy[1]=y[0]-ybar*0.02*d.XR(2); IPLm(2,xx,yy);
	    xx[0]=x[0]; xx[1]=x[0]+xbar; yy[0]=y[0]; yy[1]=y[0]-ybar*0.02*d.XR(2); IPLm(2,xx,yy);
	    xx[0]=x[1]; xx[1]=x[1]-xbar; yy[0]=y[1]; yy[1]=y[1]+ybar*0.02*d.XR(2); IPLm(2,xx,yy);
	    xx[0]=x[1]; xx[1]=x[1]+xbar; yy[0]=y[1]; yy[1]=y[1]+ybar*0.02*d.XR(2); IPLm(2,xx,yy);
	  }
	}
	else {           // ********** PICKING UP MODE **********
	  dist=TAlgo::Point2Line(px,py,x,y);
	  if(dist < dist_min){
	    dist_min=dist; ip_picked=ipl;
	  }
	}
      }
      else if (DrOpt[4]>0 &&
	       d.IType!=29 && d.IType!=32 /* Not for pixel det. yet */) {

	// *************** DRAWING PLANES from OTHER PROJ ***************

	//float range = d.Siz(1) > d.Siz(2) ? d.Siz(1) : d.Siz(2);
	float range = d.Range;
	x[0]=x[1]=d.X(0);
	double Ca = cos(Angle), Sa = sin(Angle);
	y[0]=Ca*d.X(1)+Sa*d.X(2)-range/2.;
	y[1]=Ca*d.X(1)+Sa*d.X(2)+range/2.;
	//y[0]=d.Uorig - d.Pitch/2.;
	//y[1]=d.Uorig - d.Pitch/2. + d.Pitch*d.Nwires;
	if(ifl == 0) { // drawing
	  if (p.IFlag!=0){ // if the plane is NOT OFF
	    ISPLCI(108); IPLm(2,x,y);  DrawString(" . "+d.Name, max(x[0],x[1]), max(y[0],y[1]), 108, 1.E15);
	  }
	}
	else { // picking up
	  dist=TAlgo::Point2Line(px,py,x,y);
	  if(dist < dist_min){
	    dist_min=dist; ip_picked=ipl;
	  }
	}
      }

    } // End of loop over planes in zone
  }  // End loop over zones

  if(ifl == 0) {
    ISPLCI(1);
    IUWK(0,1);
  }

  return ip_picked;

}
