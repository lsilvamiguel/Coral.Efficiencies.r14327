#include <cmath>
#include "Geometry.h"
#include "GeomPlaneMumega.h"
#include "GeomPlanePMumega.h"
#include "GeomPlaneMwpc.h"
#include "GeomPlaneGEM.h"
#include "GeomPlanePGEM.h"
#include "GeomPlaneSili.h"
#include "GeomPlaneScifiJ.h"
#include "GeomPlaneScifiG.h"
#include "GeomPlaneMuonWallA.h"
#include "TMath.h"

#include <fstream>
#if __GNUC__ > 3 || __GNUC__ == 3
# include <sstream>
#else
# include <strstream>
#endif // __GNUC__
#include <vector>
#include <map>
#include <iostream>


const float Geometry::fScale = 10.;

Geometry::Geometry(const char *geofilename) : fFilename(geofilename) {

  ifstream infile(geofilename);
  if(!infile) {
    std::cerr<<"WARNING Geometry::Geometry. Cannot open "<<geofilename<<std::endl;
    throw "Geometry::Geometry, cannot open file";
  }  
  else { 

    std::cout<<std::endl;
    std::cout<<"Geometry found : "<<geofilename<<std::endl;
    std::cout<<"********************************************************************"<<std::endl;

    std::string line;

    std::map<int,TMatrix*> rotmats;

    // reading file line by line
    while(1) {

      if (infile.bad() || infile.eof())
        break;

      getline(infile, line);

      if(line.empty()) continue;

#if __GNUC__ > 3 || __GNUC__ == 3
      std::istringstream lin(line);  
#else
      istrstream lin(line);  
#endif
      
      std::string tag;
      lin>>tag;
      if(tag == "rot") {
	TMatrix *rotM=new TMatrix(3,3);
	int index; lin>>index;
	for( int i=0; i<3; i++ ) { 
	  for( int j=0; j<3; j++ ) { 
	    lin >> (rotM->operator()(i,j));
	  }
	}
	//Geant2Coool(rotM);
	rotmats.insert(std::pair<int, TMatrix*>(index,rotM));
      }
      else if(tag == "det") {

	int    id;        lin >>id;           // detector number
	char tbname[10];    lin >>tbname;       // det. TB name
	char   name[10];    lin >> name;        // detector name
	int    unit;      lin >> unit;        // detector number in station 
	int    type;      lin >> type;        // detector type
	double rdLen;     lin >> rdLen;       // radiation length
	double xsiz;      lin >> xsiz; // detector size (cm)
	double ysiz;      lin >> ysiz; 
	double zsiz;      lin >> zsiz; 
	//Geant2Coool(xsiz,ysiz,zsiz);
	double xcm;       lin >> xcm;  // detector centre (MRS) (cm)
	double ycm;       lin >> ycm; 
	double zcm;       lin >> zcm; 
	//Geant2Coool(xcm,ycm,zcm);
	int    rotMNum;   lin >> rotMNum;     // rotation matrix number
	double wirD;      lin >> wirD;        // 1st wire offset (mm)
	//	Geant2Coool(wirD);
	double ang;       lin >> ang;         // angle of wires in DRS (degrees)
	int    nWir;      lin >> nWir;        // number of wires
	double wirP;      lin >> wirP;        // wires pitch
	//Geant2Coool(wirP);
	double eff;       lin >> eff;         // detector efficiency
	double bkg;       lin >> bkg;         // detector background
	double tGate;     lin >> tGate;       // detector time gate

	// We will not use DRS, only wire angle.
	
	std::map<int, TMatrix*>::iterator mymat= rotmats.find(rotMNum);
	if(mymat != rotmats.end()) {
	  ang += DRSAngle(mymat->second);
	} 
	else {
	  std::cerr<<"WARNING Geometry::Geometry. No rotation matrix with label "<<rotMNum<<" check "<<fFilename<<std::endl;
	}

	std::map<std::string, GeomPlane*>::iterator mydet = fAllPlanes.find(tbname);
	if(mydet != fAllPlanes.end()) {
	  try {
	    mydet->second->AddSubPlane(id,tbname,nWir,xcm,ycm,zcm,
				       xsiz,ysiz,zsiz,ang,wirP); 
	  }
	  catch(const char* errmsg) {  
	    std::cerr<<errmsg<<std::endl;
	    throw;
	  }
	}
	else {
	  GeomPlane *plane=0;
	  if(strncmp(tbname,"MM",2)==0) {
	    plane = new GeomPlaneMumega(id,tbname,nWir,xcm,ycm,zcm,
					xsiz,ysiz,zsiz,ang,wirP);
	  }
	  else if(strncmp(tbname,"MP",2)==0) {
	    plane = new GeomPlanePMumega(id,tbname,nWir,xcm,ycm,zcm,
					 xsiz,ysiz,zsiz,ang,wirP);
	  }
	  else if(strncmp(tbname,"PB",2)==0 || 
		  strncmp(tbname,"PS",2)==0 ||
		  strncmp(tbname,"PA",2)==0)
	    plane = new GeomPlaneMwpc(id,tbname,nWir,xcm,ycm,zcm,
				      xsiz,ysiz,zsiz,ang,wirP);
	  else if(strncmp(tbname,"GM",2)==0)
	    plane = new GeomPlaneGEM(id,tbname,nWir,xcm,ycm,zcm,
				     xsiz,ysiz,zsiz,ang,wirP);
	  else if(strncmp(tbname,"GP",2)==0)
	    plane = new GeomPlanePGEM(id,tbname,nWir,xcm,ycm,zcm,
				     xsiz,ysiz,zsiz,ang,wirP);
	  else if(strncmp(tbname,"SI",2)==0)
	    plane = new GeomPlaneSili(id,tbname,nWir,xcm,ycm,zcm,
				      xsiz,ysiz,zsiz,ang,wirP);
	  else if(strncmp(tbname,"FI",2)==0) {
	    char fin_str[3]; 
	    for (int i=0; i<2; i++) fin_str[i]=tbname[i+2];
                    fin_str[2] = 0;
	    int FI_num = atoi(fin_str);
	    if ( FI_num<=4 || ( FI_num==6 && tbname[4]=='V' ) || (FI_num == 15) )
	      plane = new GeomPlaneScifiJ(id,tbname,nWir,xcm,ycm,zcm,
					  xsiz,ysiz,zsiz,ang,wirP);
	    else
	      plane = new GeomPlaneScifiG(id,tbname,nWir,xcm,ycm,zcm,
					  xsiz,ysiz,zsiz,ang,wirP);
	  }
	  else 
	    plane = new GeomPlane(id,tbname,nWir,xcm,ycm,zcm,
				  xsiz,ysiz,zsiz,ang,wirP);
	  
	  fAllPlanes[tbname]=plane;
	}
      }
      else if(tag == "dead"){
      }
    }
    infile.close();
  }
}

GeomPlane* Geometry::GetPlane(const char* tbname) {
  
  std::map<std::string,GeomPlane*>::const_iterator mytracker = fAllPlanes.find(tbname);
  if(mytracker != fAllPlanes.end())
    return mytracker->second;
  else {
    std::cerr<<"WARNING Geometry::GetPlane. "<<tbname<<" not found in "<<fFilename<<std::endl;
    return 0;
  }
}


void Geometry::Geant2Coool(double &x, double &y, double &z) {
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical

  double tmp = x;
  x          = y*fScale; 
  y          = z*fScale;
  z          = tmp*fScale;  
}

void Geometry::Geant2Coool(double &x) {
  // cm to mm 
  x *= fScale;
}

void Geometry::Geant2Coool(TMatrix &m) {
  // cm to mm 
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical
  TMatrix r(3,3);
  double set[] = { 0, 1, 0, 0, 0, 1, 1, 0, 0 };
  for( int i=0; i<9; i++ ) r( i/3, i%3) = set[i];

  TMatrix rt(TMatrix::kTransposed, r);
  TMatrix mtemp(r);
  mtemp *= m;
  mtemp *= rt;
  m = mtemp;
}

double Geometry::DRSAngle(const TMatrix *m) {
  double theta=asin(m->operator()(2,1))*180/TMath::Pi();
  return theta;
}





