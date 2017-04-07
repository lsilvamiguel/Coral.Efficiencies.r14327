#include "MManGeometry.h"

#include <iostream>
#include <fstream>
#include <strstream>
#include <vector>
#include <map>


const float MManGeometry::fScale = 10.;

MManGeometry::MManGeometry(const char *geofilename) : fFilename(geofilename) {

  ifstream in(geofilename);
  if(!in) {
    throw "MManGeometry::MManGeometry, cannot open file";
  }  
  else { 
    std::cerr<<"********************************************************************"<<std::endl;
    std::cerr<<"MManGeometry found : "<<geofilename<<std::endl;
    std::cerr<<"********************************************************************"<<std::endl;

    const int linesize=1000;
    char s[linesize];
    int pos=0;    

    std::map<int,TMatrix*> rotmats;

    while(1) {
      // reading one line
      in.seekg(pos);
      in.get(s,linesize);
      if(!in) {
	std::cerr<<"End of file"<<std::endl;
	break;
      }
      pos = (int)in.tellg() + 1;
      if(std::string(s).empty()) continue;
      std::istrstream lin(s);  
      
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
	  std::cerr<<"WARNING MManGeometry::MManGeometry. No rotation matrix with label "<<rotMNum<<" check "<<fFilename<<std::endl;
	}

	std::map<std::string, Tracker*>::iterator mydet = fAllTrackers.find(tbname);
	if(mydet != fAllTrackers.end()) {
	  try {
	    mydet->second->AddSubDetector(id,tbname,nWir,xcm,ycm,zcm,
					  xsiz,ysiz,zsiz,ang,wirP,wirP); 
	  }
	  catch(const char* errmsg) {  
	    std::cerr<<errmsg<<std::endl;
	    throw;
	  }
	}
	else {
	  Tracker *tracker = new Tracker(id,tbname,nWir,xcm,ycm,zcm,
					 xsiz,ysiz,zsiz,ang,wirP,wirP);	
	  fAllTrackers[tbname]=tracker;
	}
      }
      else if(tag == "dead"){
      }
    }
    in.close();
  }
}

const Tracker* MManGeometry::GetTracker(const char* tbname) const {
  
  std::map<std::string,Tracker*>::const_iterator mytracker = fAllTrackers.find(tbname);
  if(mytracker != fAllTrackers.end())
    return mytracker->second;
  else {
    std::cerr<<"WARNING MManGeometry::GetTracker. "<<tbname<<" not found in "<<fFilename<<std::endl;
    return 0;
  }
}


void MManGeometry::Geant2Coool(double &x, double &y, double &z) {
  // From GeantRS : X // beam, Z vertical
  // To MRS :       Z // beam, Y verical

  double tmp = x;
  x          = y*fScale; 
  y          = z*fScale;
  z          = tmp*fScale;  
}

void MManGeometry::Geant2Coool(double &x) {
  // cm to mm 
  x *= fScale;
}

void MManGeometry::Geant2Coool(TMatrix &m) {
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

double MManGeometry::DRSAngle(const TMatrix *m) {
  double theta=asin(m->operator()(2,1))*180/TMath::Pi();
  return theta;
}





