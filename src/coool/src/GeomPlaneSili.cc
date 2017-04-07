#include "GeomPlaneSili.h"
#include <queue>
#include <cassert>
#include <iostream>
#include "math.h"


//ClassImp(GeomPlaneSili);

std::set<CCluster1>& GeomPlaneSili::Clusterize(const std::map<int,CDigit1>& digits) {
 
  if (firstread) {
    Resolution1 = 0.0015;
    Resolution2 = 0.0015;
    Resolution3 = 0.0025;
    char parFileName[100];
    char tps[10][100];
//     sprintf(parFileName, "/afs/e18.ph.tum.de/compass/analysis/primakoff/timecalib/analysis_coool/configs/ResolPar/Resolution_%s.par", this->GetNameC());  
//     std::cout << "Looking for "<<parFileName;
//     std::ifstream resolpar(parFileName);
//     if (resolpar) {
//       resolpar >> Resolution1 >> tps[0];
//       resolpar >> Resolution2 >> tps[1];
//       resolpar >> Resolution3 >> tps[2];
//       std::cout << " \n 1Str: "<< Resolution1
// 	   << " 2Str: "<< Resolution2
// 	   << " 3Str: "<< Resolution3 << std::endl;
//     }
//     else{std::cout<<" ...not found"<<std::endl;}
//     resolpar.close();
    firstread=false;
  } 
  
  fClusters.clear();
  
  
  for (std::map<int,CDigit1>::const_iterator id = digits.begin(); 
       id!=digits.end(); ) {
    
    const CDigit1& digit = id->second;
    std::vector<double> data;
    for (int i=0; i<3; i++) data.push_back(digit.dt[i]);
    
    //if (!strcmp(this->GetName(),"SI01X1__"))
    //std::cout.form("%3d. strip amplitudes: %5.1f %5.1f %5.1f\n", 
    //      digit.ch, data[0], data[1], data[2]);
    
    int firstchan = digit.ch;
    int lastchan = firstchan;
    int ndigits=1;             // Strip multiplicity
    int delta=0;
    double sum_ampl     = data[2];
    double sum_ampl_pos = data[2]*Wire2Pos(firstchan);
    do {
      id++;
      if (id!=digits.end()) {
        const CDigit1& htmp = id->second;
        delta = htmp.ch - lastchan;
	double s0 = data[0];
	double s1 = data[1];
	double s2 = data[2];
	double h0 = htmp.dt[0];
	double h1 = htmp.dt[1];
	double h2 = htmp.dt[2];
	//
	// coral: 
	//    if((SiIter==SiDigList.end())||(SiIter->wire!=lastwire+1) 
	//       || ! ( SiIter->a0+SiIter->a1+SiIter->a2 < 8 ||
	//	      amp0+amp1+amp2                   < 8 ||
	//	      ( fabs(SiIter->a0/SiIter->a2-amp0/amp2)<0.5 &&
	//		fabs(SiIter->a1/SiIter->a2-amp1/amp2)<0.4 ) )
	//
        if (delta <= 1   &&
	    ( h0+h1+h2 < 8 || 
	      s0+s1+s2 < 8 ||
	      (fabs(h0/h2-s0/s2)<0.5 && fabs(h1/h2-s1/s2)<0.4)
	      )
	    ) {
	  ndigits++;
	  for (int i=0; i<3; i++) data[i] += htmp.dt[i];
	  sum_ampl     += h2;
	  sum_ampl_pos += h2*Wire2Pos(htmp.ch);
	  //if (!strcmp(this->GetName(),"SI01X1__")) {
	  //std::cout.form("%3d. strip amplitudes: ", htmp.ch);
	  //std::cout.form("%5.1f %5.1f %5.1f (sum: %5.1f %5.1f %5.1f)\n", 
	  //      htmp.dt[0], htmp.dt[1], htmp.dt[2],
	  //      data[0], data[1], data[2]);
	  //}
          lastchan = htmp.ch;
        }
      }
    } while ((id!=digits.end()) && (delta<=2));
    
    // Add timing information (sample ratios)
    data.push_back(data[0]/data[2]); //a0/a2
    data.push_back(data[1]/data[2]); //a1/a2
    
    assert(lastchan>=firstchan);
    //double wire = (firstchan+lastchan)/2.;
    double pos = sum_ampl_pos/sum_ampl; //Wire2Pos(wire);
    double res = Resolution(0);
    if      (ndigits==1) res=Resolution1;
    else if (ndigits==2) res=Resolution2;
    else                 res=Resolution3;
    //std::cout << "Resolution "<<res<<std::endl;
    // this is not optimal : cluster data is data of the first digit...
    
    fClusters.insert(CCluster1(GetID(),pos,ndigits,res,data));    
    
  }
  return fClusters;
}



