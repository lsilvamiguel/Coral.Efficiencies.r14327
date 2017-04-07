#include "GeomPlaneScifiJ.h"
#include <queue>
#include "math.h"

//ClassImp(GeomPlaneScifiJ);

std::set<CCluster1>& GeomPlaneScifiJ::Clusterize(const std::map<int,CDigit1>& digits) {

  fClusters.clear();

  for (std::map<int,CDigit1>::const_iterator id = digits.begin(); 
       id!=digits.end(); ) {
    
    const CDigit1& digit = id->second;
    std::vector<double> data;
    for (unsigned int i=0; i<digit.dt.size();i++)
      data.push_back(digit.dt[i]);

    int firstchan = digit.ch;
    int lastchan = firstchan;
    int ndigits=1;
    int delta=0;
    do {
      id++;
      if (id!=digits.end()) {
        const CDigit1& htmp = id->second;
        delta = htmp.ch - lastchan;
        if (delta <= 1) {
	  ndigits++;
          lastchan = htmp.ch;
        }
      }
    } while ((id!=digits.end()) && (delta<=2));
    
    assert(lastchan>=firstchan);
    double wire = (firstchan+lastchan)/2.;
    double pos = Wire2Pos(wire);
    double res = Resolution(wire);

    // this is not optimal : cluster data is data of the first digit...
    fClusters.insert(CCluster1(GetID(),pos,ndigits,res,data));    
  }
  return fClusters;
}







