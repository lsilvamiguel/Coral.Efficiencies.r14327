#include "GeomPlaneMuonWallA.h"
#include <queue>
#include "math.h"
#include <cstdlib>

//ClassImp(GeomPlaneMuonWallA);

inline  bool isCluster(const CDigit1& c1, const CDigit1& c2) {
  int ch1 = int(c1.ch);
  int ch2 = int(c2.ch);
  int k1  = int(c1.dt[1]);
  int k2  = int(c1.dt[1]);
  if(abs(ch2 - ch1) < 2 && abs(k2 - k1) < 2) return 1;
  return 0;
}

std::set<CCluster1>& GeomPlaneMuonWallA::Clusterize(const std::map<int,CDigit1>& digits) {
//   static int ccc=0;

  fClusters.clear();
  if(digits.empty()) {
    return fClusters;
  }  
  typedef std::map<int,CDigit1>::const_iterator MIT;
  MIT st = digits.begin();
  float s1 = st->second.ch;
  float s2 = 1.;
  MIT prev = st;
  MIT h=++st;

  while(1) {
    if(h == digits.end()) {
      --h;
      std::vector<double> data;
      data.push_back(h->second.dt[0]); // time of last hit in cluster
      data.push_back(h->second.dt[1]); // same for signature
      double wire = s1/s2;
      double pos = Wire2Pos(wire);
      double res = Resolution(wire);
      fClusters.insert(CCluster1(GetID(),pos,int(s2),res,data));
      break;
    }
    else {
      if(isCluster(prev->second, h->second)) {
	s1 += h->second.ch;
	s2++;
      }
      else {
	std::vector<double> data;
	data.push_back(h->second.dt[0]); // time of last hit in cluster
	data.push_back(h->second.dt[1]); // same for signature
	double wire = s1/s2;
	double pos = Wire2Pos(wire);
	double res = Resolution(wire);
	fClusters.insert(CCluster1(GetID(),pos,int(s2),res,data));    
	s1 = h->second.ch;
	s2 = 1.;
      }
      prev = h;
    }
    h++;
  }
  
  return fClusters;
}

double GeomPlaneMuonWallA::Wire2Pos(double w) const {
  const int ss=8;
  int nmod = int(w) / ss;  
  return double(nmod)*(double(ss)+0.5) + double(int(w) % ss) - GetHalfSize();
}



