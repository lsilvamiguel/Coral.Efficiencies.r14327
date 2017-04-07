#include "GeomPlaneMumega.h"
#include <queue>
#include "math.h"

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))



//ClassImp(GeomPlaneMumega);

const float GeomPlaneMumega::fF1_TICK = 128.e-9; // in ms

std::set<CCluster1>& GeomPlaneMumega::Clusterize(const std::map<int,CDigit1>& digits) {

  fClusters.clear();

  for (std::map<int,CDigit1>::const_iterator ii = digits.begin(); ii!=digits.end(); ) {

    std::priority_queue<CDigit1,std::vector<CDigit1>,GeomPlaneMumega::lessTot> digitqueue;

    const CDigit1& digit=ii->second;
    //    if (digit.isClusterized()) continue;
    //    digit.setClusterized();
    int lastchan = digit.ch;
    int delta = 0;
    digitqueue.push(digit);
    
    do {
      ii++;
      if (ii!=digits.end()) {
        const CDigit1& htmp = ii->second;
        delta = htmp.ch - lastchan;
        if (delta <= 2) {
          lastchan = htmp.ch;
          digitqueue.push(htmp);
        }
      }
    } while ((ii!=digits.end()) && (delta<=2));

    assert(! digitqueue.empty());
    double numchan = 0, denomchan = 0;
    int chanmin = 0xffff, chanmax = -1;
    float cltot = digitqueue.top().dt[1];
#define uM_WEIGHTED_ClTIME
    // uM_WEIGHTED_ClTIME's found to improve time resolution
    // by 5.4% on 11968 (i.e. low current, SM1 off), eb 3.
    // Delay is 6.5 F1 units.
#ifdef uM_WEIGHTED_ClTIME
    double numtime = 0, denomtime = 0;
#else
    float cltime = digitqueue.top().time;
#endif
    float clamp = 0;
    int nchans = digitqueue.size();
    for (register int jj = 0; jj < 5; jj++) {
      if (digitqueue.empty()) break;
      const CDigit1& digitq=digitqueue.top();
#define uM_EXPONENTIATED_ToT 
#ifdef uM_EXPONENTIATED_ToT
      float norm = 60e-6/fF1_TICK;
      float weight = exp((float)digitq.dt[1]/norm);
      clamp += weight;
#else
      float weight = digitq.tot;
      clamp += digitq.tot;
#endif
      chanmin = min(chanmin,digitq.ch);
      chanmax = max(chanmax,digitq.ch);
      numchan += weight*digitq.ch;
      denomchan += weight;
#ifdef uM_WEIGHTED_ClTIME
      if (jj<2) {
	numtime += weight*digitq.dt[0];
	denomtime += weight;
      }
#endif
      digitqueue.pop();
    }

    double pos = -1;
    assert (denomchan != 0);
    double posch = numchan/denomchan;
    pos = Wire2Pos(posch);
    double res = Resolution(posch);
#ifdef uM_WEIGHTED_ClTIME
    float cltime = numtime/denomtime;
#endif
    
    //int nchans = abs(chanmax - chanmin) + 1;

    std::vector<double> data; 
    data.push_back(cltime);
    data.push_back(cltot);
    data.push_back(clamp);
    fClusters.insert(CCluster1(GetID(),pos,nchans,res,data));
  }
  return fClusters;
}

double GeomPlaneMumega::Resolution(double wire) const {
  return Pitch(wire)/sqrt12;
}

