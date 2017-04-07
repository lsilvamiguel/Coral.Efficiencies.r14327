// $Id: TDisplayTrackInfo.cc 13148 2011-12-28 16:55:25Z kbicker $

/*!
  Print information about track with specified Id
*/

#include <iostream>  
#include <iterator>  
#include <stdio.h>
#include "TSetup.h"
#include "TConstants.h"
#include "TDisplay.h"
#include "TEv.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/algorithm>
#else
# include <algorithm>  
#endif

using namespace std;

/*
  Version specific to "lattice" alternative.
  Changes w/ respect to "traffic/TDisplayTrackInfo":
   i) Draw QN-parameterized helix, when
   - helix is available,
   - track is single segment.
  ii) Re-built set of proj. from hits list.
*/

// Object function for use in TTrack containers
class EqId {
private:
  unsigned int id;
public:
  bool operator () (const TTrack& t) { 
    return (id == t.Id);
  }
  EqId(int i):id(i){};
};

void TDisplay::TrackInfo(int id)
{

  // ******************** PRINT TTrack w/ ID=="id" ********************

  TEv &ev = TEv::Ref(); const TSetup &setup = TSetup::Ref();
  int imod = DrOpt[6];

  //            ***** RETRIEVE TTrack w/ ID=="id" *****
  list<TTrack>::const_iterator it;  
  it = find_if(ev.lTrack().begin(),ev.lTrack().end(), EqId(id));
  if (it==ev.lTrack().end()) return; // Nothing found
  const TTrack &t = *it;

  t.Print(0);
  if (imod==1)  // ***** LONG PRINT MODE *****
    t.Print(2);  // Print TTrack's hits list

  //                   ***** HELICES *****
  t.Hfirst.Print("At first point"); t.Hlast. Print("At last  point");
  if (t.NGroups()==1 && (t.IFit&0x8)) {
    printf("QN (#hits=%d,chi2=%.3f)",t.NDics,t.Chi2aux/(t.NDics-5));
    t.Haux.Print();
  }

  //               ***** #HITS and RELATED INFO *****
  printf("Type = 0x%x, VSAT = 0x%x, NHits = %3u",t.Type,t.Scifi,t.NHits);
  if      (t.NDFs>t.NHits) printf(" (NDF=%3u)",t.NDFs);
  else if (t.NDFs<t.NHits) printf(" (NDF=%3u!!)\a\a",t.NDFs);
  if (ev.IsMC()) printf(" (%3u - from the same tr. %3d)",t.NHsame,t.IKsame);
  printf(".  Planes: from %3u to %3u (%3d)\n",
	 t.lPlnRef.front(), t.lPlnRef.back(),
	 (t.lPlnRef.back()-t.lPlnRef.front()+1) );
  if (imod==1) {
    cout<<"Used projections: ";
    const_cast<TTrack&>(t).UpdateProjs();
    copy(t.sProj.begin(), t.sProj.end(), ostream_iterator<int>(cout, "  ")); cout<<endl;
  }

  //             ***** CHI2, %X0, TIME *****
  printf("X/X0 = %.1f",t.RadLenFraction());
  if (t.HasShower) printf(" - HasShower %d\n",t.HasShower);
  else printf("\n");
  int nFree = t.Hfirst.with_mom() ? 5 : 4;
  if (t.NDFs-nFree>0)
    printf("Chi2/(NDF-%d)  = %.3f   ",nFree,t.Chi2tot/(t.NDFs-nFree));
  else
    printf("Chi2/(NDF      = %.3f   ",t.Chi2tot/t.NDFs);
  if (t.SigmaTime>0)
    printf(" - Track time = %8.4f +- %5.2f(err) +- %5.1f(disp) [ns]",
	   t.MeanTime, t.SigmaTime, t.DispTime);
  cout<<endl;
  if (t.Associate>=0) printf("Continued into track ID %d\n",t.Associate);

  int ikin = t.IKine; // ***** CORRESPONDIND MC TRACK *****
  if (ikin>=0) MCTrackInfo(ikin);

}
