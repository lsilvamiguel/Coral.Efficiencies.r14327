/*!

  Activated by setting option ReMode [11] to 3

*/

#include "TEv.h"
#include "CsZone.h"
#include "CsErrLog.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "Coral.h"
#include "TDisplay.h"
#include "TEv.h"
#include "TOpt.h"
#include "TSetup.h"
#include "TTrack.h"
#include "THit.h"
#include "TPlane.h"
#include "TConstants.h"
#include "TAlgo.h"
#include "higz.h"
#include "TWatches.h"
#include "CsHistograms.h"
#include "CA.h"

using namespace std;

extern "C"{
  float prob_(const float & chi2, const int & ndeg);
}

class TSegment;

class TNeighbour{
public:
  TSegment* segment;
  ftype chi2;
};

class TSegment:public TCATrack{
 public:
  int hallo;  // is segment constructed using beam direction (==1) 
              // or direction from target (==0)
  int MinFinalNumberOfHits;
  // == information about base space point == 
  TCAHit *base_hits[3];
  ftype base_x;                 // space coordinates & covariance matrix
  ftype base_R[6];
  ftype base_chi2;
  list<TNeighbour> NeighLeft;
  list<TNeighbour> NeighRight;
  int level;
  int found;
  int ghost;
  TSegment(){}
  int Create(TCAHit *hit1, TCAHit *hit2, int is_hallo);
  void Construct(ftype x1, ftype x2);
  TCAHit *find_nearest_hit(TLayer *layer);
};

int TSegment::Create(TCAHit *hit1, TCAHit *hit2, int is_hallo){

  if(hit1==NULL || hit2==NULL) myexit("constructing space point using NULL hits");

  Hits.clear();
  hallo = is_hallo;

  base_hits[1] = hit1;
  base_hits[2] = hit2;

  add_hit(hit1);
  add_hit(hit2);

  TLayer &layer1 = *(hit1->Layer);
  TLayer &layer2 = *(hit2->Layer);
  base_x = hit1->x;
  count_sp( hit1->x, hit1->u, layer1.Ca, layer1.Sa,
	    hit2->x, hit2->u, layer2.Ca, layer2.Sa,
	    base_R, base_x, is_hallo		      );
  base_chi2 = 0;
  return (layer1.active_zone(base_x, base_R)&&layer2.active_zone(base_x, base_R) );
}

void TSegment::Construct(ftype x1, ftype x2){

  NeighLeft.clear();
  NeighRight.clear();

  ftype sigma2 = sqr(TOpt::CAOptD[6]);

  if(x1>x2) swap(x1,x2);

  vx = x1;
  Vx = x2;

  TCAHit *h1 = base_hits[1];
  TCAHit *h2 = base_hits[2];
  TCAHit *h3 = NULL, *h4 = NULL;
  int p1 = h1->Layer->projection;
  int p2 = h2->Layer->projection;
  int p3 = -1;
  list<TCAHit*>::iterator ih=Hits.begin();
  for(; ih!= Hits.end(); ih++){
    TCAHit *hit = *ih;
    p3 = hit->Layer->projection;
    if(p3!=p1&&p3!=p2){
      h3 = hit;      
      break;
    }
  }
  if(h3!=NULL){
    for(; ih!= Hits.end(); ih++){
      TCAHit *hit = *ih;
      int p = hit->Layer->projection;
      if(p!=p1&&p!=p2&&p!=p3){
	h4 = hit;
	break;
      }
    }
  }
  if(h4!=NULL){
    TLayer &l1 = *(h1->Layer);
    TLayer &l2 = *(h2->Layer);
    TLayer &l3 = *(h3->Layer);
    TLayer &l4 = *(h4->Layer);
    count_segment(h1->x, h1->u, l1.Ca, l1.Sa, sqr(h1->sigma_u),
		  h2->x, h2->u, l2.Ca, l2.Sa, sqr(h2->sigma_u),
		  h3->x, h3->u, l3.Ca, l3.Sa, sqr(h3->sigma_u),
		  h4->x, h4->u, l4.Ca, l4.Sa, sqr(h4->sigma_u),
		  vG, vR, vx, VG, VR, Vx, sigma2
		  );
  }else if(h3!=NULL){
    TLayer &l1 = *(h1->Layer);
    TLayer &l2 = *(h2->Layer);
    TLayer &l3 = *(h3->Layer);
    ftype u1 = h3->u;
    ftype x1 = h3->x;
    ftype u0 = base_R[1]*l3.Ca + base_R[2]*l3.Sa;
    ftype v0 =-base_R[1]*l3.Sa + base_R[2]*l3.Ca;
    ftype x0 = base_x;
    ftype v1 = hallo ?v0 :v0/x0*x1;
    if(v1>l3.range_v)  v1 = l3.range_v;
    if(v1<-l3.range_v) v1 =-l3.range_v;
    count_segment(h1->x, h1->u, l1.Ca, l1.Sa, sqr(h1->sigma_u),
		  h2->x, h2->u, l2.Ca, l2.Sa, sqr(h2->sigma_u),
		  x1,   u1,  l3.Ca,  l3.Sa,  sqr(h3->sigma_u),
		  x1,   v1, -l3.Sa,  l3.Ca,  TOpt::CAOptD[7]*sqr(h3->sigma_u),
		  vG, vR, vx, VG, VR, Vx, sigma2
		  );
  }
  base_chi2/=(NumberOfHits()-2);
}


TCAHit *TSegment::find_nearest_hit(TLayer *layer){
  TCAHit *hit = NULL;
  if(layer==NULL) return hit;
  ftype dx = layer->X0 - base_x;
  ftype y = base_R[1] + dx*base_R[3];
  ftype z = base_R[2] + dx*base_R[4];
  ftype u =  y*layer->Ca+z*layer->Sa;
  ftype v = -y*layer->Sa+z*layer->Ca;
  if(fabs(u)>layer->range_u || fabs(v)>layer->range_v || sqr(u)+sqr(v)<layer->holl2) return hit;
  hit = layer->find_nearest_hit(u);
  if(hit!=NULL){
    ftype cut = fabs(dx)*(hallo ?layer->SEGMENT_CUT_HALLO :layer->SEGMENT_CUT_TARGET);
    ftype delta = fabs(hit->u-u);
    if(delta>cut) hit = NULL;
    else base_chi2+=sqr(delta)/sqr(hit->sigma_u);
  }
  return hit;
}

class TCATrackCandidate:public TCATrack{
 public:
  list<TSegment*> Segments;
  int found;
  int ghost;
  int InitialNumberOfHits;
  int MinFinalNumberOfHits;
  TCAHit *find_nearest_hit(TLayer *layer, ftype coeff);
};


TCAHit *TCATrackCandidate::find_nearest_hit(TLayer *layer, ftype coeff){
  TCAHit *hit = NULL;
  if(layer==NULL) return hit;
  ftype G[6][6], R[6], G1[6][6], R1[6];
  PropagateQK(layer->X0, R, G);
  if(!layer->active_zone(R[1],R[2])) return hit;
  count_rot(G, R, G1,R1, layer->Ca, layer->Sa);
  hit = layer->find_nearest_hit(R1[1]);
  if(hit!=NULL){
    ftype cut = 3.5*(coeff*sqrt(G1[1][1])+hit->sigma_u);
    if(fabs(hit->u-R1[1])>cut) hit = NULL;
  }
  return hit;
}

void find_branch(TNeighbour variant[10], int level, int max, TCATrackCandidate *&best_variant){
  //cout<<level<<", ";
  TSegment &segment = *(variant[level].segment);
  if(level!=1){
    for(list<TNeighbour>::iterator i=segment.NeighLeft.begin(); i!=segment.NeighLeft.end();i++){
      TSegment &s = *((*i).segment);
      if(s.found||s.ghost) continue;
      if(s.level>level-1) myexit("incorrect segment level");
      if(s.level==level-1){
	variant[level-1] = *i;
	find_branch(variant,level-1,max, best_variant);
      }
    }
  }else{
    TCATrackCandidate *t = new TCATrackCandidate;
    t->Copy(*(variant[max].segment));
    t->Segments.push_front(variant[max].segment);
    t->MinFinalNumberOfHits = (variant[max].segment)->MinFinalNumberOfHits;
    for(int i=max-1; i>0; i--){
      TSegment *s = variant[i].segment;
      t->Segments.push_front(s);
      copy(s->Hits.begin(),s->Hits.end(), inserter(t->Hits,t->Hits.begin()));
    }
    t->InitialNumberOfHits = t->NumberOfHits();
    t->Fit();
    //t.FitQK();
    /*{
      TEv& event = TEv::Ref();
     TTrack T; // create new track
      int IGROUP=0;
      T.Type=1<<IGROUP;
      for(list<TCAHit*>::iterator i=t.Hits.begin(); i!= t.Hits.end(); i++){
	TCAHit *hit = *i;
	T.AddHit(event.vHit(hit->origin_index));  // add hit to the track
      }

      mirror(t.vG);
      mirror(t.VG);
      TMtx Gf(5,5), Rf(5), Ge(5,5), Re(5);
      for(int j=1; j<=5; j++){
	Rf(j)=t.vR[j];
	Re(j)=t.VR[j];
	for(int i=j; i<=5; i++){
	  Gf(i,j) = t.vG[i][j];
	  Ge(i,j) = t.VG[i][j];
	}
      }
      for(int i=1; i<=5; i++){
	Gf(i,i) = fabs(Gf(i,i));
	Ge(i,i) = fabs(Ge(i,i));
      }

      double x = t.vx;
      T.Hfirst.Set(x,Rf,Gf);
      x = t.Vx;
      T.Hlast.Set(x,Re,Ge);
      T.Chi2tot = t.chi2*t.NumberOfHits();
      T.QuickKF(1,1);
      T.QuickKF(-1,1);
      t.chi2 = T.Chi2tot/t.NumberOfHits();
    }
    */
    ftype chi2 = t->chi2;
    if(variant[1].segment->Hits.front()->Layer->type==21){
      if(max>1){
	chi2=0;
	TSegment *s1 = variant[1].segment;
	TSegment *s2 = variant[2].segment;
	ftype x1 = s1->base_x;
	ftype x2 = s2->base_x;
	ftype y1 = s1->base_R[1];
	ftype y2 = s2->base_R[1];
	ftype z1 = s1->base_R[2];
	ftype z2 = s2->base_R[2];
	//ftype ay = (y2-y1)/(x2-x1);
	//ftype az = (z2-z1)/(x2-x1);
	ftype ay = 0;
	ftype az = 0;
	y1 = (y1+y2)/2;
	z1 = (z1+z2)/2;
	for(list<TCAHit*>::iterator ih=t->Hits.begin(); ih!= t->Hits.end(); ih++){
	  ftype x0 = (*ih)->Layer->X0;
	  ftype u0 = (*ih)->u;
	  ftype y = y1+ay*(x0-x1);
	  ftype z = z1+az*(x0-x1);
	  ftype u = (*ih)->Layer->Ca*y + (*ih)->Layer->Sa*z;
	  chi2+=sqr(u-u0);
	}
	chi2/=t->NumberOfHits();
      }
    }
    if(max==1) chi2 = variant[1].segment->base_chi2;
    //cout<<" ch="<<chi2<<" ";
    t->chi2 = chi2;
    if(t->chi2<best_variant->chi2){
      swap(best_variant,t);
    }
    delete t;
  }
}



ftype count_chi2(TSegment const& s1, TSegment const& s2){
  ftype chi2 = 0;
  /*
  ftype ay, az, by, bz;
  ftype G[4][4];
  ftype sum  = s1.Vx + s2.vx;
  ftype sum2 = sqr(sum);

  ay = 0.5*(s1.VR[3] + s2.vR[3]); 
  az = 0.5*(s1.VR[4] + s2.vR[4]);
  by = 0.5*(s1.VR[1] + s2.vR[1]) - 0.5*ay*sum;
  bz = 0.5*(s1.VR[2] + s2.vR[2]) - 0.5*az*sum;
  G[2][2] = 0.25*(s1.VG[2][2] + s2.vG[2][2]);
  G[3][3] = 0.25*(s1.VG[3][3] + s2.vG[3][3]);
  G[0][0] = 0.25*(s1.VG[0][0] + s2.vG[0][0]) + 0.25*sum2*G[2][2];
  G[1][1] = 0.25*(s1.VG[1][1] + s2.vG[1][1]) + 0.25*sum2*G[3][3];
  chi2 += (ay - s1.VR[3])*(ay - s1.VR[3])/(G[2][2] + s1.VG[2][2]);
  chi2 += (ay - s2.vR[3])*(ay - s2.vR[3])/(G[2][2] + s2.vG[2][2]);
  chi2 += (az - s1.VR[4])*(az - s1.VR[4])/(G[3][3] + s1.VG[3][3]);
  chi2 += (az - s2.vR[4])*(az - s2.vR[4])/(G[3][3] + s2.vG[3][3]);
  chi2 += (by + ay*s1.Vx - s1.VR[1])*(by + ay*s1.Vx - s1.VR[1])/(G[0][0] + s1.VG[0][0]);
  chi2 += (by + ay*s2.vx - s2.vR[1])*(by + ay*s2.vx - s2.vR[1])/(G[0][0] + s2.vG[0][0]);
  chi2 += (bz + az*s1.Vx - s1.VR[2])*(bz + az*s1.Vx - s1.VR[2])/(G[1][1] + s1.VG[1][1]);
  chi2 += (bz + az*s2.vx - s2.vR[2])*(bz + az*s2.vx - s2.vR[2])/(G[1][1] + s2.vG[1][1]);
  
  TCATrack t;
  t.Copy(s1);
  copy(s2.Hits.begin(),s2.Hits.end(), inserter(t.Hits,t.Hits.end()));
  t.Fit();
  chi2 = t.chi2;
  */
  return chi2;
}

void TEv::PrePattern3(CsZone* zone)
{
  static int nevents_ = 0;
  nevents_++;
  //if(nevents_ == 79) return;
  static TWatches watches;

  const TSetup& setup = TSetup::Ref();
  TEv& event = TEv::Ref();
  TDisplay& display = TDisplay::Ref();

  int IGROUP = 0; // we are interested only in detector group IGROUP

  static int last_event = -1;
  int current_event = int(event.ptrEvt()->getEventNumberInRun());
  if(zone==NULL||current_event==last_event) return;
  if(setup.Group2Zone(IGROUP) != zone) return;            
  if(last_event!=current_event) last_event = current_event;

  watches.Start(1);
 
  TLayer Layers[40];
  int NumberOfLayers = 0;
  TSuperLayer SuperLayers[10];
  int NumberOfSuperLayers = 0;
  vector<TSegment*> Segments;
  list<TCATrackCandidate*> tracks;

  watches.Start(2);
  //cout<<"fill super layers\n";
  Fill_Layers(Layers, NumberOfLayers);
  Fill_Super_Layers(SuperLayers, NumberOfSuperLayers, Layers, NumberOfLayers);
  watches.Stop(2);

  //   ================================================================
  //   =                                                              =
  //   =      == LOOK FOR TRACK SEGMENTS IN SUPERLAYERS ==            =
  //   =                                                              =
  //   ================================================================
  watches.Start(3);
  //cout<<"look for segments\n";

  for(int isl = 1; isl<=NumberOfSuperLayers; isl++){

    TSuperLayer &SL = SuperLayers[isl];
    ftype xv = SL.Layers[1]->X0;
    ftype xV = SL.Layers[SL.NumberOfLayers]->X0;
    SL.f_seg = Segments.size();
    SL.l_seg = Segments.size()-1;
    //if(isl==2||isl==5) continue;
    for(int i=1; i<=SL.NumberOfLayers; i++){
      for(int j=1; j<=SL.Layers[i]->NumberOfHits; j++){
	SL.Layers[i]->Hits[j].num4 = 0;
      }
    }
    
    // construction of Ref-segments from crossection of layers CrossedLayers[1] & [2]   
    
    TLayer *layer1 = SL.Layers[SL.CrossedLayers[1]];
    TLayer *layer2 = SL.Layers[SL.CrossedLayers[2]];

    ftype &S1 = layer1->Sa;
    ftype &C1 = layer1->Ca;
    ftype &S2 = layer2->Sa;
    ftype &C2 = layer2->Ca;
    ftype Det_inv = 1/(C1*S2 - S1*C2);


    for(int i=1; i<=layer1->NumberOfHits; i++){
      for(int j=1; j<=layer2->NumberOfHits; j++){
	TCAHit &h1 = layer1->Hits[i];
	TCAHit &h2 = layer2->Hits[j];
	ftype &u1 = h1.u;
	ftype &u2 = h2.u;
	ftype y = (u1*S2-u2*S1)*Det_inv;
	ftype z =-(u1*C2-u2*C1)*Det_inv;
	ftype r2 = sqr(y)+sqr(z);
	if(r2<layer1->holl2 || r2<layer2->holl2   ) continue;
	ftype v1 = -y*S1+z*C1;
	ftype v2 = -y*S2+z*C2;
	if(fabs(v1)>layer1->range_v || fabs(v2)>layer2->range_v ) continue;
	
	for(int construct_hallo=0; construct_hallo<=1; construct_hallo++){
	  if(isl==2||isl==5) construct_hallo=1;
	  TSegment *segment = new TSegment;
	  int ok = segment->Create(&h1, &h2, construct_hallo);
	  if(!ok){
	    delete segment;
	    continue;
	  }
	  for(int il=1; il<=SL.NumberOfLayers; il++){
	    TLayer *layer = SL.Layers[il];
	    if(layer==layer1 || layer==layer2) continue;
	    TCAHit *hit = segment->find_nearest_hit(layer);
	    if(hit!=NULL) segment->add_hit(hit);
	  }
	  if(segment->NumberOfHits()>= SL.RefSegmentNumberOfHits){
	    //cout<<"\nseg "<<Segments.size()<<":";
	    for(list<TCAHit*>::iterator ih=segment->Hits.begin(); ih!= segment->Hits.end(); ih++){
	      (*ih)->num4++;
	      //cout<<(*ih)->Layer->index<<",";
	    }
	    segment->Construct(xv,xV);
	    segment->MinFinalNumberOfHits = SL.RefSegmentNumberOfHits;
	    Segments.push_back(segment);
	  }else{
	    delete segment;
	  }
	}
      }
    }
    
    // construction of Extra-segments from crossection of layers CrossedLayers[i] & [i+1]   
    // at less 1 hit must be unused by Ref-segments
       
    if(isl!=2&&isl!=5){
      for(int pair=1; pair<=SL.NumberOfCrossedLayers; pair++){

	TLayer *layer1 = SL.Layers[SL.CrossedLayers[2*pair-1]];
	TLayer *layer2 = SL.Layers[SL.CrossedLayers[2*pair]];
        ftype &S1 = layer1->Sa;
	ftype &C1 = layer1->Ca;
	ftype &S2 = layer2->Sa;
	ftype &C2 = layer2->Ca;
	ftype Det_inv = 1/(C1*S2 - S1*C2);
  
	for(int i=1; i<=layer1->NumberOfHits; i++){
	  for(int j=1; j<=layer2->NumberOfHits; j++){
	    TCAHit &h1 = layer1->Hits[i];
	    TCAHit &h2 = layer2->Hits[j];
	    if((h1.num4>0)&&(h2.num4>0)) continue;
	    ftype &u1 = h1.u;
	    ftype &u2 = h2.u;
	    ftype y = (u1*S2-u2*S1)*Det_inv;
	    ftype z =-(u1*C2-u2*C1)*Det_inv;
	    ftype r2 = sqr(y)+sqr(z);
	    if(r2<layer1->holl2 || r2<layer2->holl2   ) continue;
	    ftype v1 = -y*S1+z*C1;
	    ftype v2 = -y*S2+z*C2;
	    if(fabs(v1)>layer1->range_v || fabs(v2)>layer2->range_v ) continue;
	    for(int construct_hallo=0; construct_hallo<=1; construct_hallo++){
	      TSegment *segment = new TSegment;
	      int ok = segment->Create(&(layer1->Hits[i]), &(layer2->Hits[j]), construct_hallo);
	      if(!ok){
		delete segment;
		continue;
	      }
	      for(int il=1; il<=SL.NumberOfLayers; il++){
		TLayer *layer = SL.Layers[il];
		if(layer==layer1 || layer==layer2) continue;
		TCAHit *hit = segment->find_nearest_hit(layer);
		if(hit!=NULL) segment->add_hit(hit);
	      }
	      if(segment->NumberOfHits()>= SL.ExtraSegmentNumberOfHits){ 
		/*
		cout<<"\nseg "<<Segments.size()<<":";
		for(list<TCAHit*>::iterator ih=segment->Hits.begin(); ih!= segment->Hits.end(); ih++){
		  cout<<(*ih)->Layer->index<<",";
		}
		*/
		segment->Construct(xv,xV);
		segment->MinFinalNumberOfHits = SL.ExtraSegmentNumberOfHits;
		if(segment->NumberOfHits()>= SL.RefSegmentNumberOfHits){
		  for(list<TCAHit*>::iterator ih=segment->Hits.begin(); ih!= segment->Hits.end(); ih++){
		    (*ih)->num4++;
		  }
		  construct_hallo = 2;
		}
		Segments.push_back(segment);
	      }else{
		delete segment;		
	      }
	    }
	  }
	}
      }  // pair
    }
   
    SL.l_seg = Segments.size()-1;
    //cout<<"SL "<<isl<<": NSegments="<<Segments.size()<<", f="<<SL.f_seg<<", l="<<SL.l_seg<<"\n";
  }//isl
  watches.Stop(3); 

  // ===========================================================================
  // =                                                                         =
  // =                       == LOOK FOR TRACKS ==                             =
  // =                                                                         =
  // ===========================================================================

  if(TOpt::CAOpt[16]==0){
    for(int i=1; i<=NumberOfLayers; i++){
      for(int j=1; j<=Layers[i].NumberOfHits; j++){
	Layers[i].Hits[j].used = 0;
      }
    }

    watches.Start(8); 
    for(int isl = 1; isl<=NumberOfSuperLayers; isl++){
      TSuperLayer &SL = SuperLayers[isl];
      ftype x1 = SL.Layers[1]->X0;
      ftype x2 = SL.Layers[SL.NumberOfLayers]->X0;
      for(int iseg = SL.f_seg; iseg<=SL.l_seg; iseg++){
	TSegment &segment = *(Segments[iseg]);
	segment.CountPQ();
	for(int k=1; k<=TOpt::CAOpt[14]; k++) segment.Fit(x1,x2);
	segment.found = 0;
	segment.ghost = 0;
	segment.NeighLeft.clear();
	segment.NeighRight.clear();
      }
    }
    watches.Stop(8); 
    //cout<<"Look for neighbours:\n";
    watches.Start(4); 
    
    for(int isl = 1; isl<=NumberOfSuperLayers; isl++){
      if(isl==5) continue;
      TSuperLayer &SLi = SuperLayers[isl];
      for(int iseg = SLi.f_seg; iseg<=SLi.l_seg; iseg++){
	TSegment &segment = *(Segments[iseg]);
	// look for two crossed SuperLayers along beam
	int NumNeighbours = 0;
	int jsl = isl+1;
	int jlast = NumberOfSuperLayers;
	int skip = 1;
	if(isl==2){
	  jsl = 5;
	  jlast = 5;
	  skip = 0;
	}
	for(; jsl<=jlast && NumNeighbours<TOpt::CAOpt[17]; jsl++){
	  if((jsl==2||jsl==5)&&skip) continue;
	  TSuperLayer &SL = SuperLayers[jsl];
	  TLayer &lf = *SL.Layers[1];
	  TLayer &ll = *SL.Layers[SL.NumberOfLayers];
	  ftype R[6], G[6][6];
	  segment.Propagate(lf.X0,R,G);
	  ftype yf = R[1];
	  ftype zf = R[2];
	  ftype cutf = 3.5*TOpt::CAOptD[18]*sqrt(G[1][1]+G[2][2]);
	  segment.Propagate(ll.X0,R,G);
	  ftype yl = R[1];
	  ftype zl = R[2];
	  ftype cutl = 3.5*TOpt::CAOptD[18]*sqrt(G[1][1]+G[2][2]);	
	  // look for segments;
	  if(!lf.active_zone(yf,zf,cutf)&&!ll.active_zone(yl,zl,cutl)) continue; // not neighbouring SL
	  int is_neigh = 0;
	  for(int jseg = SL.f_seg; jseg<=SL.l_seg; jseg++){
	    TSegment &s = *(Segments[jseg]);
	    ftype d1 = sqrt(sqr(s.vR[1]-yf)+sqr(s.vR[2]-zf));
	    ftype d2 = sqrt(sqr(s.VR[1]-yl)+sqr(s.VR[2]-zl));
	    ftype s_cutf = 3.5*sqrt(s.vG[1][1]+s.vG[2][2]);
	    ftype s_cutl = 3.5*sqrt(s.VG[1][1]+s.VG[2][2]);
	    if( d1<cutf+s_cutf /*&& d2<cutl+s_cutl*/){
	      s.Propagate((*SLi.Layers[SLi.NumberOfLayers]).X0,R,G);
	      ftype y = R[1];
	      ftype z = R[2];
	      ftype cut = 3.5*TOpt::CAOptD[18]*sqrt(G[1][1]+G[2][2]);	
	      ftype cut1 = 3.5*sqrt(segment.VG[1][1]+segment.VG[2][2]);	
	      ftype d = sqrt(sqr(segment.VR[1]-y)+sqr(segment.VR[2]-z));
	      if( d<cut+cut1 ){      
		is_neigh = 1;
		TNeighbour right, left;
		right.segment = &s;
		left.segment = &segment;
		left.chi2 = right.chi2 = count_chi2(segment, s);
		segment.NeighRight.push_back(right);
		s.NeighLeft.push_back(left);
	      }
	    }
	  }
	  if(is_neigh || lf.active_zone(yf,zf,-cutf)&&ll.active_zone(yl,zl,-cutl)) NumNeighbours++;
	}
      }
    }
    watches.Stop(4); 
    
    do{
      watches.Start(5); 
      //cout<<"Count levels:\n";
      int num=0;
      int max = 0; // max_level of possible track
      for(int isl = 1; isl<=NumberOfSuperLayers; isl++){
	TSuperLayer &SL = SuperLayers[isl];
	for(int iseg = SL.f_seg; iseg<=SL.l_seg; iseg++){
	  TSegment &segment = *(Segments[iseg]);
	  if(segment.found||segment.ghost) continue;
	  num++;
	  int max_level = 0;
	  for(list<TNeighbour>::iterator i=segment.NeighLeft.begin(); i!=segment.NeighLeft.end(); i++){
	    TSegment &s = *((*i).segment);
	    if(s.found||s.ghost) continue;
	    if(s.level>max_level) max_level = s.level;
	  }
	  segment.level = max_level+1;
	  if(max<segment.level) max = segment.level;  
	}
      }
      //cout <<"free segments: "<<num<<"\n";

      watches.Stop(5); 
      //cout<<"Max = "<<max<<"\n";
      if(max<1) break;
      watches.Start(6);
      list<TCATrackCandidate*> best_variants;
      for(vector<TSegment*>::iterator iseg = Segments.begin(); iseg!=Segments.end(); iseg++){
	TSegment &segment = **iseg;
	if(segment.found||segment.ghost) continue;
	if(segment.level>max) myexit("level>max");
	if(segment.level<max) continue;
	TCATrackCandidate *t = new TCATrackCandidate;
	t->chi2 = 100000000;
	TNeighbour variant[10];
	variant[max].segment = &segment;
	variant[max].chi2 = 0;
	find_branch(variant, max, max, t);
	best_variants.push_back(t);
      }
      do{
	ftype best_chi2 = 100000000;
	list<TCATrackCandidate*>::iterator best_variant = best_variants.end(), variant=best_variants.begin();
	for(; variant!=best_variants.end(); variant++){
	  if((*variant)->chi2 < best_chi2){
	    int ghost = 0;
	    for(list<TSegment*>::iterator is=(*variant)->Segments.begin(); 
		is!=(*variant)->Segments.end(); is++){
	      TSegment *s = *is;
	      if(s->found||s->ghost){
		ghost = 1;
		break;
	      }
	    }
	    if(ghost) (*variant)->chi2 = 100000000;
	    else{
	      best_chi2 = (*variant)->chi2;
	      best_variant = variant;
	    }
	  }
	}
	if(best_variant==best_variants.end()) break;
	for(list<TSegment*>::iterator is=(*best_variant)->Segments.begin(); 
	    is!=(*best_variant)->Segments.end(); is++){
	  (*is)->found = 1;
	}
	(*best_variant)->chi2 = 100000000;
	TCATrackCandidate *t = *best_variant;
	best_variants.erase(best_variant);
	// gather hits
	watches.Start(10); 

	t->FitQK();
	if(TOpt::CAOpt[15]>0){
	  vector<int> layers_to_check;
	  for(int i=1; i<t->Hits.front()->Layer->index; i++) layers_to_check.push_back(i);
	  for(list<TCAHit*>::iterator ih=t->Hits.begin(); ih!= t->Hits.end(); ih++){
	    int lfirst = (*ih)->Layer->index+1;
	    list<TCAHit*>::iterator jh = ih;
	    jh++;
	    if(jh==t->Hits.end()) break;
	    int llast  = (*jh)->Layer->index-1;
	    for(int i=lfirst; i<=llast; i++) layers_to_check.push_back(i);
	  }
	  for(int i=t->Hits.back()->Layer->index+1; i<=NumberOfLayers; i++) layers_to_check.push_back(i);   
	  for(vector<int>::iterator il=layers_to_check.begin(); il!=layers_to_check.end(); il++){
	    ftype c = TOpt::CAOptD[23];
	    if(*il<t->Hits.front()->Layer->index||*il>t->Hits.back()->Layer->index) c = TOpt::CAOptD[24];
	    TCAHit *hit = t->find_nearest_hit(&Layers[*il],c);
	    if(hit!=NULL){
	      t->add_hit(hit);
	    }
	  }	
	  t->FitQK();
	}
	watches.Stop(10);
	list<TCAHit*>::iterator ih;
	for(ih=t->Hits.begin(); ih!= t->Hits.end(); ih++){
	  (*ih)->used = 1;
	}
	list<ftype>::iterator cross;
	for(ih=t->Hits.begin(), cross=t->crosspoints.begin(); 
	    ih!= t->Hits.end()&&cross!=t->crosspoints.end(); ih++, cross++){
	  ftype U = *cross;
	  if((*ih)->mirror.size()!=1) continue;
	  ftype u1 = event.vHit((*ih)->origin_index).U;
	  ftype u2 = event.vHit((*ih)->mirror.front()).U;
	  ftype sigma  = event.vHit((*ih)->origin_index).SigU;
	  if(fabs(u2-U)<fabs(u1 -U)){
	    (*ih)->origin_index = (*ih)->mirror.front();
	  }
	}
	tracks.push_back(t);
	for(vector<TSegment*>::iterator iseg = Segments.begin(); iseg!=Segments.end(); iseg++){
	  TSegment &s = **iseg;
	  if(s.found||s.ghost) continue;
	  for(list<TCAHit*>::iterator i=s.Hits.begin(); i!= s.Hits.end(); i++){
	    if((*i)->used){
	      s.ghost = 1;
	      break;
	    }
	  } 
	}
      }while(1);
      list<TCATrackCandidate*>::iterator best_variant = best_variants.end(), variant=best_variants.begin();
      for(list<TCATrackCandidate*>::iterator i=best_variants.begin(); i!=best_variants.end(); i++){
	delete *i;
      }
      watches.Stop(6); 
    }while(1);
    
  }else{
    for(vector<TSegment*>::iterator iseg = Segments.begin(); iseg!=Segments.end(); iseg++){
      TCATrackCandidate *t = new TCATrackCandidate;
      t->Copy(**iseg);
      tracks.push_back(t);
    }	
  }

  for(vector<TSegment*>::iterator iseg = Segments.begin(); iseg!=Segments.end(); iseg++) 
    delete *iseg;
  Segments.clear();
  for(list<TCATrackCandidate*>::iterator it=tracks.begin(); it!=tracks.end(); it++){
    (*it)->ghost = 0;
    (*it)->found = 1;
  }

  // final fit
  watches.Start(9);  
  for(list<TCATrackCandidate*>::iterator it=tracks.begin(); it!=tracks.end(); it++){
    TCATrackCandidate &t = **it;
    if(!t.found||t.ghost) continue;
    //t.Fit();
    //t.FitQK();
    t.CountPQ();
  }
  watches.Stop(9); 

  // == write found tracks to extern data arrays ==
  watches.Start(7); 
   
  //cout<<"\n\nNumber of incoming tracks: "<< event.listTrack.size(); 

  
  for(list<TCATrackCandidate*>::iterator it=tracks.begin(); it!=tracks.end(); it++){
    TCATrackCandidate &t = **it;
    if(!t.found||t.ghost) continue;
    TTrack T; // create new track
    T.Type=1<<IGROUP;
    for(list<TCAHit*>::iterator i=t.Hits.begin(); i!= t.Hits.end(); i++){
      TCAHit *hit = *i;
      T.AddHit(event.vHit(hit->origin_index));  // add hit to the track
    }

    mirror(t.vG);
    mirror(t.VG);
    TMtx Gf(5,5), Rf(5), Ge(5,5), Re(5);
    for(int j=1; j<=5; j++){
      Rf(j)=t.vR[j];
      Re(j)=t.VR[j];
      for(int i=j; i<=5; i++){
	Gf(i,j) = t.vG[i][j];
	Ge(i,j) = t.VG[i][j];
      }
    }
    for(int i=1; i<=5; i++){
      Gf(i,i) = fabs(Gf(i,i));
      Ge(i,i) = fabs(Ge(i,i));
    }

    double x = t.vx;
    T.Hfirst.Set(x,Rf,Gf);
    x = t.Vx;
    T.Hlast.Set(x,Re,Ge);
    T.Chi2tot = t.chi2*t.NumberOfHits();
    //T.QuickKF(1,1);
    //T.QuickKF(-1,1);
    event.listTrack.push_back(T); // store the track
  }
  //cout<<"\nNumber of outgoing tracks: "<< event.listTrack.size()<<"\n";

  // some final operations with found tracks
  list<TTrack>::iterator it;
  for(it =  event.listTrack.begin(); it != event.listTrack.end(); it++){
    (*it).FindKine(); // in case of MC data, find corresponding Kine track
  }

  for(list<TCATrackCandidate*>::iterator it=tracks.begin(); it!=tracks.end(); it++) delete *it;

  // =========================================

  watches.Stop(7); 
 

  static int nevs=0;
  nevs++;
  watches.Stop(1);
  ftype c = 1000.0/nevs;
  //cout<<endl;
  //cout<<"total time in CA        = "<<setw(10)<<c*watches.SumTime( 1)<<"  msec/ev"<<endl;
  //cout<<"filling super layers    = "<<setw(10)<<c*watches.SumTime( 2)<<"  msec/ev"<<endl;
  //cout<<"construct segments      = "<<setw(10)<<c*watches.SumTime( 3)<<"  msec/ev"<<endl;
  //cout<<"fit segments            = "<<setw(10)<<c*watches.SumTime( 8)<<"  msec/ev"<<endl;
  //cout<<"finding neighbours      = "<<setw(10)<<c*watches.SumTime( 4)<<"  msec/ev"<<endl;
  //cout<<"count levels            = "<<setw(10)<<c*watches.SumTime( 5)<<"  msec/ev"<<endl;
  //cout<<"find best track         = "<<setw(10)<<c*watches.SumTime( 6)<<"  msec/ev"<<endl;
  //cout<<"fit + podbor hitov      = "<<setw(10)<<c*watches.SumTime(10)<<"  msec/ev"<<endl;
  //cout<<"export tracks           = "<<setw(10)<<c*watches.SumTime( 7)<<"  msec/ev"<<endl;
  //cout<<"QK-fit of final tracks  = "<<setw(10)<<c*watches.SumTime( 9)<<"  msec/ev"<<endl;
 
 return;
}

