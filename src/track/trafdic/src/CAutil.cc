#include <iostream>
#include <math.h>
#include "Coral.h"
#include "TEv.h"
#include "TOpt.h"
#include "TSetup.h"
#include "THit.h"
#include "TPlane.h"
#include "CA.h"

using namespace std;

void myexit(const char *c){
  cout<<"\nCA : "<<c<<"\n";
  assert(false);
}

void fill_zero(ftype G[6][6]){
  for(int i=0; i<6; i++)
    for(int j=0; j<6; j++)
      G[i][j] = 0;    
}

void mirror(ftype G[6][6]){
  for(int i=1; i<=5; i++)
    for(int j=i+1; j<=5; j++)
      G[j][i] = G[i][j];
}

void copy_matrix(ftype G[6][6], ftype G1[6][6]){
  for(int i=1; i<=5; i++)
    for(int j=1; j<=5; j++)
      G1[i][j] = G[i][j];
}

void copy_vector(ftype R[6], ftype R1[6]){
  for(int i=1; i<=5; i++) R1[i] = R[i];
}

void count_rot(const ftype G[6][6], const ftype R[6], ftype G1[6][6], ftype R1[6], ftype C, ftype S){

  ftype C2 = C*C;
  ftype S2 = S*S;
  ftype SC = S*C;

  G1[1][1] =  C2*G[1][1] +    2*SC*G[1][2] + S2*G[2][2];
  G1[1][2] = -SC*G[1][1] + (C2-S2)*G[1][2] + SC*G[2][2];
  G1[2][2] =  S2*G[1][1] -    2*SC*G[1][2] + C2*G[2][2];

  G1[1][3] =  C2*G[1][3] + SC*G[2][3] + SC*G[1][4] + S2*G[2][4];
  G1[2][3] = -SC*G[1][3] + C2*G[2][3] - S2*G[1][4] + SC*G[2][4];
  G1[1][4] = -SC*G[1][3] - S2*G[2][3] + C2*G[1][4] + SC*G[2][4];
  G1[2][4] =  S2*G[1][3] - SC*G[2][3] - SC*G[1][4] + C2*G[2][4];

  G1[3][3] =  C2*G[3][3] +    2*SC*G[3][4] + S2*G[4][4];
  G1[3][4] = -SC*G[3][3] + (C2-S2)*G[3][4] + SC*G[4][4];
  G1[4][4] =  S2*G[3][3] -    2*SC*G[3][4] + C2*G[4][4];

  //mirror(G1);
 
  R1[1] =  C*R[1] + S*R[2];
  R1[2] = -S*R[1] + C*R[2];
  R1[3] =  C*R[3] + S*R[4];
  R1[4] = -S*R[3] + C*R[4];
}

void count_shift(const ftype G[6][6], const ftype R[6], ftype G1[6][6], ftype R1[6], ftype dx){

  R1[1] = R[1] + dx*R[3];
  R1[2] = R[2] + dx*R[4];
  R1[3] = R[3];
  R1[4] = R[4];
  
  G1[1][1] = G[1][1] + dx*(2*G[1][3] +dx*G[3][3]);
  G1[2][2] = G[2][2] + dx*(2*G[2][4] +dx*G[4][4]);
  G1[1][2] = G[1][2] + dx*(G[2][3] + G[1][4] + dx*G[3][4]);
  G1[2][3] = G[2][3] + dx*(G[3][4]);
  G1[1][3] = G[1][3] + dx*(G[3][3]);
  G1[1][4] = G[1][4] + dx*(G[3][4]);
  G1[2][4] = G[2][4] + dx*(G[4][4]);
  G1[3][3] = G[3][3] + dx*(0); //q33
  G1[3][4] = G[3][4] + dx*(0); //q34
  G1[4][4] = G[4][4] + dx*(0); //q44

  //mirror(G1);
}

void count_sp(ftype x1, ftype u1, ftype C1, ftype S1,
	      ftype x2, ftype u2, ftype C2, ftype S2,
	      ftype R[6], ftype x0, int hallo ){

  ftype det = 1/(C1*S2 - S1*C2);
  u1*=det;
  u2*=det;
  if(hallo){
    R[1] = u1*S2-u2*S1;
    R[2] = u2*C1-u1*C2;
    R[3] = 0;
    R[4] = 0;
  }else{
    u1/=x1;
    u2/=x2;
    ftype ty = u1*S2-u2*S1;
    ftype tz = u2*C1-u1*C2;
    R[1] = ty*x0;
    R[2] = tz*x0;
    R[3] = ty;
    R[4] = tz;
  }
  R[5] = 0;
}

void count_segment(ftype x1, ftype u1, ftype C1, ftype S1, ftype sigma1, // == sigma1^2
		   ftype x2, ftype u2, ftype C2, ftype S2, ftype sigma2, // == sigma2^2
		   ftype x3, ftype u3, ftype C3, ftype S3, ftype sigma3, // == sigma3^2
		   ftype x4, ftype u4, ftype C4, ftype S4, ftype sigma4, // == sigma4^2
		   ftype G[6][6], ftype R[6], ftype x0, 
		   ftype G1[6][6], ftype R1[6], ftype x01, 
		   ftype sigma_ms){

  ftype x = x1;
  ftype X = x4;
  ftype Xxi2 = 1/sqr(X-x);
  ftype u = x1!=0? u1*x/x1 :u1;                           //  (x,u1,u2) coordinates
  ftype v = x2!=0? u2*x/x2 :u2;
  sigma1 = x1!=0? sigma1*sqr(x/x1) :sigma1;
  sigma2 = x2!=0? sigma2*sqr(x/x2) :sigma2;
  ftype det = 1/(C1*S2 - S1*C2);
  ftype det2 = sqr(det)*Xxi2;
  ftype y = (u*S2-v*S1)*det;
  ftype z =-(u*C2-v*C1)*det;              // (x,y,z) coordinates
  ftype Gx11,Gx12,Gx22;
  Gx11 = (sigma2*S1*S1+sigma1*S2*S2)*det2;
  Gx12 =-(sigma2*S1*C1+sigma1*S2*C2)*det2;
  Gx22 = (sigma2*C1*C1+sigma1*C2*C2)*det2;

  ftype U = x3!=0 ?u3*X/x3 :u3;                          //  (x,u1,u2) coordinates
  ftype V = x4!=0 ?u4*X/x4 :u4;
  sigma3 = x3!=0? sigma3*sqr(X/x3) :sigma3;
  sigma4 = x4!=0? sigma4*sqr(X/x4) :sigma4;
  ftype Det = 1/(C3*S4 - S3*C4);
  ftype Det2 = sqr(Det)*Xxi2;
  ftype Y = (U*S4-V*S3)*Det;
  ftype Z =-(U*C4-V*C3)*Det;              // (x,y,z) coordinates
  ftype GX11,GX12,GX22;
  GX11 = (sigma4*S3*S3+sigma3*S4*S4)*Det2;
  GX12 =-(sigma4*S3*C3+sigma3*S4*C4)*Det2;
  GX22 = (sigma4*C3*C3+sigma3*C4*C4)*Det2;

  if(x>X){
    swap(x,X);
    swap(y,Y);
    swap(z,Z);
    swap(Gx11,GX11);
    swap(Gx12,GX12);
    swap(Gx22,GX22);
  }
  ftype Xxi = 1/(X-x);
  y*=Xxi;
  z*=Xxi;
  Y*=Xxi;
  Z*=Xxi;
  
  sigma_ms*=2;
  ftype ms = sigma_ms*sqr(X-x);
  ftype x0mX = x0-X;
  ftype x0mX2 = sqr(x0mX);
  ftype x0mx = x0-x;
  ftype x0mx2 = sqr(x0mx);

  R[1] = Y*x0mx - y*x0mX;
  R[2] = Z*x0mx - z*x0mX;
  R[3] = Y-y;
  R[4] = Z-z;
  R[5] = 0;

  G[1][1] = Gx11*x0mX2 + GX11*x0mx2 + ms;
  G[1][2] = Gx12*x0mX2 + GX12*x0mx2 + ms;
  G[2][2] = Gx22*x0mX2 + GX22*x0mx2 + ms;
  G[1][3] = Gx11*x0mX  + GX11*x0mx;
  G[2][4] = Gx22*x0mX  + GX22*x0mx;
  G[1][4] = Gx12*x0mX  + GX12*x0mx;
  G[2][3] = G[1][4];
  G[3][3] = Gx11 + GX11 + sigma_ms;
  G[3][4] = Gx12 + GX12;
  G[4][4] = Gx22 + GX22 + sigma_ms;
  G[5][5] = 1;
  G[1][5] = G[2][5] = G[3][5] = G[4][5] = 0;

  ftype x01mX = x01-X;
  ftype x01mX2 = sqr(x01mX);
  ftype x01mx = x01-x;
  ftype x01mx2 = sqr(x01mx);

  R1[1] = Y*x01mx - y*x01mX;
  R1[2] = Z*x01mx - z*x01mX;
  R1[3] = R[3];
  R1[4] = R[4];
  R1[5] = 0;

  G1[1][1] = Gx11*x01mX2 + GX11*x01mx2 + ms;
  G1[1][2] = Gx12*x01mX2 + GX12*x01mx2 + ms;
  G1[2][2] = Gx22*x01mX2 + GX22*x01mx2 + ms;
  G1[1][3] = Gx11*x01mX  + GX11*x01mx;
  G1[2][4] = Gx22*x01mX  + GX22*x01mx;
  G1[1][4] = Gx12*x01mX  + GX12*x01mx;
  G1[2][3] =  G1[1][4];
  G1[3][3] =  G[3][3];
  G1[3][4] =  G[3][4];
  G1[4][4] =  G[4][4];
  G1[5][5] = 1;
  G1[1][5] = G1[2][5] = G1[3][5] = G1[4][5] = 0;

}


TCAHit *TLayer::find_nearest_hit(ftype u){
  TCAHit *hit=NULL;
  ftype min_dist = 10000.0;
  for(int ih=1; ih<=NumberOfHits; ih++){
    TCAHit &h = Hits[ih];
    if(h.used) continue;
    ftype dist = fabs(u - h.u);
    if(dist<min_dist){
      hit = &h;
      min_dist = dist;
    }
  } 
  return hit;
}

void TCATrack::Copy(const TCATrack& t)
{
  vx = t.vx;
  Vx = t.Vx;
  for(int i=1; i<=5; i++){
    vR[i] = t.vR[i];
    VR[i] = t.VR[i];
    for(int j=1; j<=5; j++){
      vG[i][j] = t.vG[i][j];
      VG[i][j] = t.VG[i][j];
    }
  }
  chi2 = t.chi2;
  pinv = t.pinv;
  Hits.clear();
  copy(t.Hits.begin(),t.Hits.end(), inserter(Hits,Hits.begin()));
};

void Fill_Layers(TLayer Layers[], int &NumberOfLayers){

  const TSetup& setup = TSetup::Ref();

  TEv& event = TEv::Ref();

  int IGROUP = 0; // we are interested only in detector group IGROUP 

  NumberOfLayers=0;

  for(int iplane=setup.vIplFirst()[IGROUP]; iplane<=setup.vIplLast()[IGROUP]; iplane++){ 
    const TPlane  &plane = setup.vPlane(iplane);
    const TDetect &detector = setup.vDetect(plane.IDetRef);   
    int plane_type = detector.IType;   
    int plane_proj = plane.IProj;
    NumberOfLayers++;
    TLayer &L = Layers[NumberOfLayers];
    L.index = NumberOfLayers;
    L.global_plane_index = iplane;
    L.type = plane_type;
    L.projection = plane_proj;
    L.resol = fabs(detector.Resol);
    L.pitch = fabs(detector.Pitch);
    L.range_u = fabs(detector.Range)/2;
    L.range_v = fabs(detector.Range)/2;
    L.size_x = fabs(detector.Siz(0));
    L.X0 = detector.X(0);
    L.Sa = detector.Sa;
    L.Ca = detector.Ca;
    if(plane_type==2){
      L.SEGMENT_CUT_TARGET = TOpt::CAOptD[8]/1000;
      L.SEGMENT_CUT_HALLO = TOpt::CAOptD[9]/1000;
      L.holl = 2.5;
    }else if(plane_type==11){
      L.SEGMENT_CUT_TARGET = TOpt::CAOptD[10]/1000;
      L.SEGMENT_CUT_HALLO = TOpt::CAOptD[11]/1000;
      L.holl = 19.0;
    }else{
      L.SEGMENT_CUT_TARGET = TOpt::CAOptD[12]/1000;
      L.SEGMENT_CUT_HALLO = TOpt::CAOptD[13]/1000;
      L.holl = 0.0;
    }
    L.holl2 = sqr(L.holl);
    L.NumberOfHits = 0;
    if(TOpt::CAOpt[12]==0&&plane_type==21) continue;
    vector<int> used_hits;
    for(vector<int>::const_iterator ih=plane.vHitRef().begin(); ih!=plane.vHitRef().end(); ih++){
      used_hits.push_back(0);
    }
    for(vector<int>::const_iterator ih=plane.vHitRef().begin(); ih!=plane.vHitRef().end(); ih++){
      const THit &HIT = event.vHit(*ih);
      if(HIT.IPlane!=iplane) continue; // skip hits from other planes
      if(TOpt::CAOpt[13]==0&&HIT.IKine<0) continue; // skip mirror hits
      if(TOpt::CAOpt[11]==0&&HIT.IOrig!=0) continue; // skip secondary hits
      if(used_hits[*ih-plane.vHitRef().front()] == 1) continue;
      {            // === simulation of detector unefficiensy ===
	int i=1;
	float eff[10];
	ranlux_(eff,&i);
	if(eff[0]>TOpt::CAOptD[0]) continue;
      }
      L.NumberOfHits++;
      TCAHit &hit = L.Hits[L.NumberOfHits];
      hit.origin_index = *ih;
      hit.Layer = &L;
      hit.x = detector.X(0);
      hit.u = HIT.U;
      hit.sigma_u = HIT.SigU;
      hit.v = HIT.V;
      hit.sigma_v = HIT.SigV;  
      hit.used = 0;
      hit.mirror.clear();
      used_hits[*ih-plane.vHitRef().front()] = 1;
      if(TOpt::CAOpt[13]>0&&detector.Kind==1&&HIT.sDigits().size()==1){
	const TDigit& d = *(HIT.sDigits().begin()); 
	if(d.vDigInfo.empty()) continue;
	int wire_num = d.IWire;

	for(vector<int>::const_iterator jh=plane.vHitRef().begin(); jh!=plane.vHitRef().end(); jh++){
	  const THit &H = event.vHit(*jh);
	  if(H.sDigits().size()!=1||*ih==*jh) continue;
	  if(used_hits[*jh-plane.vHitRef().front()] == 1) continue;
	  const TDigit& d = *(H.sDigits().begin()); 
	  if(d.vDigInfo.empty()) continue;
	  if(wire_num==d.IWire){
	    hit.sigma_u += fabs(H.U-hit.u)/3.5/2;
	    hit.u = (H.U+hit.u)/2;
	    hit.mirror.push_back(*jh);
	    used_hits[*jh-plane.vHitRef().front()] = 1;	    
	    break;
	  }
	}
      }
    }
  }
}

void Fill_Super_Layers(TSuperLayer SuperLayers[], int &NumberOfSuperLayers, TLayer Layers[], int NumberOfLayers){

  if(NumberOfLayers<34) myexit("Incorrect number of layers");
  NumberOfSuperLayers=7;
 
  SuperLayers[1].NumberOfLayers = 4;
  SuperLayers[1].Layers[1] = &Layers[1];
  SuperLayers[1].Layers[2] = &Layers[2];
  SuperLayers[1].Layers[3] = &Layers[3];
  SuperLayers[1].Layers[4] = &Layers[4];
  SuperLayers[1].NumberOfCrossedLayers = 2;
  SuperLayers[1].CrossedLayers[1] = 1;
  SuperLayers[1].CrossedLayers[2] = 2;
  SuperLayers[1].CrossedLayers[3] = 3;
  SuperLayers[1].CrossedLayers[4] = 4;
  SuperLayers[1].RefSegmentNumberOfHits = 4;
  SuperLayers[1].ExtraSegmentNumberOfHits = 3;

  SuperLayers[2].NumberOfLayers = 3;
  SuperLayers[2].Layers[1] = &Layers[5];
  SuperLayers[2].Layers[2] = &Layers[6];
  SuperLayers[2].Layers[3] = &Layers[7];
  SuperLayers[2].NumberOfCrossedLayers = 2;
  SuperLayers[2].CrossedLayers[1] = 1;
  SuperLayers[2].CrossedLayers[2] = 2;
  SuperLayers[2].CrossedLayers[3] = 2;
  SuperLayers[2].CrossedLayers[4] = 3;
  SuperLayers[2].RefSegmentNumberOfHits = 3;
  SuperLayers[2].ExtraSegmentNumberOfHits = 3;

  SuperLayers[3].NumberOfLayers = 8;
  SuperLayers[3].Layers[1] = &Layers[8];
  SuperLayers[3].Layers[2] = &Layers[9];
  SuperLayers[3].Layers[3] = &Layers[10];
  SuperLayers[3].Layers[4] = &Layers[11];
  SuperLayers[3].Layers[5] = &Layers[12];
  SuperLayers[3].Layers[6] = &Layers[13];
  SuperLayers[3].Layers[7] = &Layers[14];
  SuperLayers[3].Layers[8] = &Layers[15];
  SuperLayers[3].NumberOfCrossedLayers = 2;
  SuperLayers[3].CrossedLayers[1] = 1;
  SuperLayers[3].CrossedLayers[2] = 3;
  SuperLayers[3].CrossedLayers[3] = 2;
  SuperLayers[3].CrossedLayers[4] = 4;
  SuperLayers[3].CrossedLayers[5] = 5;
  SuperLayers[3].CrossedLayers[6] = 7;
  SuperLayers[3].RefSegmentNumberOfHits = 6;
  SuperLayers[3].ExtraSegmentNumberOfHits = 5;

  SuperLayers[4].NumberOfLayers = 4;
  SuperLayers[4].Layers[1] = &Layers[16];
  SuperLayers[4].Layers[2] = &Layers[17];
  SuperLayers[4].Layers[3] = &Layers[18];
  SuperLayers[4].Layers[4] = &Layers[19];
  SuperLayers[4].NumberOfCrossedLayers = 2;
  SuperLayers[4].CrossedLayers[1] = 1;
  SuperLayers[4].CrossedLayers[2] = 2;
  SuperLayers[4].CrossedLayers[3] = 3;
  SuperLayers[4].CrossedLayers[4] = 4;
  SuperLayers[4].RefSegmentNumberOfHits = 4;
  SuperLayers[4].ExtraSegmentNumberOfHits = 3;

  SuperLayers[5].NumberOfLayers = 3;
  SuperLayers[5].Layers[1] = &Layers[20];
  SuperLayers[5].Layers[2] = &Layers[21];
  SuperLayers[5].Layers[3] = &Layers[22];
  SuperLayers[5].NumberOfCrossedLayers = 2;
  SuperLayers[5].CrossedLayers[1] = 1;
  SuperLayers[5].CrossedLayers[2] = 2;
  SuperLayers[5].CrossedLayers[3] = 2;
  SuperLayers[5].CrossedLayers[4] = 3;
  SuperLayers[5].RefSegmentNumberOfHits = 3;
  SuperLayers[5].ExtraSegmentNumberOfHits = 3;

  SuperLayers[6].NumberOfLayers = 4;
  SuperLayers[6].Layers[1] = &Layers[23];
  SuperLayers[6].Layers[2] = &Layers[24];
  SuperLayers[6].Layers[3] = &Layers[25];
  SuperLayers[6].Layers[4] = &Layers[26];
  SuperLayers[6].NumberOfCrossedLayers = 2;
  SuperLayers[6].CrossedLayers[1] = 1;
  SuperLayers[6].CrossedLayers[2] = 2;
  SuperLayers[6].CrossedLayers[3] = 3;
  SuperLayers[6].CrossedLayers[4] = 4;
  SuperLayers[6].RefSegmentNumberOfHits = 4;
  SuperLayers[6].ExtraSegmentNumberOfHits = 3;

  SuperLayers[7].NumberOfLayers = 8;
  SuperLayers[7].Layers[1] = &Layers[27];
  SuperLayers[7].Layers[2] = &Layers[28];
  SuperLayers[7].Layers[3] = &Layers[29];
  SuperLayers[7].Layers[4] = &Layers[30];
  SuperLayers[7].Layers[5] = &Layers[31];
  SuperLayers[7].Layers[6] = &Layers[32];
  SuperLayers[7].Layers[7] = &Layers[33];
  SuperLayers[7].Layers[8] = &Layers[34];
  SuperLayers[7].NumberOfCrossedLayers = 2;
  SuperLayers[7].CrossedLayers[1] = 1;
  SuperLayers[7].CrossedLayers[2] = 3;
  SuperLayers[7].CrossedLayers[3] = 2;
  SuperLayers[7].CrossedLayers[4] = 4;
  SuperLayers[7].CrossedLayers[5] = 5;
  SuperLayers[7].CrossedLayers[6] = 7;
  SuperLayers[7].RefSegmentNumberOfHits = 6;
  SuperLayers[7].ExtraSegmentNumberOfHits = 5;

}

void TCATrack::add_hit(TCAHit *hit){
  if(hit==NULL) myexit("trying to add NULL hit to track");
  /*
  list<TCAHit*>::iterator ih=Hits.begin();
  while(ih!=Hits.end() && (*ih)->Layer->index<hit->Layer->index) ih++;
  Hits.insert(ih,hit);
  */

  list<TCAHit*>::iterator ih;
  for(ih=Hits.begin(); ih!= Hits.end(); ih++){
    if((*ih)->Layer->index>=hit->Layer->index) break;
  }
  Hits.insert(ih,hit);

}
