#include "TEv.h"
#include "CsHistograms.h"
#include "partab.h"
/*
  Simple V0 (Lambda0, K0) search routine,
  as final crosscheck of tracking
*/

void TEv::SearchV0()
{
  static bool first(true);
  bool hist = true;
  static CsHist1D* h[50];

  if(first && hist){
    first=false;
    CsHistograms::SetCurrentPath("/Traffic/SearchV0");
    h[0]  = new CsHist1D("h0","Zmin",  200, -1000, 1000);
    h[1]  = new CsHist1D("h1","CDA",   50, 0, 10);
    h[2]  = new CsHist1D("h2","Mass - K0 mass",       50, -0.5,  0.5);
    h[3]  = new CsHist1D("h3","Mass - Lambda mass",   50, -0.5,  0.5);
    h[4]  = new CsHist1D("h4","V0 X",   50, -25,  25);
    h[5]  = new CsHist1D("h5","V0 Y",   50, -25,  25);
    h[6]  = new CsHist1D("h6","V0 Z",   50,  0,  1000);

    CsHistograms::SetCurrentPath("/");
  }

  std::list<TTrack>::iterator itp,itn;
  for(itp = listTrack.begin(); itp != listTrack.end(); itp++){ // loop over positive TTracks
    if( (*itp).Hfirst(5) <= 0) continue;
    THlx Hp = (*itp).Hfirst;
    for(itn = listTrack.begin(); itn != listTrack.end(); itn++){ // loop over negative TTracks
      if( (*itn).Hfirst(5) >= 0) continue;
      THlx Hn = (*itn).Hfirst;
      if( ! Hp.FindCDA(Hn) ) continue; // didn't found

     // fill histograms
      if(hist) h[0]->Fill(Hp(0));
      if(hist) h[1]->Fill(Hp.Dist(Hn));

      //cuts
      if(Hp(0) < 50 || Hp.Dist(Hn) > 1.0) continue;

      THlx H0; double chi2; 
      Hp.Update(Hn, H0, chi2);
      if(hist) h[4]->Fill(H0(1));
      if(hist) h[5]->Fill(H0(2));
      if(hist) h[6]->Fill(H0(0));

      float px1,px2, py1,py2, pz1,pz2, e1,e2, m1,m2, p1,p2, E_2, P_2, mass;

      float Mpi   = Id2mass[9];
      float Mpr   = Id2mass[14];
      float Mk    = Id2mass[16];
      float Mlamb = Id2mass[18];

      p1  =  fabs(1./Hp(5));
      px1 =  Hp.DirCos(1)*p1;
      py1 =  Hp.DirCos(2)*p1;
      pz1 =  Hp.DirCos(3)*p1;

      p2  =  fabs(1./Hn(5));
      px2 =  Hn.DirCos(1)*p2;
      py2 =  Hn.DirCos(2)*p2;
      pz2 =  Hn.DirCos(3)*p2;

      P_2 = (px1+px2)*(px1+px2) + (py1+py2)*(py1+py2) +  (pz1+pz2)*(pz1+pz2);

      //cout<<"p1 = "<<p1<<" p2 = "<<p2<<" P = "<<sqrt(P_2)<<endl;
      // K0
      e1  = sqrt( p1*p1 + Mpi*Mpi );
      e2  = sqrt( p2*p2 + Mpi*Mpi );
      E_2 = (e1 + e2) * (e1 + e2);


      mass = sqrt(E_2-P_2);
      h[2]->Fill(mass-Mk);
      //cout<<"E = "<<sqrt(E_2)<<" mass = "<<mass<<endl;

      // Lambda0
      e1  = sqrt( p1*p1 + Mpr*Mpr );
      e2  = sqrt( p2*p2 + Mpi*Mpi ) ;
      E_2 = (e1 + e2) * (e1 + e2);

      mass = sqrt(E_2-P_2);
      h[3]->Fill(mass-Mlamb);
      //cout<<"E = "<<sqrt(E_2)<<" mass = "<<mass<<endl;

    }
  }
 

}

