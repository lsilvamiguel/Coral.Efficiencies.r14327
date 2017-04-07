// $Id: CsBmFiHod.cc,v 1.16 2010/01/28 12:51:24 tnagel Exp $

/*!
   \file CsBmFiHod.cc
   \brief Compass Beam SciFb reconstruction Class.
   \author  G. Khaustov
   \version $Revision: 1.16 $
   \date    $Date: 2010/01/28 12:51:24 $
*/
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdio>

#include "CsTypes.h"
#include "CsEvent.h"
#include "CsDetector.h"
#include "CsGeom.h"
#include "CsErrLog.h"
#include "CsOpt.h"
#include "CsBmFiHod.h"
#include "CsRegistry.h"
#include "CsInit.h"

using namespace std;

extern "C" float prob_(const float &hisq,const int &nf );
const int CsBmFiHod::nhod=NBMHODS;
const int CsBmFiHod::ntrackmx=BMHODTRKMX;
const string CsBmFiHod::hodnames[NBMHODS]= 
                          { "ScFbHod1X", "ScFbHod1Y", "ScFbHod2X", "ScFbHod2Y"};
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
CsBmFiHod:: CsBmFiHod(void)
{
        int ipl, ll, i, j, n;
        string tag, key, str;
        CsOpt* opt = CsOpt::Instance();
//
        tag="BmFiHod";
        key="Printlev";
        if( opt->getOpt( tag, key, Printlev ) ) ; else Printlev=0;
//
//       Reconstruction Procedure parameters
//
        for(i=0; i<10; i++) selparam[i]=0.;
        key="MinTmWin";
        if( opt->getOpt( tag, key, selparam[0]) ) ; else selparam[0]=-130.;
        key="MaxTmWin";
        if( opt->getOpt( tag, key, selparam[1]) ) ; else selparam[1]=130.;
//
        key="MxTmDiff";
        if( opt->getOpt( tag, key, selparam[2]) ) ; else selparam[2]=2.6;
//
        key="MxTmChiq";
        if( opt->getOpt( tag, key, selparam[3]) ) ; else selparam[3]=23.;
//
        key="TrckDiff";
        if( opt->getOpt( tag, key, selparam[4]) ) ; else selparam[4]=20.;
//
        key="Mn2htWin";
        if( opt->getOpt( tag, key, selparam[6]) ) ; else selparam[6]=-5.2;
        key="Mx2htWin";
        if( opt->getOpt( tag, key, selparam[7]) ) ; else selparam[7]=2.60;
//
        key="MnTrgWin";
        if( opt->getOpt( tag, key, selparam[8]) ) ; else selparam[8]=-2.;
        key="MxTrgWin";
        if( opt->getOpt( tag, key, selparam[9]) ) ; else selparam[9]=2.;
//
        key="MxNmbHit";
        if( opt->getOpt( tag, key, MxNmbHit) ) ; else MxNmbHit=INPBUFSIZE;
        if(MxNmbHit>INPBUFSIZE) MxNmbHit = INPBUFSIZE;
//
        if(Printlev>0) {
	cout << "CsBmFiHod:: Time window   " << selparam[0] << " - " 
                                             << selparam[1] <<endl;
	cout << "CsBmFiHod:: Max hit time deviation from the mean track time   " 
             << selparam[2] << endl;
	cout << "CsBmFiHod:: Max time Chiq   " << selparam[3] << endl;
	cout << "CsBmFiHod:: Min. TDC deadtime   " << selparam[4] << endl;
	cout << "CsBmFiHod:: Double hit time window   " << selparam[6] << " - "
                                                        << selparam[7] << endl;
	cout << "CsBmFiHod:: Trigger time window   " << selparam[8] << " - "
                                                        << selparam[9] << endl;
	cout << "CsBmFiHod:: Max number of hits/plane   " << MxNmbHit << endl;
        }
	if(CsOpt::Instance()->getOpt( "BeamRecons", "useTRAFFIC", useTRAFFIC));
	else useTRAFFIC=0;
//
	
	if(useTRAFFIC==0) bookhist();       
	CsRegistry reg;
	reg.SORRegistration(this);
//
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBmFiHod:: bhodrecons(void)
{
     int onehitlb, ipl, i, j, n;
//
//    decode and select data
//
      ntrack=0;
      if(Printlev>1)  rawevprint(cout);
      hist1[100]->Fill(0.);
      int declbl=decode();
      if(declbl==10) return;    // reject events with counts > 50000                  
      hist1[100]->Fill(1.);
      if(declbl==0) return;     // reject very big events
      hist1[100]->Fill(2.);
//
//     Look at Hodoscope efficiecy in trigger time window
//
      double t0,t1,t2,t3;
      bool w0,w1,w2,w3; 
      int n0,n1,n2,n3; 
      double trigwin=(selparam[9]-selparam[8])/2.+0.5;
      bool h0 = selonehit(0, trigwin, w0, n0, t0);
      bool h1 = selonehit(1, trigwin, w1, n1, t1);
      if(h0&&h1&&fabs(t0-t1)<3.)  { bool goodST1=(n0>29)&&(n0<45)&&(n1>10)&&(n1<40);
                    eff(t0, t1, 2, 120, goodST1);
                    eff(t0, t1, 3, 130, goodST1);
      }
      bool h2 = selonehit(2, trigwin, w2, n2, t2);
      bool h3 = selonehit(3, trigwin, w3, n3, t3);
      if(h2&&h3&&fabs(t2-t3)<3.)  { bool goodST2=(n2>25)&&(n2<50)&&(n3>60)&&(n3<80);
                    eff(t2, t3, 0, 100, goodST2);
                    eff(t2, t3, 1, 110, goodST2);
      }
//
      if(nfired<nhod) return;
      hist1[100]->Fill(3.);
//
//      fill some hists
//
      for(ipl=0; ipl<nhod; ipl++) {
          for(i=0; i<khit[ipl]; i++) {
              hist1[13+ipl]->Fill(thit[i][ipl]);
              hist2[17+ipl]->Fill((double)ihit[i][ipl], thit[i][ipl]);
          }
      }
//
//  some correlation histograms for 1 hit/hod events in narrow time  window
//
      if(w0&&w1&&w2&&w3) {
             hist2[21]->Fill(t0,t1);
             hist2[22]->Fill(t2,t3);
             hist2[23]->Fill(t0,t2);
             hist2[24]->Fill(t1,t3);
             if(h2&&h3) hist1[25]->Fill(t0-t1);
             if(h0&&h1) hist1[26]->Fill(t2-t3);
             hist1[27]->Fill(t0-t2);
             hist1[28]->Fill(t1-t3);
      }
//
//   look on crosscurrent and dead time
//
      double dt, dn, cut; 
      for(ipl=0; ipl<nhod; ipl++) {
            if(khit[ipl]<=1) continue;
            cut=3.*tmresol[ipl];
            for(i=0; i<khit[ipl]-1; i++) {
                for(n=i+1; n<khit[ipl]; n++) {
                    dt = thit[i][ipl]-thit[n][ipl];
                    dn = fabs(double(ihit[i][ipl]-ihit[n][ipl]));
                    if(dn==0) hist1[29+ipl]->Fill(dt);
                    if((dn==1)&&(fabs(dt)<cut)) {
                          hist1[33+ipl]->Fill(ihit[i][ipl]);
                          hist1[33+ipl]->Fill(ihit[n][ipl]);
                    }
                }
            }
      }
//
//---------------------------------------------------------------------------
//
      onehitlb=khit[0]*khit[1]*khit[2]*khit[3];
      if(onehitlb==1)    onehit();
      else multyhit();
      if(Printlev>1)  bmhodprnt(cout);       // print reconstructed event
//
//---------------------------------------------------------------------------
//       Fill output hists
//
        hist1[55]->Fill(ntrack);
        if(!ntrack) { 
                       hist1[56]->Fill(ntrack);
                       return;
        }
//
        int nhit, nbeamtrk=0, ibeamtrk=0;
        double x[4], pr;
        for(n=0; n<ntrack; n++) { 
           pr= prob_(HDtchiq[n],3);
           hist1[57]->Fill(HDtime[n]);
           hist1[58]->Fill(HDtchiq[n]);
           hist1[59]->Fill(pr);
           for(nhit=4, ipl=0; ipl<nhod; ipl++) {
              x[ipl] = (double)HDhit[0][n][ipl];
              if(HDhit[2][n][ipl]==1) { 
                     hist1[90+ipl]->Fill(HDtime[n]-HDhittm[0][n][ipl]);
              }
              else  nhit++;
//
//   extra hit distribution
//
              for(i=0;i<khit[ipl]; i++){
                  if((ihit[i][ipl]>=0)&&(HDhit[2][n][ipl]==1)) {
                          dt=HDtime[n]-thit[i][ipl];
                          dn=ihit[i][ipl]-x[ipl];
                          hist1[41+ipl]->Fill(dt);
                          if((-13<dt)&&(dt<3.)) hist1[45+ipl]->Fill(dn);
                          hist2[37+ipl]->Fill(dt,dn);
                  }
              }
           }
           hist1[60]->Fill((double)nhit);
//
//   tracks in trigger window
//
           if((selparam[8]<HDtime[n])&&(HDtime[n]<selparam[9])) {
                  nbeamtrk++;
                  ibeamtrk=n;
                  for(ipl=0; ipl<nhod; ipl++) hist1[65+ipl]->Fill(x[ipl]);
                  hist2[70]->Fill(x[0], x[2]);
                  hist2[71]->Fill(x[1], x[3]);
                  if(h0&&h1) { 
                       hist1[72]->Fill(HDtchiq[n]);
                       hist1[73]->Fill(pr);
                  }
           }
           else {
                  for(ipl=0; ipl<nhod; ipl++) hist1[61+ipl]->Fill(x[ipl]);
                  hist2[84]->Fill(x[0], x[2]);
                  hist2[85]->Fill(x[1], x[3]);
           }
        }
//
        hist1[56]->Fill(double(nbeamtrk));
        if(nbeamtrk)  hist1[100]->Fill(5.);
//
//      good beam track
//
        if(nbeamtrk==1) {
                   hist1[100]->Fill(6.);
                   if(!HDlabel[ibeamtrk])  hist1[100]->Fill(7.);
        }
//
//   correlation between tracks
//
      if(ntrack>1){ 
         for(n=0; n<ntrack-1; n++)  {
            for(i=n+1; i<ntrack; i++){
                  dt=HDtime[n]-HDtime[i];
                  if(!(HDlabel[n]+HDlabel[i]))  hist1[77]->Fill(dt);
                  hist1[78]->Fill(dt);
                  for( ipl=0; ipl<nhod; ipl++) {
                          dn=HDhit[0][n][ipl]-HDhit[0][i][ipl];
                          hist2[94+ipl]->Fill(dt,dn);
                  }
            } 
         }
      }       
//
//   noise for time separated tracks     
//
        if(ntrack<2) return;
        for(n=0; n<ntrack; n++)  {
            double dist=1000;
            for(i=0; i<ntrack; i++){
              if(i==n) continue;
                 if(fabs(HDtime[i]-HDtime[n])< dist) dist =fabs(HDtime[i]-HDtime[n]);
            }
        if(dist>2.6) hist1[160]->Fill(HDlabel[n]);
        if(dist>7.5) hist1[161]->Fill(HDlabel[n]);
        if(dist>13.) hist1[162]->Fill(HDlabel[n]);
        if(dist>26.) hist1[163]->Fill(HDlabel[n]);
        }

//       
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
bool CsBmFiHod::selonehit(int ipl, double trigwin, bool& onehit, int& nh, double& t)
{
//      double selwin=200.*.13;                  // selection for low intensity
      double selwin=70.*.13;         // selection for high intensity
      onehit=false;
      if(khit[ipl]<1) return false;
      int nwin=0, nwinb=0;
      for(int i=0; i<khit[ipl]; i++) { 
          if(fabs(thit[i][ipl])<selwin) {
             nwin++;
             if(nwin>1) return false;
             t=thit[i][ipl];
             nh=ihit[i][ipl];
             if(fabs(thit[i][ipl])<trigwin)  nwinb++;
          }
      }
      if(nwin==1) onehit=true;;
      if((nwin!=1)||(nwinb!=1)) return false;          
      return true;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBmFiHod:: eff(double t1, double t2, int ipl3, int nh, bool goodhod)
{
   int i, j, n3=0, inmb[INPBUFSIZE];
   double delt, htime[INPBUFSIZE];
//   double beamwin=70.*.13;                  // selection for low intensity
   double beamwin=40.*.13;                    // selection for high intensity
   double tmean=(t1+t2)/2.;
//
   hist1[nh+1]->Fill(0.);
   if(goodhod) hist1[nh+2]->Fill(0.);
   if(khit[ipl3]) {
      for(i=0; i<khit[ipl3];i++) {
            delt=tmean-thit[i][ipl3];
            if(fabs(delt)<beamwin){
                   htime[n3]=thit[i][ipl3];
                   inmb[n3]=ihit[i][ipl3];
                   n3++; 
            }
      }
   }
   if(n3==0) {
            hist1[nh+1]->Fill(1.);
            if(goodhod) hist1[nh+2]->Fill(1.);
  }
   else if(n3==1) {
            hist1[nh+1]->Fill(2.);
            if(goodhod) hist1[nh+2]->Fill(2.);
            hist1[nh+3]->Fill(tmean-htime[0]);
   }
   else if((n3==2)&&(abs(inmb[0]-inmb[1])<2)) {
            hist1[nh+1]->Fill(3.); 
            if(goodhod) hist1[nh+2]->Fill(3.);
            hist1[nh+4]->Fill(tmean-htime[0]);
            hist1[nh+4]->Fill(tmean-htime[1]);
            hist1[nh+5]->Fill(inmb[0]); 
            hist1[nh+5]->Fill(inmb[1]);
   }
   else { 
            hist1[nh+1]->Fill(4.);
            if(goodhod) hist1[nh+2]->Fill(4.);
            for(i=0;i<n3;i++) {
                       hist1[nh+6]->Fill(tmean-htime[i]);
                       hist1[nh+7]->Fill(inmb[i]);
            }
   }
//
}
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBmFiHod::bookhist()
{
// Book histograms
//
  static bool hfirst = true;
  char titl[80], name[40], a, b;
  int nh = 1; 
  if(hfirst){
    hfirst = false;
    string pathname =  "/BeamRecons/BmFiHod";
    CsHistograms::SetCurrentPath(pathname);
//
      for(int i=0; i< nhod; i++){
        if(i<2) a = '1'; else a = '2';
        if((i==0)||(i==2)) b = 'X'; else b = 'Y';
//
        sprintf(titl,"Number of hits, ScFb%c%c, raw events",a,b);
        sprintf(name,"ScFb%4.4i",1+i);
        hist1[1+i] = new CsHist1D(name,titl, 90, 0, 90);
//
        sprintf(name,"ScFb%4.4i",5+i);
        sprintf(titl,"Fired channels, ScFb%c%c, raw events",a,b);
        hist1[5+i] = new CsHist1D(name,titl, 96, 0, 96);
//
        sprintf(name,"ScFb%4.4i",9+i);
        sprintf(titl,"Hit time, ScFb%c%c, raw events",a,b);
        hist1[9+i] = new CsHist1D(name,titl, 200, -130. , 390.);
//
        sprintf(name,"ScFb%4.4i",13+i);
        sprintf(titl,"Hit time, ScFb%c%c, selected events",a,b);
        hist1[13+i] = new CsHist1D(name,titl, 100, -10., 10.);
//
        sprintf(name,"ScFb%4.4i",17+i);
        sprintf(titl,"Hit time vs channel number, ScFb%c%c",a,b);
        hist2[17+i] =  new CsHist2D(name, titl, 96, 0, 96, 40, -13, 13);
//
        sprintf(name,"ScFb%4.4i",29+i);
        sprintf(titl,"ScFb%c%c, dead time",a,b);
        hist1[29+i] = new CsHist1D(name, titl, 150, -150., 0.);
//
        sprintf(name,"ScFb%4.4i",33+i);
        sprintf(titl,"Hit distrib., delt.chan=1, delt.time<tcut, ScFb%c%c, multyhit events",a,b);
        hist1[33+i] =  new CsHist1D(name, titl, 100, 0, 100);
//
        sprintf(name,"ScFb%4.4i",37+i);
        hist2[37+i] = new CsHist2D(name, "Hit chnl. diff. vs hit time diff.", 50, -10, 10,60,-30.,30.);
        sprintf(name,"ScFb%4.4i",41+i);
        hist1[41+i] = new CsHist1D(name, "Track time - ad. hit time", 200, -13, 13);
        sprintf(name,"ScFb%4.4i",45+i);
        hist1[45+i] = new CsHist1D(name, "Track coord - ad. hit coord, -13<dt<3", 60, -30, 30);
//
        sprintf(name,"ScFb%4.4i",61+i);
        sprintf(titl,"Beam profile out of trig. wind.,  SciFb%c%c",a,b);
        hist1[61+i] =  new CsHist1D(name, titl, 96, 0, 96);
//
        sprintf(name,"ScFb%4.4i",65+i);
        sprintf(titl,"Beam profile, ScFb%c%c, trigger window",a,b);
        hist1[65+i] =  new CsHist1D(name, titl, 96, 0, 96);
//
        sprintf(name,"ScFb%4.4i",90+i);
        sprintf(titl,"ScFb%c%c, track hittime-tracktime, 1hit/track ev.",a,b);
        hist1[90+i] =  new CsHist1D(name, titl, 100, -5, 5);
//
        sprintf(name,"ScFb%4.4i",94+i);
        hist2[94+i] = new CsHist2D(name, "Track coord. diff. vs track time diff.", 50, -10, 10,60,-30.,30.);
//
//
        nh=100+10*i;
        sprintf(name,"ScFb%4.4i",nh+1);
        sprintf(titl,"ScFb%c%c, number of hits,one part.sel.",a,b);
        hist1[nh+1] = new CsHist1D(name, titl, 10, 0., 10.);
        sprintf(name,"ScFb%4.4i",nh+2);
        sprintf(titl,"ScFb%c%c, number of hits,1part.sel.,hod.cut",a,b);
        hist1[nh+2] = new CsHist1D(name, titl, 10, 0., 10.);
//
        sprintf(name,"ScFb%4.4i",nh+3);
        sprintf(titl,"ScFb%c%c, hit time, nhit=1, 1part.sel.",a,b);
        hist1[nh+3] = new CsHist1D(name, titl, 100, -10, 10);
//
        sprintf(name,"ScFb%4.4i",nh+4);
        sprintf(titl,"ScFb%c%c, hit time, 2 adjecent hits, 1part.sel.",a,b);
        hist1[nh+4] = new CsHist1D(name, titl, 100, -10., 10.);
        sprintf(name,"ScFb%4.4i",nh+5);
        sprintf(titl,"ScFb%c%c, fired channels,2 adjecent hits, 1part.sel.",a,b);
        hist1[nh+5] = new CsHist1D(name, titl, 96, 0., 96.);
//
        sprintf(name,"ScFb%4.4i",nh+6);
        sprintf(titl,"ScFb%c%c, hit time,too much hits, 1part.sel.",a,b);
        hist1[nh+6] = new CsHist1D(name, titl, 100, -10., 10.);
        sprintf(name,"ScFb%4.4i",nh+7);
        sprintf(titl,"ScFb%c%c, profile,nhit>1, 1part.sel.",a,b);
        hist1[nh+7] = new CsHist1D(name, titl, 96, 0., 96.);
     }
//
//
        nh=21;
        sprintf(name,"ScFb%4.4i",nh++);
        hist2[21] = new CsHist2D(name, "Hittime ScFb1x vs hittime ScFb1y, 1hit ev.sel.",
                                                           50, -5, 5, 50, -5, 5);
        sprintf(name,"ScFb%4.4i",nh++);
        hist2[22] = new CsHist2D(name, "Hittime ScFb2x vs hittime ScFb2y, 1hit ev.sel.",
                                                          50, -5, 5, 50, -5, 5);
        sprintf(name,"ScFb%4.4i",nh++);
        hist2[23] = new CsHist2D(name, "Hittime ScFb1x vs hittime ScFb2x, 1hit ev.sel.",
                                                           50, -5, 5, 50, -5, 5);
        sprintf(name,"ScFb%4.4i",nh++);
        hist2[24] = new CsHist2D(name, "Hittime ScFb1y vs hittime ScFb2y, 1hit ev.sel.",
                                                          50, -5, 5, 50, -5, 5);
//
        sprintf(name,"ScFb%4.4i",nh++);
        hist1[25]  = new CsHist1D(name, "ScFb1x-ScFb1y hit time diff., 1hit events", 50, -10, 10);
        sprintf(name,"ScFb%4.4i",nh++);
        hist1[26]  = new CsHist1D(name, "ScFb2x-ScFb2y hit time diff., 1hit events", 50, -10, 10);
        sprintf(name,"ScFb%4.4i",nh++);
        hist1[27]  = new CsHist1D(name, "ScFb1x-ScFb2x hit time diff., 1hit events", 50, -10, 10);
        sprintf(name,"ScFb%4.4i",nh++);
        hist1[28]  = new CsHist1D(name, "ScFb1y-ScFb2y hit time diff., 1hit events", 50, -10, 10);
//
        sprintf(name,"ScFb%4.4i",55);
        hist1[55] = new CsHist1D(name, "Number of tracks found", 20, 0., 20.);
        sprintf(name,"ScFb%4.4i",56);
        hist1[56] = new CsHist1D(name, "Ntracks in trigger window", 10, 0., 10.);
        sprintf(name,"ScFb%4.4i",57);
        hist1[57] = new CsHist1D(name, "Track time", 100, -10., 10.);
        sprintf(name,"ScFb%4.4i",58);
        hist1[58] = new CsHist1D(name, "Track chiq", 100, 0, 50);
        sprintf(name,"ScFb%4.4i",59);
        hist1[59] = new CsHist1D(name, "Probability", 100, 0., 1.);
        sprintf(name,"ScFb%4.4i",60);
        hist1[60] = new CsHist1D(name, "Number of hits/track", 10, 0, 10);
//
        sprintf(name,"ScFb%4.4i",70);
        hist2[70] = new CsHist2D(name, "ScFb2X vs ScFb1X, 1 part.sel., in tig wind.",
                                                           48, 0, 96, 48, 0, 96);
        sprintf(name,"ScFb%4.4i",71);
        hist2[71] = new CsHist2D(name, "ScFb2Y vs ScFb1Y, 1part.sel., in tig wind.",
                                                           48, 0, 96, 48, 0, 96);
        sprintf(name,"ScFb%4.4i",84);
        hist2[84] = new CsHist2D(name, "ScFb2X vs ScFb1X, out of trig. wind.",
                                                           48, 0, 96, 48, 0, 96);
        sprintf(name,"ScFb%4.4i",85);
        hist2[85] = new CsHist2D(name, "ScFb2Y vs ScFb1Y, out of trig. wind.",
                                                           48, 0, 96, 48, 0, 96);
//
        sprintf(name,"ScFb%4.4i",72);
        hist1[72] = new CsHist1D(name, "Track chiq,1part.sel.", 100, 0, 50);
        sprintf(name,"ScFb%4.4i",73);
        hist1[73] = new CsHist1D(name, "Track time probability,1part.sel.", 100, 0., 1.);
//
        sprintf(name,"ScFb%4.4i",77);
        hist1[77] = new CsHist1D(name, "Track time differencies", 300, -50, 50);
        sprintf(name,"ScFb%4.4i",78);
        hist1[78] = new CsHist1D(name, "Track time differencies,dirty ev.", 100, -10, 10);
        sprintf(name,"ScFb%4.4i",100);
        hist1[100] = new CsHist1D(name, "Efficiency", 10, 0., 10.);
//
        sprintf(name,"ScFb%4.4i",160);
        hist1[160] = new CsHist1D(name, "Number of additional hits, dt>2.6ns", 10, 0., 10.);
        sprintf(name,"ScFb%4.4i",161);
        hist1[161] = new CsHist1D(name, "Number of additional hits, dt>7.5ns", 10, 0., 10.);
        sprintf(name,"ScFb%4.4i",162);
        hist1[162] = new CsHist1D(name, "Number of additional hits, dt>13ns", 10, 0., 10.);
        sprintf(name,"ScFb%4.4i",163);
        hist1[163] = new CsHist1D(name, "Number of additional hits, dt>26ns", 10, 0., 10.);
//
     CsHistograms::SetCurrentPath("/");
  }
} 
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//int CsBmFiHod:: decode(void)
////
////   get, decode, calibrate and select BmFiHod data         
////   (takes raw date from CsDigit)
////
//{
//  int ipl, nhit, nchnl, rawhit[nhod];
//  double time, tm;
//  int label = 1;
////
//  nfired=0;
//  for(ipl=0; ipl<nhod; ipl++) { 
//     khit[ipl]=0; rawhit[ipl]=0; 
////
//     if(!Id[ipl]) continue;
//     list<CsDigit*>digit = Id[ipl]->getMyDigits();
//     if(!digit.size()) continue;  
//     list<CsDigit*>::iterator Idig;
//     for( Idig = digit.begin(); Idig != digit.end(); Idig ++ ) {
//        CsDigit* digWB = *Idig;
////
//	   nchnl = digWB->getAddress()-1;        // channels count   from 0 !!
//    	   time  = digWB->getDatum();
////    	   time  = digWB->getDatum()+25.;
////       
//           if((nchnl>=hodsize[ipl])||(nchnl<0)) {
//                 ostrstream message;
//                 message << " ScFbHod:: Wrong channel number:"<< nchnl
//                         << "   hodoscope: "<< hodnames[ipl] << endl;
//                 CsErrLog::Instance()->mes(elInfo, message.str());                  
//                 continue;
//           }
//           if(time>50000) return (10);
//           tm=(time-BHODt0[ipl][nchnl])*TDCslope;
////
////        some hists for raw data
//// 
//          rawhit[ipl]++;        
//           hist1[5+ipl]->Fill((double)nchnl);
//           hist1[9+ipl]->Fill(tm);
////
////   time cuts   
////
//           if((selparam[0]<tm)&&(tm<selparam[1])) {
//               nhit=khit[ipl];
//               if(nhit>=MxNmbHit) {label=0; continue;}
//               ihit[nhit][ipl]=nchnl;
//               ihitsv[nhit][ipl]=nchnl;
//               thit[nhit][ipl]=tm;
//               thitsv[nhit][ipl]=tm;
//               khit[ipl]++;
//           }
//      }  
//      if(khit[ipl]) nfired++;
//  }
//  ntothits=khit[0]+khit[1]+khit[2]+khit[3];
//  for(ipl=0; ipl<nhod; ipl++) hist1[1+ipl]->Fill((double)rawhit[ipl]);
//  return label;
//}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
int CsBmFiHod:: decode(void)
//
//   get, decode, calibrate and select BmFiHod data         
//   take data from CsCluster
//
{
  int ipl, nhit, nchnl, rawhit[nhod];
  double time, tm;
//
   int label = 1;
   nfired=0;
   for(ipl=0; ipl<nhod; ipl++) { 
     khit[ipl]=0; rawhit[ipl]=0; 
   }
//
   list<CsCluster*> clust = CsEvent::Instance()->getClusters();
   if(clust.empty()) return label;
//
   list<CsCluster*>::iterator nextclust;
   for( nextclust = clust.begin(); nextclust != clust.end(); nextclust++ ) {
      list<CsDetector*> det = (*nextclust)->getDetsList();
      if(det.empty()) continue;
      CsDetector*  idet=*(det.begin());         // It should be only 1 detector type here
          for(ipl=0; ipl<nhod; ipl++){
              if((idet)!=Id[ipl]) continue; 
              list<CsDigit*> digits = (*nextclust)->getDigitsList();
              if( digits.size()!=1 ) { 
   		    ostringstream message;
                    message << " ScFbHod:: Bad cluster - wrong number of digits: "<< digits.size()
                         << "   hodoscope: "<< hodnames[ipl] << endl;
                    CsErrLog::Instance()->mes(elInfo, message.str());                  
                    break;
              }
              nchnl = (*digits.begin())->getAddress();        // channels count   from 0 !!
//              (*nextclust)->getTime(time);  
    	      time  = (*digits.begin())->getDatum();            // just for safety reason
//       
              if((nchnl>=hodsize[ipl])||(nchnl<0)) {
		 ostringstream message;
                 message << " ScFbHod:: Wrong channel number:"<< nchnl
                         << "   hodoscope: "<< hodnames[ipl] << endl;
                 CsErrLog::Instance()->mes(elInfo, message.str());                  
                 break;
              }
              if(time>50000) return (10);                       // temporary kill such events
              tm=(time-BHODt0[ipl][nchnl])*TDCslope;            // time in ns
//
//        some hists for raw data
// 
              rawhit[ipl]++;        
              hist1[5+ipl]->Fill((double)nchnl);
              hist1[9+ipl]->Fill(tm);
//
//   time cuts   
//
              if((selparam[0]<tm)&&(tm<selparam[1])) {
                  nhit=khit[ipl];
                  if(nhit>=MxNmbHit) {label=0; break;}
                  ihit[nhit][ipl]=nchnl;
                  ihitsv[nhit][ipl]=nchnl;
                  thit[nhit][ipl]=tm;
                  thitsv[nhit][ipl]=tm;
                  inpclust[nhit][ipl]=*nextclust;
                  khit[ipl]++;
              }
      }  
  }
//
  ntothits=khit[0]+khit[1]+khit[2]+khit[3];
  for(ipl=0; ipl<nhod; ipl++) {
           hist1[1+ipl]->Fill((double)rawhit[ipl]);
           if(khit[ipl]) nfired++;
  }
  return label;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBmFiHod:: onehit(void)
//
//      four hit track, one hit/plane,
//
{
        double delt1, delt2, delt3, delt4;
        double tmean, tchiq;
//
        tmean=(thit[0][0]+thit[0][1]+thit[0][2]+thit[0][3])/4.;
        delt1=fabs(tmean-thit[0][0]);
        if(delt1>selparam[2]) return;         // time cuts
        delt2=fabs(tmean-thit[0][1]);
        if(delt2>selparam[2]) return;
        delt3=fabs(tmean-thit[0][2]);
        if(delt3>selparam[2]) return;
        delt4=fabs(tmean-thit[0][3]);
        if(delt4>selparam[2]) return;
//
        tchiq=pow(delt1,2.)/disp[0]+pow(delt2,2.)/disp[1]+
              pow(delt3,2.)/disp[2]+pow(delt4,2.)/disp[3];
        if(tchiq>selparam[3]) return;        // time chiq cut
//
//    save track
//
        ntrack++;
             HDtchiq[0]=tchiq;
             HDtime[0]=tmean;
             HDlabel[0]=0;
             for(int ipl=0; ipl<nhod; ipl++){
                 HDhittm[0][0][ipl]=thit[0][ipl];
                 HDhit[0][0][ipl]=ihit[0][ipl];
                 HDhit[2][0][ipl]=1;
                 outclust[0][0][ipl]=inpclust[0][ipl];
                 ihit[0][ipl]=-100;
             }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBmFiHod:: multyhit(void)
//
//      BEAM hodoscope track reconstruction
//
{
        int indx[4], i, j, k, l, ipl, ll, nn;
        int ktrack, seltrk[NCOMBMX];
        int trkhits[NCOMBMX][4], ifl[NCOMBMX];
	double trkchiq[NCOMBMX], trktime[NCOMBMX];
        double tmean, tchiq, chimin;
        double delt1, delt2, delt3, delt4;
//
//      four hit tracks
//
        for(i=0; i<khit[0]; i++) {
            indx[0]=i;
            for(j=0; j<khit[1]; j++){
                indx[1]=j;
                for(k=0; k<khit[2]; k++){
                    indx[2]=k;
                    for(l=0; l<khit[3]; l++) {
                       tmean=(thit[i][0]+thit[j][1]+thit[k][2]+thit[l][3])/4.;
                       delt1=fabs(tmean-thit[i][0]);
                       if(delt1>selparam[2]) continue;         // time cuts
                       delt2=fabs(tmean-thit[j][1]);
                       if(delt2>selparam[2]) continue;
                       delt3=fabs(tmean-thit[k][2]);
                       if(delt3>selparam[2]) continue;
                       delt4=fabs(tmean-thit[l][3]);
                       if(delt4>selparam[2]) continue;
//
                       tchiq=pow(delt1,2.)/disp[0]+pow(delt2,2.)/disp[1]+
                             pow(delt3,2.)/disp[2]+pow(delt4,2.)/disp[3];
                       if(tchiq>selparam[3]) continue;        // time chiq cut
//                       cout << "tchiq " << tchiq << endl;
                       indx[3]=l;
//
//    save candidates
//
                         if(ntrack==NCOMBMX) {ntrack=0; return;}
                         for(ipl=0;ipl<nhod;ipl++) trkhits[ntrack][ipl]=indx[ipl]; 
                         trkchiq[ntrack] = tchiq;
                         trktime[ntrack] = tmean;
                         ifl[ntrack]=1;
                         ntrack++;
                     }
                 }
             }
         }
//          cout << "ntrack  " << ntrack<< endl;
  if(ntrack<1) return;
  if(ntrack==1){                    // save track, return
             HDtchiq[0]=trkchiq[0];
             HDtime[0]=trktime[0];
             HDlabel[0]=0;
             for(ipl=0; ipl<nhod; ipl++){
                 k=trkhits[0][ipl];
                 HDhittm[0][0][ipl]=thit[k][ipl];
                 HDhit[0][0][ipl]=ihit[k][ipl];
                 HDhit[2][0][ipl]=1;
                 outclust[0][0][ipl]=inpclust[k][ipl];
                 ihit[k][ipl]=-100;
             }
  } else {
//
//   select tracks with best chiq, remove tracks with common hits
//      (common hits are not allowed on this stage)
//     
         int itrack=0;
         ktrack=ntrack;         
         while(ktrack) {
             for(ll=0, chimin=10e7, i=0; i<ntrack; i++) {
               if((ifl[i]==1)&&trkchiq[i]<chimin) 
                             {chimin = trkchiq[i]; ll = i;}  //track with best chiq
             }
             ifl[ll]=2; seltrk[itrack++]=ll;
             if(ktrack==1) break;
             for(i=0; i<ntrack; i++){                  // remove tracks with common points
                if((ifl[i]!=1)||
                   (fabs(trktime[i]-trktime[ll])>selparam[4])) continue;
                    for(ipl=0; ipl<nhod; ipl++) {
                       if(trkhits[ll][ipl]==trkhits[i][ipl]) 
                          {ifl[i]=0; --ktrack; break;}
                   }
             }
             --ktrack;
        }
//
//    save tracks found and mark used hits
//
       ntrack=itrack;
       if(ntrack>BMHODTRKMX){ntrack=0; return;}
       for(i=0; i<ntrack; i++){
             k=seltrk[i];
             HDtchiq[i]=trkchiq[k];
             HDtime[i]=trktime[k];
             HDlabel[i]=0;
             for(ipl=0; ipl<nhod; ipl++){
                 l=trkhits[k][ipl];
                 HDhittm[0][i][ipl]=thit[l][ipl];
                 HDhit[0][i][ipl]=ihit[l][ipl];
                 HDhit[2][i][ipl]=1;
                 outclust[0][i][ipl]=inpclust[l][ipl];
                 ihit[l][ipl]=-100;
             }
       }
   }
//         
//   find more hits belonging to the selected tracks 
//
           int trackhit, newhit, morehit;
           double hittime; 
           for(i=0; i<ntrack; i++) {
                for(ipl=0; ipl<nhod; ipl++) {
                    trackhit=HDhit[0][i][ipl];
                    morehit=-10;
                    for(j=0; j<khit[ipl]; j++) {
                        newhit=ihit[j][ipl];
                        if((newhit<0)||(abs(newhit-trackhit)!=1)) continue;
                        hittime=thit[j][ipl];
                        delt1=HDtime[i]-hittime;
                        if((delt1<selparam[6])||(delt1>selparam[7])) continue;
                        delt1=fabs(delt1);  
                        if(morehit!=-10) {                           // is new candidate better?
                           if(fabs(HDtime[i]-thit[morehit][ipl])<delt1) continue;
                        }
//
                        for(l=0; l<ntrack; l++) {
                            if(l==i) continue; 
                            if(abs(newhit-HDhit[0][l][ipl])!=1)continue;
                            delt2=HDtime[l]-hittime;
                            if((delt2<selparam[6])||(delt2>selparam[7])) continue;
                            if(fabs(delt2)<delt1) goto nexthit;
		        }
                        morehit=j;                                // OK, there is a candidate 
       nexthit:         continue;
                    }
                    if(morehit!=-10) {                               //add new hit to list
                        HDhit[1][i][ipl]=ihit[morehit][ipl];
                        HDhit[2][i][ipl]=2;
                        HDhittm[1][i][ipl]=thit[morehit][ipl];
                        outclust[1][i][ipl]=inpclust[morehit][ipl];
                        ihit[morehit][ipl]=-100;
                    }
               }
          }
//         
//   if there are more free hits in the track time window? 
//
           double tracktime; 
           for(i=0; i<ntrack; i++) {
                HDlabel[i]=0;    
                tracktime=HDtime[i];
                for(ipl=0; ipl<nhod; ipl++) {
                    for(j=0; j<khit[ipl]; j++) {
                        k=ihitsv[j][ipl];
                        if(k==HDhit[0][i][ipl]) continue;
                        if((HDhit[2][i][ipl]==2)&&(k==HDhit[1][i][ipl])) continue;
                        delt1=tracktime-thitsv[j][ipl];
                        if(fabs(delt1)>3.*tmresol[ipl]) continue;
//
                        if((ihit[j][ipl]>=0)&&(abs(k-HDhit[0][i][ipl]) ==1 )) continue;
                        HDlabel[i]++;     // very bad
                    }
               }
          }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBmFiHod:: getBmFiHodtrack(const int n, double &x, double & y, double & z,
            double& dxdz,   double& dydz, int& nhits, int& tothits, bool& trig,
            double &mntime, double &tmchiq) const
{
   double xx[NBMHODS], z0[NBMHODS];
   int i;
//
          if((n<0)||(n>=ntrack)) return;
//          for(i=0; i<nhod; i++){
//              if(HDhit[2][n][i]==1) temp = HDhit[0][n][i];
//              else {
//                   temp = (HDhit[0][n][i]+HDhit[1][n][i])/2.;
//              }
//              xx[i]=(temp+1.)*Id[i]->getWirP() + Id[i]->getWirD();    //   x = nhit*pitch+x0
//          }
          for(i=0; i<nhod; i++){
              z0[i]=(outclust[0][n][i])->getW();
              if(HDhit[2][n][i]==1) xx[i] = (outclust[0][n][i])->getU();
              else        xx[i] = ((outclust[0][n][i])->getU()+(outclust[1][n][i])->getU())/2.;
          }
          z=z0[2];
          x=xx[2];
          dxdz=(x-xx[0])/(z-z0[0]);
          dydz=(xx[3]-xx[1])/(z0[3]-z0[1]);
          y=xx[3]+ dydz*(z-z0[3]);
          mntime=HDtime[n];
          tmchiq=HDtchiq[n];
          nhits= HDlabel[n];
          tothits = ntothits;
          trig=false;           
          if((selparam[8]<mntime)&&(mntime<selparam[9])) trig=true;
          return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBmFiHod:: getBmHodmoments(const int n, double& m0, double& m1) const
{
   int i;
//
    m0=0.;
    m1=0.;
          if((n<0)||(n>=ntrack)) return;
          for(i=0; i<nhod; i++){
              m1+= HDhittm[0][n][i]/tmresol[i];
              m0+= 1./tmresol[i];
          }
          return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
double CsBmFiHod:: getNewChi2(const int n, const double mntime) const
{
int ipl; 
double chiq;
//
      if((n<0)||(n>=ntrack)) return 0.;
      for(chiq=0., ipl=0; ipl<nhod; ipl++){
            chiq+=pow((mntime-HDhittm[0][n][ipl]),2.)/disp[ipl];
      }
     return chiq;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBmFiHod:: getBmErrors(double& dx, double& dy, double& ddxdz, double& ddydz) const
{
      double a=sqrt(12.);
      double pitch0=Id[0]->getWirP();
      double pitch2=Id[2]->getWirP();
      dx=pitch2/a;
      double dz=Id[2]->getZcm()-Id[0]->getZcm();
      ddxdz=sqrt(pitch0*pitch0+pitch2*pitch2)/dz/a;
      double pitch1=Id[1]->getWirP();
      double pitch3=Id[3]->getWirP();
      dy=pitch3/a;
      dz=Id[3]->getZcm()-Id[1]->getZcm();
      ddydz=sqrt(pitch1*pitch1+pitch3*pitch3)/dz/a;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
list<CsCluster*> CsBmFiHod::  getClusters(const int n) const
{
      list<CsCluster*> myclusters;
      if((n<0)||(n>=ntrack)) return myclusters;
      for(int ipl=0; ipl<nhod; ipl++){
            myclusters.push_back(outclust[0][n][ipl]);
            if(HDhit[2][n][ipl]==2)
                     myclusters.push_back(outclust[1][n][ipl]);
      }
     return myclusters;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool  CsBmFiHod:: readcalib(void)
//
//       calibration reading
//       provisional version, should be replaced, when CDB will be OK.      
//
{
      int i, ipl;
      for(ipl=0; ipl<nhod; ipl++) {
             for(i=0; i<hodsize[ipl]; i++) {
                                     BHODt0[ipl][i]=0.;
             }
      }
      TDCslope = 1.;
      if(CsInit::Instance()->IsAMonteCarloJob())  return true; 
//
//     Temporary. Should be replaced in a future by DB reading
//
      TDCslope = 0.1297;
      CsOpt* opt = CsOpt::Instance();
      int ll, j;
      string filename, answer;
      if( !(opt->getOpt( "BmFiHod", "ReadCDB", answer ))||(answer!="YES")) {
          cout << "CsBmFiHod:: T0 calibration will read from files \n";
          for(ipl=0; ipl<nhod; ipl++) {
             if( !(opt->getOpt(hodnames[ipl], "T0calib", filename))) {
                   CsErrLog::Instance()->mes(elError, 
                                          "Error to find T0 calibration for " + hodnames[ipl]);
             }
             else   readT0 ( hodsize[ipl], filename, &BHODt0[ipl][0]); 
          }    
      } 
      else {
//
//       T0 DB read
//
           CsTime runTime=CsEvent::Instance()->getEventTime();
//           CsTime runTime(2000, 1, 5, 0, 0, 0);           // For test only
//           static CsTOFT0DbReader* pTOF =0; 
//           if(pTOF==0) {  cout << "READ"<< endl;pTOF = new CsTOFT0DbReader("scf");}
//           CsTOFT0constants tConsts; 
//           for(ipl=0; ipl<nhod; ipl++){
//               if(pTOF->readplane(hodnames[ipl], tConsts, runTime)&&
//                     (hodsize[ipl]==tConsts.getLength(hodnames[ipl]))){
//                     double* pntr=tConsts.getCoeff(hodnames[ipl]);
//                     for(i=0; i<hodsize[ipl]; i++) BHODt0[ipl][i]=*pntr++;
//               } else  CsErrLog::Instance()->mes(elError, 
//                             "Beam SciFb t0 calibration read error");
//           }
//           pTOF->commit();
       }
//
      if(Printlev>0){ 
         cout <<endl << "Beam SciFb TDC t0s" << endl;
         prntclb(cout, BHODt0);
      }
        return true;
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
void CsBmFiHod::prntclb(ostream &out, double calib[][96]) const
{
      int ipl=0, j;
      out << resetiosflags(ios::scientific);
      out << setiosflags(ios::fixed|ios::showpoint);
      out << setprecision(2); 
      out << setfill(' ');
        for(j=0; j<hodsize[0]; j++ ) {
            out << setw(3)<< j;
            for(ipl=0; ipl<nhod; ipl++) out << setw(12) << calib[ipl][j];
            out << endl;
        }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBmFiHod:: rawevprint(ostream &out) const
//
//   print raw BMS event         
//
{
  int ipl, j, nchnl;
  double time;
// 
  out << resetiosflags(ios::scientific);
  out << setiosflags(ios::fixed|ios::showpoint);
  out << setfill(' ');
  out << setprecision(1) ;
  out << "===================== raw  BEAM ScFB event ====================" << endl;
//
  list<CsDigit*>::iterator Idig;
  for(ipl=0; ipl<nhod; ipl++) {
     if(!Id[ipl]) continue;
     list<CsDigit*>digit = Id[ipl]->getMyDigits();
     out << hodnames[ipl] << ":  total number of hits: " << digit.size() 
           << "   Hits(channel/time):"<< endl;
     if(!digit.size()) continue;
     j=0;  
     for( Idig = digit.begin(); Idig !=digit.end(); Idig ++ ) {
          CsDigit* digWB = *Idig;
	  nchnl = digWB->getAddress()-1;
    	  time  = digWB->getDatum();
          time=(time-BHODt0[ipl][nchnl])*TDCslope;
          out << setw(10) << nchnl 
               << setw(8) << time;
          j++; 
          if(j==5) { j=0; out << endl;}
     }
     if(j!=0) out << endl;
  } 
  out << "--------------------- end raw  BEAM ScFB event ----------------" << endl;
  return;
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool   CsBmFiHod:: init(void)
//
//       start of job/rub initilization
//       Retrive IDs and Hod. parameters from CsDetector
//
{
  bool initOK=true; 
//
        int i, unit;
        string name;
        CsOpt* opt = CsOpt::Instance();
        list <CsDetector*>  det = CsGeom::Instance()->getDetectors();
        list<CsDetector*>::iterator idet;
        list<CsZone*> detzone; 
        list<CsZone*>::iterator izone,jzone;
        myzones.clear();
//
        for(i=0; i<nhod; i++){
             Id[i]=0;  hodsize[i]=0;
             if( opt->getOpt( hodnames[i], "Unit", unit) &&
                 opt->getOpt( hodnames[i], "Name", name)) {
                   for( idet = det.begin(); idet != det.end(); idet++ ) {
                   CsDetector* detpnt = *idet;
  	              if( detpnt ) {
                            if((detpnt->getName()==name)&&(detpnt->getUnit()==unit)) {
                               Id[i]=detpnt;
                               hodsize[i] = detpnt->getNWir();
//
                               tmresol[i]  = detpnt->getTDCResol();         // set resolution
                               if(!CsInit::Instance()->IsAMonteCarloJob()) { 
                                    tmresol[0]=3.75*0.13;
                                    tmresol[1]=3.75*0.13;
                                    tmresol[2]=3.6*0.13;
                                    tmresol[3]=3.5*0.13;
                               }
                               disp[i]=tmresol[i]*tmresol[i];
//
                               detzone=detpnt->getMyZones();                // make list of beam ScFb zones
                               if(!detzone.empty()){
                                    for( izone = detzone.begin(); izone != detzone.end(); izone++ ) {
                                            CsZone* newzone = *izone;
                                            for( jzone = myzones.begin(); jzone != myzones.end(); jzone++ ) {
                                                  if(newzone==*jzone) goto nextzone;
                                            }
                                            myzones.push_back(*izone);
                      nextzone:             continue;
                                    }
                                }
                            }
                      }
                   }  
             }
             if(!Id[i]) initOK=false;
       }
       return initOK;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsBmFiHod:: sor()
{
  if(useTRAFFIC==0){
    if(!init()) {
      CsErrLog::Instance()->mes(elFatal, 
				"Beam ScIFb initilization Error -> check Coral option file");
    }
    if(!readcalib()) 
      CsErrLog::Instance()->mes(elError,"Beam ScFb calibration reading error");
  }
  return true;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CsBmFiHod:: ~CsBmFiHod(){ }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void CsBmFiHod:: bmhodprnt(ostream &out) const
{
      int ipl,i;
//
        out << "CSBmFiHod:  " << ntrack << " BmFiHod track(s) found" << endl;
        if(!ntrack) return;        
        cout << "                           H1X                     H1Y" <<
                      "                     H2X                     H2Y"<< endl;
//
        out << "## lbl    Time TmChiq nhit/hit(s)/hittimes   nhit/hit(s)/hittimes   " <<
                "nhit/hit(s)/hittimes   nhit/hit(s)/hittimes"<< endl;
	out << resetiosflags(ios::scientific);
        out << setiosflags(ios::fixed|ios::showpoint);
//
        for(i=0;i<ntrack; i++) {
          out << setprecision(2) << setfill(' ');
          out <<setw(2)<< i << setw(4)<< HDlabel[i]<<setw(8)<< HDtime[i]
               <<setw(7)<< HDtchiq[i];
            for(ipl=0; ipl<nhod; ipl++) {
                out <<"   "<< HDhit[2][i][ipl]<<" ";
                if(HDhit[2][i][ipl]==1) 
                     out <<setw(3)<< HDhit[0][i][ipl]<<"    ";
                else out <<setw(3)<< HDhit[0][i][ipl]
                          <<setw(3)<< HDhit[1][i][ipl]<<" ";
                out << setprecision(1);
                if(HDhit[2][i][ipl]==1)
                     out <<setw(5)<<HDhittm[0][i][ipl]<<"      ";
                else out <<setw(5)<< HDhittm[0][i][ipl]
                          <<setw(5)<< HDhittm[1][i][ipl]<<" ";
            }
          out << endl;
        }
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool CsBmFiHod:: readT0 ( int nchanread, string filename, double* coeff) 
//
//         Read T0 TOF Coeff from a file
//
{
  int nchanels;
  bool status = CsOpt::expand( filename );
  ifstream in(filename.c_str());
  if(in.fail() || status == false ) {
    CsErrLog::Instance()->mes( elFatal, 
                                " Error during opening of input file " + filename);
  }
//
    char buf[256]; string tmp;
    in.getline(buf,256);
    in >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> nchanels ;
    in.getline(buf,256);
    in.getline(buf,256);
    in.getline(buf,256);
//
    if(nchanels!=nchanread) {
            ostringstream mes;
            mes << " Error during reading " << filename
                << ": exitst "<< nchanels << " values, wanted  " << nchanread<< endl;
            CsErrLog::Instance()->mes( elError,  mes.str()); 
            return false;
    }
    int nch; double temp;
    for(int i=0; i<nchanels; i++){ 
        in >> nch >> temp >> coeff[i] >> temp >> temp >> temp;
        if (!in.good()){ 
            ostringstream mes;
            mes << " Error during reading " << filename
                << "  -> only " << i << " channels are read " << endl;
            CsErrLog::Instance()->mes( elError,  mes.str()); 
            return false;  
        }
    }      
    in.close(); 
    return true;
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
