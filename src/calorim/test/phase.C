#include <cstdio>

#include "../MyPulse.h"

gStyle->SetOptStat(0);

TCanvas *c1 = new TCanvas("c1","c1", 1000, 1000);

double pd = 0.01;

TH1D *phase_cor = new TH1D("phase_cor","", 1/pd, 0, 1);
TH1D *phase_corr = new TH1D("phase_corr","", 1/pd, 0, 1);
TH2D *cfdHMAX = new TH2D("cfdHMAX","",1/pd, 0, 1, 400, -0.04, 0.01);
TH2D *phVSrph = new TH2D("phVSrph","", 1/pd, 0, 1, 1/pd, 0, 1 );

TH1D *phase_diff1 = new TH1D("phase_diff1","", 1/pd, 0, 1);
TH1D *phase_diff2 = new TH1D("phase_diff2","", 1/pd, 0, 1);

void phase() {
  
  gStyle->SetOptStat(0);
  

  c1->Divide(2,2);

  phase_cor->Clear();
  phase_cor->SetTitle("Real Corrections; phase [sample]; correction");
  phase_corr->Clear();
  phase_corr->SetTitle("Corrections; reconstructed phase [sample]; correction");
  phVSrph->Clear();
  phVSrph->SetTitle("Phase VS Reconstructed Phase (CFD); phase [sample]; reconstructed phase [sample]");
  cfdHMAX->Clear();
  cfdHMAX->SetTitle("cfd phase - hmax phase; real phase [sample]; cfd phase - hmax phase [sample]");
  phase_diff1->Clear();
  phase_diff1->SetTitle("( real phase - rec phase ) % 1; real phase");
  phase_diff2->Clear();
  phase_diff2->SetTitle("( real phase - rec phase ) % 1; real phase");

  MyPulse *p = new MyPulse("spline.bin");
  
  TF1 *F = new TF1("F",p,0,32,2,"MyPulse");
  F->SetParameters(0,1);
  F->SetNpx(10000);

  double res_list[32];
  double start = -4.;
  int i = 0;
  int max = 0;

  
  int front = -1;
  int back = 2;

  double isum = 0;
  for( int f = front; f <= back; f++)
    isum += F->Eval(f);

  for( double pos = start; pos < start + 1.; pos += pd) { 
    
    double phase = -1.;
    double phase2 = -1.;
    
    for( i=0; i < 4 - pos; i++) { 
      res_list[i] = F->Eval(pos + i); 
      if( res_list[i] > res_list[max] ) 
        max = i ;
    } 


    double sum = 0;
    for( int f = front; f <= back; f++)
      sum += res_list[max+f];
    
    for( int j = 0; j <= max; j++) { 
      if( res_list[j] >=res_list[max] /2. ) {
        double total = res_list[j] - res_list[j-1];
        double first_half = res_list[max] / 2. - res_list[j-1];
        double second_half =  - res_list[max] / 2. + res_list[j];
        phase = first_half / total; 
        break ; 
      } 
    } 
    {
      double last_result = res_list[1];
      for( int j = 2; j <= i; j++) { 
        double this_result = res_list[j] - 2 * res_list[j-2];
        if ( last_result > 0 && this_result <= 0. ) {
          phase2 = last_result / (last_result -  this_result);
          break;
        }
        last_result = this_result;
      }
    }
    phase_corr->Fill(1 - pos + start ,  1 - res_list[max]);
    phase_cor->Fill(1 - pos + start,  1 - sum/isum);
    phVSrph->Fill(1 - pos + start, phase2);
    cfdHMAX->Fill(1 - pos + start, phase2-phase);
    
    phase_diff1->Fill( 1 - pos + start,   1 - pos + start - phase );
    phase_diff2->Fill( 1 - pos + start,  1 - pos + start - phase2 );
    printf("%f %f %f %f\n", 1 - pos + start, phase, phase2,  1 - res_list[max]); 
  }
  
  c1->cd(1);
  phase_corr->SetLineColor(2);
  phase_cor->Draw();
  phase_corr->Draw("same");
  c1->cd(2);
  phVSrph->Draw("colz");
  c1->cd(3);
  cfdHMAX->Draw("colz");
  c1->cd(4);
  phase_diff2->SetLineColor(2);
  phase_diff1->Draw();
  phase_diff2->Draw("same");
  

};
