#include <cstdio>
#include <unistd.h>

#define MONITOR

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include "TStyle.h"
#include "TSystem.h"

#ifdef MONITOR
#include "../MyPulse.h"
#include "TF1.h"


MyPulse *mypulse[4][3];
TF1 *f[4][3];
 
TCanvas *c1 = NULL;
#endif

using namespace std;


TCanvas *c = NULL;
TH2D *amplcor = NULL;
TH2D *amplVScor = NULL;
TH2D *tVSchi2 = NULL;
TH1D *chi2 = NULL;
TH1D *chi2prob = NULL;
TH2D *amplVSchi2 = NULL;
TH1D *hsample = NULL;
TH2D *cutpos = NULL;

void fit(int monitor = 0, double mon_cut = 500., int id_ = -1, int plane_ = -1,double cut_low = 0., double cut_high = 0.003, char *name="/tmp/out2.bin", int analyse = -1, unsigned int skip = 0) {
 
  gStyle->SetOptStat(0);

#ifdef MONITOR
  static bool first = true;
  if( first ) {
    first = false;
    for( int id = 0; id < 4; ++id) {
      for( int plane = 0; plane < 3; ++plane) {
        mypulse[id][plane] = NULL;
      }
    }
    
    mypulse[0][0] = new MyPulse("EC01P00.par");
    mypulse[0][1] = new MyPulse("EC01P01.par");
    mypulse[0][2] = new MyPulse("EC01P02.par");
    mypulse[1][0] = new MyPulse("EC02.reg1.para");
    mypulse[1][1] = new MyPulse("EC02.reg2.para");
    mypulse[1][2] = new MyPulse("EC02.reg3.para");
    
    for( int id = 0; id < 4; ++id) {
      for( int plane = 0; plane < 3; ++plane) {
        if( mypulse[id][plane] ) {
          char file_name[20];
          sprintf( file_name, "f_%i_%i", id, plane);
          f[id][plane] = new TF1(file_name,mypulse[id][plane],0.,32.,3,"MyPulse");
          f[id][plane]->SetNpx(1000);
          printf( "Registered function 0x%x for %i:%i\n", f[id][plane], id+1, plane);
        }
      }
    }
  }
 
  gStyle->SetOptStat(0);

  if ( monitor && !c1)
    c1 = new TCanvas("c1","c1",1200,1200);
#endif

  if( !c ) {
    c = new TCanvas("c","c",1800,1200);
    c->Divide(2,2);
  }
  if( !amplcor )
    amplcor = new TH2D("amplcor","Amplitude Correction; ampl/sample[max] - 1", 1000, -.4, .4, 64*48+1, -1.5,64*48-.5 );
  else
    amplcor->Reset();
  
  if( !amplVScor )
    amplVScor = new TH2D("amplVScor", "Amplitude Vs Correction; ampl;  ampl/sample[max] - 1", 256, -.5, 4095.5,  1000, -.4, .4);
  else
    amplVScor->Reset();
  
  if( !tVSchi2 )
    tVSchi2 = new TH2D("tVSchi2","time VS chi2; time; chi2/ndf",
                       500, 0, 32*12.8,500, 0, 0.1*100*100);
  else
    tVSchi2->Reset();

  if( !chi2 )
    chi2 = new TH1D("chi2","Chi2 / NDF; chi2/ndf/ampl", 10000, 0, 100);
  else
    chi2->Reset();

  if( !chi2prob )
    chi2prob = new TH1D("chi2prob","Chi2Prob; chi2prob", 50, 0, 1);
  else
    chi2prob->Reset();

  if( !amplVSchi2 )
    amplVSchi2 = new TH2D("anmplchi2","ampl VS Chi2 / NDF; ampl; chi2/ndf/ampl*200", 256, -.5, 4095.5, 1000, 0, 30);
  else
    amplVSchi2->Reset();
  
  if( !hsample )
    hsample = new TH1D("hsample","; sample; ampl", 32 ,-.5, 31.5);
  else
    hsample->Reset();
  
  if( !cutpos )
    cutpos = new TH2D("cutpos","Cutposition; row; column", 65, -.5, 64.5, 49, -.5,48.5);

  c->cd(1);
  cutpos->Draw("colz");
  c->cd(2);
  //c->cd(2)->SetLogz();
  amplVScor->Draw("colz");
  c->cd(3);
  //c->cd(3)->SetLogy();
  chi2->Draw();
  c->cd(4);
  //c->cd(4)->SetLogz();
  amplVSchi2->Draw("colz");
  //c->cd(5);
  //tVSchi2->Draw("colz");

  FILE* file= fopen(name,"r");
  
  if( !file ) {
    fprintf( stderr, "Cannot open file '%s'\n", name);
    return;
  }

  double corr, ampl, chisq, ndf, prob, time, base[2], fbase;
  int id, plane, x , y;
  unsigned int size, run(0), event(0), last_run(0), last_event(0);
  int cut_event = 0;

  int count = 0;
  int count2 = 0;

  //amplcor->Reset(); 
  //amplVScor->Reset();
  //chi2->Reset();
  //amplVSchi2-Reset();

  fseek(file, 0, SEEK_END);
  int fsize = ftell(file);
  fseek(file, 0, SEEK_SET);
  int fstart = ftell(file);
  int fnstep = 0.1 * (count2+1) * (fsize - fstart);
  printf("Start processing!\n");

  unsigned int count_tot(0), count_used(0), count_tot_event(0), count_cut_event(0);

  while( fread(&run, sizeof(unsigned int), 1, file) == 1 && 
         fread(&event, sizeof(unsigned int), 1, file) == 1 &&
         fread(&id, sizeof(int), 1, file) == 1 &&
         fread(&plane, sizeof(int), 1, file) == 1 &&
         fread(&x, sizeof(int), 1, file) == 1 &&
         fread(&y, sizeof(int), 1, file) == 1 &&
         fread(&corr, sizeof(double), 1, file) == 1 &&
         fread(&ampl, sizeof(double), 1, file) == 1 &&
         fread(&chisq, sizeof(double), 1 , file) == 1 &&
         fread(&ndf, sizeof(double), 1 , file) == 1 &&
         fread(&prob, sizeof(double), 1 , file) == 1 &&
         fread(&time, sizeof(double), 1 , file) == 1 &&
         fread(base, sizeof(double), 2, file) == 2 &&
         fread(&fbase, sizeof(double), 1, file) == 1 &&
         fread(&size, sizeof(unsigned int), 1, file) == 1) {

    if ( last_run != run || last_event != event ) {
      if( last_run ) {
        count_tot_event++;
        if ( cut_event )
          count_cut_event++;
      }
      cut_event = 0;
      last_run = run;
      last_event = event;
    }
    
    const unsigned int Size = size;
    assert(size == 32);
    unsigned short sample[Size];
    if( fread(sample, sizeof(unsigned short), size, file) != size ) 
      continue;
    

    if ( skip ) { 
      if ( count <= skip) {
        if( count == 1 ) printf("Skipping %u ." , skip);
        if( !(count%(skip/10)) ) printf( "%2.1f\%", float(count)*100./(float)skip );
        if( count == skip ) {
          fstart = ftell(file);
          fnstep = 0.1 * (count2+1) * (fsize - fstart);
        }
        fflush(stdout);
        continue;
      }
    }
    if ( (skip && count == skip+1) || (!count && !skip ) ) {
      printf("Analyzing ");
      if ( analyse == -1 ) {
        printf( "%4f MiB .", (fsize - fstart) / 1024. / 1024. );
      } else {
        printf(".");
      }
    }
    if( analyse == -1 ) {
      int pos = ftell(file) - fstart;
      if( pos>=fnstep ) { 
        count2++;
        printf( "%2.1f%%", (float)count2*10. );
        fnstep = 0.1 * (count2+1) * (fsize - fstart);
       }
    } else {
      if( !((count-skip)%(analyse/10)) ) printf( "%2.1f%%", (count-skip)*100/analyse);
      if( count > skip + analyse  ) {
        printf( "\n" );
        break;
      }
    }
    if( !(count%50000) ) {
     for( int a = 1; a < 7; a++) {
        c->cd(a)->Modified();
        c->cd(a)->Update();
      }
    }
    if( !(count%10000) ) {
      if( count ) {
        c->cd(2)->SetLogz();
        c->cd(3)->SetLogy();
        c->cd(4)->SetLogz();
      }
      cout << ".";
    }
    cout.flush();
    if ( ndf && ampl && (ampl-0.2)) {
      if( ( id_ == -1 || id == id_ ) && 
          ( plane_ == -1 || plane_ == plane ) &&
          ( ampl >= 25 ) ) {
        count_tot++;
        if( ((chisq/ndf/ampl*200) <= cut_high*ampl+cut_low ) ) {
          count_used++;
          chi2->Fill(chisq/ndf/ampl*200/ampl*200);
          chi2prob->Fill(prob);
          amplVSchi2->Fill(ampl,chisq/ndf/ampl*200);
          int ID = x*48+y;
          if ( ID < 0 ) ID = -1;
          amplcor->Fill(corr, ID);
          amplVScor->Fill(ampl, corr);
          tVSchi2->Fill(time, chisq/ndf);
        } else {
          cutpos->Fill(x,y);
          cut_event++;
#ifdef MONITOR
          if ( monitor ) {
            hsample->Reset();
            for( int i = 0; i < Size; i++ ) {
              hsample->Fill(i, sample[i] - base[i%2] - fbase);
            }
            c1->cd();
            hsample->Draw();
            if( f[id-1][plane] ) {
              f[id-1][plane]->SetParameters((time+27.8)/12.8,ampl,-fbase);
              f[id-1][plane]->Draw("same");
            }
            c1->Modified();
            c1->Update();
            if( monitor == 1 ) {
              char *name;
              asprintf(&name, "/tmp/sig_%f_%f.pdf", (chisq/ndf/ampl*200), ampl );
              c1->SaveAs(name);
              free( name );
            }else {
              for (int j = 0; j < monitor; j++) {
                gSystem->ProcessEvents();
                usleep(100);
              }
            }
          }
#endif     
        }   
      }
    }
    count++;
    gSystem->ProcessEvents();
  }

  cout << endl;
  if( count ) {
    c->cd(2)->SetLogz();
    c->cd(3)->SetLogy();
    c->cd(4)->SetLogz();
  }

  c->cd(1);
  cutpos->Draw("colz");
  c->cd(2);
  //c->cd(2)->SetLogz();
  amplVScor->Draw("colz");
  c->cd(3);
  //c->cd(3)->SetLogy();
  chi2->Draw();
  c->cd(4);
  //c->cd(4)->SetLogz();
  amplVSchi2->Draw("colz");
  //c->cd(5);
  //tVSchi2->Draw("colz");

  printf( "total = %u, used = %u, cut = %2.1f%%\n", count_tot, count_used, (float)(count_tot - count_used) / (float)(count_tot) * 100.);
  printf( "total events = %u, cut events = %u (%2.1f%%)\n", count_tot_event, count_cut_event, (float)count_cut_event / (float)count_tot_event * 100. );
};
