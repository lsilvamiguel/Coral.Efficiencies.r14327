/*
--------------------------------------------------------------------------

 Implementation of classes to store clusters from GEM
 detectors

 Author: Bernhard Ketzer   27/11/2002   v0
                           20/05/2008   v1

--------------------------------------------------------------------------
*/

// Class declaration
#include "CsGEMCluster.h"

// C++ headers
#include <cmath>
#include <cstdlib>

//-----------------------------------------------------------------------------
// CsGEMCluster constructor
//-----------------------------------------------------------------------------
CsGEMCluster::CsGEMCluster(const std::list<CsGEMHit*> &_hits, bool _flags,
                           const CsGEMPlanePar* _par) :
    fClusterParams(_par)
{
  // Save reference to hits in cluster
  fHits = _hits;

  // does the cluster contain flagged channels
  fContainsFlaggedChannels = _flags;

  // calculate cluster amplitude from hits
  CalcAmps();

  // calculate cluster position from hits
  CalcCoG();

  // calculate cluster time
  CalcTime();
}

//-----------------------------------------------------------------------------
// AddHit method to add new hit to list of hits
//-----------------------------------------------------------------------------
bool CsGEMCluster::AddHit(CsGEMHit* _hit)
{
    bool insert = false;

    // Add if list of hits in cluster empty
    if (fHits.empty()) insert=true;

    // Check if new hit should be added to existing list of hits
    else {

        // Get list of active neighbours for new hit
        const std::list<CsGEMChan*>& neighbours = _hit->GetChan()->GetActiveNeighbours();
        for (std::list<CsGEMChan*>::const_iterator itchan=neighbours.begin(); itchan!=neighbours.end(); ++itchan) {

            // Check if any hit in cluster is a neighbour of the new hit
            for (std::list<CsGEMHit*>::const_iterator ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
                if ((*itchan) == (*ithit)->GetChan()) {
                    insert = true;
//                    std::cout<<"CsGEMCluster::AddHit() : hit "
//                             <<_hit->GetChan()->GetId()->GetDetectorChannel()
//                             <<" neighbour of "
//                             <<(*ithit)->GetChan()->GetId()->GetDetectorChannel()
//                             <<std::endl;

                // Check if there is a gap between neighbours
                    if (abs((*itchan)->GetId()->GetDetectorChannel()-(*ithit)->GetChan()->GetId()->GetDetectorChannel())>1)
                        fContainsFlaggedChannels = true;
                }
            }
        }
    }

    // Add hit to cluster
    if (insert) {
        fHits.push_back(_hit);
    }

    return insert;
}

void CsGEMCluster::CalcAmps() {
    std::vector<float> sum(3,0.);
    double sumt = 0.;
    double wsumxt = 0.;
    double ssumsq = 0.;
    double ssum;

    // calculate position by centre-of-gravity method
    // Loop over hits
    std::list<CsGEMHit*>::iterator ithit;
    for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {

        // Detector channel, hemisphere, sigma of hit
        float sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();

        // get the fraction of the hit amplitude that should be assigne to this cluster
        float frac = 1.;
        if (fClusterParams->GetShareHits() == 2)
            frac = 1. / (*ithit)->GetNrClusters();

        // Sums
        for (int i=0; i<3; i++) {
            sum[i] += (*ithit)->GetAmp()[i]*frac;   // amplitude sum for each of the samples
        }
        sumt   += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma);
        wsumxt += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma) * (*ithit)->GetXTalk();
        ssumsq += sigma*sigma;
    }

    // Cluster noise
    ssum = sqrt(ssumsq/(double)fHits.size());

    // cluster cross talk
    wsumxt /= sumt;

    fAmp     = sum;
    fNoise   = ssum;
    fXTalk   = wsumxt;

    return;
}

//-----------------------------------------------------------------------------
// Center of Gravity clustering method
//-----------------------------------------------------------------------------
void CsGEMCluster::CalcCoG() {
    std::vector<float> sum(3,0.);
    float sigma;
    int detchan;
    double sumt = 0.;
    double wsumt = 0.;
    double hem;
    double whemt = 0.;

    // Loop over hits
    std::list<CsGEMHit*>::iterator ithit;
    for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {

        // Detector channel, hemisphere, sigma of hit
        detchan = (*ithit)->GetChan()->GetId()->GetDetectorChannel();
        hem = (*ithit)->GetChan()->GetId()->GetHemisphere();
        sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();

        // get the fraction of the hit amplitude that should be assigne to this cluster
        float frac = 1.;
        if (fClusterParams->GetShareHits() == 2)
            frac = 1. / (*ithit)->GetNrClusters();

        sumt  += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma);
        wsumt += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma) * detchan;
        whemt += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma) * hem;
    }

    // Center of gravity
    double center = wsumt / sumt;

    // Position error: depends on cluster size
    double dcenter(1./sqrt(12.));
    if (fHits.size()>fClusterParams->GetClusSizeRes().size()) {
        if (fClusterParams->GetClusSizeRes().size()>0)
            dcenter = fClusterParams->GetClusSizeRes()[fClusterParams->GetClusSizeRes().size()-1];
    } else
        dcenter = fClusterParams->GetClusSizeRes()[fHits.size()-1];

    // Weighted hemisphere
    whemt = whemt / sumt;

    // Update cluster properties
    fPos = center;
    fPosErr = dcenter;
    fHemisphere = whemt;

    return;
}

//-----------------------------------------------------------------------------
// Calculate cluster time from weighted mean of hit times
//-----------------------------------------------------------------------------
void CsGEMCluster::CalcTime() {

    fHasTime = false;
    fTime = 0.;
    fTimeErr = 0.;

    double thit, ethit;

    // Loop over hits
    std::list<CsGEMHit*>::iterator ithit;
    for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
        if ( (*ithit)->GetTime(thit, ethit) ) {
            fHasTime  = true;
            fTime    += thit / (ethit*ethit);
            fTimeErr +=  1.  / (ethit*ethit);
        }
    }
    if (fHasTime) {
        fTime    = fTime / fTimeErr;
        fTimeErr = sqrt(1./fTimeErr);
    }

    return;
}

//-----------------------------------------------------------------------------
// Get the cluster time
//-----------------------------------------------------------------------------
bool CsGEMCluster::GetTime(double &_time, double &_etime) const {
    _time  = fTime;
    _etime = fTimeErr;

    return fHasTime;
}

//-----------------------------------------------------------------------------
// Comparison of 2 clusters
//-----------------------------------------------------------------------------
bool CsGEMCluster::operator< (const CsGEMCluster &cluster) const {
    return (fPos < cluster.GetPosition());
}

/*
//-----------------------------------------------------------------------------
// Fit cluster function
//-----------------------------------------------------------------------------
void CsGEMCluster::Fit(int _sample, float _thr)
{
  // Hit iterator
  std::list<CsGEMHit*>::const_iterator ithit;

  // Get number of hits in group
  int nhit = fHits.size();

  // Local arrays
  float x[nhit], y1[nhit], y2[nhit], y3[nhit], ey[nhit];
  float a[nhit+2], b[nhit+2], eb[nhit+2];
  float zero[nhit+2];
  float eyavg=0.;
  float hem(0.), hemw(0.);
  for (int i=0; i<nhit+2; i++){
    a[i] = 0.;
    b[i] = 0.;
    eb[i] = 0.;
    zero[i] = 0.;
  }
  int i;
  for (ithit=hits.begin(), i=0; ithit!=hits.end(); ithit++, i++) {
    x[i] = a[i+1] = (*ithit)->GetStripNumber();
    y1[i] = (*ithit)->GetAmp1();
    y2[i] = (*ithit)->GetAmp2();
    y3[i] = (*ithit)->GetAmp3();
    if (_sample==3) b[i+1] = (*ithit)->GetAmp3();
    else if (_sample==2) b[i+1] = (*ithit)->GetAmp2();
    else b[i+1] = (*ithit)->GetAmp1();
    eb[i+1] = ey[i] = cals[i]->GetSigma();
    eyavg += ey[i];
    hem += b[i+1] * (*ithit)->GetHemisphere();
    hemw += b[i+1];
  }

  // Normalize hemisphere
  hem /= hemw;

  // Average sigma of strips in group
  eyavg = eyavg/nhit;

  // Add two elements with zero amplitude on either side of the cluster
  a[0] = a[1]-3.;  // go left 3 strips
  b[0] = 0.; eb[0] = eyavg;
  a[nhit+1] = a[nhit]+3.;  // go right 3 strips
  b[nhit+1] = 0.; eb[nhit+1] = eyavg;

  //Calculate amplitude difference to previous hit strip
  float delta[nhit+3];
  delta[0]=0.;
  delta[nhit+2]=0.;
  for (int i=1; i<nhit+2; i++) {
    delta[i] = b[i] - b[i-1];
  }

  // Look for local maxima
  int j=0;
  int imax[nhit];
  for (int i=1; i<nhit+1; i++) {
    if (delta[i]>0. && delta[i+1]<0. &&
        (delta[i-1]>0. || delta[i+2]<0.)
        //((delta[i-1]>0. && delta[i+2]<0.) ||
         //(delta[i-1]>0. && delta[i+1]<-1.*sqrt(eb[i]*eb[i]+eb[i+1]*eb[i+1])) ||
         //(delta[i+2]<0. && delta[i]>sqrt(eb[i]*eb[i]+eb[i-1]*eb[i-1])))
        ) {
      imax[j] = i;
      //      std::cout<<"Local max found at "<<a[imax[j]]<<std::endl;
      j++;
    }
  }
  int nmax = j;

  // Maximum number of Gauss functions to fit
  // int ng = int(((float)nhit+2.)/3.);
  int ng = nmax;
  if (ng > 1) hem = 0.;
  if (ng > 8) ng=8;
  if (ng > int(((float)nhit+2.)/3.)) ng=int(((float)nhit+2.)/3.);
  if (ng == 0) ng=1;  // fit at least one gaussian

  // Create Graph
  TGraphErrors *cluster = new TGraphErrors(nhit+2,a,b,zero,eb);

  // Fit gauss function around each local maximum
  int npar = ng*3;
  double par[npar];
  TF1 *gauss[ng];
  char fcnstr[80];
  sprintf(fcnstr, "gaus(0)");
  TString totstr = fcnstr;
  for (int i=0; i<ng; i++) {
    if (ng==1) {
      //      std::cout<<"Fitting gauss from "<<a[0]<<" to "<< a[nhit+1]<<std::endl;
      gauss[i] = new TF1("g","gaus",a[0],a[nhit+1]);
    }
    else {
      //      std::cout<<"Fitting gauss from "<<a[imax[i]-1]<<" to "<< a[imax[i]+1]<<std::endl;
      gauss[i] = new TF1("g","gaus",a[imax[i]-1],a[imax[i]+1]);
    }
    cluster->Fit(gauss[i],"QR0+");
    gauss[i]->GetParameters(&par[3*i]);
    if (par[3*i+2]<0.)  par[3*i+2] *= -1.;   // sigma negative
    if (i==0) continue;
    sprintf(fcnstr, "+gaus(%i)", 3*i);
    totstr += fcnstr;
  }

  // Fit sum of all gauss functions
  //  std::cout << "Fitting "<<totstr<<" from "<<a[0]<<" to "<<a[nhit+1]<<std::endl;
  TF1 *tot = new TF1("tot", (const char *)totstr, a[0], a[nhit+1]);
  tot->SetParameters(par);
  for (int i=0; i<ng; i++) {
    //    tot->SetParLimits(3*i, 0., 1000.);          // amplitude
    tot->SetParLimits(3*i+1, -1500., 1500.);   // center
    tot->SetParLimits(3*i+2, 0.5, 100.);       // sigma
  }
  cluster->Fit(tot,"BQR0+");

  // Cluster parameters
  for (int i=0; i<ng; i++) {
    float amp = tot->GetParameter(3*i);
    float center = tot->GetParameter(3*i+1);
    float sigma = tot->GetParameter(3*i+2);
    float noise = tot->GetParError(3*i);
    float dcenter = tot->GetParError(3*i+1);
    dcenter = sqrt(pow(dcenter,2)+0.0833); // fit^2 + 1/sqrt(12)^2
    float size = 2*sqrt(2*sigma*sigma*log(amp/(thrstrip*eyavg)));
    //    float sum = amp*sigma*sqrt(2*M_PI);
    float sum = amp*sigma*sqrt(2*3.1415);

    // Add cluster
    if (sum>thrclus*noise) {
      //      std::cout<<"Cluster FIT "<<center<<" "<<sum<<" "<<size<<std::endl;

      // Get 3 amplitudes of non-flagged hit closest to cluster center
      int iclose = -1;
      for (int i=0; i<nhit; i++) {
        if (x[i] < center) continue;
        bool check = false;
        if ( i != 0) check = ((x[i]-center)<(center-x[i-1]));
        if ( (x[i]==center) || check ) {
          iclose = i;
          break;
        }
        else {
          iclose = i-1;
          break;
        }
      }

      // Hit found close to cluster center
      if ( iclose != -1 && fabs(x[iclose]-center)<3. )  {
        double y;
//        std::cout << "Closest strip " << x[iclose] << std::endl;
        if (_sample==3) y = y3[iclose];
        else if (_sample==2) y = y2[iclose];
        else y = y1[iclose];
        double sum1 = y1[iclose] / y * sum;
        double sum2 = y2[iclose] / y * sum;
        double sum3 = y3[iclose] / y * sum;
//        std::cout << sum1 << " " << sum2 << " " << sum3 << std::endl;

        // Calculate time of arrival
        double time = getGEMtime(sum1, sum2, sum3, timecals);
        double dtime = 12.; // 12ns measured time resolution
        if (time == 1000) dtime = 1.E+10;

        // Add cluster
        AddCluster(center,dcenter,sum1,sum2,sum3,noise,size,time,dtime,hits,hem);
      }
      //      std::cout<<"------------------------------------"<<std::endl;
    }

    // Clean up
    gauss[i]->Delete();
  }

//   // Plot
//   TCanvas *fitc = new TCanvas("fitc","fitc",400,400);
//   fitc->cd();
//   cluster->Draw("AP");
//   tot->Draw("SAME");
//   TLine *l = new TLine(a[0], thrstrip*eyavg, a[nhit+1], thrstrip*eyavg);
//   l->Draw();
//   fitc->Modified();
//   fitc->Update();
//   int dum; cin>>dum;

  tot->Delete();
  cluster->Delete();
}
*/
