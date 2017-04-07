
#include <iostream>
#include <queue>

#include "GeomPlaneGEM.h"
#include "math.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TString.h"
#include "ChipAPV.h"

// in case no calibration is loaded we need a way to find the source ID of this detector
#if USE_DATABASE==1
#else
#include "DetID.h"
#endif

//ClassImp(GeomPlaneGEM);

GeomPlaneGEM::GeomPlaneGEM(int id,const char* name,int nwir,
			   double x,double y,double z,
			   double dx,double dy,double dz,
			   double ang,double pitch):   
  GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch)   {
  fDefFlag  =   1;
  fDefPed   = 750.;
  fDefSigma =   4.5;

  fPlane = new CsGEMPlane(name);
} 

GeomPlaneGEM::~GeomPlaneGEM() {
  delete fPlane;
}

void GeomPlaneGEM::ResetPlane() {
  fPlane->Clear();
}

bool GeomPlaneGEM::CalcTime(const CDigit1 &digit, double &time) {
  std::vector<float> amps;
  amps.push_back(digit.dt[0]);
  amps.push_back(digit.dt[1]);
  amps.push_back(digit.dt[2]);

  std::map<CsGEMChanId,CsGEMChan*>::const_iterator chit = fPlane->GetChannels().find(CsGEMChanId(digit.ch, (int) digit.dt[9]));
  if ( chit == fPlane->GetChannels().end() ) {
    std::cout << "Trying to calculate a time for non-existing channel" << std::endl;
    return false;
  }

  float sigma = chit->second->GetCal()->GetPedSigma();

  double timeerr;

  return fPlane->GetTimeCals()->CalcTime(amps, sigma, time, timeerr);
}

std::set<CCluster1>& GeomPlaneGEM::Clusterize(const std::map<int,CDigit1>& digits) { 
  fClusters.clear();

  // fill single hits into plane
  // loop over hits in plane (they are sorted by strip number (map))
  for (std::map<int,CDigit1>::const_iterator id = digits.begin(); id!=digits.end(); id++) {
    const CDigit1& digit = id->second;

    fPlane->AddHit(digit.ch, (int) digit.dt[9], digit.dt[0], digit.dt[1], digit.dt[2]);
  } // for(map....)

  // hits were already filled start clustering
  // call clusterize method
  fPlane->Clusterize();

  // read clusters from plane and add to cluster list
  std::list<CsGEMCluster*>::const_iterator cluit;
  for (cluit=fPlane->GetClusters().begin(); cluit!=fPlane->GetClusters().end(); cluit++) {
    std::vector<double> data;
    std::vector<double> error;

    // calculate cluster position, hardcode spatial resolution to 70um
    double pos = Wire2Pos((*cluit)->GetPosition());
    double res = 0.007;

    // save the cluster amplitudes
    data.push_back((*cluit)->GetAmp()[0]); error.push_back((*cluit)->GetNoise());
    data.push_back((*cluit)->GetAmp()[1]); error.push_back((*cluit)->GetNoise());
    data.push_back((*cluit)->GetAmp()[2]); error.push_back((*cluit)->GetNoise());

    // save strip position and hemisphere
    data.push_back ((*cluit)->GetPosition());
    error.push_back((*cluit)->GetHemisphere());

    // save cluster time if it could be calculated
    double time, etime;
    if ((*cluit)->GetTime(time, etime)) {
      data.push_back(time); error.push_back(etime);
    }

    fClusters.insert( CCluster1( GetID(), pos, (*cluit)->GetSize(), res, data, error ) );
  }
 
  return fClusters;
}

#if USE_DATABASE==1
void GeomPlaneGEM::SetCalibrations(const std::map<unsigned int, PlaneGEM::APVcalib> &c, const CS::Chip::Maps *maps, const std::vector<float> &ct) {
  // handle pedestal calibration first, time calibration in the ending
  // delete old calibrations
  fPlane->ClearChannels();

  // need GEM.xml map to sort calibration data
  typedef CS::Chip::Maps::const_iterator m_it;
  
  // assume same source id for all apvs/channels from one GEM
  std::set<uint16> srcIDs;
  maps->GetSrcIDs( CS::DetID(GetName()), srcIDs );
  if ( srcIDs.size() != 1 )
    std::cerr << GetName() << ": One plane must be connected to only one source ID!" << std::endl;

  unsigned int srcId = *srcIDs.begin();
  m_it m_low = maps->lower_bound(uint64(srcId)<<48);
  for( m_it cc=m_low; ((*cc).first>>48)==srcId; cc++ ) {
    const CS::ChipAPV::Digit *m = dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
    if ( m==NULL ) {
      std::cerr << "ChipAPV wrong map." << std::endl;
      continue;
    }
    
    // check that detector in mapping is really this detector
    if (strcmp(m->GetDetID().GetName().c_str(), GetName()))
      continue;
    
    int   flag  = fDefFlag;
    float ped   = fDefPed;
    float sigma = fDefSigma;

    // test if there are valid pedestals for current channel
    if ( c.count(m->GetChip()) && c.find(m->GetChip())->second.channel.size()>m->GetChipChannel() ) {
      const PlaneGEM::APVchannel a_ch = c.find(m->GetChip())->second.channel[m->GetChipChannel()];
      flag  = a_ch.flag;
      ped   = a_ch.ped;
      sigma = a_ch.sigma;
    } else { // use default and print a warning
      static bool printed=false;
      if (!printed)
        std::cerr << "Broken pedestals/mapping for " << GetName() << std::endl;
      printed=true;
    }

    fPlane->AddChan(m->GetChannel(), m->GetChanPos(),
                   flag, ped, sigma,
                   m->GetChip(), m->GetChipChannel());
  }

  // time calibration from here on
  std::stringbuf buf;
  std::iostream  stream(&buf);

  // write the constants into a buffer
  for (unsigned int i=0; i<ct.size(); i++)
    stream << ct[i] << " ";
  stream.seekg(0, std::ios_base::beg);

  // read calibration from buffer
  CsGEMTimeCals cals; cals.Clear();
  stream >> cals;
  fPlane->SetTimeCals(cals);

  // if time calibration is not valid, use the standard
  if ( !fPlane->GetTimeCals()->IsValid() ) {
    CsGEMTimeCalOld *cal02 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
    CsGEMTimeCalOld *cal12 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
    fPlane->SetTimeCals(CsGEMTimeCals(cal02, cal12));
  }
}
#else
void GeomPlaneGEM::SetCalibrations(const CS::Chip::Maps *maps) {
  // no database, no real calibrations can be loaded, load default ones
  // handle pedestal calibration first, time calibration in the ending
  // delete old calibrations
  fPlane->ClearChannels();

  // need GEM.xml map to sort calibration data
  typedef CS::Chip::Maps::const_iterator m_it;
  
  // assume same source id for all apvs/channels from one GEM
  std::set<uint16> srcIDs;
  maps->GetSrcIDs( CS::DetID(GetName()), srcIDs );
  if ( srcIDs.size() != 1 )
    std::cerr << GetName() << ": One plane must be connected to only one source ID!" << std::endl;;
  unsigned int srcId = *srcIDs.begin();
  m_it m_low = maps->lower_bound(uint64(srcId)<<48);
  for( m_it cc=m_low; ((*cc).first>>48)==srcId; cc++ ) {
    const CS::ChipAPV::Digit *m = dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
    if ( m==NULL ) {
      std::cerr << "ChipAPV wrong map." << std::endl;
      continue;
    }
    
    // check that detector in mapping is really this detector
    if (strcmp(m->GetDetID().GetName().c_str(), GetName()))
      continue;
    
    fPlane->AddChan(m->GetChannel(), m->GetChanPos(),
                   fDefFlag, fDefPed, fDefSigma,
                   m->GetChip(), m->GetChipChannel());
  }

  // time calibration from here on
  CsGEMTimeCalOld *cal02 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
  CsGEMTimeCalOld *cal12 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
  fPlane->SetTimeCals(CsGEMTimeCals(cal02, cal12));
}
#endif // USE_DATABASE==1

void GeomPlaneGEM::SetParameters(float ts, float tc, int s) {
  fPlane->GetPar()->SetThrHit(ts);
  fPlane->GetPar()->SetThrClus(tc);
  fPlane->GetPar()->SetSample(s);
};



