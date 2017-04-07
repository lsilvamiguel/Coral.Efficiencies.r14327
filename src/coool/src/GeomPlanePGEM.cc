
#include <iostream>
#include <queue>

#include "GeomPlanePGEM.h"
#include "math.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TString.h"
#include "ChipAPV.h"

GeomPlanePGEM::GeomPlanePGEM(int id,const char* name,int nwir,
			   double x,double y,double z,
			   double dx,double dy,double dz,
			   double ang,double pitch):   
  GeomPlaneGEM(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch)   {
  fDefFlag  =   1;
  fDefPed   = 750.;
  fDefSigma =   4.5;

  fNwirV  = 0;
  fPitchV = 0.;
  
  fPlane = new CsGEMPlane(GetName());
  fPlanePix = 0;
} 

GeomPlanePGEM::~GeomPlanePGEM() {
  if (fPlane)
    delete fPlane;
  if (fPlanePix)
    delete fPlanePix;
}

void GeomPlanePGEM::ResetPlane() {
  if (fPlane) {
    GeomPlaneGEM::ResetPlane();
  } else {
    fPlanePix->Clear();
  }
}

bool GeomPlanePGEM::CalcTime(const CDigit1 &digit, double &time) {
  if (fPlane) {
    return GeomPlaneGEM::CalcTime(digit, time);
  } else {
    std::vector<float> amps;
    amps.push_back(digit.dt[2]);
    amps.push_back(digit.dt[3]);
    amps.push_back(digit.dt[4]);

    std::map<CsGEMChanId,CsGEMChan*>::const_iterator chit = fPlanePix->GetChannels().find(CsGEMChanId(digit.ch, (int) digit.dt[0], (int) digit.dt[1]));

    if ( chit == fPlanePix->GetChannels().end() ) {
      std::cout << "Trying to calculate a time for non-existing channel" << std::endl;
      return false;
    }

    float sigma = chit->second->GetCal()->GetPedSigma();

    double timeerr;
          
    return fPlanePix->GetTimeCals()->CalcTime(amps, sigma, time, timeerr);
  }
}
            
std::set<CCluster1>& GeomPlanePGEM::Clusterize(const std::map<int,CDigit1>& digits) {
  fClusters.clear();

  if (fPlane) {
    GeomPlaneGEM::Clusterize(digits);
  } else if (fPlanePix) {
    // fill single hits into plane
    // loop over hits in plane (they are sorted by strip number (map))
    for (std::map<int,CDigit1>::const_iterator id = digits.begin(); id!=digits.end(); id++) {
      const CDigit1& digit = id->second;

      fPlanePix->AddHit(digit.ch, (int) digit.dt[0], (int) digit.dt[1], digit.dt[2], digit.dt[3], digit.dt[4]);
    } // for(map....)

    // hits were already filled start clustering
    // call clusterize method
    fPlanePix->Clusterize();

    // read clusters from plane and add to cluster list
    std::list<CsPixelGEMCluster*>::const_iterator cluit;
    for (cluit=fPlanePix->GetClusters().begin(); cluit!=fPlanePix->GetClusters().end(); cluit++) {
      std::vector<double> data;
      std::vector<double> error;

      // calculate cluster position, hardcode spatial resolution to 90um
      double posX, posY;
      Pad2Pos((*cluit)->GetPositionX(), (*cluit)->GetPositionY(), posX, posY);
      double res = 0.009;

      // save the cluster amplitudes
      data.push_back((*cluit)->GetAmp()[0]); error.push_back((*cluit)->GetNoise());
      data.push_back((*cluit)->GetAmp()[1]); error.push_back((*cluit)->GetNoise());
      data.push_back((*cluit)->GetAmp()[2]); error.push_back((*cluit)->GetNoise());

      // save strip position and hemisphere
      data.push_back ((*cluit)->GetPositionX());
      error.push_back( posX );
      data.push_back ((*cluit)->GetPositionY());
      error.push_back( posY );

      // save cluster time if it could be calculated
      double time, etime;
      if ((*cluit)->GetTime(time, etime)) {
        data.push_back(time); error.push_back(etime);
      }

      fClusters.insert( CCluster1( GetID(), posX, (*cluit)->GetSize(), res, data, error ) );
    }
  }

  return fClusters;
}

void GeomPlanePGEM::AddSubPlane (int id, const char *name, int nwires,
                                 double x, double y, double z,
                                 double dx, double dy, double dz,
                                 double angle, double inpitch) {
  // check the validity of the sub plane adding
  if ( strcmp(GetName(), name) )
    throw "GeomPlanePGEM::AddSubPlane : different names !";
      
  if ( fabs(fabs(angle-GetAngle())-90.) > 0.001 )
    throw "GeomPlanePGEM::AddSubPlane : two PixelGEM subplanes have to be perpendicular!";

  // this is a pixel plane, so delete strip plane and create pixel plane
  delete fPlane; fPlane = 0;
  fPlanePix = new CsPixelGEMPlane(GetName());

  // set number of wires and pitch for the second projection
  fNwirV  = nwires;
  fPitchV = inpitch;
}

#if USE_DATABASE==1
void GeomPlanePGEM::SetCalibrations(const std::map<unsigned int, PlaneGEM::APVcalib> &c, const CS::Chip::Maps *maps, const std::vector<float> &ct) {
  if (fPlane) {
    GeomPlaneGEM::SetCalibrations(c, maps, ct);
  } else if (fPlanePix) {
    // handle pedestal calibration first, time calibration in the ending
    // delete old calibrations
    fPlanePix->ClearChannels();

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
      // cc should now in any case be a CS::ChipAPV::Digit
      // (probably also a CS::ChipAPV::DigitPixel, but we do not care yet)
      const CS::ChipAPV::Digit *mt = dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
      if ( mt==NULL ) {
        std::cerr << "ChipAPV wrong map." << std::endl;
        continue;
      }

      // check that detector in mapping is really this detector
      if (strcmp(mt->GetDetID().GetName().c_str(), GetName()))
        continue;

      // now we can be more sure that this is a DigitPixel
      const CS::ChipAPV::DigitPixel *m = dynamic_cast<const CS::ChipAPV::DigitPixel*>(mt);
      if ( m==NULL ) {
        std::cerr << "ChipAPV wrong map. CS::ChipAPV::Digit is not a CS::ChipAPV::DigitPixel." << std::endl;
        continue;
      }

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
 
      std::pair<int,int> xy=pixgem::detch2xy(m->GetChannel());
      int pixX = xy.first;  if (m->GetDetOrientation() < 0) pixX=31-pixX;
      int pixY = xy.second;
      fPlanePix->AddChan(m->GetChannel(), pixX, pixY,
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
    fPlanePix->SetTimeCals(cals);

    // if time calibration is not valid, use the standard
    if ( !fPlanePix->GetTimeCals()->IsValid() ) {
      CsGEMTimeCalOld *cal02 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
      CsGEMTimeCalOld *cal12 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
      fPlanePix->SetTimeCals(CsGEMTimeCals(cal02, cal12));
    }
  }
}
#else
void GeomPlanePGEM::SetCalibrations(const CS::Chip::Maps *maps) {
  // no database, no real calibrations can be loaded, load default ones
  if (fPlane) {
    GeomPlaneGEM::SetCalibrations(maps);
  } else if (fPlanePix) {
    // handle pedestal calibration first, time calibration in the ending
    // delete old calibrations
    fPlanePix->ClearChannels();

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
      // cc should now in any case be a CS::ChipAPV::Digit
      // (probably also a CS::ChipAPV::DigitPixel, but we do not care yet)
      const CS::ChipAPV::Digit *mt = dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
      if ( mt==NULL ) {
        std::cerr << "ChipAPV wrong map." << std::endl;
        continue;
      }
    
      // check that detector in mapping is really this detector
      if (strcmp(mt->GetDetID().GetName().c_str(), GetName()))
        continue;

      // now we can be more sure that this is a DigitPixel
      const CS::ChipAPV::DigitPixel *m = dynamic_cast<const CS::ChipAPV::DigitPixel*>(mt);
      if ( m==NULL ) {
        std::cerr << "ChipAPV wrong map. CS::ChipAPV::Digit is not a CS::ChipAPV::DigitPixel." << std::endl;
        continue;
      }
 
      std::pair<int,int> xy=pixgem::detch2xy(m->GetChannel());
      int pixX = xy.first;  if (m->GetDetOrientation() < 0) pixX=31-pixX;
      int pixY = xy.second;
    
      fPlanePix->AddChan(m->GetChannel(), pixX, pixY,
                         fDefFlag, fDefPed, fDefSigma,
                         m->GetChip(), m->GetChipChannel());
    }

    // time calibration from here on
    CsGEMTimeCalOld *cal02 = new CsGEMTimeCalOld(0.2,1.1,-40.,22.,1.65);
    CsGEMTimeCalOld *cal12 = new CsGEMTimeCalOld(0.3,1.0,0.,27.,1.23);
    fPlanePix->SetTimeCals(CsGEMTimeCals(cal02, cal12));
  }
}
#endif

void GeomPlanePGEM::SetParameters(float ts, float tc, int s) {
  if (fPlane) {
    GeomPlaneGEM::SetParameters(ts, tc, s);
  } else if (fPlanePix) {
    fPlanePix->GetPar()->SetThrHit(ts);
    fPlanePix->GetPar()->SetThrClus(tc);
    fPlanePix->GetPar()->SetSample(s);
  }
}

void GeomPlanePGEM::Pad2Pos(double padX, double padY, double &posX, double &posY) {
  double ctrX = GetNWires()/2. - .5;
  posX = (padX-ctrX) * GetPitch();

  double ctrY = GetNWiresV()/2. - .5;
  posY = (padY-ctrY) * GetPitchV();
}

