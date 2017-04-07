
#include <iostream>
#include <queue>

#include "GeomPlanePMumega.h"
#include "math.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TString.h"
#include "ChipAPV.h"



GeomPlanePMumega::GeomPlanePMumega(int id,const char* name,int nwir,
			           double x,double y,double z,
			           double dx,double dy,double dz,
			           double ang,double pitch):   
  GeomPlane(id,name,nwir,x,y,z,dx,dy,dz,ang,pitch), fPlane(0), fPlanePix(0)  {
  fDefFlag  =   1;
  fDefPed   = 750.;
  fDefSigma =   4.5;

  fNwirV  = 0;
  fPitchV = 0.;
  
  fPlane = 0;
  fPlanePix = 0;
  fPlanePixMM = 0;

  if (name[4] == 'P') fPlanePix = new CsPixelMumegaPlane(GetName());
  else if (name[4] == 'M') fPlanePixMM = new CsRectPixelMumegaPlane(GetName());
  else fPlane = new CsMumegaPlane(GetName());

  if (fPlanePixMM) fPlanePixMM->SetPar(0); //  SimpleClustering

}



GeomPlanePMumega::~GeomPlanePMumega() {
  if (fPlane)    delete fPlane;
  if (fPlanePix) delete fPlanePix;
  if (fPlanePixMM) delete fPlanePixMM;
}



void GeomPlanePMumega::ResetPlane() {
  if (fPlane)    fPlane->Clear();
  if (fPlanePix) fPlanePix->Clear();
  if (fPlanePixMM) fPlanePixMM->Clear();
}



bool GeomPlanePMumega::CalcTime(const CDigit1 &digit, double &time) {
  std::vector<float> amps;
  double timeerr = 0;
  float sigma = 0;

  if (fPlane) {
    amps.push_back(digit.dt[0]);
    amps.push_back(digit.dt[1]);
    amps.push_back(digit.dt[2]);
  }
  if (fPlanePix || fPlanePixMM) {
    amps.push_back(digit.dt[2]);
    amps.push_back(digit.dt[3]);
    amps.push_back(digit.dt[4]);
  }

  std::map<CsMumegaChanId,CsMumegaChan*>::const_iterator chit;
  if (fPlane) {
    chit = fPlane->GetChannels().find(CsMumegaChanId(digit.ch, (int) digit.dt[9]));
    if ( chit == fPlane->GetChannels().end() ) {
      std::cout << "Trying to calculate a time for non-existing channel" << std::endl;
      return false;
    }

    sigma = chit->second->GetCal()->GetPedSigma();
    return fPlane->GetTimeCals()->CalcTime(amps, sigma, time, timeerr);
  }

  if (fPlanePix) {
    chit = fPlanePix->GetChannels().find(CsMumegaChanId(digit.ch, (int) digit.dt[0], (int) digit.dt[1]));
    if ( chit == fPlanePix->GetChannels().end() ) {
      std::cout << "Trying to calculate a time for non-existing channel" << std::endl;
      return false;
    }

    sigma = chit->second->GetCal()->GetPedSigma();
    return fPlanePix->GetTimeCals()->CalcTime(amps, sigma, time, timeerr);
  }

  if (fPlanePixMM) {
    chit = fPlanePixMM->GetChannels().find(CsMumegaChanId(digit.ch, (int) digit.dt[11], (float) digit.dt[0], (float) digit.dt[1]));
    if ( chit == fPlanePixMM->GetChannels().end() ) {
      std::cout << "Trying to calculate a time for non-existing channel" << std::endl;
      return false;
    }

    sigma = chit->second->GetCal()->GetPedSigma();
    return fPlanePixMM->GetTimeCals()->CalcTime(amps, sigma, time, timeerr);
  }

  return false;
}



std::set<CCluster1>& GeomPlanePMumega::Clusterize(const std::map<int,CDigit1>& digits) {
  fClusters.clear();

  // fill single hits into plane
  // loop over hits in plane (they are sorted by strip number (map))
  for (std::map<int,CDigit1>::const_iterator id = digits.begin(); id!=digits.end(); id++) {
    const CDigit1& digit = id->second;

    if (fPlane) fPlane->AddHit(digit.ch, (int) digit.dt[9], digit.dt[0], digit.dt[1], digit.dt[2]);
    if (fPlanePix) fPlanePix->AddHit(digit.ch, (int) digit.dt[0], (int) digit.dt[1], digit.dt[2], digit.dt[3], digit.dt[4]);
    if (fPlanePixMM) fPlanePixMM->AddHit(digit.ch, (int) digit.dt[11], digit.dt[0], digit.dt[1], digit.dt[2], digit.dt[3], digit.dt[4]);
  } // for(map....)

  // hits were already filled start clustering
  // call clusterize method
  if (fPlane) fPlane->Clusterize();
  if (fPlanePix) fPlanePix->Clusterize();
  if (fPlanePixMM) fPlanePixMM->Clusterize();


  // read clusters from plane and add to cluster list
  if (fPlane) {
    std::list<CsMumegaCluster*>::const_iterator cluit;
    double res = 0.007;
    double pos;
    std::vector<double> data;
    std::vector<double> error;
    for (cluit=fPlane->GetClusters().begin(); cluit!=fPlane->GetClusters().end(); cluit++) {
      data.clear();
      error.clear();

      // calculate cluster position, hardcode spatial resolution to 70um
      pos = Wire2Pos((*cluit)->GetPosition());

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
  }

  if (fPlanePix) {
    std::list<CsPixelMumegaCluster*>::const_iterator cluit;
    double posX, posY;
    double res = 0.009;
    double time, etime;
    std::vector<double> data;
    std::vector<double> error;
    for (cluit=fPlanePix->GetClusters().begin(); cluit!=fPlanePix->GetClusters().end(); cluit++) {
      data.clear();
      error.clear();

      // calculate cluster position, hardcode spatial resolution to 90um
      Pad2Pos((*cluit)->GetPositionX(), (*cluit)->GetPositionY(), posX, posY);

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
      if ((*cluit)->GetTime(time, etime)) {
        data.push_back(time); error.push_back(etime);
      }

      fClusters.insert( CCluster1( GetID(), posX, (*cluit)->GetSize(), res, data, error ) );
    }
  }

  if (fPlanePixMM) {
    std::list<CsRectPixelMumegaCluster*>::const_iterator cluit;
    double posX, posY;
    double res = 0.009;
    double time, etime;
    std::vector<double> data;
    std::vector<double> error;
    for (cluit=fPlanePixMM->GetClusters().begin(); cluit!=fPlanePixMM->GetClusters().end(); cluit++) {
      data.clear();
      error.clear();

      // calculate cluster position, hardcode spatial resolution to 90um
      posX = (*cluit)->GetPositionX();
      posY = (*cluit)->GetPositionY();

      // save the cluster amplitudes
      data.push_back((*cluit)->GetAmp()[0]); error.push_back((*cluit)->GetNoise());
      data.push_back((*cluit)->GetAmp()[1]); error.push_back((*cluit)->GetNoise());
      data.push_back((*cluit)->GetAmp()[2]); error.push_back((*cluit)->GetNoise());

      // save strip position and hemisphere

#warning GeomPlanePMumega::Clusterize: data clusters to improve
      data.push_back ( posX );
      error.push_back( posX );
      data.push_back ( posY );
      error.push_back( posY );

      // save cluster time if it could be calculated
      if ((*cluit)->GetTime(time, etime)) {
        data.push_back(time); error.push_back(etime);
      }

      fClusters.insert( CCluster1( GetID(), posX, (*cluit)->GetSize(), res, data, error ) );
    }
  }

  return fClusters;
}



void GeomPlanePMumega::AddSubPlane (int id, const char *name, Int_t nwires,
                                    double x, double y, double z,
                                    double dx, double dy, double dz,
                                    double angle, double inpitch) {
  // check the validity of the sub plane adding
  if ( strcmp(GetName(), name) )
    throw "GeomPlanePMumega::AddSubPlane : different names !";

  if (fPlane) {  // strips
    GeomPlane::AddSubPlane(id, name, nwires, x, y, z, dx, dy, dz, angle, inpitch);
  } else {  // pixels PGEM or PMM
    if ( fabs(fabs(angle-GetAngle())-90.) > 0.001 )
      throw "GeomPlanePMumega::AddSubPlane : two PixelMumega subplanes have to be perpendicular!";

    // set number of wires and pitch for the second projection
    fNwirV  = nwires;
    fPitchV = inpitch;
  }
}



#if USE_DATABASE==1
void GeomPlanePMumega::SetCalibrations(const std::map<unsigned int, PlanePMumega::APVcalib> &c, const CS::Chip::Maps *maps, const std::vector<float> &ct) {

  // handle pedestal calibration first, time calibration in the ending
  // delete old calibrations
  if (fPlane) fPlane->ClearChannels();
  if (fPlanePix) fPlanePix->ClearChannels();
  if (fPlanePixMM) fPlanePixMM->ClearChannels();

  // need PixelMM.xml map to sort calibration data
  typedef CS::Chip::Maps::const_iterator m_it;

  // assume same source id for all apvs/channels from one PixelMM plane (probably false later..)
  std::set<uint16> srcIDs;
  maps->GetSrcIDs( CS::DetID(GetName()), srcIDs );
  if ( srcIDs.size() != 1 )
    std::cerr << GetName() << ": One plane must be connected to only one source ID!" << std::endl;

  unsigned int srcId = *srcIDs.begin();
  m_it m_low = maps->lower_bound(uint64(srcId)<<48);

  register int   flag  = 0;
  register float ped   = 0;
  register float sigma = 0;
  for( m_it cc=m_low; ((*cc).first>>48)==srcId; cc++ ) {
    // cc should now in any case be a CS::ChipAPV::Digit
    // (probably also a CS::ChipAPV::DigitPixel, but we do not care yet)
    const CS::ChipAPV::Digit *m = dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
    if ( m==NULL ) {
      std::cerr << "GeomPlanePMumega::SetCalibrations: ChipAPV wrong map,digit is not a CS::ChipAPV::Digit\n";
      continue;
    }
    
    // check that detector in mapping is really this detector
    if (strcmp(m->GetDetID().GetName().c_str(), GetName()))
      continue;

    // now we can be more sure that this is a DigitPixel
    const CS::ChipAPV::DigitPixel *mp = dynamic_cast<const CS::ChipAPV::DigitPixel*>(m);
    if ( fPlanePix && (mp==NULL) ) {
      std::cerr << "GeomPlanePMumega::SetCalibrations: ChipAPV wrong map, CS::ChipAPV::Digit is not a CS::ChipAPV::DigitPixel." << std::endl;
      continue;
    }
    
    // same for DigitPixelMM
    const CS::ChipAPV::DigitPixelMM *mpmm = dynamic_cast<const CS::ChipAPV::DigitPixelMM*>(m);
    if ( fPlanePixMM && (mpmm==NULL) ) {
      std::cerr << "GeomPlanePMumega::SetCalibrations: ChipAPV wrong map, CS::ChipAPV::Digit is not a CS::ChipAPV::DigitPixelMM." << std::endl;
      continue;
    }

    flag  = fDefFlag;
    ped   = fDefPed;
    sigma = fDefSigma;

    const CS::ChipAPV::DataID &d = reinterpret_cast<const CS::ChipAPV::DataID &>(m->GetDataID());
    unsigned int adrid = d.u.s.adc_id*16+d.u.s.chip_id;
    // test if there are valid pedestals for current channel
//     if ( c.count(m->GetChip()) && c.find(m->GetChip())->second.channel.size()>m->GetChipChannel() ) {
//       const PlanePMumega::APVchannel a_ch = c.find(m->GetChip())->second.channel[m->GetChipChannel()];
    if ( c.count(adrid) && c.find(adrid)->second.channel.size()>m->GetChipChannel() ) {
      const PlanePMumega::APVchannel a_ch = c.find(adrid)->second.channel[m->GetChipChannel()];
      flag  = a_ch.flag;
      ped   = a_ch.ped;
      sigma = a_ch.sigma;
    } else { // use default and print a warning
      static bool printed=false;
      if (!printed)
        std::cerr << "Broken mapping or pedestals in calib for "<<GetName()<<" srcId "<<srcId
          <<" adr_id "<<adrid<<" chip "<<m->GetChip()<<" chipchannel "<<m->GetChipChannel()<< std::endl;
//       printed=true;
    }

    if (fPlane) fPlane->AddChan(m->GetChannel(), m->GetChanPos(),
                                flag, ped, sigma,
                                m->GetChip(), m->GetChipChannel());

    if (fPlanePix) {
      std::pair<int,int> xy=pixgem::detch2xy(mp->GetChannel());
      register int pixX = xy.first;  if (mp->GetDetOrientation() < 0) pixX=31-pixX;
      register int pixY = xy.second;
      fPlanePix->AddChan(mp->GetChannel(), pixX, pixY,
                         flag, ped, sigma,
                         mp->GetChip(), mp->GetChipChannel());
    }

    if (fPlanePixMM) {
      int pixelversion = mpmm->GetPixelMMversion();
      int conn = mpmm->GetConnNb();
      int chan = mpmm->GetChannel();
      register int chandet = chan + (pixmm::Instance(pixelversion)->GetNConn(conn))*MAX_APV_NBCH;
      int pixnb = pixmm::Instance(pixelversion)->GetPixNb(conn,chan);
      register float pixX = pixmm::Instance(pixelversion)->GetXPix(pixnb);
      register float pixY = pixmm::Instance(pixelversion)->GetYPix(pixnb);
      if ( mpmm->GetDetOrientation() < 0 ) {
        pixX = -pixX;
      }
// cerr<<GetName()<<": chandet "<<chandet<<" conn "<<conn<<" chan "<<chan<<" pixnb "<<pixnb<<" pixX "<<pixX<<" pixY "<<pixY<<" getchip "<<mpmm->GetChip()<<" chipchannel "<<mpmm->GetChipChannel()<<endl;
      fPlanePixMM->AddChan(chandet, pixnb, pixX, pixY,
                           flag, ped, sigma,
                           mpmm->GetChip(), mpmm->GetChipChannel());
    }

  }

  // time calibration from here on
  std::stringbuf buf;
  std::iostream  stream(&buf);

  // write the constants into a buffer
  for (register unsigned int i=0; i<ct.size(); i++)
    stream << ct[i] << " ";
  stream.seekg(0, std::ios_base::beg);

  // read calibration from buffer
  CsMumegaTimeCals cals; cals.Clear();
  stream >> cals;
  bool validfg = false;   
  if (fPlane) { fPlane->SetTimeCals(cals); validfg = fPlane->GetTimeCals()->IsValid(); }
  if (fPlanePix) { fPlanePix->SetTimeCals(cals); validfg = fPlanePix->GetTimeCals()->IsValid(); }
  if (fPlanePixMM) { fPlanePixMM->SetTimeCals(cals); validfg = fPlanePixMM->GetTimeCals()->IsValid(); }

  // if time calibration is not valid, use the standard
  if ( !validfg ) {
    CsMumegaTimeCalOld *cal02 = new CsMumegaTimeCalOld(0.2,1.1,-40.,22.,1.65);
    CsMumegaTimeCalOld *cal12 = new CsMumegaTimeCalOld(0.3,1.0,0.,27.,1.23);
    if (fPlane) fPlane->SetTimeCals(CsMumegaTimeCals(cal02, cal12));
    if (fPlanePix) fPlanePix->SetTimeCals(CsMumegaTimeCals(cal02, cal12));
    if (fPlanePixMM) fPlanePixMM->SetTimeCals(CsMumegaTimeCals(cal02, cal12));
  }
}

#else

void GeomPlanePMumega::SetCalibrations(const CS::Chip::Maps *maps) {
  // no database, no real calibrations can be loaded, load default ones
  // handle pedestal calibration first, time calibration in the ending
  // delete old calibrations
  if (fPlane) fPlane->ClearChannels();
  if (fPlanePix) fPlanePix->ClearChannels();
  if (fPlanePixMM) fPlanePixMM->ClearChannels();

  // need PixelMM.xml map to sort calibration data
  typedef CS::Chip::Maps::const_iterator m_it;
  
  // assume same source id for all apvs/channels from one PixelMM
  std::set<uint16> srcIDs;
  maps->GetSrcIDs( CS::DetID(GetName()), srcIDs );
  if ( srcIDs.size() != 1 )
    std::cerr << GetName() << ": One plane must be connected to only one source ID!" << std::endl;;

  unsigned int srcId = *srcIDs.begin();
  m_it m_low = maps->lower_bound(uint64(srcId)<<48);

  for( m_it cc=m_low; ((*cc).first>>48)==srcId; cc++ ) {
    // cc should now in any case be a CS::ChipAPV::Digit
    // (probably also a CS::ChipAPV::DigitPixel, but we do not care yet)
    register const CS::ChipAPV::Digit *m = dynamic_cast<const CS::ChipAPV::Digit*>((*cc).second);
    if ( m==NULL ) {
      std::cerr << "ChipAPV wrong map.\n";
      continue;
    }
    
    // check that detector in mapping is really this detector
    if (strcmp(m->GetDetID().GetName().c_str(), GetName()))
      continue;

    // now we can be more sure that this is a DigitPixel
    register const CS::ChipAPV::DigitPixel *mp = dynamic_cast<const CS::ChipAPV::DigitPixel*>(mt);
    if ( fPlanePix && (mp==NULL) ) {
      std::cerr << "ChipAPV wrong map. CS::ChipAPV::Digit is not a CS::ChipAPV::DigitPixel." << std::endl;
      continue;
    }


    if (fPlane) fPlane->AddChan(m->GetChannel(), m->GetChanPos(),
                                fDefFlag, fDefPed, fDefSigma,
                                m->GetChip(), m->GetChipChannel());

    if (fPlanePix) {
      std::pair<int,int> xy=pixgem::detch2xy(mp->GetChannel());
      register int pixX = xy.first;  if (mp->GetDetOrientation() < 0) pixX=31-pixX;
      register int pixY = xy.second;

      fPlanePix->AddChan(mp->GetChannel(), pixX, pixY,
                         fDefFlag, fDefPed, fDefSigma,
                         mp->GetChip(), mp->GetChipChannel());
    }

    if (fPlanePixMM) {
      int conn = mpmm->GetConnNb();
      int chan = mpmm->GetChannel();
      register int chandet = chan + conn*MAX_APV_NBCH;
      int pixnb = pixmm::GetPixNb(conn,chan);
      register float pixX = pixmm::GetXPix(pixnb);
      register float pixY = pixmm::GetYPix(pixnb);
      if ( mpmm->GetDetOrientation() < 0 ) {
        pixX = -pixX;
      }
      fPlanePixMM->AddChan(chandet, pixnb, pixX, pixY,
                           flag, ped, sigma,
                           mpmm->GetChip(), mpmm->GetChipChannel());
    }
  }

  // time calibration from here on
  CsMumegaTimeCalOld *cal02 = new CsMumegaTimeCalOld(0.2,1.1,-40.,22.,1.65);
  CsMumegaTimeCalOld *cal12 = new CsMumegaTimeCalOld(0.3,1.0,0.,27.,1.23);
  if (fPlane) fPlane->SetTimeCals(CsMumegaTimeCalOld(cal02, cal12));
  if (fPlanePix) fPlanePix->SetTimeCals(CsMumegaTimeCals(cal02, cal12));
  if (fPlanePixMM) fPlanePix->SetTimeCals(CsMumegaTimeCals(cal02, cal12));
}
#endif



void GeomPlanePMumega::SetParameters(float ts, float tc, int s) {
  if (fPlane) {
    fPlane->GetPar()->SetThrHit(ts);
    fPlane->GetPar()->SetThrClus(tc);
    fPlane->GetPar()->SetSample(s);
  }
  if (fPlanePix) {
    fPlanePix->GetPar()->SetThrHit(ts);
    fPlanePix->GetPar()->SetThrClus(tc);
    fPlanePix->GetPar()->SetSample(s);
  }
  if (fPlanePixMM) {
    fPlanePixMM->GetPar()->SetThrHit(ts);
    fPlanePixMM->GetPar()->SetThrClus(tc);
    fPlanePixMM->GetPar()->SetSample(s);
  }
}



void GeomPlanePMumega::Pad2Pos(double padX, double padY, double &posX, double &posY) {
  register double ctrX = GetNWires()/2. - .5;
  posX = (padX-ctrX) * GetPitch();

  register double ctrY = GetNWiresV()/2. - .5;
  posY = (padY-ctrY) * GetPitchV();
}

