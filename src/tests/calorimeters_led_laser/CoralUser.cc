#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "Reco/Calorimeter.h"
#include "Reco/CellDataRaw.h"

#include "DaqDataDecoding/DaqEvent.h"

#include "Coral.h"
#include "CsDet.h"
#include "CsOpt.h"
#include "CsCalorimeter.h"

#include "CoralUser.h"
#include "RunDataBase.h"
#include "StoreLED.h"
#include "StoreDeltaADC.h"
#include "DevMonitorFEM.h"

using namespace std;
using namespace Reco;
using namespace MN;

int run = 0;
void ReadDB( void );
void WriteDB( void );
bool ProcessLEDEvent( bool off_spill_event );
bool ProcessPhysicsEvent( void );

RunDataBase * data_base_write;
RunDataBase * data_base_read;
StoreDeltaADC *store_delta_adc;

unsigned int events_counter;
unsigned int calib_events_counter;
unsigned int physic_events_counter;
unsigned int led_events_counter;
unsigned int led_in_spill_events_counter;
bool            print_event_info;


// Initialization of DB
void CoralUserInit() {

  store_delta_adc=NULL;

  events_counter = 0;
  calib_events_counter = 0;
  physic_events_counter = 0;
  led_events_counter = 0;
  led_in_spill_events_counter = 0;

  data_base_write = NULL;
  data_base_read = NULL;
  string path_calib_db_write("DB_COMPASS_NEW");
  string path_calib_db_read("DB_COMPASS_OLD");
  CsOpt *opt = CsOpt::Instance(); // Instantiate the Option Interpreter
  //   ******************** DETECTOR LOCAL CALIB DB ********************
  std::string tag = ""; 
  std::string  key = "calib_db_read path";
  std::string str;
  if( opt->getOpt( tag, key, str ) ) path_calib_db_read = str;
  key = "calib_db_write path";
  if( opt->getOpt( tag, key, str ) ) path_calib_db_write = str;
  
  try
  {
    if( data_base_write == NULL )
    {
      data_base_write = new RunDataBase(path_calib_db_write.c_str(),
                                    RunDataBase::WRITE|RunDataBase::READ);
    }
    if( data_base_read == NULL )
    {
      data_base_read  = new RunDataBase(path_calib_db_read.c_str(),
                                    RunDataBase::WRITE|RunDataBase::READ);
    }
  }
  catch( std::exception &e )
  {
    cerr << " DB Init problems " << endl;
    exit(1);
  }
  catch( ... )
  {
    cerr << " DB Init problems " << endl;
    exit(1);
  }

   print_event_info = false; 
   tag = "";  key = "print event info"; int n;
   if( opt->getOpt( tag, key, n ) ) {
      print_event_info = true;
   }
   cout <<" print_event_info " << print_event_info << endl;
   vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();
   for(unsigned ci = 0; ci < calorimeters.size(); ci++) 
   {
     StoreLED *sl = new StoreLED(calorimeters[ci]);
     calorimeters[ci]->SetStoreLED(sl);

     std::string value;
     if( opt->getOpt( calorimeters[ci]->GetName(), "FEM_IS_IN_USE", value ) )
     {
       string value_nocase = value;
       bool fem_is_configured = false;
       for( size_t i=0; i< value.length(); i++ )
       value_nocase[i] = toupper(value[i]);
       if( value_nocase == "NO" )
         fem_is_configured= false;
       else
         fem_is_configured= true;
       value.clear();
       if( fem_is_configured )
       {
         DevMonitorFEM * fem = new DevMonitorFEM();
         calorimeters[ci]->SetFEM(fem);
         calorimeters[ci]->GetFEM2Modify()->Init( calorimeters[ci] );
         calorimeters[ci]->GetFEM2Modify()->InitReadOut();
       }  
     }
   }
   
   store_delta_adc = new StoreDeltaADC();
}

// --------------------------------------------------
// This function is called after new event reconstruction.
// Put you event analysis here...
void CoralUserEvent() {
//   bool debug = true;
  bool debug = false;
  debug = print_event_info;
//   static int counter = 0;
//   cout <<" CoralUserEvent **********************************  CoralUserEvent() cnt = " << counter << endl;
//   counter++;      
  events_counter++;
  static bool first = true;
  CsEvent* e = Coral::Instance()->getEvent();
  if(first) {
          run = e->getRunNumber();
          ReadDB();
          first = false;
  }
  const CS::DaqEvent& ddde = e->getDaqEvent();
//  Detectors cleaning is done in coral?? Need to clarify. For fem and probabbly other coral independent stuff we need to do it here or at the end of events processing.
// Ne ice!
//   if( fem_ != NULL ) fem_->Clear();

  if( debug ) cout <<" *******************************************************                Event Type " << ddde.GetHeader().GetEventType() << endl;
  if( debug ) cout <<" More about Times " << 
                                                 " Time first " << ddde.GetTime().first << endl;
// Is empty since ..? quite a while                                                           cout << " Time second " << ddde.GetTime().second  <<
//  At this point is not ready ... which is probabbly  can be changed             cout << " TTimeInSpill " << ddde.GetTT().GetTimeInSpill()  << endl;

  if(ddde.GetHeader().GetEventType() == CS::DaqEvent::CALIBRATION_EVENT) {
    calib_events_counter++;
    if( debug ) cout <<" This is calibration event " << endl;
	// LED event selection
    unsigned mask =0xF8000000;
    unsigned iand =  ddde.GetHeader().GetTypeAttributes()[0]&mask;
//    if( debug ) cout <<" This is calibration event  iand = " << iand << endl;
    if( debug ) printf(" This is calibration event %x  \n", iand);

    bool led_event_off_spill = ( iand == 0x68000000 );
    bool led_event_in_spill = ( iand == 0x60000000 );
    
    if( led_event_off_spill ) {
      led_events_counter++;
      bool led_event_processed = ProcessLEDEvent(led_event_off_spill);
      if( debug ) cout <<" This is OFF SPILL LED event  iand = " << iand <<
                                   " burst # " << ddde.GetHeader().GetBurstNumber() <<  
                                      " event # in burst " << ddde.GetHeader().GetEventNumberInBurst()  << endl;
    } else if( led_event_in_spill ) {
    // Process in spill LED/LASER calibartion event
      led_in_spill_events_counter++;
       bool led_event_in_spill_processed = ProcessLEDEvent(led_event_off_spill);
      if( debug ) cout <<" This is IN SPILL LED event  iand = " << iand <<
                                   " burst # " << ddde.GetHeader().GetBurstNumber() <<  
                                      " event # in burst " << ddde.GetHeader().GetEventNumberInBurst()  << endl;
    } else {                                                                                                   // What is the attribute ??  0xf80000000
    // Process calibartion but not LED/LASER calibartion event
      if( debug ) cout <<" This is NOT LED event  iand = " << iand <<
                                   " burst # " << ddde.GetHeader().GetBurstNumber() <<  
                                      " event # in burst " << ddde.GetHeader().GetEventNumberInBurst()  << endl;
    }
   if( debug ) cout <<" ProcessCalibrationEvents after DecodeEvent still ... about Times " << 
                                                 " Time first " << e->getDaqEvent().GetTime().first <<
                                                           " Time second " << e->getDaqEvent().GetTime().second  <<
                                                                                             "  Just TimeInSpill " << e->getTimeInSpill()  << 
                                                                                                  " TTimeInSpill " << e->getDaqEvent().GetTT().GetTimeInSpill()  <<
                                                                                                     " TTCSPhaseTime " << e->getDaqEvent().GetTT().GetPhaseTCS()  << endl;
  } else {
    physic_events_counter++;  
    bool physics_event_processed = ProcessPhysicsEvent();
     if( debug ) cout <<" ProcessPhysicsEvent Now after DecodeEvent still ... about Times " << 
                                                 " Time first " << e->getDaqEvent().GetTime().first <<
                                                           " Time second " << e->getDaqEvent().GetTime().second  <<
                                                                                             "  Just TimeInSpill " << e->getTimeInSpill()  << 
                                                                                                  " TTimeInSpill " << e->getDaqEvent().GetTT().GetTimeInSpill()  <<
                                                                                                     " TTCSPhaseTime " << e->getDaqEvent().GetTT().GetPhaseTCS()  << endl;
    if( debug ) cout <<" This is PHYSIC event  " << endl;
    // Process Physics event
  }
}

////////////////////////////////////////////////////////////////////////////////

bool ProcessLEDEvent( bool led_event_off_spill ) {
  bool debug = false;
//   bool debug = true;
//   static int led_counter = 0;
//   cout <<" ProcessLEDEvent **********  ProcessLEDEvent() cnt = " << led_counter << endl;
//   led_counter++;      
  vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();
  CsEvent* event = Coral::Instance()->getEvent();
  bool decoded = false;
  try {
    decoded = CsInit::Instance()->getDaqEventsManager().DecodeEvent();
  }
  catch(std::exception &e) {
    cerr << "ProcessLEDEvent: exception:\n" << e.what() << "\n";
    return( false ); 
  } 
  catch(std::string &e) {
    cerr << "ProcessLEDEvent: exception:\n" << e << "\n";
    return( false ); 
  } 
  catch(char *e) {
    cerr << "ProcessLEDEvent: exception:\n" << e << "\n";
    return( false ); 
  } 
  catch(...) {
    cerr << "ProcessLEDEvent: unknown exception\n";
    return( false );
  }
  if( debug ) cout <<"ProcessLEDEvent Now after DecodeEvent still ... about Times " << 
                                                 " Time first " << event->getDaqEvent().GetTime().first <<
                                                           " Time second " << event->getDaqEvent().GetTime().second  <<
                                                                                          "  Just TimeInSpill " << event->getTimeInSpill()  << 
                                                                                             " TTimeInSpill " << event->getDaqEvent().GetTT().GetTimeInSpill()  << 
                                                                                                 " TTCSPhaseTime " << event->getDaqEvent().GetTT().GetPhaseTCS()  << endl;

  double tcsphase = 40; //TCSPhaseTime();
  Coral::Instance()->getEvent()->setTCSPhaseTime(tcsphase);
// Coral time in spill stay not initialized in my case. This is bad.  
  
// TODO for in_spill LED events TCSPhaseTime probably make sence

  for(unsigned ci = 0; ci < calorimeters.size(); ci++) {
    calorimeters[ci]->SetTCSPhase(tcsphase);
    calorimeters[ci]->SetEventIDInfo();
    calorimeters[ci]->DecodeChipDigitsLEDEvent(event->getChipDigits());
    if( calorimeters[ci]->GetFEM2Modify() != NULL )
    {
       calorimeters[ci]->GetFEM2Modify()->Clear();
       calorimeters[ci]->GetFEM2Modify()->DecodeChipDigits(event->getChipDigits()); 
       bool all_fem_corrections_are_ok = true;
       for ( std::vector<Reco::CellDataRaw>::iterator it= calorimeters[ci]->GetSignals2Modify().begin(); it!= calorimeters[ci]->GetSignals2Modify().end(); it++)
       {
         unsigned icell = it->GetCellIdx();
         std::pair <bool,double> fem_corr = calorimeters[ci]->GetFEM()->GetCorrectionFactor(icell);
         if( !fem_corr.first )
         {
           all_fem_corrections_are_ok = false; 
           continue;
         }  
         if( fem_corr.second < 0.001 ) 
         {
           all_fem_corrections_are_ok = false;
           std::cerr <<" ApplyCorrectionsMonitorFEM " << calorimeters[ci]->GetName() <<" cell " << calorimeters[ci]->GetCellName(icell) << " Crazy FEM corrections " <<  fem_corr.second <<" Skipped !! " << std::endl;
           continue;
         } 
         double e =  it->GetAmplitude();
         it->SetAmplitude(e/ fem_corr.second );
       }
       if( all_fem_corrections_are_ok ) calorimeters[ci]->GetFEM2Modify()->TurnedOutGood4Corrections();
       calorimeters[ci]->GetFEM2Modify()->CalculateInternalCrossNormalization();
       calorimeters[ci]->GetFEM2Modify()->MonitorInSpill();
    }

    if( led_event_off_spill )
      calorimeters[ci]->ProcLED();

    calorimeters[ci]->GetStoreLED2Modify()->AddEvent( led_event_off_spill );
  }
//  cout <<" ProcessLEDEvent result " << decoded << endl;
//   return decoded;
  return true;
}

////////////////////////////////////////////////////////////////////////////////

bool ProcessPhysicsEvent() {
  bool debug = false;
//   bool debug = true;
  CsEvent* event = Coral::Instance()->getEvent();
  bool decoded = false;
  try {
    decoded = CsInit::Instance()->getDaqEventsManager().DecodeEvent();
  }
  catch(std::exception &e) {
    cerr << "ProcessPhysicsEvent: exception:\n" << e.what() << "\n";
    return( false ); 
  } 
  catch(std::string &e) {
    cerr << "ProcessPhysicsEvent: exception:\n" << e << "\n";
    return( false ); 
  } 
  catch(char *e) {
    cerr << "ProcessPhysicsEvent: exception:\n" << e << "\n";
    return( false ); 
  } 
  catch(...) {
    cerr << "ProcessPhysicsEvent: unknown exception\n";
    return( false );
  }
//   cout <<" ProcessPhysicsEvent DecodeEvent result " << decoded << endl;
// Actually as on 10.10.11 Coral ( CsEvent.cc ) ignore output of DecodeEvent function, .. we also ..  
  if( debug ) cout <<" ProcessPhysicsEvent Now after DecodeEvent still ... about Times " << 
                                                 " Time first " << event->getDaqEvent().GetTime().first <<
                                                           " Time second " << event->getDaqEvent().GetTime().second  <<
                                                                                             "  Just TimeInSpill " << event->getTimeInSpill()  << 
                                                                                                  " TTimeInSpill " << event->getDaqEvent().GetTT().GetTimeInSpill()  <<
                                                                                                     " TTCSPhaseTime " << event->getDaqEvent().GetTT().GetPhaseTCS()  << endl;


  FillDeltaSADC(event->getChipDigits());
  return true;
}

////////////////////////////////////////////////////////////////////////////////

void  FillDeltaSADC( const CS::Chip::Digits& digits)
{
//   bool debug = true;
  bool debug = false;
  if( debug ) cout << " FillDeltaSADC " << endl;
  vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();
  typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
  std::pair<m_it,m_it> m_range;
  int previous_src_id = -1;
  for(unsigned ci = 0; ci < calorimeters.size(); ci++) 
  {
    if( debug ) cout <<" Get Digits " << calorimeters[ci]->GetName() << endl;
    if( calorimeters[ci]->GetName() == "EC01P1__" )
    {
      m_range = digits.equal_range(  std::string("EC01P00") );
      for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
      {
        const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(d_it->second);
        if( ds == NULL ) continue;
        int src_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetSourceID();
        int port_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetPort();
        int chip_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetChip();
        if( src_id != previous_src_id )
        {
          if( debug ) cout <<" New SrcID " << src_id << endl;
          previous_src_id = src_id;
        }
        const std::vector<uint16>& samples = ds->GetSamples();
        std::pair<int,int> minimax = store_delta_adc->MiniMaxSADC(samples);
        size_t delta = minimax.second - minimax.first;
        store_delta_adc->UpdateDelta( std::pair<size_t,size_t>(src_id,chip_id), delta );
        if( debug ) cout<< " SrcID " <<  src_id  <<  " Port " <<  port_id  << " Chip " <<  chip_id  <<  " delta " <<  delta  << endl;
     }
      m_range = digits.equal_range(   std::string("EC01P01") );
      for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
      {
        const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(d_it->second);
        if( ds == NULL ) continue;
        int src_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetSourceID();
        int port_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetPort();
        int chip_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetChip();
        if( src_id != previous_src_id )
        {
          if( debug ) cout <<" New SrcID " << src_id << endl;
          previous_src_id = src_id;
        }
        const std::vector<uint16>& samples = ds->GetSamples();
        std::pair<int,int> minimax = store_delta_adc->MiniMaxSADC(samples);
        size_t delta = minimax.second - minimax.first;
        store_delta_adc->UpdateDelta( std::pair<size_t,size_t>(src_id,chip_id), delta );
        if( debug ) cout<< " SrcID " <<  src_id  <<  " Port " <<  port_id  <<  " Chip " <<  chip_id  << " delta " <<  delta  << endl;
      }
      m_range = digits.equal_range(   std::string("EC01P02") );
      for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
      {
        const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(d_it->second);
        if( ds == NULL ) continue;
        int src_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetSourceID();
        int port_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetPort();
        int chip_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetChip();
        if( src_id != previous_src_id )
        {
          if( debug ) cout <<" New SrcID " << src_id << endl;
          previous_src_id = src_id;
        }
        const std::vector<uint16>& samples = ds->GetSamples();
        std::pair<int,int> minimax = store_delta_adc->MiniMaxSADC(samples);
        size_t delta = minimax.second - minimax.first;
        store_delta_adc->UpdateDelta( std::pair<size_t,size_t>(src_id,chip_id), delta );
        if( debug ) cout<< " SrcID " <<  src_id  <<  " Port " <<  port_id  << " Chip " <<  chip_id  <<  " delta " <<  delta  << endl;
      }
    }
    else
    {

      m_range = digits.equal_range(  calorimeters[ci]->GetName() );
      for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
      {
        const CS::ChipSADC::Digit *ds = dynamic_cast<const CS::ChipSADC::Digit *>(d_it->second);
        if( ds == NULL )
        {
           const CS::ChipADC::Digit *dadc = dynamic_cast<const CS::ChipADC::Digit *>(d_it->second);
           if( dadc == NULL ) continue;
           int32 amp = dadc->GetAmplitude();
            int src_id = static_cast<CS::ChipADC::DataID>(dadc->GetDataID()).u.s.src_id;
            int geo_id = static_cast<CS::ChipADC::DataID>(dadc->GetDataID()).u.s.geo_id;
           store_delta_adc->UpdateDelta( std::pair<size_t,size_t>(src_id,geo_id), (size_t)amp );
           if( debug ) cout << " SrcID " <<  src_id  <<  " GeoID " <<  geo_id  <<  " ADC amp " <<  (size_t)amp  << endl;
           continue;
        }

        int src_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetSourceID();
        int port_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetPort();
        int chip_id = static_cast<CS::ChipSADC::DataID>(ds->GetDataID()).GetChip();
        if( src_id != previous_src_id )
        {
          if( debug ) cout <<" New SrcID " << src_id << endl;
          previous_src_id = src_id;
        }
        const std::vector<uint16>& samples = ds->GetSamples();
        std::pair<int,int> minimax = store_delta_adc->MiniMaxSADC(samples);
        size_t delta = minimax.second - minimax.first;
        if( calorimeters[ci]->IsItMSADCSrcID(src_id) )
          store_delta_adc->UpdateDelta( std::pair<size_t,size_t>(src_id,port_id), delta );
        else  
          store_delta_adc->UpdateDelta( std::pair<size_t,size_t>(src_id,chip_id), delta );
        if( debug ) cout<< " SrcID " <<  src_id  <<  " Port " <<  port_id  <<  " Chip " <<  chip_id  << " delta " <<  delta  << endl;
//        if( calorimeters[ci]->GetName() == "HC02P1__")  cout << " delta HC2 " <<  delta  <<" SrcID " << src_id <<   endl;
//        if( calorimeters[ci]->GetName() == "HC01P1__")  cout << " delta HC1 " <<  delta  <<" SrcID " << src_id <<  endl;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

namespace {

  class ForFileReading
  {
    public:
      enum {N=10000000};
      friend istream & operator >> (istream &,ForFileReading &c);

      char s[N];

  };

  istream & operator >> (istream &in,ForFileReading &c)
  {
    in.read(c.s,ForFileReading::N);
    if( in.gcount()>=ForFileReading::N )
    {
      cerr << " Fatal! Setup::ReadCalib(): too small internal buffer! needed " << in.gcount() <<
                                          " but only " << ForFileReading::N << " available " << endl;
      exit(1);                                        
    }


    c.s[in.gcount()]=0;
    in.clear();
    return in;
  }
}

// --------------------------------------------------
// Database files reading
void ReadDB( void ) {
  bool debug = false;
  if( debug ) cout << " reading DB for run " << run << endl;
  vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();
  for( unsigned i=0; i < calorimeters.size(); i++)
  {
    const std::vector< std::string > &calib_tags = calorimeters[i]->GetInputInfoNames();
    if( debug ) cout << " For Calorimeter " << calorimeters[i]->GetName() << 
                                       " Input list size = " << calib_tags.size() << endl;
    for( unsigned itc=0; itc < calib_tags.size(); itc++)
    {
      if( debug ) cout << " Try to read " << calib_tags[itc] << endl;
      
      ForFileReading *ar = new ForFileReading;
      ForFileReading &r = *ar;
      try
      {
// int CsCalorimeter::InputFEMInfo( size_t when, const string &s)
// {
//   if( fem_ !=NULL ) return fem_->InputFEMInfo( when, s);
//   return 0;
// }
        if( data_base_read->Read(string(calorimeters[i]->GetName()),calib_tags[itc],r,run) )
        {
          string s(r.s); 
          int ierr_code = 0;
          if( calib_tags[itc] == string("_TimeInSpillLED") )
          {
            ierr_code = calorimeters[i]->GetStoreLED2Modify()->InputTimeInSpillLED( s);
          }  
          else if( calib_tags[itc] == string("_FEMInfo" ) )
          {
            ierr_code =  calorimeters[i]->GetFEM2Modify()->InputFEMInfo(Reco::Calorimeter::OLD, s);
          }  
          else if( calib_tags[itc] == string("_FEMInfoRef" ) )
          {
            ierr_code =  calorimeters[i]->GetFEM2Modify()->InputFEMInfo(Reco::Calorimeter::PRIM, s);
          }
          else  
            ierr_code = calorimeters[i]->InputAnyCalibInfo( calib_tags[itc], s);

          if( ierr_code == 0 )
          {
            if( debug ) cout << " Reading of " << calib_tags[itc] << " is OK !!! " << endl;
          }
	  else
	  {
            cerr << " Problems with InputAnyCalibInfo " << calib_tags[itc] << " for " << 
	                                                    calorimeters[i]->GetName() << endl;
 	  }
        }
      }
      catch( const std::exception &e )
      {
        cerr << " ReadDB(): " << calorimeters[i]->GetName() << 
                                         " error in reading CALIB for run  " << run << endl;
        cerr << e.what() << "\n";
      }
      catch( const char *s )
      {
        cerr << s << "\n";
      }
      catch( ... )
      {
        cerr << "Setup::ReadCalib(): Unknown exception.\n";
      }
      delete ar;

    }
    calorimeters[i]->UpdateInternalAfterNewSettings();
  }

}
// --------------------------------------------------
// Database files saving
void WriteDB() {
//   bool debug = true;
  bool debug = false;
  if( debug ) cout <<" WriteDB debug " << endl;
  string s("");
  int run_start = run;
  int run_finish = run;
  vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();
  for( unsigned i=0; i < calorimeters.size(); i++)
  {
    if( debug ) cout <<" WriteDB Output " << calorimeters[i]->GetName() << endl;

    const std::vector< std::string > &calib_tags = calorimeters[i]->GetOutputInfoNames();
    for( unsigned itc=0; itc < calib_tags.size(); itc++)
    {
      s.clear();
      int ierr_code = 0;
      std::string comment("");
      if( calib_tags[itc] == "_LED" || calib_tags[itc] == "_LEDinSpills" || calib_tags[itc] == "_LEDinSpillsXY" )
      {
        if( calorimeters[i]->GetFEM() != NULL )
        {
           comment = calorimeters[i]->GetFEM()->GetSummaryComment();
        }
      }
           
      if( calib_tags[itc] == string("_TimeInSpillLED") )
      {
        ierr_code = calorimeters[i]->GetStoreLED()->OutputTimeInSpillLED( s);
      }  
      else if( calib_tags[itc] == string("_LEDinSpills") )
      {    
        if( !calorimeters[i]->XYInitialized() ) 
          ierr_code = calorimeters[i]->GetStoreLED2Modify()->OutputLEDInfoInSpills( s,comment );
        else  
          ierr_code = calorimeters[i]->GetStoreLED2Modify()->OutputLEDInfoInSpillsXY( s,comment );
      }  
      else if( calib_tags[itc] == string("_LEDinSpillsXY") )
      {
        ierr_code = calorimeters[i]->GetStoreLED2Modify()->OutputLEDInfoInSpillsXY( s,comment );
      }  
      else if(  calib_tags[itc] == "_FEMInfo" )
      {
        if( calorimeters[i]->GetFEM() != NULL )
          ierr_code = calorimeters[i]->GetFEM()->OutputFEMInfo(s);
      }  
      else if(  calib_tags[itc] == "_FEMInfoInSpills" )
      {
        if( calorimeters[i]->GetFEM() != NULL )
          ierr_code = calorimeters[i]->GetFEM()->OutputFEMInfoInSpills(s);
      }    
      else
      {
        ierr_code = calorimeters[i]->OutputAnyCalibInfo( calib_tags[itc], s,comment);
      }
      
      if( s.empty() )
      {
        cerr <<" Empty Output From " << calorimeters[i]->GetName() <<" tag " <<  calib_tags[itc] << endl;
        ierr_code = 1;
      }  

      if( ierr_code == 0 )
      {
        data_base_write->Write(string(calorimeters[i]->GetName()),calib_tags[itc],s,run_start,run_finish);
      }
    }
  }
  if( store_delta_adc != NULL )
  {
     s.clear();
     int ierr_code = store_delta_adc->OutputDeltaADC(s);
      if( ierr_code == 0 )
      {
        data_base_write->Write(string("DAQ"),string("_DeltaADC"),s,run_start,run_finish);
      }
  }
}
// --------------------------------------------------
void PrintStatistic() {
  cout << " Processed events " << events_counter <<" Calib.ev " << calib_events_counter <<
                               " LED.ev " << led_events_counter <<" LEDinSpill.ev " << led_in_spill_events_counter <<" PHYS.ev " << physic_events_counter << endl;
}
// --------------------------------------------------
// Database files saving
void CoralUserEnd() {
  
  PrintStatistic();
  vector<CsCalorimeter*> &calorimeters = CsGeom::Instance()->getCalorimeters();
  for( unsigned i=0; i < calorimeters.size(); i++)
  {
    calorimeters[i]->GetStoreLED2Modify()->Process();
    calorimeters[i]->EndOfJob();
  }
  WriteDB();
}
