#include <csignal>
#include <iostream>
#include <fstream>
#include <strstream>
#include <exception>

#include "DaqEvent.h"
#include "DaqEquipment.h"
#include "ChipF1.h"
#include "ChipADC.h"

#include "Reco_config.h"
#include "GUICalorimeter.h"
#include "CalorimeterGAMS.h"
#include "DataBase.h"

#if USE_Qt
#include <qapplication.h>
#include <qwidget.h>
#endif

#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
TROOT ROOT("","");

//#include "Arguments.h"
#include "monitoring/monitor.h"

////////////////////////////////////////////////////////////////////////////////

// Signals handler (to catch Ctrl-C during program run)

static bool flag_end=false;

void operation_system_signal(int n)
{
  cerr.form("\n"
            "============================================\n"
            "=== The program has received signal: %3d ===\n"
            "============================================\n\n",n);
  if( flag_end )
  {
    cerr << "Forcing exit.\n\n";
    exit(1);
  }
  else
    flag_end = true;
}

////////////////////////////////////////////////////////////////////////////////

using namespace CS;
using namespace Reco;

int main(int argc,char **argv)
{
  #if USE_Qt
    QApplication app_Qt(argc,argv );
  #endif
  
    TApplication app_ROOT("gui",&argc,argv);

  //--------------------------------------------------------------------------
  // Set signal handler
  (void) signal(0,      (void(*)(int))operation_system_signal); /* 0 */
  (void) signal(SIGHUP, (void(*)(int))operation_system_signal); /* 1 */
  (void) signal(SIGINT, (void(*)(int))operation_system_signal); /* 2 */
  (void) signal(SIGQUIT,(void(*)(int))operation_system_signal); /* 3 */
  (void) signal(SIGALRM,(void(*)(int))operation_system_signal); /* 4 */
  (void) signal(SIGTERM,(void(*)(int))operation_system_signal); /* 5 */
  //--------------------------------------------------------------------------

  if( sizeof(DaqEquipment::DaqSLink)==3 )
  {
    cerr << "The problem!  sizeof(DaqSLink)=" << sizeof(DaqEquipment::DaqSLink) << "  (must equal to 3)\n";
    return 1;
  }

  static ofstream log("log");
  clog=log;

//  Arguments arg(argc,argv);
//  bool flag_print_event = arg.receivedArgument("--print-event");
//  bool flag_print_equipment = arg.receivedArgument("--print-equipment");

  if( argc!=2 )
  {
    cerr << "Bad input arguments. Usage:\n"
         << "  " << argv[0] << " <data_source>\n";
    return 1;
  }

  DataBase db("DB",DataBase::WRITE|DataBase::READ|DataBase::CREATE);

  CalorimeterGAMS GAMS("GAMS",2000);
//  GAMS.GetInfoFromDataBase(db);

  #if USE_Qt
  GAMS.CreateGUI();
  GAMS.GetGUI().show();
  #endif

//   new TCanvas("c","cccccccccccccccc",444,333); // Open ROOT canvas

  int status;
  
  status = monitorSetDataSource(argv[1]);
  if( status )
  {
    cerr << "monitorSetDataSource: " << monitorDecodeError(status) << "\n";
    exit(1);
  }

  bool flag_print_event=0, flag_print_equipment=0;
  
  char *monitoring_policy[] =
  {
//     "All_events",              "yes",
//     "Start_of_run",            "yes",
//     "Start_of_run_files",      "yes",
//     "Start_of_burst",          "yes",
//     "Calibration_event",       "yes",
     "Physics_event",           "yes",
//     "Event_format_error",      "yes",
//     "End_of_burst",            "yes",
//     "End_of_run_files",        "yes",
//     "End_of_link",             "yes",
//     "End_of_run",              "yes",
    0
  };
  
  status = monitorDeclareTable(monitoring_policy);
  if( status )
  {
    cerr << "monitorDeclareTable: " << monitorDecodeError(status) << "\n";
    exit(1);
  }

  while(!flag_end)
  {
    #if USE_Qt
      app_Qt.processEvents();
    #endif
    
      gSystem->ProcessEvents();

    void *event_buffer=NULL;

    status = monitorGetEventDynamic(&event_buffer);
    if( status )
      cerr << "monitorGetEventDynamic: " << monitorDecodeError(status) << "\n";

    if( event_buffer==NULL )
      break;

    try
    {
      vector<CellDataRaw> cal_data;

      CS::DaqEvent event(event_buffer);
      if( flag_print_event )
        event.Print();
      event.ReadEquipment();

      for( vector<DaqEquipment*>::iterator it=event.GetEquipment().begin(); it!=event.GetEquipment().end(); it++ )
      {
        if( flag_print_equipment )
          (**it).Print();

        try
        {
          void ADC(const DaqEquipment &eq, vector<CellDataRaw> &cal_data, const vector<CellDataRaw> *PED=NULL);
          ADC(**it,cal_data);
        }
        catch( const std::exception &e )
        {
          cerr << "Exception:\n" << e.what() << endl;
        }
        catch( const char *s )
        {
          cerr << "Exception:\n" << s << endl;
        }
        catch( ... )
        {
          cerr << __FILE__ << " in " << __FUNCTION__ << " at " << __LINE__ << "Unknown exception.\n";
        }
      }
      
    }
    catch( const std::exception &e )
    {
      cerr << "Exception:\n" << e.what() << endl;
    }
    catch( const char *s )
    {
      cerr << "Exception:\n" << s << endl;
    }
    catch( ... )
    {
      cerr << __FILE__ << " in " << __FUNCTION__ << " at " << __LINE__  << "Unknown exception.\n";
    }

    free(event_buffer);
  }

  DaqEquipment::PrintStatistics();

  GAMS.PutInfoToDataBase(db);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void ADC(const DaqEquipment &eq, vector<CellDataRaw> &cal_data, const vector<CellDataRaw> *ped)
{
  const size_t N_hist=1024;
  static TH1F *h[N_hist];
  static bool init=true;
  if( init )
  {
    init=false;
    for( size_t i=0; i<N_hist; i++ )
    {
      strstream name;
      name.form("h%5.5d",i);
      h[i] = new TH1F(name.str(),string(string("Title: ")+name.str()).c_str(),1000,0,5000000);
    }
  }

  const ChipADC *chip = dynamic_cast<const ChipADC*>(&eq);
  if( chip==NULL )
    return;

  while(true)
  {
    const uint32 *d=NULL;
    try
    {
      d=eq.GetNextData();
    }
//     catch( ChipADC::FatalError &e )
//     {
//       // Fatal error in this equipment, I will pass it to a caller method.
//       throw;
//     }
    catch( std::exception &e )
    {
      // This data line is incorrect.
      // We will print error message (with problem description) and
      // will go to next line
      cerr << e.what() << endl;
      continue;
    }
    catch( ... )
    {
      // Unknown error. I will pass it to a caller method.
      throw;
    }

    if( d==NULL )
      break;
  
    ChipADC::HeaderTrailer &dd = *(ChipADC::HeaderTrailer*)d;
    if( dd.is_data )
    {
      ChipADC::Data &data = *(ChipADC::Data*)d;

      int ped_value=0;
      if( ped!=NULL )
      {
        register int address = chip->GetDataSourceID(*d);
        for( size_t i=0; i<ped->size(); i++ )
          if( (*ped)[i].GetAddress()==address )
          {
            ped_value = (int) (*ped)[i].GetAmplitude();
            break;
          }
      }
      cal_data.push_back(CellDataRaw(chip->GetDataSourceID(*d),data.GetData()-ped_value));

//      if( chip->GetGeoID()>=0 && chip->GetGeoID()<N_hist )
//        h[chip->GetGeoID()] -> Fill(data);

//       if( chip->GetGeoID()==-1 )
//         cout << "ADC data: geoID is unknown\n";
    }
    else
      if( !dd.is_trailer )
      {
//        cout << "geoID = " << chip->GetGeoID() << "\n";
//        cout << "ADC: header\n";
      }
      else
      {
//        cout << "ADC: trailer\n";
      }
  }
}

////////////////////////////////////////////////////////////////////////////////
