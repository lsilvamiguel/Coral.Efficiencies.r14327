#include <fstream>
#include <cstdlib>
#include <csignal>

#include "DaqEventsManager.h"
#include "DaqEvent.h"
#include "ChipGassiplex.h"

using namespace CS;

// -----------------------------------------------------------------------------

void create_ChipGassiplex(ChipGassiplex &cg,uint16 srcID,uint16 event_number);

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

// -----------------------------------------------------------------------------

int main(int argc,char **argv)
{
  //--------------------------------------------------------------------------
  // Set signal handler to allow user to abort the program.
  (void) signal(0,      (void(*)(int))operation_system_signal); /* 0 */
  (void) signal(SIGHUP, (void(*)(int))operation_system_signal); /* 1 */
  (void) signal(SIGINT, (void(*)(int))operation_system_signal); /* 2 */
  (void) signal(SIGQUIT,(void(*)(int))operation_system_signal); /* 3 */
  (void) signal(SIGALRM,(void(*)(int))operation_system_signal); /* 4 */
  (void) signal(SIGTERM,(void(*)(int))operation_system_signal); /* 5 */
  //--------------------------------------------------------------------------

  try
  {
    switch( argc )
    {
      case 1:
      {
        uint16 srcID=1022, map_ver=1;
        Chip::Maps maps;
        maps.insert(pair<int,const Chip::Map*>(srcID,new ChipGassiplex::Map("0 RI01P 1226",map_ver,"runs=\"0\" chambers=\"7 6 5 4 3 2 1 0\"")));

        size_t event_number=11111;
        DaqEvent fake_event;
        ChipGassiplex cg;
        create_ChipGassiplex(cg,srcID,event_number);

        cg.Print(cout,"ChipGassiplex: ");
        //cg.Dump();

        Chip::Digits digits;
        cg.Decode(maps,digits);
        digits.Print();

        fake_event.Add(cg);
        fake_event.ReadChips();

        fake_event.Print(cout,"Event: ");
        break;
      }
      
      case 4:
      {
        size_t event_size = atoi(argv[1]);

        DaqEventsManager manager("DaqDataDecoding");
        manager.AddDataSource(argv[2]);
        ofstream f(argv[3]);
        if( !f.is_open() )
          throw Exception("Can not open file \"%s\"",argv[3]);

        try
        {
          while(!flag_end)
          {
            CS::DaqEvent event_orig = manager.GetEvent();
            CS::DaqEvent event(event_orig.GetBuffer(),true);

            ChipGassiplex cg;
            create_ChipGassiplex(cg,1111,event.GetEventNumberInRun());

            while( event.GetLength()<event_size )
              event.Add(cg);
            
            f.write(event.GetBuffer(),event.GetLength());
          }
          
        }
        catch( CS::DaqEvent::ExceptionEndOfStream )
        {
          cout << "End of data\n";
        }
        
        break;
      }
      
      default:
        cerr.form("Usage:\n  %s\n  %s <event-size-min>  <file-in.dat> <file-out.dat>\n\n",
                   argv[0],argv[0]);
    }
  }
  catch( const std::exception &e )
  {
    cerr << "exception:\n" << e.what() << "\n";
  }
  catch( const char * s )
  {
    cerr << "exception:\n" << s << "\n";
  }
  catch( ... )
  {
    cerr << "Oops, unknown exception!\n";
  }
}

void create_ChipGassiplex(ChipGassiplex &cg,uint16 srcID,uint16 event_number)
{
  cg.SetSourceID(srcID);
  size_t geoID=1222;
  const bool HEADER=false, TRAILER=true;

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddData(ChipGassiplex::Data(8,1,11,504,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,505,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,506,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,508,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,509,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,510,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,534,1011));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddData(ChipGassiplex::Data(8,1,11,504,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,432,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,216,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,  0,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,288,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,360,1011));
  cg.AddData(ChipGassiplex::Data(8,1,11,324,1011));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));

  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,HEADER,false));
  cg.AddHeaderTrailer(ChipCol::HeaderTrailer(event_number,geoID,TRAILER ,false));
}
