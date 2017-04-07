#include <cstdio>
#include <iostream>             // C++ I/O
#include <fstream>              // C++ files I/O
#include <csignal>
#include <cstdlib>
#include "utils.h"
#include "DaqEvent.h"           // DaqEvent class interface
#include "Chip.h"               // Interface to data from COMPASS chips
#include "ChipF1.h"
#include "ChipADC.h"
#include "ChipAPV.h"

using namespace std;

// -----------------------------------------------------------------------------

static bool flag_end=false;

void decode(CS::Chip&);

// -----------------------------------------------------------------------------

void operation_system_signal(int n)
{
  cerr << endl
       << "============================================" << endl
       << "=== The program has received signal: %3d ===" << endl
       << "============================================" << endl << endl;
  if( flag_end )
  {
    cerr << "Forcing exit.\n\n";
    exit(1);
  }
  else
    flag_end = true;
}

// -----------------------------------------------------------------------------

// Program start point
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

  try // This is try-catch block for catching C++ exceptions
  {
    // Test input arguments...
    if( argc!=2 )
      throw "Usage:  <EXE-file> <data file name>";
    
    // Open file
    ifstream f(argv[1]);
    if( !f.is_open() )
      throw "Can not open file";
    
    // Loop over all DATE events in the stream
    while(!flag_end)
      CS::DaqEvent(f).ReadChips(decode);
  }
  catch( CS::DaqEvent::ExceptionEndOfStream )
  {
    // This is the normal exit from the loop
    cout << "End of file.\n";
  }

  // Something is wrong...
  // Print error message and finish the program.

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

  CS::Chip::PrintStatistics();

  return 0;
}

// -----------------------------------------------------------------------------

void decode_ChipAPV(const CS::ChipAPV& chip);
void decode_ChipADC(const CS::ChipADC& chip);
void decode_ChipF1 (const CS::ChipF1 & chip);

void decode(CS::Chip& chip)
{
  chip.Print();
  
  const CS::ChipF1  *chip_F1  = dynamic_cast<const CS::ChipF1 *>(&chip);
  if( chip_F1!=NULL )
    decode_ChipF1(*chip_F1);

  const CS::ChipADC *chip_ADC = dynamic_cast<const CS::ChipADC*>(&chip);
  if( chip_ADC!=NULL )
    decode_ChipADC(*chip_ADC);

  const CS::ChipAPV *chip_APV = dynamic_cast<const CS::ChipAPV*>(&chip);
  if( chip_APV!=NULL )
    decode_ChipAPV(*chip_APV);
}

// -----------------------------------------------------------------------------

void decode_ChipF1(const CS::ChipF1& chip)
{
}

// -----------------------------------------------------------------------------

void decode_ChipADC(const CS::ChipADC& chip)
{
}

// -----------------------------------------------------------------------------

void decode_ChipAPV(const CS::ChipAPV& chip)
{
}

// -----------------------------------------------------------------------------
