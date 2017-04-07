#include <CsPlatforms.h>
#include <CsTypes.h>

#include <CsSTD.h>

#include "CsDateReader.h"
#include "DaqDataDecoding/DaqEvent.h"

int main(int argc, char *argv[]) {

  cout << "dateReader " << CS_PLATFORM << " version" << endl;

  int karg = argc-1;

  if(karg!=1) {
    cout << "usage: dateReader -option" << endl;
    cout << "       -verbose (full dump of each event)" << endl;
    cout << "       -scan    (quick dump of each dump)" << endl;
    cout << "       -summary (summary table only)" << endl;
    cout << "Note: link the date file as date.dat" << endl; 
    exit(1);
  }

//  char opt[100];
//  int iret = sprintf(opt,"%s",argv[1]);

  string option = argv[1];

  CsDateReader* myInput = CsDateReader::Instance();

  int ia = myInput->decoderAttach("date.dat");

  if(ia != 1) {
    cerr << "Cannot open date.dat" << endl;
    exit(1);
  }

  const int mbuf(50000);
  const int mbuf4 = mbuf*sizeof(uint8);
  uint8 *ibuf = new uint8[mbuf];

  int kevt(0);
  int kread(0);

  int eventTypeSize[11] = {0,0,0,0,0,0,0,0,0,0,0};
  int eventTypeSummary[11] = {0,0,0,0,0,0,0,0,0,0,0};
  string eventTypeLabel[11];
  eventTypeLabel[ 0] = "Not implemented (should be 0)";
  eventTypeLabel[ 1] = "Start Of Run";
  eventTypeLabel[ 2] = "End Of Run";
  eventTypeLabel[ 3] = "Start Of Run file";
  eventTypeLabel[ 4] = "End Of Run file";
  eventTypeLabel[ 5] = "Start Of Burst";
  eventTypeLabel[ 6] = "End Of Burst";
  eventTypeLabel[ 7] = "Physics Event";
  eventTypeLabel[ 8] = "Calibration event";
  eventTypeLabel[ 9] = "End of Link";
  eventTypeLabel[10] = "Date Error";

  unsigned long minTime(UINT_MAX);
  unsigned long maxTime(0);

  try {
    
    while(myInput->decoderGetBuffer()==1) {
      
      try {

	kread++;
	
	CS::DaqEvent myH(myInput->buffer());
	
	int i = myH.GetType();
	if(i>10) i=0;
	if(i<=0) i=0;
	eventTypeSummary[i]++;
	eventTypeSize[i] += myH.GetLength();
	
	unsigned long it = myH.GetTime().first;
	if (it>maxTime) { maxTime = it; }
	if (it<minTime) { minTime = it; }

	if (option != "-summary" ) {
	  cout << " Event " << kevt << " length " << 
	    myH.GetLength() << " bytes" << endl;
	 
	}

	if (option == "-verbose") { 
	  
	  myH.Print();
	  
	  const int nPerPage(10);
	  
	  int jmod(nPerPage);
	  ibuf = myInput->buffer();
	  int iline(0);
	  cout << setfill('0') << setw(4) << iline << ":";
	  int id;
	  for (id=0;id<myH.GetLength();id++) {
	    cout << setfill(' ') << setw(8) << ((int32) ibuf[id]);
	    if(id%jmod==nPerPage-1) { 
	      iline += sizeof(uint8); 
	      cout << endl; 
	      cout << setfill('0') << setw(4) << iline << ":";
	    }
	  }
	  if(id%jmod!=nPerPage-1) { cout << endl; }
	}
      
	kevt++;
      }
      catch( const CS::Exception &e )
	{
	  cerr << "exception:\n" << e.what() << "\n";
	}
    }
  }
  catch( CS::DaqEvent::ExceptionEndOfStream )
    {
      // This is the normal exit from the loop
      cout << "End of file.\n";
    }
  catch( const char * s )
    {
      cerr << "exception:\n" << s << "\n";
    }
  catch (const std::exception &e )
    {
      cerr << "exception:\n" << e.what() << "\n";
      return 1;  // return failure
    }
  catch( ... )
    {
      cerr << "Oops, unknown exception!\n";
    }
  
  cout << "dateReader: decoderStatus " << myInput->decoderGetBuffer() << endl; 
  cout << "dateReader: " << kevt << " events read" << endl; 
  
  if (option == "-verbose" || option == "-summary") { 
    
    cout << "dateReader: Event summary " << endl;
    cout << "Input events from file" << " = " << kread << endl;
    for (int i=0;i<11;i++) {
      cout << eventTypeLabel[i] << " = " << eventTypeSummary[i] << endl;
    }
    cout << "Mean sizes " << endl;
    for (int i=0;i<11;i++) {
      int kbs = eventTypeSize[i];
      if(eventTypeSummary[i]>0) {
	kbs = kbs/eventTypeSummary[i];
	cout << eventTypeLabel[i] << " = " << kbs << " Bytes" << endl;
      }
    }
    cout << "Min/Max time " << endl;
    cout << minTime << " " << maxTime << endl;

    
    unsigned int bs = 0;
    for(int i=0;i<11;i++) {
      bs += eventTypeSize[i];
    }
    cout << "Total size read " << bs << " Bytes"<< endl;

  }
}
