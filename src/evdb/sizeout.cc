#include <CsPlatforms.h>
#include <CsTypes.h>

#include <CsSTD.h>

#include "CsDateReader.h"
#include "DaqDataDecoding/DaqEvent.h"

#include <fcntl.h>
#include <shift.h>

int main(int argc, char *argv[]) {

  cout << "sizeout " << CS_PLATFORM << " version" << endl;

  if(argc!=3) {
    printf("usage: eventout.exe mSize outFile\n"); 
    printf("mSize int number of kB\n"); 
    exit(1);
  }

  int mSize = atoi(argv[1]);
  char errString[100];
  char outFile[100];
  
  int iret = sprintf(outFile,"%s",argv[2]);
 
  CsDateReader* myInput = CsDateReader::Instance();

  int ia = myInput->decoderAttach("date.dat");

  if(ia != 1) {
    cerr << "Cannot open date.dat" << endl;
    exit(1);
  }

  int fh = rfio_open(outFile,O_WRONLY | O_CREAT | O_TRUNC, 0770);
  if( fh == -1 ) {
    //     if( fh == NULL ) {
    iret = sprintf(errString,"Error opening %s (fh=%d)",outFile,fh);
    sperror(errString);
    exit(1);
  }

  const int mbuf(50000);
  const int mbuf4 = mbuf*sizeof(uint8);
  uint8 *ibuf = new uint8[mbuf];

  int kwrit(0);
  int kread(0);

  int newSize = mSize * 1000; // kB conversion

  try {
    
    while(myInput->decoderGetBuffer()==1) {
      
      try {

	kread++;
	
	CS::DaqEvent myH(myInput->buffer(),true); //true=means can be modified
	
	int bl = myH.GetLength();
	const CS::uint32* buf;

	if(bl<newSize) {
	  myH.SetLength(newSize);
	  bl = myH.GetLength();
	}

	buf = myH.GetBuffer();
	iret = rfio_write(fh,(char*)buf, bl );

	kwrit++;

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

  rfio_close (fh);
  
  cout << "sizeout: " << kread << " events read" << endl; 
  cout << "sizeout: " << kwrit << " events written to " << outFile << endl; 
  
}

