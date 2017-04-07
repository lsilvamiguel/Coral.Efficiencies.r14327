#include <errno.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>

#include <cstdio>
#include <iostream>

#include "coral_config.h"
#if USE_RFIO
#  include <shift.h>
#else
#  define RFIO_READOPT 0
#  define RFIO_STREAM 0
#  define rfiosetopt(a,b,c) 0
#  define rfio_serror() ""
   int serrno = 0;
   int rfio_errno = 0;
#endif

#include "FlatFile.h"
#include "OracleErrorCode.h"

using namespace std;

FlatFile::FlatFile()
  : _fh(-1), _size(0), _openErrorRetrySleepTime(300),
    _EINTR_OpenRetries(5), _EBUSY_OpenRetries(1)
{
  // 100 hours of retry
  _ENOSPC_OpenRetries = 100*3600/_openErrorRetrySleepTime;
  _offset = 0;

  int vNeg = RFIO_STREAM; // Negotiation flag: RFIO_STREAM for _V3
  if ( rfiosetopt(RFIO_READOPT, &vNeg, sizeof(vNeg)) != 0 )
    cerr << "FlatFile::Create(): rfiosetopt failed." << endl;
}

FlatFile::~FlatFile()
{
  if( isOpen() )
    {
      if(!Close())
	{
	  throw RFIO_FILE_CLOSE_ERROR;
	}
    }
}

string	FlatFile::Path() const
{
  string path = _basedir + "/" + _fname;
  string::size_type idx = path.find("//", 0, strlen("//"));
  while(idx!=string::npos)
    {
      path.erase(idx, 1);  // remove spurious "/"
      idx = path.find("//", 0, strlen("//"));
    }
  return path;
}

bool FlatFile::Create(const string &dir, const string &name )
{
  _basedir = dir;
  _fname = name;
  cout<<"Create new file: " << Path() << endl;
  return _open(Path().c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0664);
}

bool FlatFile::Open(const string &dir, const string &name) 
{
  _basedir = dir;
  _fname = name;
  return _open(Path().c_str(), O_RDONLY, 0);
}

bool FlatFile::Read(void *buf, FileOffset length, FileOffset pos)
{
  if(!Seek(pos)) return false;
  if( read(_fh, buf, length) != (int)length )
    {
      cerr <<  "FlatFile::Read(): Error reading from " << Path() << " offset: " 
	   << (unsigned)pos << "\n" << rfio_serror() 
	   << "\nserrno = " << serrno << ", errno=" << errno << endl;
      return false;
    }
  _offset += length;
  return true;
}

bool FlatFile::Append(void *buffer, const FileOffset length )
{
  if( !isOpen() )
    {
      cerr << "FlatFile::Append(): Attempt to write to a file that is not open." << endl;
      return false;
    }
  if( _size + length > _maxFileSize )
    {
      cerr << "Error writing to " << Path() << ", maximum file size exceeded." << endl;
      return false;
    }
  int	sleepCounter = 0;
  FileOffset iret = ::write(_fh, buffer, length);
  if( iret != length )
    {
      // RFIO errors are final, no retry possible
      cerr << "Error writing to " <<  Path() 
	   << "\n" << rfio_serror() 
	   << "\nserrno = " << serrno << ", errno=" << errno << endl;
      return false;
    }
  _size += length;
  _offset = _size;
  return true;
}

bool FlatFile::Close()
{
  if( !isOpen() )
    {
      cerr << "Attempt to close file " << Path() << " which is not open" << endl;
      return true;
    }
  int rc = close( _fh );
  if(rc == 0)
    {
      cout << "File: " << Path() << " is closed successfully." << endl;
      _fh = -1;
      return true;
    }
  cerr << "FlatFile::Close(): Error of closing file " << Path() << "\n"
       << rfio_serror() 
       << "\nserrno = " << serrno << ", errno=" << errno << endl;
  return false;
}

bool FlatFile::Discard()
{
  if(Close())
    {
      cout << "Cleaning up file " << Path() << endl;
      if( unlink(const_cast<char*>(Path().c_str()) ))
	{
	  cerr << "FlatFile::Discard(): cannot remove file " << Path() << "\n"
	       << rfio_serror() 
	       << "\nserrno = " << serrno << ", errno=" << errno << endl;
	  return false;
	}
    }
  return true;
}

bool FlatFile::_open(const char *file, unsigned flags, unsigned bits)
{
  rfio_errno = serrno = errno = 0;
  int  sleep_ENOSPC_count = 0;
  int  sleep_EINTR_count = 0;
  int  sleep_EBUSY_count = 0;
  bool retry = false;
  char local_file[1024];
  strcpy(local_file,file);
  char* physical_path = &local_file[0];
  while(1)
    {
      _fh = open(physical_path, flags, bits);
      if( _fh >= 0 )
	{
	  cout << "done!" << endl;
	  break;
	}
      int error =  rfio_errno? rfio_errno : ( serrno? serrno : errno );
      if( error == EINTR )
	{
	  cout << "FlatFile::_open(): EINTR while opening file: " << Path() << endl;
	  retry = (sleep_EINTR_count++ < _EINTR_OpenRetries);
	}
      else if( error == EBUSY )
	{
	  cout << "FlatFile::_open(): CASTOR file busy while opening: " 
	       << Path() << endl;
	  retry = (sleep_EBUSY_count++ < _EBUSY_OpenRetries);
	}
      else if( error == ENOSPC )
	{
	  cout << "FlatFile::_open(): CASTOR: NO SPACE while opening: " 
	       << Path() << endl;
	  retry = (sleep_ENOSPC_count++ < _ENOSPC_OpenRetries);	    
	}
      if(!retry)
	{
	  perror("FlatFile::_open()");
	  cerr << "FlatFile::_open(): Error opening " << Path() << "\n"
	       << rfio_serror() 
	       << "\nserrno = " << serrno << ", errno=" << errno << endl;
	  return false;
	}
      cout << "FlatFile::_open(): Sleeping for " << _openErrorRetrySleepTime 
	   << "s." << endl;
      sleep(_openErrorRetrySleepTime);
      cout << "FlatFile::_open(): Retrying open" << endl;
    }
  return true;
}

bool FlatFile::Read(void* buf, FileOffset length)
{
  return Read(buf,length,_offset);
}

bool FlatFile::Write(void* buf, FileOffset length)
{
  return Append(buf,length);
}

bool FlatFile::Seek(FileOffset pos)
{
  if(_offset == pos) return true;
  if( lseek(_fh, pos, SEEK_SET) < 0 )
    {
      cerr << "FlatFile::Seek(): Error positioning file " << Path()
	   <<" to offset " << (unsigned)pos << "\n" << rfio_serror()
	   << "\nserrno = " << serrno << ", errno=" << errno << endl;
      return false;
    }
  _offset = pos;
  return true;
}
