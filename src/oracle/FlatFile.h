#ifndef __FLAT_FILE_H__
#define __FLAT_FILE_H__

#include <string>
#include <unistd.h>

typedef off_t	FileOffset;
typedef int	FlatFileID;


class FlatFile
{
 public:
  static const FileOffset  _maxFileSize = 1024U*1024*1024*2-1;

  FlatFile();
  virtual ~FlatFile();

  std::string	Path() const;
  bool          isOpen() const { return _fh>0; }

  // Create a new empty file (overwrite old)
  bool Create(const std::string &dir, const std::string &fname );
  bool Open(const std::string &dir, const std::string &name );
  bool Read(void *buf, FileOffset length, FileOffset pos);
  bool Read(void* buf, FileOffset length);
  bool Write(void* buf, FileOffset length);
  bool Seek(FileOffset pos);
  bool Append(void *buffer, const FileOffset length );
  bool Close();
  bool Discard();  // close and delete the file
  FileOffset getFileOffset() const { return _offset; }
  std::string getDir() const { return _basedir; }
  std::string getFileName() const { return _fname; }

 protected:
  // Castor open with loop on EBUSY
  bool _open(const char *file, unsigned flags, unsigned bits);

  std::string	_basedir;
  std::string	_fname;
  int		_fh;
  FileOffset	_size;
  FileOffset    _offset;

  int		_openErrorRetrySleepTime;
  int		_EINTR_OpenRetries;
  int		_EBUSY_OpenRetries;
  int		_ENOSPC_OpenRetries;
};

#endif
