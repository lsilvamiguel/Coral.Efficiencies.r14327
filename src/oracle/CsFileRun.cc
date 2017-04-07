#if USE_FileStore

#include "CsFileRun.h"
#include "CsDstObject.h"
#include "CsDst_04.h"

using namespace std;

CsFileRun::CsFileRun(int runnum) : CsDstRun()
{
  runNumber = runnum;
  dstEvent = NULL;
  curChunkIdx = 0;
}

CsFileRun::~CsFileRun()
{
  dstFile.close();
  delete dstEvent;
}

bool CsFileRun::Init(std::vector<std::string>& file_names, std::string& chunkName)
{
  for(unsigned i = 0; i < file_names.size(); i++)
    {
      if(chunkName == "")
	fileNames.push_back(file_names[i]);
      else if(strstr(file_names[i].c_str(),chunkName.c_str()))
	fileNames.push_back(file_names[i]);
    }
  if(fileNames.size() == 0)
    {
      cerr << "CsFileRun::Init: No any chunk found for the run " << runNumber << endl;
      return false;
    }
  if(!OpenFile(fileNames[curChunkIdx]))
    {
      cerr << "CsFileRun::Init: Cannot open file: " << fileNames[curChunkIdx] << endl;
      return false;
    }
  uint32 chunkType = ERROR_TYPE;
  if(!dstFile.read((char*)&chunkType,sizeof(chunkType)))
    {
      cerr << "CsFileRun::Init: File reading error." << endl;
      return false;
    }
  dstFile.seek(0);
  switch (chunkType)
    {
    case DST_CHUNK_TYPE_04:
      dstChunk = new CsDstChunk_04;
      ((CsDstChunk_04*)dstChunk)->Init(chunkName.c_str());
      break;
    default:
      cerr << "CsFileRun::Init: Incorrect chunk type: " << chunkType << endl;
      return false;
    }
  if(!dstChunk->Load(dstFile))
    {
      cerr << "CsFileRun::Init: Chunk Information reading error" << endl;
      return false;
    }
  runLog = dstChunk->GetTBNames();
  return true;
}

bool CsFileRun::NextEvent()
{
  delete dstEvent;
  dstEvent = NULL;
  if(!dstFile.isOpen()) return false;
  dstEvent = dstChunk->GetNextEvent(dstFile);
  if(!dstEvent) return false;
  return true;
}

bool CsFileRun::SelectNextChunk(int, const std::string&)
{
  curChunkIdx++;

  if(curChunkIdx >= fileNames.size())
    {
      return false;
    }
  
  delete dstChunk;
  dstChunk = NULL;
  dstFile.close();

  if(!OpenFile(fileNames[curChunkIdx]))
    {
      cerr << "CsFileRun::SelectNextChunk: Cannot open file: " << fileNames[curChunkIdx] << endl;
      return false;
    }
  uint32 chunkType = ERROR_TYPE;
  if(!dstFile.read((char*)&chunkType,sizeof(chunkType)))
    {
      cerr << "CsFileRun::SelectNextChunk: File reading error." << endl;
      return false;
    }
  dstFile.seek(0);
  switch (chunkType)
    {
    case DST_CHUNK_TYPE_04:
      dstChunk = new CsDstChunk_04;
      ((CsDstChunk_04*)dstChunk)->Init("");
      break;
    default:
      cerr << "CsFileRun::Init: Incorrect chunk type: " << chunkType << endl;
      return false;
    }
  if(!dstChunk->Load(dstFile))
    {
      cerr << "CsFileRun::SelectNextChunk: Chunk Information reading error" << endl;
      return false;
    }
  return true;
}

bool CsFileRun::OpenFile(std::string& path)
{
  string dir;
  string name;
  string::size_type delim = path.find_last_of("/");
  if(delim == string::npos)
    {
      cerr << "CsFileRun::OpenFile: incorrect path: " << path << endl;
      return false;
    }
  name = path.substr(delim+1);
  dir  = path.substr(0,delim+1);
  if(!dstFile.open(dir,name))
    {
      cerr << "CsFileRun::OpenFile: cannot open file: " << dir << name << endl;
      return false;
    }
  return true;
}

#endif //#if USE_FileStore
