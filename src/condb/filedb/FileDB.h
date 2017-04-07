#ifndef FileDB_h
#define FileDB_h

#include <sstream>
#include <string>

#include "CDB.h"
#include "Reco/DataBase.h"

class FileDB : public CDB{
 public:

  FileDB(std::string dbPath);

  void read(const std::string &folder, std::string &data, Time timePoint, const char* keyword=0);

  Reco::DataBase * getRecoDataBase(void) { return d;}

 protected:
  Reco::DataBase *d;

  class Calib {
  public:
    int ch;
    std::vector<float> data;
    Calib() : ch(0){}
    friend std::istream& operator>>(std::istream& in,FileDB::Calib &c) {
      char line[256];
      float tmp;
      in.getline(line,256,'\n');
      std::istringstream is(line);
      is>>c.ch;
      is>>tmp;
      do{
	c.data.push_back(tmp);
	is>>tmp;
      }while(!is.eof());
      return in;
    }
  };

  class String {
  public:
    std::string data;
    String(){}
    friend std::istream& operator>>(std::istream& in,FileDB::String &c) {
     char line[1000];
      in.getline(line,1000,EOF);
      c.data = line;
      if( in.eof() ) in.clear();
      return in;
    }

  };

};

#endif


