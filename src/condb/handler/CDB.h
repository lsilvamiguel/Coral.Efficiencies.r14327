#ifndef CDB_h
#define CDB_h

class CDB {

 public:
  typedef std::pair<time_t,unsigned> Time;

 public:
  virtual ~CDB(void){};

          void read(const std::string &folder, std::string &data, Time tPoint, const std::string &keyword) {
	    read(folder, data, tPoint, keyword.c_str());
	  }

  virtual void read(const std::string &folder, std::string &data, Time tPoint, const char* keyword=0) = 0;

  virtual bool   ConnectDB(void)        { return true; }
  virtual void   DisconnectDB(void)     {}
  virtual bool   isConnected(void)      { return true; }
};

#endif
