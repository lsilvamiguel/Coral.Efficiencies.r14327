#ifndef __Containers__
#define __Containers__

#include <vector>

class CDigit {
 public:
  CDigit() {}
  virtual ~CDigit() {}
};

class CDigit1 : public CDigit {
 public:
  int ch;
  std::vector<double> dt;
  CDigit1() : CDigit() {}
  CDigit1(int channel, const std::vector<double>& data) : CDigit(), ch(channel), dt(data)
    {}
  bool operator< (const CDigit1& other) const {
    return ch<other.ch;
  }
};

class CCluster {
 public:
  CCluster() {}
  virtual ~CCluster() {}
};

class CCluster1 : public CCluster {
 public:
  int id;
  double pos;
  int size;
  double res;
  std::vector<double> dt;
  std::vector<double> er;
  CCluster1() : CCluster() {}
  CCluster1(int pid, double position, int csize, double cres, const std::vector<double>& data) : 
      CCluster(), id(pid), pos(position), size(csize),res(cres),dt(data) {}
  CCluster1(int pid, double position, int csize, double cres, const std::vector<double>& data, const
                  std::vector<double>& error) : 
      CCluster(), id(pid), pos(position), size(csize),res(cres),dt(data), er(error) {}
  bool operator< (const CCluster1& other) const {
    return pos<other.pos;
  }
};

#endif




