#ifndef __MMLibrary__
#define __MMLibrary__

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>

// this library contains some code necessary for Micromegas.

namespace MM {

  // a set of constants -------------------------------------------------
  class Constant {
  public:
    // amplitude = exp( tot/sTotnorm )
    // static const float sTotnorm; 
    static const int sRetainMaxHits; 
    static const int sRetainMaxHitsTime; 
  };


  // describes one line in the calibration files ------------------------
  class ChannelCalibExt {
  public:
    int ch;
    float t0;
    int flag;
    ChannelCalibExt() : ch(0),t0(0),flag(0) {}
    ChannelCalibExt(int ch,float t0,int flag) :ch(ch),t0(t0),flag(flag){}  
    
    friend std::istream& operator>>(std::istream& in,MM::ChannelCalibExt &c) {
    std::string line;
    std::getline(in, line);
    std::istringstream is(line);
    is >> c.ch;
    is >> c.t0;
    c.flag = 0;
    is >> c.flag;
    return in;
    }
  };
  
  // triangle cut in the tot vs t plane ----------------------------------
  //
  //           \  yes  /
  //            \     /
  //     no      \   /  no
  //              \ /
  //               x
  //          (fTpos,fApos)
  //
  //    fPp is the slope for the right line
  //    fPm                      left  line 
  //    fPm is POSITIVE !

  class ECutTriangle {
    
  private:
    double fTpos;   
    double fApos;
    double fPp;
    double fPm;
    double fHeight;
    
  public:
    ECutTriangle(double tpos, double apos, double pp, double pm, double height)
      : fTpos(tpos),fApos(apos),fPp(pp),fPm(pm), fHeight(height) {}
    
    ECutTriangle(const std::vector<double>& pars)
      : fTpos(pars[0]),fApos(pars[1]),fPp(pars[2]),fPm(pars[3]), 
      fHeight(pars[4]) {}
    
    bool IsInside(double t, double a) {
      return (a > (fTpos+t*fPp) )  && (a > (fTpos-t*fPm) );
    } 

    void Dump() {
      std::cout<<"ECutTriangle "<<fTpos<<" "<<fApos<<" "<<fPp<<" "<<fPm<<" "<<fHeight<<std::endl;
    }
    
    friend std::ostream& operator<<(std::ostream& out,ECutTriangle& cut) {
      if(!out) return out;
      out<<cut.fTpos<<"\t"<<cut.fApos<<"\t"<<cut.fPp<<"\t"<<cut.fPm<<"\t"<<cut.fHeight;
      return out;
    }
    
    friend std::istream& operator>>(std::istream& in,ECutTriangle& cut) {
      if(!in) return in;
      in>>cut.fTpos>>cut.fApos>>cut.fPp>>cut.fPm>>cut.fHeight;
      return in;
    }
  };


  // digit class ---------------------------------------------------------
  class Digit {
  public:
    int ch;       // digit channel
    double t;     // digit time
    double a;     // digit amplitude
    double tot;   // digit tot
    
    static int avstotmode; // which mode to calculate a vs Tot ?

    enum avstot_modes {STANDARD=0, REBOURGEARD, REBOURGEARD2, NEW};

    Digit(int ch, double t, double tot, bool largePitch);

    static void SetAvsTot(unsigned mode); 

    friend std::ostream& operator<<(std::ostream& out, Digit& d);
    
    ~Digit() {}
  };  


  // cluster class ---------------------------------------------------------
  class Cluster {
    
  protected:
    std::vector<Digit> fDigits;
    
  public:   
    double wire;   // uCluster (in chan num)
    double t;      // cluster time
    double tot;    // cluster tot
    double a;      // amplitude (sum of amplitudes of all digits in the cluster)
    int    s;      // cluster size
    
    Cluster(double wire, double t, double tot, int s);
    Cluster(const Cluster& c);
    Cluster(const std::vector<Digit>& v);
    virtual ~Cluster();
    
    const std::vector<Digit>& GetDigits() const {return fDigits;}

    void SetDigits(const std::vector<Digit>& v) {fDigits=v;}

    //void AddDigit(int ch, double time, double tot, bool largePitch) {
    //  fDigits.push_back( Digit(ch,time,tot,largePitch) );
    //}
    
    // splits the cluster and returns the products, including this
    std::vector<Cluster*> Split();
    
    // recalculates cluster parameters, looking at the digits.
    void Update();
    
    void Dump();
  };
}

#endif
