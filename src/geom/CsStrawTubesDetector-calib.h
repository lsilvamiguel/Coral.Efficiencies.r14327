#ifndef CsStrawTubesDetector_Calib_inc
#define CsStrawTubesDetector_Calib_inc

#include <string>
#include <map>
#include <vector>

namespace CS
{

class CsStrawTubesDetector_Calib
{
  public:
                        CsStrawTubesDetector_Calib  (void) {}
    virtual            ~CsStrawTubesDetector_Calib  () {}
    virtual void        read                        (const std::string &st) {}
};

class CsStrawTubesDetector_CalibXray: public CsStrawTubesDetector_Calib
{
  public:
    void read(const std::string &st);

    /// Return the X-ray correction in [mm]
    float xray_correction(int channel,float y) const;
    void make_pictures(const std::string &stdc);

  public:
    std::vector<float> spacers_pos;      // coordinates of spacers
    std::map<int,std::vector<float> >  xray;  // Channel number --> List of spacer corrections.
};

class CsStrawTubesDetector_CalibT0: public CsStrawTubesDetector_Calib
{
  public:

    class ChanPos
    {
      public:
        ChanPos (int c=-1,int p=0) : chan(c), pos(p) {}
        bool operator() (const ChanPos &cp1,const ChanPos &cp2 ) const
        {
            if(cp1.chan!=cp2.chan )
                return cp1.chan<cp2.chan;
            else
                return cp1.pos < cp2.pos;
        }
        
        int chan,pos;
    };

  public:

    void read(const std::string &st);

    /// Return channel T0, or throw an exception
    float GetT0(int channel,int pos=0) const;

  public:

    std::map<ChanPos,float,ChanPos>  T0s;  // Channel number + pos ==>> T0
};

}

#endif //CsStrawTubesDetector_Calib_inc
