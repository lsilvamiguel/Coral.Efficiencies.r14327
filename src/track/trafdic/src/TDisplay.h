// $Id: TDisplay.h 14069 2015-09-17 20:44:46Z lsilva $

#ifndef TDisplay_h
#define TDisplay_h

/*!
  \class TDisplay
  \brief Event Display
  
  COMPASS spectrometer Event Display
  for tracking algorithms development and debugging
  
  \warning Only one instance of this class is allowed
  \author Sergei.Gerassimov@cern.ch
*/

#include <cassert>
#include "CsSTD.h"
#include "TSetup.h"
#include "TDetect.h"
#include "TPlane.h"
#include "TTrack.h"
#ifdef COMPASS_USE_OSPACE_STD
# include <ospace/std/stack>
#else
# include <stack>
#endif

class TDisplay {

public:

  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  TDisplay();   //!< Constructor
  ~TDisplay();  //!< Destructor

  //! Returns pointer to this object
  static TDisplay* Ptr();

  //! Returns reference to this object
  static TDisplay& Ref();

  //! Main drawing function
  int Draw(int key);

  //! Method for misceleneouse debugging drawings
  void UserDraw();
  
  //! get flag to draw Runge-Kutta extrapolation trajectories
  bool Rkutraj() { return(rkutraj); }

  //! get flag to draw extrapolation trajectorie by material map call coordinates
  bool Mmaptraj() { return(mmaptraj); }

  //! set flag to draw Runge-Kutta extrapolation trajectories
  void SetRkutraj(bool flg) {rkutraj = flg;}

  //! set flag to draw extrapolation trajectorie by material map call coordinates
  void SetMmaptraj(bool flg) {mmaptraj = flg;}

  //! draw point X[3] with specified color index
  void point(double x[], int color);

private:

  static TDisplay* address;

  bool Init();               //!< Display initialization. (Called from the constructor)

  // data members
  int   Proj;     // current projection number
  float AngProj;  // current projection angle (Warning: rounded value)
  int   DrOpt[20];// Misc. drawing options
  float Xvf[2];   // current X view field
  float Yvf[2];   // current Y view field
  float XvfN[2];  // current normalized X view field
  float YvfN[2];  // current normalized Y view field
  float Xmenu[2]; // Menu left upper corner
  float Wmenu[2]; // Menu size (in X - per simbol, in Y - per item)
  bool  rkutraj;  // flag to draw Runge-Kutta extrapolation trajectories (rkutcoord.cc)
  bool  mmaptraj; // flag to draw extrapolation trajectorie by material map call coordinates
  int   ev_count; // counter of drawn events
  std::stack <float> stXvf0, stXvf1, stYvf0, stYvf1; // stackes for zoomed veiw fields

  //         ********** RE-ORDER RECONSTRUCTION ZONES **********
  int *iZone2Zone, *iplLasts;

  // Drell-Yan absorber
  int hasAbsorber;               // 0: no absorber, >0: else
  float xA[13], yA[13], zA[13];	 // Starting Jura upstream, clockwise

  // private methods (called by "Draw" function)
 
  void Clear();
  void DrawFrame();
  int  DrawMenu();
  void DrawModeMenu();
  void DrawMag();
  void DrawMaterialMap();
  void DrawMisc();
  void DrawText();
  void DrawString(std::string s, float x, float y, int color, float slope);
  int  DrawHits  (int flg = 0);
  int  DrawHitsMC(int flg = 0);
  int  DrawDet   (int flg = 0);
  int  DrawTracks  (int flg = 0);
  bool DrawPolyLine(bool smooth, int color,
		    int ifl, float px, float py, float &dist,
		    const TTrack &t);
  bool DrawExtrapolate(CsVertex *bestVertex, int color,
		       int ifl, float px, float py, float &dist,
		       const TTrack &t);
  int  DrawMCTracks(int flg = 0);
  void TrackInfo   (int id);
  void MCTrackInfo (int id);
  void DetInfo     (int ip);
  void HitInfo     (int ih);



  void IPLm(int n, float* xx, float* yy); //!< modified version of HIGZ IPL function 
  void IPMm(float x, float y);            //!< modified version of HIGZ IPM function
  void Distort(float& x, float& y);       //!< coodinate transformations

};

// Inline functions

inline TDisplay* TDisplay::Ptr()
{ 
  // No check on existence of the object. Normaly, TDisplay::Ref() has to be used.
  return(address);
}

inline TDisplay& TDisplay::Ref()
{
  if (address != 0) {
    return(*address);
  }
  std::cout<<"TDisplay::Ref() ==> The object of TDisplay class is not yet created or already destructed"<<std::endl;
  assert(false);
  return(*address);  // just to get rid of compiler warnings
}

#endif










