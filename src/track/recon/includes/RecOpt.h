#ifndef RecOpt_h
#define RecOpt_h

/*
  \class RecOpt
  \brief  Recon options

  Parameters, modes, cuts etc. for Recon, taken from CsOpt.
  
  Abstract class with static data members
  to provide direct, fast and global access to flags and cut
  
*/



#include "CsSTD.h"

class RecOpt {
 public:
  
  
  //  data
  
  static int     inPar [4];            //number of iterations ,min. number of hits used 
                                // for reconstruction and refit in spline parameter
  static int     McPar [3];          //MC switch
  static float   roadPar [17];       // roadwidth in y,z,theta (det.1-2)
  static float   activeRegion [7];    // active region in telescope 1 and 2
  static int   nbPlanes [16];      // min. numbers of addit. planes in tel.1 and tel.2 (th and all)
  static float ReProj  [16];        //tolerances in assoc (z and zp) for proj. 
  static float roadLengh [4];       //fraction of rad.lenght for each plane in micros,straws/SDC,sci-fi,gems
  static float chiPar  [17];       //chi parameters
  static int   deadPlanes1 [50];    // numbers of planes "killed" in reconstruction before SM1
  static int   deadPlanes2 [100];    // numbers of planes "killed" in reconstruction behind SM1
  static float ReTarget [4];      //target parameters
  static int Switch;              //switch Recon  on or off  
  static int blacklist [3];              //switch planes in reconstruction on or off
 
   // different methods  

   virtual void dum()=0;           //!< to make this class abstract
   static bool getRecOptions();       //!< get and store Recon options
   
 private:
   template <class Typ> static bool getOptArray(std::string, std::string, Typ*, int);
}; 




#ifdef INI_STATIC_HERE
int   RecOpt::inPar [4];
int   RecOpt::McPar [3] ;
float RecOpt::roadPar [17];
float RecOpt::activeRegion [7];
int   RecOpt::nbPlanes [16];
int   RecOpt::deadPlanes1 [50];
int   RecOpt::deadPlanes2 [100];
float RecOpt::ReProj  [16];
float RecOpt::roadLengh [4];
float RecOpt::chiPar  [17];
float RecOpt::ReTarget  [4];
int   RecOpt::Switch;
int   RecOpt::blacklist [3];

#undef INI_STATIC_HERE
#endif


#endif









