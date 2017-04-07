#ifndef Recon_h
#define Recon_h


/*!
  \class Recon
  \brief Recon package
  
  General interfaces class
  
  \warning Only one instance of this class is allowed
*/

#include <list>
#include <cstdlib>
#include <iostream>
#include "CsRegistry.h"
#include "CsEndOfJob.h"
#include <CsSTD.h>



class CsCluster;

class Recon : public CsEndOfJob, 
	      public CsEndOfEvent 
{
  
 public:

  Recon();                 // Constructor
  
  virtual ~Recon();        //Destructor
  
  static Recon&  ref();    //returns reference to this object

  //data members
  double magType[2];              // Magnes position
  bool   isScheme1;               // reconstruction schema flag
  bool   isMonteCarlo;            // Monte Carlo flag         
  std::list<CsCluster*> listUnusedClus;  
   
  //methods
 
  bool end();             // end of job method CsEndofEvent class

  bool eoe();             // end of event method CsEndofEvent class

  //! Increment timer by specified value (sec.) for using traffic watches
   void   countTime(float t) {time+=t;}
   void   countTime1(float t) {time1+=t;}
   void   countTime2(float t) {time2+=t;}
   void   countTime3(float t) {time3+=t;}
   void   countTime6(float t) {time6+=t;}
   void   countTime4(float t) {time4+=t;}
   void   countTime7(float t) {time7+=t;}

 private:

  static Recon* address;
  float time;
  float time1;
  float time2;
  float time3;
  float time4;
  float time6;
  float time7;
}; 


inline Recon& Recon::ref()
{
  if (address != 0) {
    return(*address);
  } else {
    std::cout<<"The object of Recon class is not yet created or already destructed"<<std::endl;
    exit(1);
  }
}

#endif //Recon_h
  
  
  









