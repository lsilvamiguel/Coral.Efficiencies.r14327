#ifndef RecCall_h
#define RecCall_h


/*!
  \class RecCall
  \brief Recon package
  
  
*/


#include <list>
#include <CsSTD.h>



class CsTrack;
class CsCluster;

class RecCall  
{
  
 public:

  RecCall(){}; 
                // Constructor
 

  virtual ~RecCall(){};        //Destructor
  
   static RecCall&  ref();    //returns reference to this object 
 
  //methods
 
  void  ReconIni();       // get input for recon_ini.F  
 
  void GetClusters(std::list<CsCluster>& listClus); //get unused clusters  
  
    
  void  RecClus(std::list<CsTrack*>& lCsTrk); //export tracks to coral
  
 private:  
  
  
 
}; 


#endif //RecCall_h
  
  
  







