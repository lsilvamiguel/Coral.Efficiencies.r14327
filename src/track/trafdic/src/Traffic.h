#ifndef Traffic_h
#define Traffic_h

/*!
  \class Traffic
  \brief Traffic package
  
  General interfaces class
  
  \warning Only one instance of this class is allowed
  \author Sergei.Gerassimov@cern.ch

*/

#include <iostream>
#include <cassert>
#include "CsSTD.h"
#include "CsRegistry.h"
#include "TWatches.h"

class Traffic:  public CsEndOfJob, 
		public CsEndOfEvent {
public:

  static unsigned int NobjCreated;
  static unsigned int NobjDestructed;

  Traffic();   //!< Constructor
  virtual ~Traffic();  //!< Destructor

  //! Returns pointer to this object

  static Traffic* Ptr();

  //! Returns reference to this object

  static Traffic& Ref();

  //! "End of job" method
  bool end();

  //! "End of event" method
  bool eoe();

  // misc. functions for (general "statistics" only)
  
  //! Timer 
  TWatches Stopwatch;


private:

  static Traffic* address;
  static int      Nevt;

  // counters (just for final statistics printout in "end()" method)
  long unsigned int nevs;
  long unsigned int ntracks;
  long unsigned int ntracks_with_P;
  long unsigned int ntracks_beam;
  long unsigned int ntracks_beam_with_P;
  long unsigned int n_selected_MC;
  long unsigned int n_reconstructed_MC;
  
  friend class TEv; 

};

// Inline functions

inline Traffic* Traffic::Ptr()
{
  // No check for existence of the object. Normaly, TEv::Ref() has to be used.
    return(address);
}

inline Traffic& Traffic::Ref()
{
  if (address != 0) {
    return(*address);
  }
  std::cout<<"The object of TRaffic class is not yet created or already destructed"<<std::endl;
  assert(false);
  return(*address); // just to get rid of compiler warnings
}
#endif

















