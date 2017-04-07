// $Id: Coral.h,v 1.13 2010/02/08 18:52:04 tnagel Exp $

/*!
   \file    Coral.h
   \brief   Compass Reconstruction and AnaLysis Package.
   \version $Revision: 1.13 $
   \author  Benigno Gobbo
   \date    $Date: 2010/02/08 18:52:04 $
*/

#ifndef Coral_h
#define Coral_h

#include "CsSTD.h"
#include "CsTypes.h"

#include "CsInit.h"
#include "CsRegistrySing.h"
#include "CsEvent.h"
#include "CsGeant3.h"
#include "CsGeom.h"

/*! \class Coral
    \brief   Compass Reconstruction and AnaLysis Package.

    The Coral User interface. This singleton must be istantiated from
    the user main program with the same argc and argv of the main 
    function. The class is mainly intended to allow the user to access
    the coral package objects and methods is a easier way.
*/

class Coral {

 public:

  /*! \fn static Coral* init( int argc, char** argv );
    \brief First singleton instantiation.

    At the moment a single argument is accepted and mandatory: the options
    file name. Only the -h option is accepted, if specified Coral prints a 
    short help message and exits.
    \param argc The same argument of main function  
    \param argv The same argument of main function
  */
  static Coral* init( int argc, char **argv );

  /*! \fn static Coral* Instance();
    \brief singleton instantiation (but first).
  */
  static Coral* Instance();
  
  /*! \fn bool getNextEvent();
    \brief Read next event and fill the related objects. Returns \c true
    of the event is read, \c false otherwise.
  */
  inline bool getNextEvent() { return( CsEvent::Instance()->getNextEvent() ); }

  /*! \fn inline void end();
    \brief This method should be called at the end of process. It 
    automatically sends a call to all registered packages.  
  */
  inline void end( ) { CsRegistrySing::Instance()->callEndMethods(); }

  /*! \fn inline CsInit* getCsInit();
    \brief Returns the pointer to the CsInit singleton object
  */
  inline CsInit* getCsInit() { return( CsInit::Instance() ); }

  /*! \fn inline CsEvent* getEvent();
    \brief Returns the pointer to the CsEvent singleton object
  */
  inline CsEvent* getEvent() { return( CsEvent::Instance() ); }

  /*! \fn inline bool isAMonteCarloEvent() const;
    \brief Returns \c true if the current event is a Monte Carlo one
  */
  inline bool isAMonteCarloEvent() const { 
    return(CsEvent::Instance()->isAMonteCarloEvent()); 
  }

  /*! \fn inline std::list<CsMCTrack*> getMCTracks();
    \brief Returns the list of pointers to the current event Monte Carlo true
    tracks.
  */
  inline std::list<CsMCTrack*> getMCTracks() { 
    return(CsEvent::Instance()->getMCTracks()); 
  }

  /*! \fn inline std::list<CsMCVertex*> getMCVertices();
    \brief Returns the list of pointers to the current event Monte Carlo true
    vertices.
  */
  inline std::list<CsMCVertex*> getMCVertices() { 
    return(CsEvent::Instance()->getMCVertices()); 
  }

  /*! \fn inline std::list<CsMCHit*> getMCHits();
    \brief Returns the list of pointers to the current event Monte Carlo true
    hits.
  */
  inline std::list<CsMCHit*> getMCHits() { 
    return(CsEvent::Instance()->getMCHits()); 
  }

  /*! \fn std::list <CsDetector*> getDetectors();
    \brief Returns the list of pointers to the CsDetector objects.
  */
  std::list <CsDetector*> getDetectors() { 
    return(CsGeom::Instance()->getDetectors()); 
  }

  /*! \fn inline std::string getMCFileName();
    \brief Returns the current Monte Carlo file name
  */
  inline std::string getMCFileName() { 
    return( CsGeant3::Instance()->getMCFileName()); 
  }

  /*! \fn inline bool isADataEvent() const;
    \brief Returns \c true if the current event is a Read Data one
  */
  inline bool isADataEvent() const { 
    return(CsEvent::Instance()->isADataEvent()); 
  }

  /*! \fn inline uint32 getRunNumber() const;
    \brief In Data Events returns the Run Number. It is 0 if MC events.
  */
  inline uint32  getRunNumber() const { 
    return(CsEvent::Instance()->getRunNumber()); 
  }

  /*! \fn inline uint32  getEventNumberInRun() const;
    \brief In Data Events returns the Event Number in Run. It is 0 if MC 
    events.
  */
  inline uint32  getEventNumberInRun() const { 
    return(CsEvent::Instance()->getEventNumberInRun()); 
  }

  /*! \fn inline uint32  getBurstNumber();
    \brief In Data Events returns the Burst Number. It is 0 if MC events.
  */
  inline uint32  getBurstNumber() { 
    return(CsEvent::Instance()->getBurstNumber());
  }

  /*! \fn inline uint32  getEventNumberInBurst();
    \brief In Data Events returns the Event Number in Burst. It is 0 if MC 
    events.
  */
  inline uint32  getEventNumberInBurst() { 
    return(CsEvent::Instance()->getEventNumberInBurst()); 
  }

  /*! \fn inline uint32  getTriggerMask();
    \brief In Data Events returns the Event Trigger Mask. It is 0 if MC 
    events.
  */
  inline uint32  getTriggerMask() { 
    return(CsEvent::Instance()->getTriggerMask()); 
  }

  /*! \fn inline uint32  getErrorCode();
    \brief In Data Events returns the Error Code. It is 0 if MC 
    events.
  */
  inline uint32  getErrorCode() { 
    return(CsEvent::Instance()->getErrorCode()); 
  }

//   /*! \fn inline uint32  getRawBuffer();
//     \brief In Data Events returns pointer to the Raw Data Buffer. It is 0 if 
//     MC events.
//   */
//   inline uint8*  getRawBuffer() { 
//     return(CsEvent::Instance()->getRawBuffer()); 
//   }

  /*! \fn inline CsTime  getEventTime();
    \brief In Data Events returns the Event Time in CsTime format (see this
    class description for details). In MC events return an empty object.
  */
  inline CsTime  getEventTime() { 
    return(CsEvent::Instance()->getEventTime()); 
  }

//   /*! \fn inline void    getRawEvent();
//     \brief Gets this Raw Data event. 
//   */
//   inline void getRawEvent() { 
//     CsEvent::Instance()->getRawEvent(); 
//   }
// 
//   /*! \fn inline std::list<CsDate> getSubEvents()
//     \brief In Data Events returns the list of Date Sub Events (see the Date
//     package or the CsDate description for details). In MC events returns 
//     an empty list.
//   */
//   inline std::list<CsDate> getSubEvents() { 
//     return( CsEvent::Instance()->getSubEvents()); 
//   }
// 
//   /*! \fn inline std::list<CsDateEquipment> getEquipments();
//     \brief In Data Events returns the list of Date Equipments (see the Date
//     package or the CsDateEquipment description for details). In MC events 
//     returns an empty list.
//   */
//   inline std::list<CsDateEquipment> getEquipments() { 
//     return( CsEvent::Instance()->getEquipments()); 
//   }
// 
//   /*! \fn inline bool sameEndianess();
//     \brief Keeps track of the relative endianess of the machine where
//     the job is running with respect to the machine where the Date file
//     where produced. Returns \c true if the endianess is the same, \c false
//     otherwise.
//   */
//   inline bool sameEndianess() { 
//     return( CsEvent::Instance()->sameEndianess()); 
//   }

  /*! \fn inline std::list<CsDigit*> getDigits();
    \brief Performs the Digit simulation from the Monte Carlo true Hits.
    Returns the list of pointers to the generated CsDigit objects.
  */
  inline std::list<CsDigit*> getDigits() { 
    return(CsEvent::Instance()->getDigits()); 
  }

  /*! \fn inline std::list<CsCluster*> getClusters();
    \brief Performs the Cluster simulation from Digits.
    Returns the list of pointers to the generated CsCluster objects.
  */
  inline std::list<CsCluster*> getClusters() { 
    return(CsEvent::Instance()->getClusters()); 
  }

 protected:

  Coral( int argc, char **argv ); //!< The Protected Singleton Constructor

 private:
  static Coral* instance_; //!< The singleton static attribute

};

#endif // Coral_h
