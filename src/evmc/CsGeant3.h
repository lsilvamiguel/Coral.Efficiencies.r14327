// $Id: CsGeant3.h,v 1.20 2003/05/23 06:50:27 benigno Exp $

/*!
   \file    CsGeant3.h
   \brief   Geant Interface Class.
   \author  Benigno Gobbo
   \version $Revision: 1.20 $
   \date    $Date: 2003/05/23 06:50:27 $

   \par History:
   19990420  Start of history.<br>
   19990422  First running stuff.<br>
   19990426  Class is a singleton now.<br>
   19990504  Modified due to CsMCHit changes.<br>
   19990522  Add Rotation (Geant ref)->(Compass ref).<br>
   19990720  Set up Wires Reference System and Compass Units.<br>
   19990726  Added fill of CsFiels class.<br>
   19990803  Added Comgeant n-tuple outputs.<br>
   19990906  Removed the geometry part.<br>
   20000512  Moved to COMGEANT 5. Also backward compatible.<br>
   20000622  Some stuffs moved from CsEvent to here.<br>
   20000824  Mickey Mouse Monte Carlo.<br>
*/

#ifndef CsGeant3_h
#define CsGeant3_h

#include "coral_config.h"

#include "CsSTD.h"
#include "CsMCHit.h"
#include "CsMCTrack.h"
#include "CsMCVertex.h"
#include "CsTime.h"
#include "CsTmpTrigger.h"
#include <CLHEP/Matrix/Matrix.h>
#include <string.h>


#if USE_TGEANT
  #include "includeTGEANT/T4Event.hh"
  #include "includeTGEANT/T4OutputBackEnd.hh"
  #include "includeTGEANT/T4OutputROOT.hh"
  #include "includeTGEANT/T4OutputASCII.hh"
#endif

#define Lund_Structures_Version (1)

typedef struct 
{ float x, y, w2, q2, u; int lst[35]; 
  float parl[30]; int genType; float cut[14]; 
  float parhfl[6], cuthfl[8]; int lsthfl[10];
} ludata;

typedef struct 
{ float x, y, w2, q2, u; float uservar[20]; int lst[40]; 
  float parl[30]; int genType; float cut[14]; 
  float parhfl[10], cuthfl[8]; int lsthfl[10];
} ludatanew;

typedef struct 
{ int k[5]; float p[5], v[5]; int lu2kine; 
} lujet;

typedef struct
{ int msel; int mselpd; int msub[500]; int kfin[2][80]; float ckin[200];
} pysubs;

typedef struct
{ int mstp[200]; float parp[200]; int msti[200]; float pari[200];
} pypars;

/*! \class CsGeant3 
  \brief Geant Interface Class.

  This class reads Monte Carlo files generated with Comgeant. CsMCHit,
  CsMCTrack and CsMCVertex objects are created from Comgeant output.
  At present only zebra output files can be read.
*/

class CsGeant3 {

public:

  /*! \fn static CsGeant3* Instance();
    \brief singleton instantiation.
  */
  static CsGeant3* Instance();

  /*! \fn bool getNextEvent( unsigned long selTrigMask = 0xffffffff );
    \brief look for next event, read it and fill objects lists
  */
  bool getNextEvent( unsigned long selTrigMask = 0xffffffff );

#if USE_TGEANT   
    /*! \fn inline bool isATGeantFile();
      \brief Returns \c true if the current event is from a TGEANT file
    */
    inline bool isATgeantFile() {
        return( _TGeantInUse );
    }
#endif
  
  /*! \fn inline bool isAnNtFile();
    \brief Returns \c true if the current event is from an NTuple file 
  */
  inline bool isAnNtFile() { return( _NTFile ); }

  /*! \fn std::string getMCFileName();
    \brief Returns the current Monte Carlo file name
  */
  std::string getMCFileName(); 

  /*! \fn CsTime getEventTime();
    \brief Not implemented at present.
  */
  inline CsTime getEventTime() { return( _time ); }

  /*! \fn std::list <CsMCTrack*>  &getMCTracks(void);
    \brief Returns the list of pointers to CsMCTrack objects
  */
  const std::list<CsMCTrack*>  &getMCTracks(void) const { return _tracks; }

  /*! \fn std::list <CsMCVertex*> &getMCVertices();
    \brief Returns the list of pointers to CsMCVertex objects
  */
  const std::list<CsMCVertex*> &getMCVertices(void) const { return _vertices; }

  /*! \fn std::list <CsMCHit*> &getMCHits(void);
    \brief Returns the list of pointers to CsMCHit objects
  */
  const std::list<CsMCHit*>    &getMCHits(void) const { return( _hits ); }

  /*! \fn int    getRunNb(void);
    \brief Returns the Run Number
  */
  int                   getRunNb(void) const { return _run; }

  /*! \fn int    getEventNb(void);
    \brief Returns the Event Number
  */
  int                   getEventNb(void) const { return _event; }

  /*! \fn int openGeantFile(const std::string* fname, int lrecl, const char *option)
    \brief This function assumes a Geant file in Zebra/FZ format
    created with the C I/O option "L". Returns 0 if OK 
    \param name The fine name
    \param lrecl The record lenght
    \param options CFOPEN function options. Defaults are: <br> 
    "I" = input mode, <br> 
    "X" = exchange mode, <br> 
    "L" = C library mode <br>
    (see Zebra doc. for details)
  */
  int  openGeantFile( const std::string* name, int lrecl=900, 
		       const char* option = "XIL");

  /*! \fn int openGeantNTFile(const std::string* fname, const int lrecl )
    \brief Opens an NT file. Returns 0 if OK.
    \param fname The file name
    \param lrecl The record lenght
  */
  int openGeantNTFile(const std::string* fname, const int lrecl=1024 );

  /*! \fn int readGeantEvent(const char *option)
    \brief Read the next event from the Geant FZ file.
    Search for the next HEAD record. If found, read VERT,KINE,HITS,OHIT and 
    DIGI records. Returns 0 if nothing read or EOF reached, 1 otherwise.
    \param option Options parameter. Possible values are:<br>
    ""  = only HEAD and VERT banks are read<br>
    "K" = read also KINE banks<br>
    "H" = read also HITS banks<br>
    "D" = read also DIGI banks<br>
    "O" = read also special Comgeant OHIT banks.
  */
  bool  readGeantEvent( const char* option = "KHDO");

  /*! \fn bool readGeantNTEvent()
    \brief Read the next event from the Geant Ntuple file.
    returns \c true if all OK
  */
  bool readGeantNTEvent();

  /*! \fn int readGeantHead()
    \brief Fill the informations from Geant3 HEAD bank.
    Return 0 if nothing read, 1 othervise.
  */
  int   readGeantHead();

  bool  readGeantLund();

  /*! \fn int readGeantNTHead()
    \brief Fill the informations from Header block of Comgeant Ntuples.
    Return 0 if nothing read, 1 othervise.
  */
  int   readGeantNTHead();

  /*! \fn int readMickeyMouseHead()
    \brief Fill the "head" informations from Mickey Mouse MC.
    Returns always one.
  */
  int   readMickeyMouseHead();

  /*! \fn int readGeantKine()
    \brief Fill MC Track and Vertex classes from Geant3 KINE and VERT banks.
    Returns 0 if nothing read, the Number of Tracks otherwise.
  */
  int   readGeantKine();

  /*! \fn int readGeantNTKine()
    \brief Fill MC Track and Vertex classes from Comgeant Ntuples.
    Returns 0 if nothing read, the Number of Tracks otherwise.
  */
  int   readGeantNTKine();

  /*! \fn int readMickeyMouseKine()
    \brief Fill MC Track and Vertex classes from Mickey Mouse internal
    Monte Carlo.
  */
  int   readMickeyMouseKine();

  /*! \fn int readGeantHits()
    \brief Fill MC Hit classe from Geant3 OHIT/HITS bank.
    At present (19990422) uses OHIT bank only.
    Implement JHIT in future? Will see...
  */
  int   readGeantHits();

  /*! \fn int readGeantNTHits()
    \brief Fill MC Hit classe from Geant3 Ntuples.
  */
  int   readGeantNTHits();

  /*! \fn int readMickeyMouseHits();
    \brief Fill MC Hit classe from internal Mickey Mouse Monte Carlo.
  */
  int   readMickeyMouseHits();

  /*! \fn void setNtuple();
    \brief Setup CWR ntuple parameters.
  */
  void setNtuple();

  //! clear the private lists
  void clear();

  /*! \fn unsigned long  getTriggerMask(void) const; 
    \brief Returns the Event Trigger Mask.  
    events. 
  */
  inline unsigned long getTriggerMask(void) const { return (_TrigMask); };
  
  //! Does loop over MC hits of mu' and checks which hodoscope is fired
  void setTriggerMask(void);
  
  // --- Lund stuff
  //! returns ludata structure (old ludata version, CG<7)
  ludata             getLudata(void) const { return( _ludata ); }
  const ludata&      getLudataRef(void) const { return( _ludata ); }
  //! returns ludatanew structure (new ludata version, GC>=7)
  ludatanew          getLudatanew(void) const { return( _ludatanew ); }
  const ludatanew&   getLudatanewRef(void) const { return( _ludatanew ); }
  //! returns lujet structures vector
  std::vector<lujet>        getLujets(void) const { return( _lujets ); }
  const std::vector<lujet>& getLujetsRef(void) const { return( _lujets ); }
  //! returns pysubs structure
  pysubs             getPysubs(void) const { return( _pysubs ); }
  const pysubs&      getPysubsRef(void) const { return( _pysubs ); }
  //! returns pypars structure
  pypars             getPypars(void) const { return( _pypars ); }
  const pypars&      getPyparsRef(void) const { return( _pypars ); }
  //! generator type: 2,5 Lepto format, 3 Aroma format, 6 new Pythia format
  int                getGenType(void) { 
    if( _useludatanew ) return( _ludatanew.genType );
    else                return( _ludata.genType );
  }
  //! true if new ludata format (ludatanew)
  bool               isLudataNew(void) const { return( _useludatanew ); }
  //! true if TLND bank available
  bool               isTlndOk(void) const { return( _tlndok ); }
  //! true if RLND bank available
  bool               isRlndOk(void) const { return( _rlndok ); }

  //! returns Comgeant version,release: e.g. 6.9, 7.0
  float              getCGVersion(void) const { return( _CGVersion ); }

  //! returns Geometry year,version: e.g. 2002.03
  float              getGeoVersion(void) const { return( _GeoVersion ); }

 protected:

  CsGeant3();                               //!< The Constructor

 private:

  static CsGeant3*    _instance;            //!< The Singleton Static pointer
  CsTime              _time;                //!< event time. Not yet used 
  std::string         _Gname;               //!< current file name
  int                 _run;                 //!< run number
  int                 _event;               //!< event number
  unsigned long       _TrigMask;            //!< trigger mask
  float               _CGVersion;           //!< comgeant version
  float               _GeoVersion;          //!< geometry version
  
  CsTmpTrigger _trig;                       //!< trigger object- added am

  std::list <CsMCTrack*>   _tracks;   //!< List of pointers to the MC True Tracks
  std::list <CsMCVertex*>  _vertices; //!< List of pointers to the MC True Vertices
  std::list <CsMCHit*>     _hits;     //!< List of pointers to the MC True Hits

  bool _NTFile;             //!< \c true if ntuple input file
  std::list<std::string*> _MCFiles;   //!< list of ptrs to the MC bin. zebra input files
  std::list<std::string*> _MCNtFiles; //!< list of ptrs to the MC ntp. zebra input files
  std::list<std::string*>::iterator _currentMCFile;   //<! current MC bin. zebra file
  std::list<std::string*>::iterator _currentMCNtFile; //<! current MC ntp. zebra file
  bool                    _mickeyfirst;
  std::vector< std::vector<float> > _mickeytracks;
  std::vector< std::vector<float> > _mickeyvertices; 

  ludata             _ludata;
  ludatanew          _ludatanew;
  std::vector<lujet> _lujets;
  pysubs             _pysubs;
  pypars             _pypars;
  bool               _useludatanew;
  bool               _tlndok, _rlndok;
  void               _clearMCstructs( void );
  void               _halfclearMCstructs( void );

  /*! \fn void  _geaRef2CsRefVec( double& x, double& y, double& z );
     \brief Rotates the axes from the Geant Reference System ( X along
     beam, Z vertical) to the Compass Reference System ( Z along beam,
     Y vertical).
     \param x x coordinate. WARNING: its value will be modified!
     \param y y coordinate. WARNING: its value will be modified!
     \param z z coordinate. WARNING: its value will be modified!
  */
  void  _geaRef2CsRefVec( double& x, double& y, double& z );

  /*! \fn void  void  GeaRef2CsRefMat( HepMatrix& a );
     \brief Rotates matrices from the Geant Reference System ( X along
     beam, Z vertical) to the Compass Reference System ( Z along beam,
     Y vertical)
     \param a the matrice to be rotated. WARNING: it will be modified!
  */
  void  _geaRef2CsRefMat( CLHEP::HepMatrix& a );

  struct _sortMCHits;   //!< To sort the list of pointers


private:
#if USE_TGEANT
// ###### TGEANT CODE START ######
    void openTgeantFile(std::string _fileName);
    void readTgeantEvent(void);
    void readTgeantLund(void);
    void readTgeantKine(void);
    void readTgeantHits(void);
    
    int  getIndexOfVertex(double t,double x,double y,double z);
    double zBeamStart;

    void createTgeantTrackingHit(T4HitData&, bool isTrigger=false);
    void createTgeantCaloHit(T4HitData&);
    void createTgeantRICHHit(T4RichData&);
    double calcTgeantTime(double _time, double* _hitPosition);

    T4Event* currentEvent;

    std::vector<CsMCVertex*> _tVertices;
    std::vector<CsMCTrack*> _tTracks;
    std::map<int,int> _tTrackId;

    T4OutputBackEnd* _outputBackEnd;
    std::vector<T4Event>* _t4Run;

// ###### TGEANT CODE END ######
#endif
    bool _TGeantInUse;
};
 
#endif // CsGeant3_h
