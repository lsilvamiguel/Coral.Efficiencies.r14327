// $Id: CsInit.h 14070 2015-09-20 01:36:58Z lsilva $

/*!
   \file    CsInit.h
   \brief   Compass Initialization Class.
   \author  Benigno Gobbo 
   \version $Revision: 14070 $
   \date    $Date: 2015-09-20 03:36:58 +0200 (Sun, 20 Sep 2015) $
*/

#ifndef CsInit_h
#define CsInit_h

#include "coral_config.h"
#include "DaqDataDecoding/DaqEventsManager.h"

#include <string>
#include <list>

class CDB;


/*! \class CsInit 
    \brief Compass Initialization Class. 

    Performs all initialization stuff. It is a singleton class and is
    Instantiated first time by the Coral interface.  It is, first of all,
    responsible for the first CsOpt singleton instantiation. It keeps in
    memory parameters kept from options file. 
*/

class CsInit {

 public:

  /*! \fn static CsInit* Instance( int argc, char **argv );
    \brief First singleton instantiation.

    At the moment a single argument is accepted and mandatory: the options
    file name. Only the \c -h option is accepted, if specified Coral prints a 
    short help message and exits.
    \param argc parameter from main function
    \param argv parameter from main function
   */
  static CsInit* Instance( int argc, char **argv );

  /*! \fn static CsInit* Instance();
    \brief singleton instantiation (but first).
  */
  static CsInit* Instance();

  /*! \fn void dumpOptions( int dumplevel );
      \brief dumps the list of strings read from input file, and eventally
      the list of unused lines and the list of lines with errors inside.
      \param dumplevel if 0 only the list of read lines is printed;
      if 1 also the list of unused lines is printes; if greather that 1
      also the list of lines with errors is printed.
  */
  void dumpOptions( int dumplevel );

  /*! \fn inline string* getDetectorTable();
    \brief Returns the detector table file name.
  */
  std::string* getDetectorTable() { return &detTable_; }

  /*! \fn inline string* getPitchTable();
    \brief Returns the variable-sized pitch table file name.
  */
  std::string* getPitchTable() { return &pitchTable_; }

  /*! \fn inline bool IsVarPitch();
    \brief Returns \c true if a file for variable-sized pitches is given. 
    Returns \c false otherwise.    
  */
  bool IsVarPitch() const  { return varPitchTable_; }

  /*! \fn map<unsigned int,char>& getTCSMasks();
    \brief Returns the variable-sized pitch table file name.
  */
  const std::map<unsigned int,char>& getTCSMasks() { return TCSMasks_; } 

  // Next tree methods are for Monte Carlo processing only.

  /*! \fn inline bool IsAMonteCarloJob();
    \brief Returns \c true if this is a Monte Carlo processing job. 
    Returns \c false otherwise.    
  */
  bool IsAMonteCarloJob()const  { return(mcJob_);} 

  /*! \fn inline list<string*> getMCFilesList();
    \brief In Monte Carlo processing jobs return the list of selected
    input binary zebra files.
  */
  const std::list<std::string*>& getMCFilesList() const { return(mcFilesPtr_);} 

  /*! \fn inline list<string*> getMCNtFilesList();
    \brief In Monte Carlo processing jobs return the list of selected
    input ntuple zebra files.
  */
  const std::list<std::string*>& getMCNtFilesList() const { return(mcNtFilesPtr_);} 

  // Next file methods are for Monte Carlo processing only.

  /*! \fn inline bool   IsADataJob()
    \brief Returns \c true if this is a Real Raw Data processing job. 
    Returns \c false otherwise.    
  */
  bool IsADataJob() const  { return(dtJob_);} 

  /*! \fn inline bool   IsAHadronJob()
    \brief Returns \c true if this is (declared as) a Hadron Data(2008:9) processing job. 
    Returns \c false otherwise.    
  */
  bool IsAHadronJob() const  { return(hadronJob_);} 

  /*! \fn inline bool   IsABCSJob()
    \brief Returns \c true if this is (declared as) a BeamCharge&Spin Data(2008:9) processing job. 
    Returns \c false otherwise.    
  */
  bool IsABCSJob() const  { return(BCSJob_);} 

  /*! \fn inline int    getYear();
    \brief Returns the Year when selected data was taken.
  */
  int   getYear() const    { return(year_);}

  /*! \fn inline string getPeriod();
    \brief Returns the SPS period when selected data was taken
  */
  const std::string& getPeriod() const  { return(period_);}

  /*! \fn inline string getDataType();
    \brief Returns the type of real data (e.g. "raw", "dst", etc.)
  */
  const std::string& getDataType() const { return(dataType_);}

  /*! \fn inline list<int>    getRuns();
    \brief Return the list of runs to be read.
  */
  const std::list<int>&    getRuns() const    { return(runs_);}

  /*! \fn inline string getContainer();
    \brief Return the DB container to be read.
  */
   const std::string& getContainer()  const     { return(container_);}

  /*! \fn inline bool isFromDB()
    \brief In real raw data job Returns \c true if data is from Data Base. 
    Returns \c false otherwise.    
  */
  bool   isFromDB() const { return(fromDB_);} 

  /*! \fn inline list<string*> getDateFilesList();
    \brief In Data processing jobs return the list of selected
    input date binary files (not from DB).
  */
  const std::list<std::string*>& getDateFilesList() { return(dtFilesPtr_);} 

  /*! \fn unsigned int getSkipEvents();
    \brief Returns the number of events to skip from first read. 
  */
  inline unsigned int getSkipEvents()const { return( skipEvents_ ); } 

  /*! \fn unsigned int getMaxEvents();
    \brief Returns the maximun number of events to read. 
  */
  inline unsigned int getMaxEvents()const { return( maxEvents_ ); }

  /*! \fn unsigned int getMaxConsecutiveSkips();
    \brief Returns the max. # of consecutive events w/ decoding error that we tolerate (and hence skip) before giving up. 
  */
  inline unsigned int getMaxConsecutiveSkips()const { return( maxConsecutiveSkips_ ); }

  /*! \fn float getMinTimeInSpill();
    \brief Returns the cut on time elapsed since BoS. 
  */
  inline float getMinTimeInSpill()const { return( minTimeInSpill_ ); }

  //! Returns the pointer to the zebra buffer space
  inline int* getZebra() { return( zebra_ ); }

  //! Returns the CORAL major release number 
  int getCoralMajorRelease();

  //! Returns the CORAL minor release number
  int getCoralMinorRelease();

  //! Returns the CORAL build number 
  int getCoralBuild();

  //! Returns CDB 
  CDB *getDB(){return _cdb;}

  //! Returns \c true if Zebra was initialized 
  inline bool initZebraDone() const{ return(initzebradone_); }

/*! \fn void initZebra(const int nzebra );
    \author B.Gobbo
    \brief Initialize Zebra. The space is allocated dynamically
    \param nzebra the size of the Zebra store in words.
*/
  void initZebra();

  //! Returns \c true if parallel reconstruction is set
  bool parallelReconstruction() const{ return( parallel_ ); }

  //! Returns \c true if reconstructed events have to be stored on DBs
  bool saveRecoEvents() const{ return( saveRecoEvents_ ); }

  //! Returns \c true if Mickey Mouse hits simulation is set
  bool mickeyHits() const{ return( _mickeyhits ); }

  //! Returns \c true if complete Mickey Mouse simulation is set
  bool mickeyAll() const{ return( _mickeyall ); }

  // undocumented
  bool mickeyOnReco() const{ return( _mickeyonreco ); }

  //! Number of events to be generated by Mickey Mouse MC
  int mickeyNumberOfEvents() const{ return( _mickeynevents ); }

  //! Returns \c true if resetting random seed @ each event is requested
  bool resetRandomSeed() const { return _resetRandomSeed; } 
  
  bool useCalibration(void) const {return use_calibration;}
  
  //! \returns \c true if any of CDB's is used, whether ConditionsDB or FileDB
  bool useCDB()          const { return cdbSwitch_!=USE_NoCDB; }
  //! \returns \c true if ConditionsDB CDB is used
  bool useConditionsDB() const { return cdbSwitch_==USE_ConditionsDB; }
  //! \returns \c true if FileDB CDB is used
  bool useFileDB()       const { return cdbSwitch_==USE_FileDB; }
  //! \returns \c true if MySQLDB CDB is used
  bool useMySQLDB()      const { return cdbSwitch_==USE_MySQLDB; }
  //! \returns the run time to use with CDB
  int  getCDBUseTime()      const { return cdbUseTime; }

  const
  CS::DaqEventsManager &getDaqEventsManager(void) const {return manager;}
  CS::DaqEventsManager &getDaqEventsManager(void)       {return manager;}

  const CS::Chip::Maps& getDaqMaps              (void) const {return getDaqEventsManager().GetMaps();}
        CS::Chip::Maps& getDaqMaps              (void)       {return getDaqEventsManager().GetMaps();}
  const CS::DaqOption& getDaqOptions(void) const {return getDaqEventsManager().GetDaqOptions();}
        CS::DaqOption& getDaqOptions(void)       {return getDaqEventsManager().GetDaqOptions();}

  int getStartOfRun(){return startOfRun_;}
  int getEndOfRun(){  return endOfRun_;}

  //! \Returns \c true if the updating the solenoid field is enabled in options.
  bool updateSolenoidField() const { return updateSolenoidField_; }

 protected:

  /*! \fn CsInit( int argc, char **argv );
    \brief The Constructor.
    At the moment a single argument is accepted and mandatory: the options
    file name. Only the -h option is accepted, if specified Coral prints a 
    short help message and exits.
    \param argc parameter from main function
    \param argv parameter from main function
  */ 
  CsInit( int argc, char **argv );

  /*! \fn void GetMagnetInfo(int Run, unsigned int ignore, int &solSign, int &dipole, int &smSign);
     \brief Retrieve magnet info for argument run#, from MySQL or File depending upon option.
     \param ignore:  0x1: ignore Target, 0x2: ignore SM1/2
     \param solSign: Sign of solenoid field (+/-1 or 0 or -2: not available).
     \param dipole:  Dipole is on/off (1 or 0 or -2: not available).
     \param solSign: Sign of SM1/2 fields (+/-1 or -2: not available).
   */
  void GetMagnetInfo(int Run, unsigned int ignore,
		     int &solSign, int &dipole, int &smSign);

  /*! \fn string GetDetectorsDat(string path, int run, int jobType);
     \brief Derive detectors.dat from argument path. If path points to a directory, the file therein most appropriate for argument run is returned, accounting for relevant magnet info if requested. Magnet info discriminates among "plus/minus/transv" extensions, based on arget magnet setting, and among "mu+/mu-" based on SM1/2 polarity, if no target field.
     \param jobType: 0: spin asymmetry, 1: hadron, 2: charge&spin asymmetry.
   */
  std::string GetDetectorsDat(std::string path, int run, int jobType = 0);

 private:

  static CsInit* instance_;    //!< The singleton pointer

  unsigned int skipEvents_;    //!< Number of events to skip
  unsigned int maxEvents_;     //!< Max number of events to read
  unsigned int maxConsecutiveSkips_;  //!< Max. # of consecutive events w/ decoding error that we tolerate (and hence skip) before giving up.
  float minTimeInSpill_;       //!< Cut on time elapsed since BoS

  CDB* _cdb;                   //!< Pointer to the CDB
  int startOfRun_;
  int endOfRun_;

  // Data Run. Needs at least a data file name. At present uses the same
  // detector table than MC.
  bool      dtJob_;                        //!< \c true if data job
  int       year_;                         //!< Year of data taking
  std::string    period_;                  //!< SPS period of data taking
  std::string    dataType_;                //!< Real Data Type (raw, dst, etc.)
  bool      fromDB_;                       //!< \c true if data from Data Base
  std::list<int> runs_;                    //!< List of runs to be processed
  std::string    container_;               //!< The DB container, if specified
  bool      saveRecoEvents_;               //!< \c true if Reco Events to be stored in DB 
  bool      parallel_;                     //!< \c true if parallel reconstruction is set
  std::list<std::string>  dtFiles_;        //!< Binary data files (DATE format)
  std::list<std::string*> dtFilesPtr_;     //!< The corresponding list of pointers

  bool      mcJob_;                        //!< \c true if Monte Carlo job
  std::list<std::string>  mcFiles_;        //!< Zebra binary Monte Carlo files
  std::list<std::string*> mcFilesPtr_;     //!< The correspondig list of pointers
  std::list<std::string>  mcNtFiles_;      //!< Zebra ntuple Monte Carlo files
  std::list<std::string*> mcNtFilesPtr_;   //!< The correspondig list of pointers

  const static int nZebra_ = 500000;       //!< Size of to the zebra buffer space
  static int zebra_[nZebra_];              //!< Pointer to the zebra buffer space

  bool      _mickeyhits;         //!< true if Mickey Mouse extrap. 
  bool      _mickeyall;          //!< true if Mickey Mouse all
  bool      _mickeyonreco;       // undocumented
  int       _mickeynevents;      //!< Number of events to be generated 
  bool      _resetRandomSeed;    //!< \c true if resetting random seed @ each event is requested

  bool      hadronJob_;      //!< \c true if (declared as) a Hadron(2008:9) job
  bool      BCSJob_;         //!< \c true if (declared as) a BCS asymmetry job

  std::string    detTable_;       //!< Detector data file name 
  bool           varPitchTable_;  //!< \c true if variable-sized pitch table found
  std::string    pitchTable_;     //!< Data file name for variable-sized pitch table

  std::map<unsigned int,char> TCSMasks_; //!< Bit mask <-> trigger type, designated by a one character type code. 

  bool      use_calibration; //!< Use calibration data

  //! \enum CDBSwitch
  //! \brief Defines the type of CDB used
  enum CDBSwitch {
    USE_NoCDB        = 0, //!< Default = No CDB
    USE_ConditionsDB = 1, //!< ConditionsDB CDB
    USE_FileDB       = 2, //!< FileDB CDB      
    USE_MySQLDB      = 3  //!< MySQLDB CDB      
  };
  CDBSwitch cdbSwitch_;      //!< Switch for CDB use
  time_t    cdbUseTime;      //!< run time to use with CDB instead of real start-of-run time

  bool initzebradone_;         //!< \c true if zebra was initialized
  
  CS::DaqEventsManager manager;

  bool updateSolenoidField_;   //!< \c true if the updating the solenoid field is enabled in options.
};

#endif // CsInit_h
