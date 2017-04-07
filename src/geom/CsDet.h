/*!
   \file    CsDet.h
   \brief   General properties of any COMPASS detector.
   \version $Revision: 1.28 $
   \author  COMPASS Software Group
   \date    $Date: 2010/01/24 16:10:40 $
*/

#ifndef CsDet_h
#define CsDet_h

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <ctime>

#include "DaqDataDecoding/DetID.h"
#include "DaqDataDecoding/Chip.h"
#include "DaqDataDecoding/Exception.h"
class CDB;

namespace CS
{
  class DaqEvent;
  class Detector;
  class Detectors;
}

class CsDigit;
class CsCluster;
class CsTrack;

/*! \brief Base abstract class of \c any COMPASS detector
*/
class CsDet
{
  // ---------------------------------------------------------------------------
  // Constructors, destructor
  // ---------------------------------------------------------------------------

  public:

    /// Destructor
    virtual            ~CsDet                   (void);

    /// Base constructor
                        CsDet                   (const CS::DetID &id, const std::string& _TBname);

  private:
    /// Copy constructor is not allowed
                        CsDet                   (const CsDet &d);

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  private:
    /// Assignment operator is not allowed
    CsDet &             operator =              (const CsDet &d);

  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

  public:

    /// Find detector with given identification
    static CsDet       *FindDetector            (const CS::DetID &id);

    /// Find detector with given name
    static CsDet       *FindDetector            (const std::string &name);

    /// \return Detector identification
    const CS::DetID         GetID                   (void) const {return id;}

    /// \return detector's name
    const char         *getName                 (void) const {return id.GetName().c_str();}

    /*! \return technical board name

       See http://wwwcompass.cern.ch/compass/tech-board/detname.html
    */
    const std::string&       GetTBName               (void) const {return TBname;}

    /// Decode the MC data for this detector
    virtual void        makeMCDecoding          (void) = 0;

    /// Clusterize the digits
    virtual void        clusterize              (void) = 0;


    /*! \brief Decode real data for given detector.

        This function will be pure virtual a bit later.
    */
    ///  Provide R/O info from DaqDataDecoding(mapping files) to Det(s) (if needed)
    virtual void SetDaqDataDecodingInfoGeneral( const std::multimap<CS::Chip::DataID,CS::Chip::Digit*> &daqmap ) {}
    virtual void        DecodeChipDigit         (const CS::Chip::Digit &digit) {}
    /// To provide specific decoding in LED event ( SADC decoding in calorimeters using pulse shape)
    virtual void        DecodeChipDigitLEDEvent  (const CS::Chip::Digit &digit) { DecodeChipDigit(digit);}

    virtual void        DecodeChipDigits          (const CS::Chip::Digits &digits);
    /// To provide specific decoding in LED event ( SADC decoding in calorimeters using pulse shape)
    virtual void        DecodeChipDigitsLEDEvent  (const CS::Chip::Digits &digits) { DecodeChipDigits(digits);}

    /*! \brief Add MonteCarlo hit to detector
        \return \b number of detectors which accepted the hit

        This method will try the MC hit to every detector from list all_detectors.
    */
    static unsigned     AddMCHitAll             (int detector_number,const void *data);

    static std::map<std::string,CsDet*> &GetAllDetectors      (void) {return all_detectors;}

    /*! \brief Add MonteCarlo hit to detector
        \return \b true if hit was added

        This will be pure virtual function in future.
    */
    virtual bool        AddMCHit                (int detector_number,const void *data);

    /// \returns list of pointers to the digits associated to this detector
    const std::list<CsDigit*> &getMyDigits           (void) const {return myDigits_;}

    /// Monitoring developmnet, please dont use this function, contact Vladimir.Kolosov@cern.ch 
    std::list<CsDigit*> &getMyDigits4MN           (void) {return myDigits_;}

    /// \c true if this detector data was arleady clusterized clusterize()
    bool                clusterized             (void) const { return( clusteringDone_ ); }

    /// set this detector as clusterized
    void                setClusteringDone       (void) { clusteringDone_ = true; }

    /// Clear the list of detector clusters.
    void                clearClusterList        (void) { clusteringDone_ = false; myClusters_.clear(); }

    /*! \brief Add a cluster to this detector.
        \param clus the cluster to be added
    */
    void                addCluster              ( CsCluster& clus ) { myClusters_.push_back( &clus ); }

    void                getAssociatedClusters   (std::list<CsCluster*> &clusters) const;

    /// sort myClusters_ according to cluster->getU();
    void                sortClusters            (void);

    /*! \brief Add Digit to this detector.
        \param digit the digit to be added
    */
    void                addDigit                (CsDigit& digit) {myDigits_.push_back(&digit);}

    /// \returns list of pointers to the clusters associated to this detector
    const std::list<CsCluster*>& getMyClusters       (void) const { return( myClusters_ ); }

    /// \returns list of pointers to the clusters associated to this detector
    std::list<CsCluster*>&       getMyClusters       (void)       { return( myClusters_ ); }

    /// update clusters related to this detector for given track. \return false if no clusters were updated
    virtual bool        updateClustersForTrack  (CsCluster &,const CsTrack &) {return false;}

    /// Clear detector. Be ready for next event.
    virtual void        Clear                   (void);

    /*! \brief  Read calibration data.
          Every detector should overwrite this method.
        Example:
        \code
          class MyDet: public CsDet
          {
            ...
            vector<int> calib_data;
            double      another_calib_data;
          }

          MyDet::ReadCalib(const tm &point)
          {
            ReadFromDataBase("MyDet-calib", calib_data,        point);
            ReadFromDataBase("MyDet-calib2",another_calib_data,point);
          }
       \endcode
    */
//     virtual void        ReadCalib               (const tm &point) {}

    virtual void readCalibration(time_t timePoint);

    /*! \brief Write data to a data base.
         Write data 'v' with name 'container' and valide time interval [start,end]
         to the data base specified in 'DataBase' option of CORAL-config file (see below).
         Example:
        \code
          vector<float> v;
          v.push_back(111);
          v.push_back(-2.2);

          tm t[2];  // Start,Finish

          // Start time
          t[0].tm_year = 2000 - 1900;     // year 2000
          t[0].tm_mon  = 4    - 1;        // month number 4 - April
          t[0].tm_mday = 25;              // day 25
          t[0].tm_hour = 0;
          t[0].tm_min  = 0;
          t[0].tm_sec  = 0;

          // Finish time
          t[1].tm_year = 2000 - 1900;
          t[1].tm_mon  = 4    - 1;
          t[1].tm_mday = 25;
          t[1].tm_hour = 23;
          t[1].tm_min  = 59;
          t[1].tm_sec  = 59;

          WriteToDataBase("MyContainer",v,t[0],t[1]);
        \endcode
    */
    virtual void readMCCalibration(time_t timePoint);

    //LS                                                                                            
    virtual void readMCEffMaps(time_t timePoint);

    template <class T>
    void                WriteToDataBase         (const std::string &container,const T& v, const tm &start, const tm &end) const;

    /*! \brief Read data from data base

         Read element  'v' with name 'container' and from the data base. The time 'point' should
         be inside the validity interval of the element 'v'.
         Example:
         \code
           int i;
           tm t;
           t.tm_year = ... // Initialise 't'
           ReadFromDataBase("read_me",i,t);
         \endcode
    */
    template <class T>
    void                ReadFromDataBase        (const std::string &container,T& v, const tm &point) const;

  // ---------------------------------------------------------------------------
  // Attributes, data
  // ---------------------------------------------------------------------------

  private:

    /// Detector's identification
    CS::DetID               id;

    std::string              TBname;

  protected:

    /// digits associated to this detector
    std::list<CsDigit*>      myDigits_;

    /// Clusters associated to this detector
    std::list<CsCluster*>    myClusters_;

    CDB *cdb_;
    std::string folderSet_;

    CDB *mccdb_;

  private:

    /// true if Detector clusterized
    bool                clusteringDone_;

    struct              sortClusters_;

    /// This is list of all detectors.
    static std::map<std::string,CsDet*> all_detectors;
};


//------------------------------------------------------------------------------

#ifndef OPERATOR_ISTR_VECTOR_T
#define OPERATOR_ISTR_VECTOR_T

template <class T>
std::istream& operator>>(std::istream &in,std::vector<T> &v)
{
  v.clear();

  do
    v.push_back(T());
  while( in>>v.back() );

  v.pop_back();  // Remove last element becase it was not read.

  if( in.eof() )
    in.clear();

  return in;
}

#endif


//------------------------------------------------------------------------------





//------------------------------------------------------------------------------

extern std::ostream& operator<<(std::ostream &o, const tm* t);

//------------------------------------------------------------------------------


#include "Reco/DataBase.h"
extern Reco::DataBase *data_base;


template <class T>
void CsDet::WriteToDataBase         (const std::string &container,const T& v, const tm &start, const tm &end) const
{
  if( data_base==NULL )
    throw CS::Exception("CsDet::WriteToDataBase():  data base is not opened.");
  data_base->Write(container,v,start,end);
}

template <class T>
void CsDet::ReadFromDataBase        (const std::string &container,T& v, const tm &point) const
{
  if( data_base==NULL )
    throw CS::Exception("CsDet::ReadFromDataBase():  data base is not opened.");
  data_base->Read(container,v,point);
}

#endif // CsDet_h
