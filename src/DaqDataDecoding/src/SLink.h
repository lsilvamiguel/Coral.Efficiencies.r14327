#ifndef CompassSoft_SLink__include
#define CompassSoft_SLink__include

#include <iostream>
#include <string>
#include "config.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

class Chip;
class DaqOption;

/*! \brief Data structure associated with S-Link multiplexer

    Some documentation was extracted from "Compass-Note 2001-8"

    \author Alexander Zvyagin
*/
class SLinkMultiplexer
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /*! \brief Default constructor.
    
        Set all attributes to zero.
    */
                        SLinkMultiplexer        (void);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// Print S-Link info
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;
    
    bool                IsSLinkMultiplexer      (void) const;

    /// \return Size of event in 32-bit words.
    unsigned            GetEventSize            (void) const {return event_size;}

    /// Set event size
    void                SetEventSize            (unsigned v) {event_size=v;}

    /*! \brief \return Equipment number of the CATCH board - distributed via VME to the CATCH during equipment initialisation

        The link source ID is the equipment
        number of the CATCH board sending the data. This number is distributed via VME to all CATCHes during equipment initialisation and stored
        locally. The source ID provides the distinctive identification for the decoding routines of the different front end types. Depending on this 10 bit
        number different encoding routines must be applied to the data to extract hardware channel addresses, time, amplitude or scaler information of
        the data stream which follows the S-Link header.
    */
    uint16              GetSourceID             (void) const {return source_id;}
    
    /// Set source ID
    void                SetSourceID             (uint16 v) {source_id=v;}

    /*! \brief \return Event type
    */
    unsigned            GetEventType            (void) const {return event_type;}
    
    /// Set event type
    void                SetEventType            (unsigned v) {event_type=v;}

    /// \return error flag
    bool                IsError                 (void) const {return error_flag;}
    
    /// Set error flag to \b true or \b false
    void                SetError                (bool v) {error_flag=v;}

    /// \return Number of the event within the spill
    unsigned            GetEventNumber          (void) const {return event_number;}
    
    /// Set event number
    void                SetEventNumber          (unsigned v) {event_number=v;}

    /// \return Number of the spill in which the event was taken
    unsigned            GetSpillNumber          (void) const {return spill_number;}
    
    /// Set spill number
    void                SetSpillNumber          (unsigned v) {spill_number=v;}

    /// \return stat
    unsigned            GetStat                 (void) const {return stat;}
    
    /// Set stat
    void                SetStat                 (unsigned v) {stat=v;}

  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    uint32              event_size     : 16;    ///< size of the event (in words ?)
    uint32              source_id      : 10;    ///< source ID of the catch board
    uint32              event_type     :  5;    ///< event type
    uint32              error_flag     :  1;    ///< this event contain an error

    uint32              event_number   : 20;    ///< event number
    uint32              spill_number   : 11;    ///< spill number
    uint32              stat           :  1;    ///< mode of the TCS receiver
};

////////////////////////////////////////////////////////////////////////////////

/*! \brief Data structure associated with S-Link

    Some documentation was extracted from "Compass-Note 2000-8"
    http://hpfr02.physik.uni-freiburg.de/projects/compass/electronics/notes/dataformat-2000-8/format.html

    \author Alexander Zvyagin
*/
class SLink: public SLinkMultiplexer
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /*! \brief Default constructor.
    
        Set all attributes to zero.
    */
                        SLink                   (void);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// Print S-Link info
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /*! \brief Equipment status word.

        The equipment status word, encoded in the least significant byte, is not defined yet. It may contain
        information like setup ok, timeout, skip event, error in header/trailer, error in data, etc...
    */
    unsigned            GetStatus               (void) const {return status;}
    
    /// Set status word.
    void                SetStatus               (unsigned v) {status=v;}

    /// The TCS error information where all bits equal to zero means no error
    unsigned            GetTCSError             (void) const {return TCS_error;}
    
    /// Set TCS error information
    void                SetTCSError             (unsigned v) {TCS_error=v;}

    /// Count of errors recognised in the event from this particular equipment.
    unsigned            GetErrorsCounter        (void) const {return errors_counter;}
    
    /// Set errors counter
    void                SetErrorsCounter        (unsigned v) {errors_counter=v;}

    /*! \brief Data format description.

        This number is used to discriminate between different readout modes, i.e. debug or sparsified readout, and different encoding
        versions, i.e. possible CATCH software upgrade in year 2001 or later. The format byte may also be used for further specifications of the
        encoding method which must be applied. Further bit 7 of the format word is used to discern the first event in a run which contains detector
        describing information and no measured data (see table).
    */
    unsigned            GetFormat               (void) const {return format;}
    
    //// Set format
    void                SetFormat               (unsigned v) {format=v;}
    
    /// \return \b true if this is first event in run.
    bool                IsFirstEventInRun       (void) const {return GetFormat()&128;}

    /// Analyse S-Link errors and put the result to the chip's errors list.
    void                AnalyseErrors           (Chip &chip,DaqOption &opt) const;
    
  //============================================================================
  // Attributes, data
  //============================================================================

  private:

    uint32              status         :  8;    ///< status info
    uint32              TCS_error      :  8;    ///< error info from TCS
    uint32              errors_counter :  8;    ///< number of error words
    uint32              format         :  8;    ///< format of the data
};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_DaqChip__include
