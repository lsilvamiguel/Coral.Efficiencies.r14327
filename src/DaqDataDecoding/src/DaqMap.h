#ifndef CompassSoft_DaqMap__include
#define CompassSoft_DaqMap__include

#include "config.h"
#include <iostream>
#include <string>
#include <vector>
#include <list>
#include "ObjectXML.h"
#include "utils.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief Map class
    
    \author Alexander Zvyagin
*/
class DaqMap: public ObjectXML
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    /// Destructor
    virtual            ~DaqMap                  (void) {}
    
    /*! \brief Base constructor
    */
                        DaqMap                  (const ObjectXML &o);

  //============================================================================
  // Methods
  //============================================================================

  public:

    /// \return Map version
    int                 GetVersion              (void) const {return version;}
    
    const std::string&       GetOptions              (void) const {return options;}
    
    /// Reset all internal variables.
    virtual void        Clear                   (void);

    /*! \return list of runs that are valid for this map

        List of runs is vector and every element of this vector is pair [first_run,second_run]
    */
    const std::vector< std::pair<uint32,uint32> > &
                        GetRunsList             (void) const {return runs_list;}

    /// Test that run is in the range.
    bool                TestRunsRange           (uint32 run) const;
        
    bool                IsOption                (const std::string &o) const {return string_match(options,"\\b"+o+"\\b");}

    /// Set runs list for this map.
    void                SetRunsList             (const std::string &runs_list);

  //============================================================================
  // Attributes
  //============================================================================

  protected:

    /// Version number
    int                 version;
    
    std::string         options;
    
    /// List of valid runs
    std::vector< std::pair<uint32,uint32> > runs_list;
};

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft_DaqMap__include
