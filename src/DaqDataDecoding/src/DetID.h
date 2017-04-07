/*!
   \file    DetID.h
   \brief   Detector identification.
   \author  Alexander Zvyagin
*/

#ifndef CompassSoft__DetID_h
#define CompassSoft__DetID_h

#include <iostream>
#include <string>

namespace CS {

/*! \brief COMPASS detector identification

    \author Alexander Zvyagin
*/
class DetID
{
  // ---------------------------------------------------------------------------
  // Constructors, destructor
  // ---------------------------------------------------------------------------

  public:
  
    /// Destructor
    virtual            ~DetID                   (void) {}

    /// Base constructor
                        DetID                   (const std::string &name_,unsigned n=unsigned(-1)) : name(name_), number(n) {}

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

  public:

    /// Test for equality
    bool                operator ==             (const DetID &id) const;

    /// Test for inequality
    bool                operator !=             (const DetID &id) const {return !(id==*this);}

    /// Print to output stream
    friend std::ostream     &operator <<             (std::ostream &o, const DetID &id) {id.Print(o);return o;}
    
                        operator int            (void) const {return int(number);}

  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

  public:

    /// Print detector's identification to output stream.
    void                Print                   (std::ostream &o=std::cout,const std::string &prefix="") const {o<<prefix<<name<<" "<<number;}

    /// \return detector's name
    const std::string       &GetName                 (void) const {return name;}
    
    /// \return detector's number
    unsigned            GetNumber               (void) const {return number;}

    void                Clear                   (void) {name="DetUnknown"; number=unsigned(-1);}

  // ---------------------------------------------------------------------------
  // Attributes, data
  // ---------------------------------------------------------------------------

  private:
  
    /// This is name of the detector
    std::string              name;
    
    unsigned            number;
};

////////////////////////////////////////////////////////////////////////////////

inline bool DetID::operator == (const DetID &id) const
{
    if( name!=id.name )
        return false;

    if( number!=unsigned(-1) && id.number!=unsigned(-1) )
        return number==id.number;
    
    return true;
}

////////////////////////////////////////////////////////////////////////////////

/// Test for inequality
inline bool operator <(const DetID &id1,const DetID &id2)
    {return id1.GetName()<id2.GetName();}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // CompassSoft__DetID_h
