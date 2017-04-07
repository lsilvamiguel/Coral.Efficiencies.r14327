/*
   Authors:
     Alexander Zvyagin   ( Alexander.Zvyagin@cern.ch)

   Copyright(C): 2000-2002  A.Zvyagin

     This library is free software; you can redistribute it and/or
     modify it under the terms of the GNU Library General Public
     License as published by the Free Software Foundation; either
     version 2 of the License, or (at your option) any later version.

     This library is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
     Library General Public License for more details.

     You should have received a copy of the GNU Library General Public
     License along with this library; if not, write to the Free
     Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef CompassSoft_Exception__include
#define CompassSoft_Exception__include

#include <exception>
#include <string>
#include <string.h>
#include <map>
#include <memory>
#include <vector>
#include <iostream>
#include <cstdarg>

namespace CS {

/*! \brief Class for throwing a C++ exceptions.

    \author Alexander Zvyagin
*/
class Exception : public std::exception
{
  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:

    /// Destructor (do nothing)
    virtual            ~Exception               (void) throw () {}

    /*! \brief base Constructor

        Example:
        \code
          if( i!=1 )
            throw Exception("Bad value 'i':  i=%d",i);
        \endcode
    */
                        Exception               (const char *name,...);

                        Exception               (const char *name,va_list &ap);

  //============================================================================
  // Static methods
  //============================================================================

  public:

    /*! \return exception's level (0 on default).
    */
    static
    int                 GetLevel                (const std::string &name) {return map__name_level.count(name)>0?map__name_level[name]:0;}
    
    /*! \brief Set exception level.
    
        If level==-10 the only exception format will be added to map__format_n,
        otherwise map__name_n also will be updated.
        
        Use this method to suppress looooong list of similar exceptions.
    */
    static
    void                SetLevel                (const std::string &name,int level) {map__name_level[name]=level;}
    
  //============================================================================
  // Methods
  //============================================================================

  public:

    /// \return exception name. Method from std::exception
    virtual const char *what                    (void) const throw () {return GetName().c_str();}

    /// \return exception format.
    const std::string  &GetFormat               (void) const {return format;}

    /// \return exception name.
    const std::string  &GetName                 (void) const {return the_name;}

    /// Print exception if it's level >0
    virtual void        Print                   (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// Print exceptions statistics
    static void         PrintStatistics         (std::ostream &o=std::cout,const std::string &prefix="",
                                                 unsigned int events_norm=0,const std::string &options="");

    /// Set maximum memory used by the Exception class.
    static void         SetMemoryMax            (size_t m) {memory_max=m;}

  protected:
  
    /// Set the exception name.
    void                SetName                 (const char *name,va_list &ap);

  //============================================================================
  // Attributes
  //============================================================================

  private:

    /// Exception's format.
    std::string        format;

    /// Exception's name.
    std::string        the_name;

    /// Map of used formats
    static std::map<const std::string,unsigned int> map__format_n;

    /// Map of used exceptions
    static std::map<const std::string,unsigned int> map__name_n;
    
    /// Map of exception's levels.
    static std::map<const std::string,         int> map__name_level;
    
    /// Maximum memory for error messages.
    static size_t                         memory_max;

    /// Current memory consumption.
    static size_t                         memory_cur;
};

} // namespace CS

#endif // CompassSoft_Exception__include
