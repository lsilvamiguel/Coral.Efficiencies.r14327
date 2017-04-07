/*!
   \file    CsHome.h
   \brief   CORAL persistence object.
   \version 0.1
   \author  Alexander Zvyagin
   \date    March 2000
*/

#ifndef CsHome_h
#define CsHome_h

#include "coral_config.h"

#include <string>

/*! \class CsHome
    \author Alexander Zvyagin
    \brief Interface to CORAL persistence object.
    
    This is abstract class with general properties of persistence object.
    The base properties of a CsHome object are:
    \arg 1. CsHome object has name
    \arg 2. CsHome object has title
    \arg 3. there is path to CsHome object (it is located somewere)
    
    \sa CsHomeNamed
*/
class CsHome
{
  // ---------------------------------------------------------------------------
  // Constructors, destructor.
  // ---------------------------------------------------------------------------

  public:
  
    /// The destructor.
    virtual            ~CsHome                  (void) {}

    /// Copy constructor.
                        CsHome                  (const CsHome &home);

    /// Base constructor.
                        CsHome                  (void) {}
//                        CsHome                  (const std::string &name="Unnamed",const std::string &title="Untitled",const std::string &path=".") {}

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

    CsHome             &operator =              (const CsHome &home);
    
  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

    /// \return an object's name.
    virtual std::string      GetName                 (void) const            = 0;
    
    /// \return an object's title.
    virtual std::string      GetTitle                (void) const            = 0;
    
    /// \return path to an object.
    virtual std::string      GetPath                 (void) const            = 0;
    
    /// Set object's name.
    virtual void        SetName                 (const std::string &name_ )  = 0;

    /// Set object's title.
    virtual void        SetTitle                (const std::string &title_)  = 0;
    
//    virtual void          Rename                (const std::string &name_,const std::string &title_,const std::string &path_) = 0;
    
//    virtual void          Rename                (const std::string &name_,const std::string &title_);

    /*! \brief Set path to object.
        Examples:
        \arg "/"           Set object in top directory.
        \arg "."           Set object in current directory.
        \arg "dir"         Set object in subdirectory "dir" from current directory.
        \arg "/dir1/dir2"  Set object in absolute directory.
        \arg ""            An object will not be connected to (and will not be stored in) any directory.
    */
    virtual void        SetPath                 (const std::string &path) = 0;
    
//     /// \return name of package responsible for object saving.
//     virtual std::string       &PersistenceType       (void) const = 0;

//     /// \return current path.
//     static std::string       GetCurrentPath          (const std::string &package);
//     
//     /// \return current path.
//     static void         SetCurrentPath          (const std::string &path,const std::string &package);
//     
//     /// \return \b true if \c name is valid.
//     static bool         CheckName               (const std::string &name);
};

////////////////////////////////////////////////////////////////////////////////

/*! \class CsHomeNamed
    \author Alexander Zvyagin
    \brief Interface to CORAL persistence object.
    
    This is concrete class of CsHome.
*/
class CsHomeNamed : public CsHome
{
  // ---------------------------------------------------------------------------
  // Constructors, destructor.
  // ---------------------------------------------------------------------------

  public:
  
    /// The destructor.
    virtual            ~CsHomeNamed             (void) {}

    /// Copy constructor.
                        CsHomeNamed             (const CsHomeNamed &home);

    /// Base constructor.
                        CsHomeNamed             (const std::string &name_="Unnamed",const std::string &title_="Untitled",const std::string &path_=".") :
                                                 name(name_), title(title_), path(path_) {}

  // ---------------------------------------------------------------------------
  // Operators
  // ---------------------------------------------------------------------------

    CsHomeNamed        &operator =              (const CsHomeNamed &home);
    
  // ---------------------------------------------------------------------------
  // Methods
  // ---------------------------------------------------------------------------

    /// \return an object's name.
    virtual std::string      GetName                 (void) const            {return name;}
    
    /// \return an object's title.
    virtual std::string      GetTitle                (void) const            {return title;}
    
    /// \return path to an object.
    virtual std::string      GetPath                 (void) const            {return path;}
    
    /// Set object's name.
    virtual void        SetName                 (const std::string &name_ )  {name=name_;}

    /// Set object's title.
    virtual void        SetTitle                (const std::string &title_)  {title=title_;}
    
//    virtual void          Rename                (const std::string &name_,const std::string &title_,const std::string &path_) = 0;
    
//    virtual void          Rename                (const std::string &name_,const std::string &title_);

    /*! \brief Set path to object.
        Examples:
        \arg "/"           Set object in top directory.
        \arg "."           Set object in current directory.
        \arg "dir"         Set object in subdirectory "dir" from current directory.
        \arg "/dir1/dir2"  Set object in absolute directory.
        \arg ""            An object will not be connected to (and will not be stored in) any directory.
    */
    virtual void        SetPath                 (const std::string &path_) {path=path_;}
    
//     /// \return name of package responsible for object saving.
//     virtual std::string       &PersistenceType       (void) const = 0;

//     /// \return current path.
//     static std::string       GetCurrentPath          (const std::string &package);
//     
//     /// \return current path.
//     static void         SetCurrentPath          (const std::string &path,const std::string &package);
// 
  private:
    
    std::string name,title,path;
};

#endif // CsHome_h
