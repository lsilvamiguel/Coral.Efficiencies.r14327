#ifndef ObjectXML___include
#define ObjectXML___include

#include <iostream>
#include <sstream>
#include <string>
#include <list>
#include <map>
#include "Exception.h"

namespace CS {

class ObjectXML
{
  //============================================================================
  // Types
  //============================================================================

  public:
  
    typedef std::pair<std::string,std::string> Attribute;

  //============================================================================
  // Constructors and destructor
  //============================================================================

  public:
    
    virtual            ~ObjectXML               (void);

                        ObjectXML               (const ObjectXML &o);

                        ObjectXML               (const std::string &_name,ObjectXML *parent=NULL);

  //============================================================================
  // Operators
  //============================================================================

  public:
  
    /// Duplicate the full tree of objects.
    ObjectXML&          operator =              (const ObjectXML &o);
    
  //============================================================================
  // Static methods
  //============================================================================

  public:

    static void         CopyAttributes          (std::list<Attribute> &attributes,const char **atts);

    static bool         IsEmpty                 (const std::string &s);

    // Parse an XML file
    static void         ParseFile               (const std::string &file,ObjectXML &objs);
    
    // Parse file or directory with the ".xml" extension
    static void         Parse                   (const std::string &file,ObjectXML &objs);

    // Parse a stream
    static int          Parse                   (std::istream &in,ObjectXML &objs);

  //============================================================================
  // Methods
  //============================================================================
  
  public:

    virtual void            Print               (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return the object name
    virtual const std::string&   GetName        (void) const {return name;}

    /// \return the object name
    virtual void            SetName             (const std::string &n) {name=n;}
    
    const std::list<ObjectXML*>& GetChildren    (void) const {return children;}

          std::list<ObjectXML*>& GetChildren    (void)       {return children;}

    void                    AddChildren         (ObjectXML *child);
    
    const std::list<Attribute>&  GetAttributes  (void) const {return attributes;}

          std::list<Attribute>&  GetAttributes  (void)       {return attributes;}

    const std::list<std::string>&     GetBody   (void) const {return body;}

          std::list<std::string>&     GetBody   (void)       {return body;}

    const std::list<std::string>&     GetComments  (void) const {return comments;}

          std::list<std::string>&     GetComments  (void)       {return comments;}

    /// \return pointer to parent ObjectXML object (may be NULL)
    const ObjectXML*        GetParent           (void) const {return parent;}

    /// \return pointer to parent ObjectXML object (may be NULL)
          ObjectXML*        GetParent           (void)       {return parent;}

    /// Set the parent ObjectXML for this object.
    void                    SetParent           (ObjectXML *p);
    
    /*! \brief Get the attribute value with name \c name
    
        \return pointer to the found attribute (and filled \b x argument) if \b name was found
        \return \b NULL if there is no attribute with this name
    */
    template<class T>
    const Attribute*        GetAttribute        (const std::string &name,T &x) const {return const_cast<ObjectXML*>(this)->GetAttribute(name,x);}

    template<class T>
          Attribute*        GetAttribute        (const std::string &name,T &x);
    
    template<class T>
    void                    SetAttribute        (const std::string &name,const T &x);

  private:

    /// Remove the object from parent's list of objects
    void                    RemoveFromParent    (void) const;

    /// Kill all children and release all memory allocated by them.
    void                    KillChildren        (void);

  //============================================================================
  // Attributes
  //============================================================================

  protected:

    /// Pointer to the parent object
    ObjectXML*               parent;

    /// List of children objects
    std::list<ObjectXML*>    children;

  private:

    /// Name of the XML object: <name attr1=v1 attr2=v2/>
    std::string              name;

    /// List of attributes
    std::list<Attribute>     attributes;

    /// The body
    std::list<std::string>   body;

    /// All comments messages
    std::list<std::string>   comments;
};

////////////////////////////////////////////////////////////////////////////////

template<class T>
ObjectXML::Attribute* ObjectXML::GetAttribute(const std::string &n,T &x)
{
  for( std::list<Attribute>::iterator it=attributes.begin(); it!=attributes.end(); it++ )
    if( it->first==n )
    {
      std::istringstream s(it->second);

      s>>x;
      if( s.fail() )
        throw Exception("ObjectXML::GetAttribute(): bad format: ObjectXML::name=\"%s\"  %s=\"%s\"",
                         GetName().c_str(),it->first.c_str(),it->second.c_str());

      return &*it;  // return pointer to the found attribute
    }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////

template<> inline
ObjectXML::Attribute* ObjectXML::GetAttribute<std::string>(const std::string &n,std::string &x)
{
    for( std::list<Attribute>::iterator it=attributes.begin(); it!=attributes.end(); it++ )
        if( it->first==n )
        {
            x=it->second;
            return &*it;  // return pointer to the found attribute
        }

  return NULL;
}

////////////////////////////////////////////////////////////////////////////////

template<class T>
void ObjectXML::SetAttribute(const std::string &n,const T &x)
{
    std::string *s=NULL;
    
    for( std::list<Attribute>::iterator it=attributes.begin(); it!=attributes.end(); it++ )
        if( it->first==n )
        {
            s = &it->second;
            break;
        }

    if( s==NULL )
    {
        attributes.push_back(Attribute(n,""));
        s = &attributes.back().second;
    }

    std::string sbuf;
    std::ostringstream ss(sbuf);

    ss << x;
    if( ss.fail() )
        throw Exception("ObjectXML::Set(): can not set attribute \"%s\" for name \"%s\"",
                         n.c_str(),GetName().c_str());
    *s = ss.str();
}

////////////////////////////////////////////////////////////////////////////////

template<> inline
void ObjectXML::SetAttribute<std::string>(const std::string &n,const std::string &x)
{
    for( std::list<Attribute>::iterator it=attributes.begin(); it!=attributes.end(); it++ )
        if( it->first==n )
        {
            it->second = x;
            return;
        }

    attributes.push_back(Attribute(n,x));
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS

#endif // ObjectXML___include
