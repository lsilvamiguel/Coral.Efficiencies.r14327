/*!
   \file    CsHBOOK.h
   \brief   Header file for internal classes of histograms package.
   \version $Revision: 1.9 $
   \author  Alexander Zvyagin
   \date    $Date: 2010/02/03 18:22:23 $
*/

#ifndef CsHBOOK_h
#define CsHBOOK_h

#include "coral_config.h"

#if USE_HBOOK

#include <cassert>
#include <cstdlib>
#include <string>
#include "CsHbookProto.h"
#include "CsHome.h"
#include "CsTypes.h"

////////////////////////////////////////////////////////////////////////////////

class CsHbookObject : public CsHomeNamed
{
  public:

                       ~CsHbookObject           (void) {}

                        CsHbookObject           (const CsHbookObject &h);

                        CsHbookObject           (const std::string &name, const std::string &title, const std::string &path=".")
                                                  : CsHomeNamed(name,title,path), hbook_title(title)
                                                {
                                                  if( std::string::npos!=name.find(" | ") )
                                                    throw "CsHbookObject::CsHbookObject() histogram name can not contain string \" | \"";
                                                
                                                  id = atoi(name.c_str()); // Attempt to extract number from name. example: "123"
                                                  if( id<=0 )
                                                  {
                                                    id = NewID();
                                                    hbook_title = GetName() + " | " + GetTitle();
                                                  }
                                                  if( hbook_title.length()>80 )
                                                    CsErrLog::Instance()->msg(elWarning,__FILE__,__LINE__,
                                                     "CsHbookObject::CreateHbookTitle: HBOOK title will be truncated: \"%s\"",
                                                      hbook_title.c_str());
                                                }

    CsHbookObject&      operator =              (const CsHbookObject &h)
                                                {
                                                  if( this!=&h )
                                                  {
                                                    id = h.id;
                                                    hbook_title = h.hbook_title;
                                                    CsHomeNamed::operator=(h);
                                                  }
                                                  return *this;
                                                }
                                                
    size_t              ID                      (void) const {return id;}
    
    const std::string&       GetHbookTitle           (void) const {return hbook_title;}
    
  protected:

    void                DecodeNameAndTitle      (void)
                                                {
                                                  size_t p=hbook_title.find(" | ");
                                                  if( p==std::string::npos )
                                                  {
                                                    char s[66];
                                                    sprintf(s,"%zu",ID());
                                                    SetName(s);
                                                    SetTitle(hbook_title);
                                                  }
                                                  else
                                                  {
                                                    SetName(hbook_title.substr(0,p));
                                                    SetTitle(hbook_title.substr(p+3));
                                                  }
                                                }
    
  private:
    
    static size_t       NewID                   (void)
                                                {
                                                  for( size_t i=1; i>0; i++ )
                                                    if( !hexist(i) )
                                                      return i;
                                                  CsErrLog::Instance()->mes(elFatal,"CsHbookObject::NewID: Can not find new ID for an HBOOK object!");
                                                  return 0;
                                                }

  protected:
  
    virtual std::string      WhoAmI                  (void) const = 0;

  private:

    size_t              id;
    std::string              hbook_title;
};

////////////////////////////////////////////////////////////////////////////////

class CsHbookHist : public CsHbookObject
{
  public:

    /// Get histogram from HBOOK
//                        CsHbookHist             (int id) : CsHbookObject(id) { Init(id); }

                        CsHbookHist             (const std::string &name, const std::string &title)
                                                : CsHbookObject(name,title) {}

    int nx, ny, nwt;
    float xmin,xmax,ymin,ymax;
    
  protected:

    std::string              WhoAmI                  (void) const {return "Histogram";}

  private:

                        CsHbookHist             (const CsHbookHist &c);

    CsHbookHist        &operator =              (const CsHbookHist &c);
    
    void                Init                    (int id)
                                                {
                                                  assert(id!=0);
                                                  char title[81];
                                                  hgive(id,title,nx,xmin,xmax,ny,ymin,ymax,nwt);
                                                  DecodeNameAndTitle();
                                                }
};

////////////////////////////////////////////////////////////////////////////////

class CsHbookNtuple : public CsHbookObject
{
  public:

    /// Get Ntuple from HBOOK
                        CsHbookNtuple           (const std::string &name, const std::string &title)
                                                : CsHbookObject(name,title), Nvar(0)                                                
                                                {}

    float32            *buffer;
    int                 Nvar;
    
  protected:

    std::string              WhoAmI                  (void) const {return "NTuple";}

  private:

                        CsHbookNtuple           (const CsHbookNtuple &c);

    CsHbookNtuple      &operator =              (const CsHbookNtuple &c);
};

////////////////////////////////////////////////////////////////////////////////

#endif  // USE_HBOOK
#endif  // CsHBOOK_h
