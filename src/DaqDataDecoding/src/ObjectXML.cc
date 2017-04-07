#include <iostream>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <expat.h>
#include <string.h>
#include "ObjectXML.h"
#include "Exception.h"

using namespace std;

namespace CS {

////////////////////////////////////////////////////////////////////////////////

ObjectXML::~ObjectXML(void)
{
    KillChildren();
    RemoveFromParent();
}

////////////////////////////////////////////////////////////////////////////////

ObjectXML::ObjectXML(const string &_name,ObjectXML* _parent) :
  parent(NULL),
  name(_name)
{
    SetParent(_parent);
}

////////////////////////////////////////////////////////////////////////////////

ObjectXML::ObjectXML(const ObjectXML &o) :
  parent(NULL)
{
    *this = o;
}

////////////////////////////////////////////////////////////////////////////////

ObjectXML& ObjectXML::operator=(const ObjectXML &o)
{
    KillChildren();

    name        = o.name;
    attributes  = o.attributes;
    body        = o.body;
    comments    = o.comments;

    SetParent(o.parent);
    
    // Now copy the tree of children
    for( list<ObjectXML*>::const_iterator it=o.children.begin();
         it!=o.children.end(); it++ )
    {
        ObjectXML *c = new ObjectXML(**it);
        c->SetParent(this);
    }

    return *this;
}

////////////////////////////////////////////////////////////////////////////////

void ObjectXML::KillChildren(void)
{
    while( !children.empty() )
        delete children.front();
    children.clear();
}

////////////////////////////////////////////////////////////////////////////////

void ObjectXML::RemoveFromParent(void) const
{
    if( parent==NULL )
        return;

    // serach for the object in parents' list of objects
    for( list<ObjectXML*>::iterator it=parent->children.begin();
         it!=parent->children.end(); it++ )
        if( this==*it )
        {
            // The child "this" was found.
            // Remove the child from parents children list

            parent->children.erase(it);
            return;
        }

    throw Exception("ObjectXML::RemoveFromParent(): Internal error. Can not find the child! Memory leak?");
}

////////////////////////////////////////////////////////////////////////////////

void ObjectXML::SetParent(ObjectXML *p)
{
    RemoveFromParent();
    parent = p;
    if( parent!=NULL )
        parent->children.push_back(this);
}

////////////////////////////////////////////////////////////////////////////////

void ObjectXML::AddChildren(ObjectXML *child)
{
    if( child==NULL )
        return;

    if( children.end()!=find(children.begin(),children.end(),child) )
        throw Exception("ObjectXML::AddChildren(): two identical children: object=\"%s\" child=\"%s\"",
                         GetName().c_str(),child->GetName().c_str());

    //children.push_back(child);
    child->SetParent(this);
}

////////////////////////////////////////////////////////////////////////////////

void ObjectXML::Print(ostream &o,const string &prefix) const
{
  o << prefix << "<" << name;

  for( list<Attribute>::const_iterator it=attributes.begin(); it!=attributes.end(); it++ )
    o << " " << it->first << "=\"" << it->second << "\"";
    
  if( GetBody().empty() && GetChildren().empty() )
  {
    o << "/>\n"; // no body, no children... Save a line of output:    <name ...... />
    return;
  }

  o << ">\n";
  
  for( list<string>::const_iterator it=GetBody().begin(); it!=GetBody().end(); it++ )
    o << prefix << (*it) << "\n";

  for( list<ObjectXML*>::const_iterator it=GetChildren().begin(); it!=GetChildren().end(); it++ )
    (*it)->Print(o,prefix+"    ");
  
  o << prefix << "</" << name << ">\n";
}

////////////////////////////////////////////////////////////////////////////////

void
ObjectXML::CopyAttributes(list<Attribute> &attributes,const char **atts)
{
  for( int i=0; atts[i]!=0; i+=2 )
    if( atts[i+1]==0 )
    {
      cerr << "XML:startElement(): an attribute is zero!!\n";
      break;
    }
    else
      attributes.push_back(Attribute(atts[i],atts[i+1]));
}

////////////////////////////////////////////////////////////////////////////////

bool
ObjectXML::IsEmpty(const string &s)
{
  char ss[s.length()+1];
  return 1!=sscanf(s.c_str(),"%s",ss);
}

////////////////////////////////////////////////////////////////////////////////

class UserData
{
  public:
    ObjectXML *me;
};

////////////////////////////////////////////////////////////////////////////////

// XML_StartElementHandler
void
startElement(void *userData, const char *name, const char **atts)
{
  try
  {
    UserData &u = *reinterpret_cast<UserData*>(userData);
    ObjectXML *o=new ObjectXML(name,u.me);
    u.me=o;
    ObjectXML::CopyAttributes(o->GetAttributes(),atts);
  }
  catch( const Exception &e )
  {
    e.Print(cerr);
  }
  catch( ... )
  {
    cerr << "startElement: unknown exception!\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

// XML_EndElementHandler
void
endElement(void *userData, const char *name)
{
  try
  {
    UserData &u = *reinterpret_cast<UserData*>(userData);

    if( u.me==NULL )
      throw Exception("XML::endElement(): objects list is empty!");

    if( u.me->GetName()!=name )
      throw Exception("XML::endElement(): Expect \"%s\" got \"%s\"",
                      u.me->GetName().c_str(),name);

    u.me = u.me->GetParent();
  }
  catch( const Exception &e )
  {
    e.Print(cerr);
  }
  catch( ... )
  {
    cerr << "endElement: unknown exception!\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

// XML_DefaultHandler
void
default_handler(void *userData,const XML_Char *s,int len)
{
  try
  {
    UserData &u = *reinterpret_cast<UserData*>(userData);

    string ss(s,len);
    if( !ObjectXML::IsEmpty(ss) )
    {
      if( u.me==NULL )
        throw Exception("XML::default_handler(): objects list is empty!");

      u.me->GetBody().push_back(ss);
    }
  }
  catch( const Exception &e )
  {
    e.Print(cerr);
  }
  catch( ... )
  {
    cerr << "default_handler: unknown exception!\n";
  }
}

////////////////////////////////////////////////////////////////////////////////

// XML_CommentHandler
void
comment_handler(void *userData, const XML_Char *data)
{
    UserData &u = *reinterpret_cast<UserData*>(userData);
    if( u.me!=NULL )
        u.me->GetComments().push_back(data);
}

////////////////////////////////////////////////////////////////////////////////

void
ObjectXML::Parse(const string &path,ObjectXML &objs)
{
  // First of all we try to open the "path"
  struct stat info;
  if( 0!=stat(path.c_str(),&info) )
    throw Exception("ObjectXML::Parse(): bad path \"%s\"",path.c_str());

  if( S_ISDIR(info.st_mode) )
  {
    // The path is a directory.
    // Scan that directory for files with the extension ".xml"

    struct dirent **namelist;
    int n = scandir(path.c_str(),&namelist,0,0);

    if (n < 0)
      throw Exception("ObjectXML::Parse(): can not scan directory \"%s\"",path.c_str());

    while(n--)
    {
      const char *name=namelist[n]->d_name;  // file name
      if( strlen(name)>4 &&
          0==strcmp(name+strlen(name)-4,".xml") )
      {
        string s=path+'/'+name;
        try
        {
          Parse(s,objs);
        }
        catch( const std::exception &e )
        {
          cerr << e.what() << "\n";
        }
        catch( const char *s )
        {
            cerr << s << "\n";
        }
        catch( ... )
        {
          cerr << "ObjectXML::Parse(): Unknown problem in the map file \"" << s << "\"\n";
        }
      }
      free(namelist[n]);
    }
    free(namelist);
  }
  else
  if( S_ISREG(info.st_mode) )
    ObjectXML::ParseFile(path.c_str(),objs);
  else
    throw Exception("ObjectXML::Parse(): not a file or directory: \"%s\"",path.c_str());
}

////////////////////////////////////////////////////////////////////////////////

void ObjectXML::ParseFile(const string &path,ObjectXML &objs)
{
    ifstream in(path.c_str());
    if( !in.is_open() )
    {
        fprintf(stderr,"ObjectXML::Parse(): can not open file \"%s\"",path.c_str());
        return;
    }

    if( Parse(in,objs) )
        printf("File \"%s\" contains error(s), see above.\n",path.c_str());
}

////////////////////////////////////////////////////////////////////////////////

int ObjectXML::Parse(istream &in,ObjectXML &objs)
{
  int errors=0;
  const unsigned buf_size=1000000;
  char *buf = new char[buf_size];

  UserData u;
  u.me=&objs;

  XML_Parser parser = XML_ParserCreate(NULL);
  XML_SetUserData                       (parser, &u);
  XML_SetElementHandler                 (parser, startElement, endElement);
//XML_SetCharacterDataHandler           (parser, characters_handler);
//XML_SetProcessingInstructionHandler   (parser, instruction_handler);
  XML_SetDefaultHandler                 (parser, default_handler);
  XML_SetCommentHandler                 (parser, comment_handler);

  try
  {
    in.read(buf,buf_size);
    if( size_t(in.gcount())==buf_size )
      throw Exception("ObjectXML::Parse(): too short buffer.");

    int done=0;
    if( !XML_Parse(parser, buf, in.gcount(), done) )
    {
      Exception("ObjectXML::Parse(): Error: %s  at line %d",
                XML_ErrorString(XML_GetErrorCode(parser)),
                XML_GetCurrentLineNumber(parser)).Print();
      errors++;
    }
  }
  catch(const Exception &e)
  {
    errors++;
    e.Print(cerr);
  }
  catch(...)
  {
    errors++;
    cerr << "ObjectXML::Parse(): unknown error!\n";
  }

  XML_ParserFree(parser);
  delete [] buf;
  return errors;
}

} // namespace CS
