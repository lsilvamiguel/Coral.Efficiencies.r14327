#include <cstdio>
#include <cassert>
#include "Exception.h"

#warning "Improve Print() style: idents for multiline messages, add IsPrintable() method"

namespace CS {

std::map<const std::string,unsigned int> Exception::map__format_n, Exception::map__name_n;
std::map<const std::string,         int> Exception::map__name_level;
size_t Exception::memory_max=10000000;
size_t Exception::memory_cur=0;

////////////////////////////////////////////////////////////////////////////////

Exception::Exception(const char *name,...)
{
    if( name!=NULL && name[0]!=0 )
    {
        va_list ap;
        va_start(ap,name);
        SetName(name,ap);
        va_end(ap);
    }
}

////////////////////////////////////////////////////////////////////////////////

Exception::Exception(const char *name,va_list &ap)
{
    if( name!=NULL && name[0]!=0 )
        SetName(name,ap);
}

////////////////////////////////////////////////////////////////////////////////

void Exception::SetName(const char *name, va_list &ap)
{
    char s[11111];

    if( memory_cur>memory_max )
    {
        memory_cur=memory_max;
        format = "-- Exception --   all %d bytes of memory were consumed! See SetMemoryMax()";
        sprintf(s,format.c_str(),memory_max);

        static bool first=true;
        if( first )
        {
            first=false;
            std::cerr << s << "\n";
            std::cerr << "\n";
            std::cerr << "XXXXXXXXXXXXXXXXXXXXXXXZXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
            std::cerr << "XXX  ATTENTION! YOU HITTED MAXIMUM MEMORY LIMIT IN ERRORS HANDLING!!!  XXX\n";
            std::cerr << "XXX            NEW ERROR MESSAGES WILL *NOT* BE REPORTED               XXX\n";
            std::cerr << "XXXXXXXXXXXXXXXXXXXXXXXZXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
      }
    }
    else
    {
        format=name;
        vsnprintf(s,sizeof(s)-1,name,ap);
        s[sizeof(s)-1]=0;
    }

    assert(strlen(s)<=sizeof(s));

    the_name=s; // save the full error message

    size_t n1=map__format_n.size(), n2=map__name_n.size(); // map sizes

    if( map__name_level.count(format)>0 && map__name_level[format]<0 )
        map__format_n [format+"    [full output is suppressed]"]++;
    else
    {
        map__format_n [format  ]++;
        map__name_n   [the_name]++;
    }

    if( n1!=map__format_n.size() )
        memory_cur += format.length()+4;

    if( n2!=map__name_n.size() )
        memory_cur += the_name.length()+4;
}

////////////////////////////////////////////////////////////////////////////////

void Exception::Print(std::ostream &o,const std::string &prefix) const
{
    if( map__name_level.count(GetFormat())>=0 && map__name_level[GetName()]>=0 )
        o << prefix << what() << "\n";
}

////////////////////////////////////////////////////////////////////////////////

void Exception::PrintStatistics(std::ostream &o,const std::string &prefix,unsigned int events_norm,const std::string &options)
{
    o << "\n";
    o << prefix << "Exception::PrintStatistics():\n\n";
  
    if( events_norm==0 )
        events_norm=(unsigned)-1; // The big number

    if( events_norm!=(unsigned)-1 )
        o << prefix << "Error messages were normalised on " << events_norm << " events\n";
  
    const char format[] = "%10.2f%%  %9d message(s): %s\n";

    if( options!="reduced" )
    {
        o << prefix << " ---- Errors full names ----\n";

        for( std::map<const std::string,unsigned int>::const_iterator it=map__name_n.begin(); it!=map__name_n.end(); it++ )
        {
            char s[strlen(format)+100+it->first.length()];
            sprintf(s,format,it->second*100./events_norm,it->second,it->first.c_str());
            o << s;
        }
        o << "\n";
    }

    o << prefix << " ---- Errors formats ----\n";

    for( std::map<const std::string,unsigned int>::const_iterator it=map__format_n.begin(); it!=map__format_n.end(); it++ )
    {
        char s[strlen(format)+100+it->first.length()];
        sprintf(s,format,it->second*100./events_norm,it->second,it->first.c_str());
        o << s;
    }
    o << "\n";

    if( events_norm!=(unsigned)-1 )
        o << prefix << "Error messages were normalised on " << events_norm << " events\n";

    //o.form("%sAbout %d bytes were used for storing all messages.\n",prefix.c_str(),memory_cur);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
