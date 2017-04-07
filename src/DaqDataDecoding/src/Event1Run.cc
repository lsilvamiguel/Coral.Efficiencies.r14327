#include <cstring>
#include <cstdio>
#include "Event1Run.h"

namespace CS {

using namespace std;

////////////////////////////////////////////////////////////////////////////////

void Event1Run::Add(unsigned int r,const void *buf)
{
    if( r!=run )
        throw "Event1Run::Add(): wrong run number.";

    const SLink &slink=*reinterpret_cast<const SLink *>(buf);
    if( !slink.IsFirstEventInRun() )
        printf("Event1Run::Add(): wrong S-Link format for the first event of run.\n");

    const CatchInfo::Header &header=*reinterpret_cast<const CatchInfo::Header *>(((char*)buf)+sizeof(SLink));
    if( catch_info.count(header.catch_number)>0 )
        printf("Event1Run::Add(): Replacing catch %d info.\n",header.catch_number);
    
    catch_info[slink.GetSourceID()] = new CatchInfo(buf);
}

////////////////////////////////////////////////////////////////////////////////

void Event1Run::Clear(void)
{
    for( map<int,CatchInfo*>::iterator it=catch_info.begin(); it!=catch_info.end(); it++ )
        delete it->second;
    catch_info.clear();
}

////////////////////////////////////////////////////////////////////////////////

void Event1Run::Print(const char *prefix) const
{
    printf("%sThis is the first-event-of-run number %d.\n",prefix,run);

    for( map<int,CatchInfo*>::const_iterator it=catch_info.begin(); it!=catch_info.end(); it++ )
    {
        it->second->Print(prefix);
        printf("%s\n",prefix);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Event1Run::Write(std::ostream &o) const
{
    unsigned int b[2]={run,catch_info.size()};
    if( !o.write(reinterpret_cast<char*>(b),sizeof(b)) )
        throw "Event1Run::Write(): failed to write to the stream";
    
    for( map<int,CatchInfo*>::const_iterator it=catch_info.begin(); it!=catch_info.end(); it++ )
        it->second->Write(o);
}

////////////////////////////////////////////////////////////////////////////////

void Event1Run::Read(std::istream &o)
{
    Clear();

    unsigned int b[2], &size=b[1];
    if( !o.read(reinterpret_cast<char*>(b),sizeof(b)) )
        throw "Event1Run::Read(): failed to read from the stream";
    run=b[0];

    for( unsigned int i=0; i<size; i++ )
    {
        CatchInfo *c=new CatchInfo;
        try
        {
            c->Read(o);
        }
        catch(...)
        {
            delete c;
            throw;
        }
        catch_info[c->GetCatch()]=c;
    }
    
    if( size!=catch_info.size() )
    {
        printf("Event1Run::Read(): expected %d catches; it was read %zu catches.\n",size,catch_info.size());
        throw "Event1Run::Read(): the data are corrupted.";
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace
