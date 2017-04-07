#include "ObjectsCounter.h"

std::map<std::string,int> ObjectsCounterMaster::counter;
std::ostream *ObjectsCounterMaster::stream=NULL;
bool ObjectsCounterMaster::init=ObjectsCounterMaster::Init();

bool ObjectsCounterMaster::Init(void)
{
    std::atexit(ObjectsCounterMaster::End);
    return true;
}

void ObjectsCounterMaster::End(void)
{
    Print();
}

void ObjectsCounterMaster::SetStream(std::ostream *o)
{
    stream=o;
}

void ObjectsCounterMaster::Print(void)
{
    if( stream==NULL )
        return;
    (*stream) << "ObjectsCounter usage: allocated objects:\n";
    for( std::map<std::string,int>::const_iterator it=counter.begin(); it!=counter.end(); it++ )
        (*stream) << it->first << " " << it->second << "\n";
}
