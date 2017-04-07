#include <cstdio>
#include "DateEquipment.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

void DateEquipment::Print(std::ostream &o,const std::string &prefix) const
{
    o << prefix;
    char s[222];
    sprintf(s,"version=%d type=%u equipmentId=%u length=%u full-length=%u",
               GetVersion(), GetType(), GetEquipmentID(),
               GetLength(),GetFullLength());
    o << s;
    
    if( GetVersion()>0xffff )
        o << " basic_elem_size=" << GetEquipment36().basic_element_size;
    
    o << '\n';
}

////////////////////////////////////////////////////////////////////////////////

} // namespace CS
