#include <cmath>
#include <cstdio>

#include "Trigger.h"
#include "Exception.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

Trigger::Trigger(const ObjectXML &xml) :
    bit(unsigned(-1)),
    precision(0),
    mt_correction(0),
    mt_offset(0),
    verbose(false)
{
    if( xml.GetName()!="Trigger" && xml.GetName()!="TriggerSettings" )
        throw Exception("Trigger::Trigger(): XML-object is not a trigger!");

    xml.GetAttribute("name",name);
    xml.GetAttribute("bit",bit);
    xml.GetAttribute("precision",precision);
    xml.GetAttribute("MT_offset",mt_offset);
    xml.GetAttribute("MT_correction",mt_correction);
    xml.GetAttribute("verbose",verbose);
}

////////////////////////////////////////////////////////////////////////////////

void Trigger::Print(const std::string &prefix) const
{
    printf("%sname=%20s   bit=%2d   precision=%4.1f   mt_correction=%7g   mt_offset=%7g\n",
             prefix.c_str(), GetName().c_str(),GetBit(),precision,mt_correction,mt_offset);
}

////////////////////////////////////////////////////////////////////////////////

} // namespace
