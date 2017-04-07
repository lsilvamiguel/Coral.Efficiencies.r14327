#include <map>
#include <cstring>
#include <cstdio>
#include "CatchInfo.h"
 
namespace CS {

using namespace std;

////////////////////////////////////////////////////////////////////////////////

CatchInfo::CatchInfo(const void *buf) :
    buffer(NULL)
{
    const SLink &slink=*reinterpret_cast<const SLink *>(buf);
    if( slink.GetEventSize()<sizeof(Header)/4 )
        throw "CatchInfo::CatchInfo(): the size of SLink is too small.";
    unsigned int size=(slink.GetEventSize()+sizeof(SLink))*4;
    buffer = new uint8[size];
    memcpy(buffer,buf,size);
}

////////////////////////////////////////////////////////////////////////////////

unsigned int CatchInfo::GetBufferSize(void) const
{
    if( buffer==NULL )
        throw "CatchInfo::GetBufferSize(): CatchInfo is empty.";

    return GetSLink().GetEventSize()*4+sizeof(SLink);
}

////////////////////////////////////////////////////////////////////////////////

unsigned int CatchInfo::GetDataSize(void) const
{
    if( buffer==NULL )
        throw "CatchInfo::GetDataSize(): CatchInfo is empty.";

    int n=GetBufferSize()-sizeof(SLink)-sizeof(Header);

    if( n<0 )
        throw "CatchInfo::GetDataSize(): negative data size!";
    
    return n;
}

////////////////////////////////////////////////////////////////////////////////

void CatchInfo::Write(std::ostream &o) const
{
    if( buffer==NULL )
        throw "CatchInfo::Write(): CatchInfo is empty.";
    if( !o.write(reinterpret_cast<char*>(buffer),GetBufferSize() ) )
        throw "CatchInfo::Write(): failed to write to a stream.";
}

////////////////////////////////////////////////////////////////////////////////

void CatchInfo::Read(std::istream &o)
{
    delete [] buffer;
    
    // First we read SLink header
    SLink slink;
    if( !o.read(reinterpret_cast<char*>(&slink),sizeof(slink)) )
        throw "CatchInfo::Read(): failed to read a S-Link from a stream.";

    // Check that the data size is reasonable
    if( slink.GetEventSize()*4<sizeof(Header) )
    {
        printf("CatchInfo::Read(): Slink size is %d words. It is too small. I need %zu words minimum\n",
                 slink.GetEventSize(),sizeof(Header)/4);
        throw "CatchInfo::Read(): too little data";
    }

    // Now we know how much memory we need
    buffer=reinterpret_cast<uint8*>(&slink);  // Set buffer pointer so we can call later GetBufferSize()
    unsigned size=GetBufferSize();
    buffer = new uint8 [size];
    memcpy(buffer,&slink,sizeof(slink));
    
    // And finally read the rest of the data
    if( !o.read(reinterpret_cast<char*>(buffer)+sizeof(slink),size-sizeof(slink)) )
        throw "CatchInfo::Read(): failed to read Header+Data from a stream.";
}

////////////////////////////////////////////////////////////////////////////////

void CatchInfo::Print(const char *prefix) const
{
    static map<int,string> sender_name;
    if( sender_name.size()==0 )
    {
        sender_name[APV]            = "GeSiCA-APV";
        sender_name[ScalerCMC]      = "Scaler-CMC";
        sender_name[ADC]            = "FI-ADC HOTLink";
        sender_name[RICH]           = "RICH HOTFiber";
        sender_name[ScalerHOTLink]  = "Scaler HOTLink";
    }

    if( buffer==NULL )
        printf("%sCatchInfo is empty.\n",prefix);
    printf("%sSender is %d %s\n",prefix,GetSender(),sender_name[GetSender()].c_str());
    GetSLink().Print(cout,prefix);
    GetHeader().Print(prefix);
}

////////////////////////////////////////////////////////////////////////////////

void CatchInfo::Header::Print(const char *prefix) const
{
    const char *on_off[]={"OFF","ON"};
    printf("%scatch S/N=%u(0x%X)\n"
           "%sram=%u; TCS FPGA %u; S-Link FPGA %u; Formetter %u; Merger %u\n"
           "%scmc1=%u(0x%X) %s; cmc2=%u(0x%X) %s; cmc3=%u(0x%X) %s; cmc4=%u(0x%X) %s\n",
            prefix, catch_number, catch_number,
            prefix, ram, TCS_FPGA_rev, SLink_FPGA_rev, formatter_rev, merger_rev,
            prefix, CMC1_number, CMC1_number, on_off[bool(CMC&1)],
                    CMC2_number, CMC2_number, on_off[bool(CMC&2)],
                    CMC3_number, CMC3_number, on_off[bool(CMC&4)],
                    CMC4_number, CMC4_number, on_off[bool(CMC&8)]);

    for( unsigned int i=0; i<16; i++ )
    {
        const PortLine &p=port_lines[i];
        printf("%s  port=%2u  geoID=%3u(0x%X)  FE serial number is %u(0x%X)\n",
                prefix, int(p.port), int(p.geoID), int(p.geoID), int(p.fe_number), int(p.fe_number));
    }
}

////////////////////////////////////////////////////////////////////////////////

} // namespace
