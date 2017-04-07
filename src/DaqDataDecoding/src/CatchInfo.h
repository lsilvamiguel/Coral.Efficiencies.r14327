#ifndef CatchInfo___include
#define CatchInfo___include

#include <iostream>
#include "config.h"
#include "SLink.h"

namespace CS {

class CatchInfo
{
  public:

    enum Sender {APV=0,ScalerCMC=3,ADC=4,RICH=5,ScalerHOTLink=7};

    class Header
    {
      public:

        class PortLine
        {
          public:
                        PortLine                (uint32 n=0) {reinterpret_cast<uint32&>(*this)=n;}
            uint32      fe_number:16;
            uint32      geoID:10;
            uint32      one_one:2;
            uint32      port:4;
        };

      public:
      
        void            Print                   (const char *prefix="") const;
      
        uint32          ram:1;
        uint32          x:1;
        uint32          zero:10;
        uint32          CMC:4;
        uint32          catch_number:16;

        uint32          TCS_FPGA_rev:8;
        uint32          SLink_FPGA_rev:8;
        uint32          formatter_rev:8;
        uint32          merger_rev:8;

        uint32          CMC2_number:16;
        uint32          CMC1_number:16;
        uint32          CMC4_number:16;
        uint32          CMC3_number:16;
        
        PortLine        port_lines[16];
    };

  public:
  
    /// Destroy the CatchInfo object
                       ~CatchInfo               (void) {delete []buffer;}

    /// Construct an empty CatchInfo object
                        CatchInfo               (void) {buffer=NULL;}

    /// Construct CatchInfo object from the existing buffer
                        CatchInfo               (const void *buf);

  private:

    /// You can not use Copy constructor
                        CatchInfo               (const CatchInfo &c);// : buffer(NULL) {*this=c;}

  private:

    /// You can not use assignment operator
    CatchInfo &         operator =              (const CatchInfo &c);
    
  public:

    /// Print the CatchInfo
    void                Print                   (const char *prefix="") const;

    /// \return S-Link object.
    const SLink &       GetSLink                (void) const {if(buffer==NULL) throw "CatchInfo::GetSLink(): CatchInfo is empty.";
                                                              return *reinterpret_cast<const SLink *>(buffer);}

    Sender              GetSender               (void) const {return Sender((GetSLink().GetFormat()>>4)&7);}

    /// \return First event of run header.
    const Header &      GetHeader               (void) const {if(buffer==NULL) throw "CatchInfo::GetHeader(): CatchInfo is empty.";
                                                              return *reinterpret_cast<const Header *>(buffer+sizeof(SLink));}

    /// \return Pointer to the raw data structure (starting from the S-Link header).
    const uint8 *       GetBuffer               (void) const {return buffer;}

    /// \return The size of buffer in bytes
    unsigned int        GetBufferSize           (void) const;

    /// \return Pointer to the raw data structure (starting from the S-Link header).
    const uint8*        GetData                 (void) const {return buffer+sizeof(SLink)+sizeof(Header);}

    unsigned int        GetDataSize             (void) const;

    /// Write to output stream
    void                Write                   (std::ostream &o) const;

    /// Read from input stream
    void                Read                    (std::istream &o);
    
    unsigned int        GetCatch                (void) const {return GetHeader().catch_number;}

  private:

    /// Pointer to the buffer with the raw data
    uint8 *             buffer;
};

} // namespace CS

#endif // CatchInfo___include
