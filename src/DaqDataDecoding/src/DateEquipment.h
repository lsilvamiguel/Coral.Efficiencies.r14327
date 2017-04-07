#ifndef CompassSoft_DateEquipment__include
#define CompassSoft_DateEquipment__include

#include <iostream>
#include <string>
#include "config.h"

namespace CS {

////////////////////////////////////////////////////////////////////////////////

/*! \brief DATE Equipment header

    This is the equipment header structure format from DATE library
*/
class DateEquipment
{
  private:

    struct EquipmentOld
    {
      uint16          headerExtLen;           ///< length in bytes of an optional extension
      uint16          type;                   ///< the equipment type identifier
      uint8           reserved;               ///< reserved byte
      uint8           rawByteAlign;           ///< length (in bytes) of the word read from hw
      uint16          equipmentId;            ///< equipment identifier
      uint32          rawDataLen;             ///< length (in bytes) of the data block (no header)
    };
    
    struct Equipment36
    {
      uint32          size;                     ///< Equipment size with the header length.
      uint32          type;
      uint32          equipmentId;
      uint32          type_attributes[3];
      uint32          basic_element_size;
    };

  private:

    /// Base constructor.
                    DateEquipment               (void);

  public:

    /*! \brief Construct the header from given buffer.

        \arg buffer is a pointer to at least 12 bytes memory region.
    */
                    DateEquipment               (uint32 const * const buffer,uint32 ver=0) : buf(buffer),version(ver) {}
  private:

    const EquipmentOld &  GetEquipmentOld(void) const {return *reinterpret_cast<const EquipmentOld*>(buf);}
    const Equipment36  &  GetEquipment36 (void) const {return *reinterpret_cast<const Equipment36 *>(buf);}

  public:

    uint32          GetVersion                  (void) const {return version;}
    
    uint32          SizeOf                      (void) const {return GetVersion()<0xffff ? sizeof(EquipmentOld): sizeof(Equipment36);}

    /// \return true if there is an S-Link information
    bool            IsSLink                     (void) const {return GetType()==64;}

    /// \return S-Link type
    uint16          GetType                     (void) const {return GetVersion()<0xffff ? GetEquipmentOld().type : GetEquipment36().type;}

    /// \return DATE equipment identification number
    uint16          GetEquipmentID              (void) const {return GetVersion()<0xffff ? GetEquipmentOld().equipmentId : GetEquipment36().equipmentId;}

    /// \return header's length in bytes, including extendent length.
    uint32          GetLength                   (void) const {return GetVersion()<0xffff ? sizeof(EquipmentOld)+GetEquipmentOld().headerExtLen : GetEquipment36().size;}

    // /// \return length (in bytes) of the data block (no header)
    //uint32          GetDataLength               (void) const {return GetVersion()<0xffff ? GetEquipmentOld().rawDataLen : throw "DateEquipment::GetDataLength(): no information";}
    
    uint32          GetFullLength               (void) const {return GetVersion()<0xffff ? GetLength()+GetEquipmentOld().rawDataLen : GetEquipment36().size;}


    /// Print DATE header information
    void            Print                       (std::ostream &o=std::cout,const std::string &prefix="") const;

    /// \return Pointer to start of data.
          uint8    *GetBuffer                   (void)       {return (      uint8*) buf;}

    /// \return Pointer to start of data.
    const uint8    *GetBuffer                   (void) const {return (const uint8*) buf;}

  private:
    
    uint32 const * const buf;
    uint32               version;
};

} // namespace CS

#endif // DateEquipment
