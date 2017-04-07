#ifndef __CS_STORE_MISC_H__
#define __CS_STORE_MISC_H__

#include <string>
#include <vector>
#include <iostream>

typedef unsigned char       uint8;
typedef unsigned int        uint32;
typedef unsigned short      uint16;
typedef short int           int16;
typedef float               float32;
typedef double              float64;

std::string itostr(int i);
std::string utostr(unsigned u);

void AddValueToVector(std::vector<uint8>& buffer, const void* value, size_t size);

class VectorBuffer : public std::streambuf
{
public:
  VectorBuffer(std::vector<uint8>& buffer, unsigned size = 0)
    {
      _M_in_beg = (char*)&buffer[0];
      _M_in_cur = _M_in_beg;
#if (__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ > 3)
      _M_in_end = _M_in_beg + (int)((size == 0) ? buffer.size() : size);
#elif __GNUC__ == 3
      _M_buf = (char*)&buffer[0];
      _M_buf_size = (size == 0) ? buffer.size() : size;
      _M_buf_unified = true;
      _M_mode = std::ios_base::in;
      _M_in_end = _M_in_beg+_M_buf_size;
#endif
    }
  ~VectorBuffer() {};
};

#endif
