#ifndef CompassSoft_utils
#define CompassSoft_utils

#include <cassert>
#include <cstdio>
#include <iostream>
#include <vector>
#include "config.h"
#include <string>

#include <sstream>

namespace CS {

////////////////////////////////////////////////////////////////////////////////

inline
unsigned long long bits_mask(unsigned int bit_start, unsigned int n)
{
    unsigned long long f=0;
    for( unsigned int k=bit_start; k<bit_start+n; k++ )
        f |= 1ULL<<k;
    return f;
}

////////////////////////////////////////////////////////////////////////////////

inline
unsigned long long bits_get_number(unsigned long long i, unsigned int bit_start, unsigned int n)
{
    return (i & bits_mask(bit_start,n))>>bit_start;
}

////////////////////////////////////////////////////////////////////////////////

inline
std::string bits_string(unsigned long long var, unsigned int bit_start=0, unsigned int length=32)
{
    assert( length <= 32 );
    char s[32+1];
    int i;
    for( i=0; i<(int)length; i++ )
    {
        int k = bit_start+length-i-1;
        s[i] = ((1ULL<<k)&var)?'1':'0';
    }
    s[i]=0;

    return s;
}

////////////////////////////////////////////////////////////////////////////////

inline
void  bits_print(std::ostream &o,unsigned long long var, unsigned int bit_start=0, unsigned int n=32, const char *format="%8x")
{
    for( int k=bit_start+n-1; k>=(int)bit_start; k-- )
        o << (((1ULL<<k)&var)?'1':'0');
    o << "  ";
    char s[200];
    snprintf(s,200,format,(int)((var & bits_mask(bit_start,n))>>bit_start));
    o << s;
}

////////////////////////////////////////////////////////////////////////////////

inline
std::vector<int> bits_from_mask(uint32 x)
{
    std::vector<int> bits;
    for( int i=0; i<32; i++ )
        if( (1UL<<i) & x )
            bits.push_back(i);
    return bits;
}

////////////////////////////////////////////////////////////////////////////////

inline
void swap(int16 *buffer,unsigned n)
{
    const int8 *end = (int8*)(buffer+n);
    for( int8 *p=(int8*)buffer; p!=end; p+=2 )
    {
        int8 tmp=p[0];
        p[0]=p[1];
        p[1]=tmp;
    }
}

////////////////////////////////////////////////////////////////////////////////

inline
void swap(uint16 *buffer,unsigned n)
{
    swap( (int16*)buffer,n );
}


////////////////////////////////////////////////////////////////////////////////

inline
void swap(int32 *buffer,unsigned n)
{
    const int8 *end = (int8*)(buffer+n);
    for( int8 *p=(int8*)buffer; p!=end; p+=4 )
    {
        int8 tmp=p[0];
        p[0]=p[3];
        p[3]=tmp;
        tmp=p[1];
        p[1]=p[2];
        p[2]=tmp;
    }
}

////////////////////////////////////////////////////////////////////////////////

inline
void swap(uint32 *buffer,unsigned n)
{
    swap( (int32*)buffer,n );
}

////////////////////////////////////////////////////////////////////////////////

bool string_match(const std::string &str, const std::string &pattern);

////////////////////////////////////////////////////////////////////////////////

}

#endif // CompassSoft_utils
