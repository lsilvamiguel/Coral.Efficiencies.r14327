/* src/config.h.  Generated automatically by configure.  */
#ifndef cool_config___include
#define cool_config___include



#ifndef USE_DATABASE
#define USE_DATABASE 1
#endif

/* #undef WORDS_BIGENDIAN */

#define SIZEOF_CHAR 1
#define SIZEOF_SHORT 2
#define SIZEOF_INT 4
#define SIZEOF_LONG 4
#define SIZEOF_LONG_LONG 8
/* #undef __CHAR_UNSIGNED__ */

#if  SIZEOF_CHAR==1
  #ifndef __CHAR_UNSIGNED__
    typedef          char         int8;
    typedef unsigned char        uint8;
  #else
    typedef   signed char         int8;
    typedef          char        uint8;
  #endif
#else
  #error Can not define types int8,uint8
#endif

#if  SIZEOF_SHORT==2
  typedef          short        int16;
  typedef unsigned short       uint16;
#elif SIZEOF_INT==2
  typedef          int          int16;
  typedef unsigned int         uint16;
#else
  #error Can not define types int16,uint16
#endif

#if  SIZEOF_SHORT==4
  typedef          short        int32;
  typedef unsigned short       uint32;
  #if USE_DATE_LIB
    #define long32              short
  #endif
#elif SIZEOF_INT==4
  typedef          int          int32;
  typedef unsigned int         uint32;
  #if USE_DATE_LIB
    #define long32                int
  #endif
#elif SIZEOF_LONG==4
  typedef          long         int32;
  typedef unsigned long        uint32;
  #if USE_DATE_LIB
    #define long32               long
  #endif
#else
  #error Can not define types int32,uint32
#endif

#if  SIZEOF_SHORT==8
  typedef          short        int64;
  typedef unsigned short       uint64;
#elif SIZEOF_INT==8
  typedef          int          int64;
  typedef unsigned int         uint64;
#elif SIZEOF_LONG==8
  typedef          long         int64;
  typedef unsigned long        uint64;
#elif SIZEOF_LONG_LONG==8
  typedef          long long    int64;
  typedef unsigned long long   uint64;
#else
  #warning Can not define types int64,uint64
#endif

#endif // cool_config___include
