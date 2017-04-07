#ifndef CS_TIMING_INFO_FORMAT_H_
#define CS_TIMING_INFO_FORMAT_H_

struct ec02_tcorr_data_t {
  unsigned int bin[5];
};
struct ec02_tcorr_header_t{
  unsigned int version;
  unsigned int run;
  unsigned int nu_1;
  unsigned int nu_2;
  time_t time;
};
struct ec02_tcorr_header2_t{
  unsigned int checksum[200];
};

#endif // CS_TIMING_INFO_FORMAT_H_

