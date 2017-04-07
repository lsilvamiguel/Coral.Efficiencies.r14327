#ifndef DigitizerSADCP_h
#define DigitizerSADCP_h

#include <vector>

class MyPulse;

class DigitizerSADCP : public DigitizerSADCBase {
  public:
    virtual ~DigitizerSADCP();
    DigitizerSADCP ( CsCalorimeter* parent, Type type );
    virtual void PrintType ();

  private:
   struct cfd_result_s {
     float time; /*! reconstructed cfd time [clk]*/
     unsigned int ampl : 12; /*! estimated amplitude [ADC channel]*/
   };
   typedef cfd_result_s cfd_result_t;

   struct cfd_option_s {
     unsigned int thr : 12; /*! Amplitude threshold to apply [ADC channel]*/
     unsigned int smin; /*! first sample to analyse */
     unsigned int smax; /*! last sample to analyse */
     unsigned int delay : 12; /*! how long to delay the signal */
     float ampl; /* amplification factor */
   };
   typedef cfd_option_s cfd_option_t;

   MyPulse                        *mypulse_;
   int                             plane_;

   const result_t &PulseFit( const std::vector<uint16> &sample, const unsigned int icell);

};

#endif
