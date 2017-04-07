#include "CsCalorimeter.h"
#include "CsDigitizerSADC.h"
#include "CsErrLog.h"

#include "Reco/CellDataRaw.h"

int CsCalorimeter::StoreSADCDigit (const int& icell, const CS::ChipSADC::Digit& digit,
                                   const bool& isLed) {
    bool debug = false;

    // Get cell position
    int32 x_digit = digit.GetX();
    int32 y_digit = digit.GetY();


    // Get samples
    const std::vector<CS::uint16>& sample = digit.GetSamples();

    if ( debug ) {
        std::cout << "CsCalorimeter::StoreSADCDigit " << GetName() << " debug." << std::endl;
        std::cout << "  x_digit = " << x_digit << " y_digit = " << y_digit << std::endl;
        std::cout << "  sample size = " << sample.size() << " integral size = " << digit.GetIntegrals().size() << std::endl;
        std::cout << "  sadc_format_version_ " << sadc_format_version_ << " sadc_decode_version_ " << sadc_decode_version_ << std::endl;
    }

    sadc_samples_.insert(std::pair<int,const std::vector<CS::uint16> * >(icell,&sample) );

//    CsDigitizerSADC *dig = GetDigitizersSADC()[icell];
   if( debug ) std::cout << "CsCalorimeter::StoreSADCDigit " << GetName() <<
                             " NDigitizersSADC " << GetDigitizersSADC().size() <<
                                " NLedDigitizersSADC " << GetLedDigitizersSADC().size() << std::endl;
   if( GetDigitizersSADC().size() != NCells() ||  GetLedDigitizersSADC().size() != NCells() )
   {
     std::cerr <<" DigitizersSADC fatal configuration problem " << std::endl;
     exit(1);
   }

   if( GetDigitizersSADC()[icell]  == NULL ||   GetLedDigitizersSADC()[icell] == NULL )
   {
     std::cerr <<" DigitizersSADC fatal configuration problem for cell " << icell <<
                   " rdig " << GetDigitizersSADC()[icell] <<" ldig " << GetLedDigitizersSADC()[icell] << std::endl;
     exit(1);
   }
   DigitizerSADCBase *dig = isLed ? GetLedDigitizersSADC()[icell] : GetDigitizersSADC()[icell];
   dig -> Fit(sample, icell);

    //  Store decoded info
    double t0=cells_info[Reco::Calorimeter::TIME][OLD][icell].GetMean();
    double signal= dig->GetSignal();
    double time = dig->GetTime();
    double ped =  dig->GetPed();
    cells_info[PED][NEW][icell].Add(ped);


    // TCSphase correction
    double TCS_T0 = 40.; // maximum of TCS phase distr.
    double tcs_cor = GetTCSPhase()-TCS_T0;
    if ( make_tcs_corrections_ ) {
        time += tcs_cor;
    }

    double time_err = sadc_clock_;

    if ( debug ) {
        std::cout << "  Official signal " << signal << " pedestal " << ped << " time " << time << std::endl;
        PrintSADCSettings();
        digit.Print();
        std::cout << dig->GetClassName() <<std::endl;
    }
    bool filterok = true;
    if( sadc_shape_filter_apply_ ) { filterok = !dig->IsNoiseByShape(); }

    // Cuts
    if ( signal > raw_amp_cut_delta_ && time-t0 > raw_time_cut_min_ && time-t0 < raw_time_cut_max_ && filterok ) {
        if ( debug )
            std::cout << "  store this information" << std::endl;

        double energy;
        if (isLed)
          energy = signal * cells_info[CALIB][OLD][icell].GetMean();
        else
          energy = SignalToEnergy(icell, signal, GetTimeInSpill());

        signals_.push_back( Reco::CellDataRaw(icell, energy,
                                              time - cells_info[Reco::Calorimeter::TIME][OLD][icell].GetMean()) );
        signals_.back().SetAmplitude(signal);
        signals_.back().SetTimeErr(time_err);
        return 1;
    }

    return 0;
}

// ****************************************************************************

void CsCalorimeter::DecodeChipSADCDigit(const CS::ChipSADC::Digit& digit, const bool& isLed) {
//    bool debug = true;
    bool debug = false;

    int32 x_digit = digit.GetX();
    int32 y_digit = digit.GetY();

    if ( debug ) {
        std::cout << "CsCalorimeter::DecodeChipSADCDigit: ChipSADC Calorimeter " << getName()
                  << " x= " << x_digit << " y= " << y_digit << std::endl;
        std::cout << " Samples size " << digit.GetSamples().size()
                  << " Integrals size " << digit.GetIntegrals().size() << std::endl;
        digit.Print();
    }

    int icell = GetCellOfColumnRow(x_digit,y_digit);
    if ( icell < 0 ) {
        if( GetName() == "HC01P1__" && x_digit == 27 && y_digit == 0 ) { // Not elegant solution to fix complains on HCAL1 pin diod
          return;
        }
        std::cerr << "ERROR CsCalorimeter::DecodeChipSADCDigit: wrong SADC digit position in " << getName()
                  << " x=" << x_digit << " y=" << y_digit << std::endl;
        digit.Print();
        return;
    }

    if ( digit.GetSamples().size() == 0 ) {
        std::cerr << "ERROR CsCalorimeter::DecodeChipSADCDigit: wrong sample size " << digit.GetSamples().size()
                  << " in " << getName() << " x=" << x_digit << " y=" << y_digit << std::endl;
        digit.Print();
        return;
    }

    StoreSADCDigit(icell, digit, isLed);
}

// ****************************************************************************

void CsCalorimeter::DecodeChipADCDigit(const CS::ChipADC::Digit& digit) {
//    bool debug = true;
    bool debug = false;

    int32 x_digit = digit.GetX();
    int32 y_digit = digit.GetY();
    int32 e_digit = digit.GetAmplitude();

    if ( debug )
        std::cout << "CsCalorimeter::DecodeChipADCDigit: ChipADC Calorimeter " << getName()
                  << " x= " << x_digit << " y= " << y_digit << " E= " << e_digit << std::endl;

    int icell = GetCellOfColumnRow(x_digit,y_digit);
    if ( icell < 0 ) {
        std::cerr << "ERROR CsCalorimeter::DecodeChipADCDigit: wrong ADC digit position in " << getName()
                  << " x=" << x_digit << " y=" << y_digit << std::endl;
        return;
    }

    // In case the overflow bit of an ADC value is set in the DATE data stream,
    // DaqDataDecoding subtracts the ADC value from 4096, so negative e_digits
    // mean the data is in overflow. Usually this should mean, that also the
    // ADC value is close to 4096, however there are some (very rare) cases
    // where this is not the case. This obviously broken data is ignored,
    // other values in the overflow are still used
    if (e_digit < 0) { // overflow bit was set
        if (e_digit > -100) { // ADC value close to 4096, this seems to be a
                              // valid overflow value
            e_digit = 4096;
            CsErrLog::msg(elInfo, __FILE__, __LINE__,
                          "%s: ADC value in overflow.", GetName().c_str());
        } else { // broken data, ignore it
            CsErrLog::msg(elWarning, __FILE__, __LINE__,
                          "%s: Overflow bit in data stream set, but ADC value not in overflow. Information from one cell ignored!", GetName().c_str());
            return;
        }
    }

    // Amplitudes from calorimeter cells
    signals_.push_back( Reco::CellDataRaw(icell) );
    signals_.back().SetAmplitude(e_digit);
    signals_.back().SetEnergy(e_digit * cells_info[CALIB][OLD][icell].GetMean());
}

// ****************************************************************************

void CsCalorimeter::DecodeChipDigit(const CS::Chip::Digit& digit) {
//    bool debug = true;
    bool debug = false;

    if ( debug )
        std::cout << "CsCalorimeter::DecodeChipDigit: Calorimeter " << getName() << " new digit for you." << std::endl;

    // first try ChipSADC
    try {
        const CS::ChipSADC::Digit& d = dynamic_cast<const CS::ChipSADC::Digit&>(digit);

        // so the digit is from a ChipSADC
        DecodeChipSADCDigit(d, false);
    } catch (std::bad_cast& e) {
        // if digit is not from ChipSADC, try ChipADC next
        try {
            const CS::ChipADC::Digit& d = dynamic_cast<const CS::ChipADC::Digit&>(digit);

            // so the digit is from a ChipADC
            DecodeChipADCDigit(d);
        } catch (std::bad_cast& e) {
            // digit neither from ChipADC nor from ChipSADC
            std::cerr << "CsCalorimeter::DecodeChipDigit: Wrong digit type."  << std::endl;
            throw CS::Exception("CsCalorimeter::DecodeChipDigitLEDEvent: Wrong digit type.");
        }
    }
}

// ****************************************************************************

void CsCalorimeter::DecodeChipDigits(const CS::Chip::Digits& digits) {
    if( skip_decoding_ ) return;

    CsDet::DecodeChipDigits(digits);

    MakeCsDigits();
    StoreRawDataStatistic();

    if ( make_sadc_histo_ )
        DecodingTestSADC();
}

// ****************************************************************************

void CsCalorimeter::DecodeChipDigitLEDEvent(const CS::Chip::Digit& digit)
{
//    bool debug = true;
    bool debug = false;

    if( debug )
        std::cout << "CsCalorimeter::DecodeChipDigitLEDEvent: Calorimeter " << getName() << " new digit for you." << std::endl;

    // first try ChipSADC
    try {
        const CS::ChipSADC::Digit& d = dynamic_cast<const CS::ChipSADC::Digit&>(digit);

        // so the digit is from a ChipSADC
        DecodeChipSADCDigit(d, true);
    } catch (std::bad_cast& e) {
        // if digit is not from ChipSADC, try ChipADC next
        try {
            const CS::ChipADC::Digit& d = dynamic_cast<const CS::ChipADC::Digit&>(digit);

            // so the digit is from a ChipADC
            DecodeChipADCDigit(d);
        } catch (std::bad_cast& e) {
            // digit neither from ChipADC nor from ChipSADC
            std::cerr << "CsCalorimeter::DecodeChipDigitLEDEvent: Wrong digit type."  << std::endl;
            throw CS::Exception("CsCalorimeter::DecodeChipDigitLEDEvent: Wrong digit type.");
        }
    }
}

// ****************************************************************************

void CsCalorimeter::DecodeChipDigitsLEDEvent(const CS::Chip::Digits& digits) {
    if( skip_led_decoding_ ) return;

    typedef std::multimap<CS::DetID,CS::Chip::Digit*>::const_iterator m_it; // iterator type
    // get all digits for the detector
    std::pair<m_it,m_it> m_range = digits.equal_range(GetTBName());

    // loop on all found digits
    for( m_it d_it=m_range.first; d_it!=m_range.second; d_it++ )
        DecodeChipDigitLEDEvent(*d_it->second);

    MakeCsDigits();
    StoreRawDataStatistic();
}

// ****************************************************************************
