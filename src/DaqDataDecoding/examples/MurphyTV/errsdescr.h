MurphyTV::ErrZoo::Err MurphyTV::ErrZoo::allerrs[MurphyTV::ErrZoo::NUM_ERRS]=
  {
    {DaqError::EVENT_MISSING_SRCID,"SourceID is missing",
     "There was an event in which a SourceID did not appear  \n"
     "in the data, but is available in the mapping file.\n"
    },
    {DaqError::EVENT_UNKNOWN_SRCID,"unknown SourceID",
     "There was an event in which a SourceID did appear  \n"
     "in the data, but no mapping file is available for it.\n"
    },
    {DaqError::EVENT_SRCID_TOO_BIG,"wrong SourceID",
     "There was an event in which a SourceID did appear  \n"
     "in the data, which is > 1024.\n"
    },   
    {DaqError::SLINKM_BAD_SIZE,"Bad Size from MUX",
     "The amount of data coming from CATCH boards connected to a SLink Multiplexer is \n"
     "not equal to the expected size from the Multiplexer SLink header.\n"
    },
    {DaqError::SLINK_ERRFLAG,"Error-flag set",
     "The error flag for the following data in an SLink header was set.\n"
    },
    {DaqError::SLINK_ERRNR,"Error counter<>0",
     "If the error count field in the SLink header was not zero, this error is reported\n"
     "once regardless of the counter's value.\n"
     "This field is set by the CATCH board in case of errors it is aware of"
    },
    {DaqError::SLINK_WRONG_EVENT_NUMBER,"Wrong Slink event number",
     "The event number in the Slink header does not match the \n" 
     "event number in the date header.\n"
    }, 
    {DaqError::TCS_FIFO,"TCS FIFO full",
     "The \"FIFO full\" bit in the \"TCS error\" field of an SLink header was set.\n"
     "This bit is set by the CATCH boards.\n"
    },
    {DaqError::TCS_SYNC,"TCS synchronization",
     "The \"synchronization error\" bit in the \"TCS error\" field of an SLink header\n"
     "was set.\n"
     "This bit is set by the CATCH boards.\n"
    },
    {DaqError::TCS_UNDEF,"TCS-other-bits",
     "There were some bits set in the \"TCS error\" field which are not defined but \n"
     "usually should be zero"
    },
    {DaqError::TDC_ERR1,"missing header",
     "A data or setup word was received instead of an expected header word.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR2,"Timeout(no data)",
     "Timeout if the CATCH did not receive any data from a frontend.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR3,"event number smaller",
     "A header or trailer word with a smaller event number than expected was received.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR4,"No CMC connected",
     "No CMC was connected.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR5,"Timeout(no trailer)",
     "Timeout if no trailer or to many data words were received.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR6,"FIFO full",
     "FIFO is full, port is off till end of burst.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR7,"event number larger",
     "A header or trailer word with a larger event number than expected was received.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR10,"HOTLink errors",
     "There was a HOTLink transmission error.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_ERR12,"Data skipped",
     "CATCH skipped data because data rate was too high.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_SERR1,"TBO (word)",
     "Trigger buffer overflow.Only for TDC data.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_SERR2,"CMC FIFO full",
     "FIFO buffer overflow on a common mezzanine card.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_SERR3,"Bit 0..23 transmission",
     "HotLink Transmission error in bit 23..0 \n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_WRGEOID,"Wrong GeographicID",
     "A data/header word with an invalid Geographic ID was received.\n"
     "This is reported by the CATCH board via an error word in the data.\n"
    },
    {DaqError::TDC_UNDEF8,"TDC#8",
     "A CATCH generated an error word with code 8 ,which is not defined and should \n"
     "not appear at all.\n"
    },
    {DaqError::TDC_UNDEF9,"TDC#9",
     "A CATCH generated an error word with code 9 ,which is not defined and should \n"
     "not appear at all.\n"
    },
    {DaqError::TDC_UNDEF11,"TDC#11",
     "A CATCH generated an error word with code 11 ,which is not defined and should \n"
     "not appear at all.\n"
    },
    {DaqError::TDC_UNDEF13,"TDC#13",
     "A CATCH generated an error word with code 13 ,which is not defined and should \n"
     "not appear at all.\n"
    },
    {DaqError::TDC_UNDEF14,"TDC#14",
     "A CATCH generated an error word with code 14 ,which is not defined and should \n"
     "not appear at all.\n"
    },
    {DaqError::TDC_UNDEF15,"TDC#15",
     "A CATCH generated an error word with code 15 ,which is not defined and should \n"
     "not appear at all.\n"
    },
    {DaqError::TDC_WRPLL,"Wrong PLL (data)",
     "The PLL locked bits were not all set in a data word coming from aTDC-CMC.\n"
     "If it was a TDC-CMC is determined from the format field of the SLink header.\n"
     "If bit 6 of this field is not set correctly, this error can also be reported\n"
     "as a consequence.\n"
    },
    {DaqError::TDC_WRPLL_H,"Wrong PLL (header)",
     "The PLL locked bits were not all set in a header word coming from aTDC-CMC.\n"
     "If it was a TDC-CMC is determined from the format field of the SLink header.\n"
     "If bit 6 of this field is not set correctly, this error can also be reported\n"
     "as a consequence.\n"
    },
    {DaqError::TDC_WRPORT,"Wrong Port#(data)",
     "In TDC data, the port number in a data word did not match the port number \n"
     "in the preceeding header word.\n"
    },
    {DaqError::TDC_WRPORT_T,"Wrong Port#(trailer)",
     "In TDC data, the port number in a trailer word did not match the port number \n"
     "in the preceeding header word.\n"
    },
    {DaqError::TDC_ERR1U,"no header (unreported)",
     "In TDC data, there was a data word without header and no preceeding error word.\n"
    },
    {DaqError::TDC_WRERR,"faked error 3/7",
     "There was an CATCH error word telling that the event number of a header/trailer\n"
     "was larger/smaller ,but the event number of the header/trailer was equal.\n"
    },
    {DaqError::TDC_ERR7U,"event# larger(unreported)",
     "In TDC data, there was a header/trailer word which had an event number larger\n"
     "than it should have according to the SLink header, but there was no preceeding\n"
     "error word.\n"
    },
    {DaqError::TDC_ERR3U,"event# smaller(unreported)",
     "In TDC data, there was a header/trailer word which had an event number smaller\n"
     "than it should have according to the SLink header, but there was no preceeding\n"
     "error word.\n"
    },
    {DaqError::TDC_TBO,"TBO(bit)",
     "In TDC data, the trigger buffer overflow bit of an header/trailer word was set.\n"
    },
    {DaqError::TDC_WRTIME,"Trigger times differ",
     "In TDC data, the trigger time in a header/trailer differed more than 1 from the\n"
     "one in the preceeding header/trailer (Rollovers are taken into account.).\n"
     "Note that a correct time after an incorrect one is also reported as an error since\n"
     "the program does not try to guess which one is correct. If the correct word after \n"
     "the incorrect one belongs to another port, errors will be reported for both ports."
    },
    {DaqError::TDC_TBO2,"TBO(bit)",
     "In TDC data, the trigger buffer overflow bit of an header word was set.\n"
    },
    {DaqError::TDC_TEMPERATURE,"High Temperature",
     "Front End temperature exceeded set threshold\n"
    },
    {DaqError::TDC_ERR20,"No Header",
     "Maximum size of TDC data received without header word.\n"
    },
    {DaqError::TDC_ERR21,"No Trailer",
     "Maximum size of TDC data received without header word.\n"
    },
    {DaqError::TDC_ERR22,"No Header or Trailer",
     "Maximum size of TDC data received without header word.\n"
    },
    {DaqError::ADC_W_GEOID_EW,"Wrong geoID (extra word)",
     "In ADC data, the optional (and not yet implemented ) extra word after an header \n"
     "contained a geographicID different from the one in the header.\n"
    },
    {DaqError::ADC_W_GEOID_TRAILER,"Wrong geoID(trailer)",
     "In ADC data, a trailer word contained a geographicID different from the one in the\n"
     "header.\n"
    },
    {DaqError::ADC_DATA_WO_HEADER,"Missing header",
     "In ADC data, a data word occured directly after a trailer word.\n"
    },
    {DaqError::ADC_NO_TRAILER,"Missing trailer",
     "In ADC data, there was a header word without a preceeding trailer word.\n"
     "(exception: the first header word).\n"
    },
    {DaqError::ADC_TRAILER_WO_HEADER,"Double Trailer",
     "In ADC data, a trailer word occured directly after another trailer word.\n"
    },
    {DaqError::ADC_WR_EV_HEADER,"Wr.event# (header)",
     "In ADC data, the event number in a header word did not match the SLink event \n"
     "number.\n"
    },
    {DaqError::ADC_WR_EV_TRAILER,"Wr.event# (trailer)",
     "In ADC data, the event number in a trailer word did not match the SLink event\n"
     "number.\n"
    },
    {DaqError::EVENT_WRONG_GEOID_ORDER,"Wrong GeoID order",
     "In ADC (or debug mode F1) data,a header occured with a Geographic ID that differs \n"
     "from the one in the header on the same position in the previous event.\n"
    },
    {DaqError::APV_H_S_E,"APV Header sample error",
     "In APV data, an header word had some of the sampling error bits set.\n"
    },
    {DaqError::APV_ADC_OUT,"Wrong ADC data size",
     "In APV data,the size of ADC data block (as taken from ADC header) is bigger than\n"
     "surrounding data block, or <=0. \n"
    },
    {DaqError::APV_LAST_TRIGGER_TICKS_FAILED,"Wrong APV trigger ticks",
     "APV measured wrong ticks for the last tigger. \n"
     "This should be reported to APV frontend experts.\n"
    },
    {DaqError::APV_H_WR_EV,"APV header event number mismatch",
     "APV header event number mismatches with SLINK ones.\n"
    },
    {DaqError::APV_APV_OUT,"Wrong APV data size",
     "The size of APV data block (as taken from APV header) is bigger than\n"
     "surrounding ADC data block, or <=0. \n"
    },
    {DaqError::APV_MISS_ADC,"Missing ADC data",
     "An ADC card that is present in the mapping file was not found in the \n"
     "data stream.\n"
    },    
    {DaqError::GASSIPLEX_BAD_MSB,"GASSIPLEX bad MSB",
     "GASSIPLEX data word did not start with %1000 \n"
    },
    {DaqError::GASSIPLEX_BIG_BORA,"GASSIPLEX big bora number",
     "GASSIPLEX data BORA number>=24 encountered. \n"
    },
    {DaqError::GASSIPLEX_BAD_CHAN,"GASSIPLEX bad channel#",
     "GASSIPLEX data channel no. %4 !=0 or channel no. >575 \n"
    },
    {DaqError::GASSIPLEX_BAD_GEOID,"GASSIPLEX bad bora/chamber ID",
     "GASSIPLEX data word bora and chamber ID do not match ADC GeoID. \n"
    },
    {DaqError::GASSIPLEX_BAD_CHANID,"GASSIPLEX bad channelID",
     "GASSIPLEX data word contained unused ChannelID. \n"
    },
    {DaqError::EVENT_CORRUPTED0,"The subevent is corrupt",
     "The data of this subevent are corrupt. \n"
     "The data of this subevent (detector) could not be decoded.\n"
    }, 
    {DaqError::EVENT_CORRUPTED1,"Bad equipment length",
     "Equipment length too small or not 4 byte aligned. \n"
    },
    {DaqError::EVENT_CORRUPTED2,"Equipment length !=SLink size",
     "Equipment length != size expected from SLink header. \n"
    },
    {DaqError::EVENT_SUBEVENT_HEADER_DIFF,"Detector header mismatch",
     "The header of the detector and the event do not match.\n"
    },
    {DaqError::EVENT_UNKNOWN,"Unknown event type",
     "Unkown event type. \n"
    },
    {DaqError::EXCEPTION,"Decoding exception",
     "The decoding of the detector data created an exception. \n"
     "in the decoding library. \n"
    },
    {DaqError::EVENT_BAD_CHIP_SIZE,"Event bad chip size",
     "Bad size of data block from a SourceID.\n"
    },    
    {DaqError::TCS_ECC,"Slink ECC checksum error",
     "The ECC checkum error bit of the TCS system is set. This \n"
     "error bit is transferred through the Slink then.\n"
     "Note that in the SLINK data the error bit is inverted, \n"
     "if it is 0 there was a TCS ECC error detected.\n"
    },
    {DaqError::SADC_WRONG_EVENT, "wrong event# SADC header",
     "The event number in the SADC header does not match the event \n"
     "number in the slink. \n"
    },
    {DaqError::SADC_BAD_H_ADC, "corrupted SADC header",
    "The SADC header is corrupted and can not be decoded. \n"
    },
    {DaqError::SADC_H_ADC_BAD_SIZE, "wrong size SADC header",
    "The blocksize in the SADC header is not correct. \n"
    },
    {DaqError::SADC_BAD_DATA, "bad SADC data word",
    "In the SADC data an invalid data word has been detected. \n"
    },
    {DaqError::SADC_BAD_INTEGRAL, "bad SADC integral word",
    "The integral word in the SADC data is not correct. \n"
    },
    {DaqError::SADC_BAD_INTEGRALS, "incon. SADC int. words",
    "The integral words in the SADC data are not consistent. \n"
    },
    {DaqError::SADC_BAD_CHAN_N, "bad SADC channel#",
    "The channel numbers in the SADC data are not in correct order or \n"
    "there are bad entries in the channel number field. \n"
    },
    {DaqError::SADC_NO_DATA_LINE, "SADC no data line",
    "The SADC data do not have a sampling data line for one\n"
    "or several sampling data channels. \n"
    },
    {DaqError::SADC_UNKNOWN_DATA, "SADC unknown data line",
    "The SADC data have a sampling data line for one  \n"
    "or several channels which is corrupted and cannot be decoded. \n"
    },
    {DaqError::SADC_SHORT_SAMPLES, "SADC too few samples",
    "The SADC data for one channel contain too few data \n"
    "data samples or some data samples are lost. \n"
    },
    {DaqError::SADC_BAD_N_SAMPLES, "SADC bad number of samples",
    "The SADC data for one channel contain a different number of\n"
    "data samples than declared in the data format.\n"
    },
    {DaqError::SADC_UNKNOWN_MODE, "SADC bad mode",
    "The SADC is running in an unknown mode\n"
    "to the decoding library.\n"
    },
    {DaqError::SADC_H_ADC_ERROR, "SADC error bit set",
    "Error bit is set in ADC header\n"
    },
    {DaqError::SADC_BAD_DATA_HEADER, "Bad data header",
    "The data header signature is bad,\n"
    "probably the data are corrupted.\n"
    },
    {DaqError::SADC_DATA_OUTSIDE, "Data format error",
    "Probably the data are corrupted.\n"
    },
    {DaqError::SADC_MODULE_MISSING, "ADC header is without data",
    "The module is not connected.\n"
    },
    {DaqError::TT_WRONG_HITS_IN_MASTER_TIME, "wrong hit # in mastertime",
    "The mastertime has a wrong number of hits.\n"
    "This is a serious problem and immediate action is needed! \n"
    "All recorded data might be unusable without mastertime information.\n"
    },
    {DaqError::TT_BIG_SIGMA, "Mastertime too big sigma",
    "The mastertime data show a too big sigma.\n"
    "This is a serious problem and immediate action is needed! \n"
    "All recorded data might be unusable with bad mastertime information.\n"
    },
    {DaqError::TT_MT_SHIFT_CHANGED, "Mastertime shift changed",
    "The mastertime data shows a new offset towards the trigger time.\n"
    "The trigger timing should be checked. The bits of the trigger mask\n"
    "affected are displayed in the value field (in decimal notation).\n"
    },
    {DaqError::WRONG_HITS_IN_TRIGGER_MASK, "Wrong trigger mask",
    "The number of hits in the trigger mask is 0 or too large.\n"
    "This is a serious problem and immediate action is needed! \n"
    "All recorded data might be unusable without trigger mask information.\n"
    },
    {DaqError::UNKNOWN_TYPE, "unknown error type",
    "An error of an unknown type happened in this detector. \n"
    "This should not happen during data taking. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_EVENT1, "GeSiCa: event number mismatch",
    "In GeSiCa data an event number mismatch between Slink \n"
    "and event header was found. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_EVENT2, "GeSiCa: header/trailer mismatch",
    "In GeSiCa data an event number mismatch between \n"
    "header and trailer was found. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_EVENT3, "GeSiCa: header/port mismatch",
    "In GeSiCa data an event number mismatch between \n"
    "header and port data was found. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_EVENT4, "GeSiCa: trailer/port mismatch",
    "In GeSiCa data an event number mismatch between \n"
    "trailer and port data was found. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_PORT_PORT, "GeSiCa: port mismatch",
    "In GeSiCa data an port number mismatch between \n"
    "port header and trailer was found. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_CHIP_PORT, "GeSiCa: chip header mismatch",
    "In GeSiCa data an port number mismatch between \n"
    "chip and port header was found. \n"
    },
    {DaqError::ChipHotGeSiCA_BAD_EVENT, "GeSiCa: bad event structure",
    "In GeSiCa data an corrupted event structure was found. \n"
    },
    {DaqError::ChipHotGeSiCA_BAD_CHIP, "GeSiCa: bad chip structure",
    "In GeSiCa data an corrupted chip data structure was \n"
    "found. \n"
    },
    {DaqError::ChipHotGeSiCA_UNEXPECTED_DT, "GeSiCa: unexpected trailer",
    "In GeSiCa data a trailer was found in the wrong place. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_CHIP, "GeSiCa: wrong chip",
    "In GeSiCa data a wrong chip number was found. \n"
    },
    {DaqError::ChipHotGeSiCA_WRONG_EVENT, "GeSiCa: wrong event number",
    "In GeSiCa data a wrong event number was found. \n"
    },
    {DaqError::ChipHotGeSiCA_PORT_TIMEOUT, "GeSiCa: port timeout",
    "In GeSiCa data a timeout for one or more ports was found. \n"
    },
    {DaqError::ChipHotGeSiCA_EVENT_PORT_LOST, "GeSiCa: full event lost",
    "The GeSiCa has skipped a full event from one port. \n"
    },
    {DaqError::ChipHotGeSiCA_DATA_PORT_LOST, "GeSiCa: port data lost",
    "The GeSiCa has truncated data from one port. \n"
    },
    {DaqError::ChipHotGeSiCA_SKIP_SPILL, "GeSiCa: spill was skipped",
    "The GeSiCa has skipped the data of a spill. \n"
    },
    {DaqError::TCS_PHASE_NO_DIGITS, "missing TCS phase",
    "The TCS phase measurement is missing. This is a serious error. \n"
    "All data from GeSiCas cannot be decoded and are lost \n." 
    },
    {DaqError::TCS_PHASE_MANY_DIGITS, "TCS phase double hits",
    "The TCS phase measurement has multiple hits. \n"
    "This is a serious error. \n"
    "All data from GeSiCas cannot be decoded and are lost \n." 
    },
    {DaqError::APV_LOCAL_TIMETAG_FAILED,
    		"APV TimeTag difference",
    		"APV ERROR: ADC time tag differs across a single GeSiCA. \n"
    },   
  }; 
