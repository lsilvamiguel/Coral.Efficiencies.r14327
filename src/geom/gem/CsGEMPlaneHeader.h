#ifndef CsGEMPlaneHeader_H
#define CsGEMPlaneHeader_H

/*
--------------------------------------------------------------------------

 Declaration of classes to store/analyze events from COMPASS GEM
 detectors

--------------------------------------------------------------------------

class CsGEMPlaneHeader

                 Header object containing event number, run number, date
                 (not used for the time being)

--------------------------------------------------------------------------

v1.0     19/11/2002    by    Bernhard Ketzer

v2.0     27/05/2008    by    Bernhard Ketzer
--------------------------------------------------------------------------
*/

//-----------------------------------------------------------------------------
// CsGEMPlaneHeader class declaration
//-----------------------------------------------------------------------------
class CsGEMPlaneHeader {
    private:
        unsigned int   fDate;          // Date
        unsigned int   fRun;           // Run number
        unsigned int   fEvtNumInRun;   // Event number in run
        unsigned int   fEvtNumInSpill; // Event number in spill

    public:
                        CsGEMPlaneHeader() : fDate(0), fRun(0), fEvtNumInRun(0), fEvtNumInSpill(0) { };

        virtual        ~CsGEMPlaneHeader() { };

        void            Set(unsigned int _d, unsigned int _r, unsigned int _er, unsigned int _es) { fDate=_d; fRun=_r; fEvtNumInRun=_er, fEvtNumInSpill=_es; }

        unsigned int    GetDate()                         const { return fDate;          }
        unsigned int    GetRun()                          const { return fRun;           }
        unsigned int    GetEvtNumInRun()                  const { return fEvtNumInRun;   }
        unsigned int    GetEvtNumInSpill()                const { return fEvtNumInSpill; }
};

#endif

