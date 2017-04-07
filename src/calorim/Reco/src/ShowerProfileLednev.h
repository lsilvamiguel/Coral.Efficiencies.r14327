/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ShowerProfileLednev.h,v $
   $Date: 2011/01/30 15:30:40 $
   $Revision: 1.3 $
*/

#ifndef ShowerProfileLednev_Reco________include
#define ShowerProfileLednev_Reco________include

#include "ShowerProfile.h"

namespace Reco {

class ShowerProfileLednev : public ShowerProfile {

    // ==========================================
    // Constructors and destructor
    // ==========================================

    public:

        /// Default constructor
                        ShowerProfileLednev                 (const Calorimeter* c, const CellType& ct, std::vector<double> params);

        /// Destructor
        virtual        ~ShowerProfileLednev                 (void) {};

    // ==========================================
    // Methods
    // ==========================================

    public:

        virtual double  GetEnergyInCell                     (double distX, double distY) const;
        virtual double  GetDeDxInCell                       (double distX, double distY) const;
        virtual double  GetDeDyInCell                       (double distX, double distY) const;

        virtual void    Print                               (std::ostream& out, const std::string& title="") const;

    private:
        double          F                                   (double distX, double distY) const;
        double          dFdX                                (double distX, double distY) const;
        double          dFdY                                (double distX, double distY) const;

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        unsigned int    countPar;   // number of parameter sets
        double          paramsA[3];
        double          paramsB[3];

};

} // namespace Reco

#endif // ShowerProfileLednev_Reco________include

