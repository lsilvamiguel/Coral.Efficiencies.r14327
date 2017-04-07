/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ShowerProfile.h,v $
   $Date: 2010/11/29 19:16:07 $
   $Revision: 1.2 $
*/

#ifndef ShowerProfile_Reco________include
#define ShowerProfile_Reco________include

#include <map>

#include "CellType.h"

namespace Reco {

class Calorimeter;

class ShowerProfile {

    // ==========================================
    // Constructors and destructor
    //
    // may not be called directly, only from
    // inheriting classes
    // ==========================================

    protected:

        /// Default constructor
                                    ShowerProfile       (const Calorimeter* c, const CellType& ct);

    public:

        /// Destructor
        virtual                    ~ShowerProfile       (void);

    // ==========================================
    // Methods
    // ==========================================

    public:

        const Calorimeter*          GetCalorimeter      () const { return calorimeter_; }
        const CellType&             GetCellType         () const { return cellType_; }

        /* \brief Predict the energy in a cell
         *
         * Predict the energy deposit in a cell, which's centre is distX and distY from the shower centre
         *
         * \param ct the cell type
         * \param distX distance of cell to shower centre along X
         * \param distY distance of cell to shower centre along Y
         */
        virtual double              GetEnergyInCell     (double distX, double distY) const = 0;

        virtual double              GetDeDxInCell       (double distX, double distY) const = 0;
        virtual double              GetDeDyInCell       (double distX, double distY) const = 0;

        /* \brief Print information on this shower profile
         */
        virtual void                Print               (std::ostream& out, const std::string& title="") const = 0;

    // ==========================================
    //  Attributes, data
    // ==========================================

    private:

        const Calorimeter*      calorimeter_;
        const CellType&         cellType_;

};

} // namespace Reco

#endif // ShowerProfile_Reco________include

