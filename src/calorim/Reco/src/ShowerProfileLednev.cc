/*
   $Source: /afs/cern.ch/user/g/gobbo/w0/coral/src/calorim/Reco/src/ShowerProfileLednev.cc,v $
   $Date: 2010/11/29 19:16:07 $
   $Revision: 1.2 $
*/

#include <cmath>

#include "ShowerProfileLednev.h"

#include "Calorimeter.h"

namespace Reco {

////////////////////////////////////////////////////////////////////////////////

ShowerProfileLednev::ShowerProfileLednev(const Calorimeter* c, const CellType& ct, std::vector<double> params) :
    ShowerProfile(c, ct), countPar(0) {
    if (params.size()!=4 && params.size()!=6) {
        std::cerr << "ShowerProfileLednev::ShowerProfileLednev: exactly four or six parameters expected, got " << params.size()
                  << " for calorimeter " << GetCalorimeter()->GetName()
                  << ", cell type " << CellType::GetHwTypeString(GetCellType().GetHwType()) << ". Fatal error!" << std::endl;
        exit(1);

        // TODO: should be fatal exception, but is caught somewhere
        throw Exception("ShowerProfileLednev::ShowerProfileLednev: exactly four or six parameters expected, got %zu for calorimeter %s, cell type %s. Fatal error!",
                        params.size(), GetCalorimeter()->GetName().c_str(), CellType::GetHwTypeString(GetCellType().GetHwType()).c_str());
    }

    if (params.size()==4) {
        countPar   = 2;
        paramsA[0] = params[0];
        paramsA[1] = params[1];
        paramsB[0] = params[2];
        paramsB[1] = params[3];
    } else if (params.size()==6) {
        countPar   = 3;
        paramsA[0] = params[0];
        paramsA[1] = params[1];
        paramsA[2] = params[2];
        paramsB[0] = params[3];
        paramsB[1] = params[4];
        paramsB[2] = params[5];
    }
}

////////////////////////////////////////////////////////////////////////////////

void ShowerProfileLednev::Print(std::ostream& out, const std::string& title) const {
    out << title;

    out << "Lednev Shower profile with " << countPar << " pairs of parameters:";

    for (unsigned int i=0; i<countPar; i++)
        out << " (" << paramsA[i] << "," << paramsB[i] << ")";

    out << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

double ShowerProfileLednev::GetEnergyInCell(double distX, double distY) const {
    double sum(0.);

    sum+=F(distX+GetCellType().GetStepX()/2., distY+GetCellType().GetStepY()/2.);
    sum-=F(distX+GetCellType().GetStepX()/2., distY-GetCellType().GetStepY()/2.);
    sum-=F(distX-GetCellType().GetStepX()/2., distY+GetCellType().GetStepY()/2.);
    sum+=F(distX-GetCellType().GetStepX()/2., distY-GetCellType().GetStepY()/2.);

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double ShowerProfileLednev::GetDeDxInCell(double distX, double distY) const {
    double sum(0.);

    sum+=dFdX(distX+GetCellType().GetStepX()/2., distY+GetCellType().GetStepY()/2.);
    sum-=dFdX(distX+GetCellType().GetStepX()/2., distY-GetCellType().GetStepY()/2.);
    sum-=dFdX(distX-GetCellType().GetStepX()/2., distY+GetCellType().GetStepY()/2.);
    sum+=dFdX(distX-GetCellType().GetStepX()/2., distY-GetCellType().GetStepY()/2.);

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double ShowerProfileLednev::GetDeDyInCell(double distX, double distY) const {
    double sum(0.);

    sum+=dFdY(distX+GetCellType().GetStepX()/2., distY+GetCellType().GetStepY()/2.);
    sum-=dFdY(distX+GetCellType().GetStepX()/2., distY-GetCellType().GetStepY()/2.);
    sum-=dFdY(distX-GetCellType().GetStepX()/2., distY+GetCellType().GetStepY()/2.);
    sum+=dFdY(distX-GetCellType().GetStepX()/2., distY-GetCellType().GetStepY()/2.);

    return sum;
}

////////////////////////////////////////////////////////////////////////////////

double ShowerProfileLednev::F(double distX, double distY) const {
    double result(0.);

    for (unsigned int i=0; i<countPar; i++)
        result += paramsA[i]
                  * (atan(distX/paramsB[i])
                     + atan(distY/paramsB[i])
                     + atan2(distX*distY, paramsB[i]*sqrt(paramsB[i]*paramsB[i]+distX*distX+distY*distY))
                    );

    result /= (2.*M_PI);

    result += 0.25;

    return result;
}

////////////////////////////////////////////////////////////////////////////////

double ShowerProfileLednev::dFdX(double distX, double distY) const {
    double result(0.);

    for (unsigned int i=0; i<countPar; i++)
        result += paramsA[i] * ( paramsB[i]/(distX*distX+paramsB[i]*paramsB[i]) *
                                                  (1. + distY/sqrt(paramsB[i]*paramsB[i]+distX*distX+distY*distY))
                                      );

    result /= (2.*M_PI);

    return result;
}

////////////////////////////////////////////////////////////////////////////////

double ShowerProfileLednev::dFdY(double distX, double distY) const {
    double result(0.);

    for (unsigned int i=0; i<countPar; i++)
        result += paramsA[i] * ( paramsB[i]/(distY*distY+paramsB[i]*paramsB[i]) *
                                                  (1. + distX/sqrt(paramsB[i]*paramsB[i]+distX*distX+distY*distY))
                                      );

    result /= (2.*M_PI);

    return result;
}

////////////////////////////////////////////////////////////////////////////////


} // namespace Reco

