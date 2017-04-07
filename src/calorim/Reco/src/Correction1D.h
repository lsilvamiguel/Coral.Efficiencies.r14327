#ifndef _correction1d_h_
#define _correction1d_h_

#include <map>
#include <string>

#include "Exception.h"

namespace Reco {

  /// class to store correction factors which are not constant (0D), but
  /// depending on an external variable, eg. the energy or the time in spill.
  class Correction1D {
  public:

    /// set a single data point
    void SetDataPoint(double variable, double correction) {
      if ( ! fmap.insert(std::make_pair(variable, correction)).second )
        throw Exception( (std::string(__func__)+" duplicate variable: %f").c_str(), variable);
    }

    /// \return number of data points that have been set so far
    int NPoints(void) const { return fmap.size(); }

    /// \return all data points
    const std::map<double, double>& GetMap() const { return fmap; }

    /// \return variable-dependent correction factor, linearly interpolated
    /// between data points
    double Interpolate(double variable) const {
      if ( fmap.size() == 0 )
        return 1.;

      // if we only have one data point, we return it independent of variable
      if ( fmap.size() == 1 )
        return fmap.begin()->second;

      std::map<double,double>::const_iterator hi = fmap.lower_bound(variable);
      // variable is larger than highest data point?
      if ( hi == fmap.end() )
        return fmap.rbegin()->second;

      // variable is equal or smaller than the lowest data point?
      if ( hi == fmap.begin() )
        return fmap.begin()->second;

      std::map<double,double>::const_iterator lo = hi;
      lo--;

      // linear interpolation
      return lo->second
        + (hi->second - lo->second) * (variable - lo->first) / (hi->first - lo->first);
    }

  private:
    std::map<double,double> fmap;

  };  // class Correction1D

}  // namespace Reco

#endif
