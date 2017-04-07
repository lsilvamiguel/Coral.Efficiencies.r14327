#ifndef __ToyCalorimeter_h__
#define __ToyCalorimeter_h__

#include <Calorimeter.h>

class ToyCalorimeter : public Reco::Calorimeter {
    public:
                    ToyCalorimeter(double cellThreshold=0.);

        virtual    ~ToyCalorimeter();
};

#endif // __ToyCalorimeter_h__
