#ifndef OPTCORR_H
#define OPTCORR_H


class CsRCOptCorr {

private:

//080922 - Added by Paolo
  static CsRCOptCorr* instance_;

  double direction[2][7];
  double position[16][2][7][7];

  bool GetGrid( int channel, double dx, double dy, int &idx, int &idy );
  bool DoInterpolation( int channel, double direction_x, double direction_y,
                        int gridx, int gridy,
                        double& position_x, double &position_y );

//080922 - Added by Paolo
  bool ExtendGrid();
  bool ExtrapG( int, double*, double*, double, double& );
  void checkNo9999();

public:

//080922 - Added by Paolo
  static CsRCOptCorr* Instance();

  CsRCOptCorr();
  bool GetPhotonPosition( int channel,
                          double direction_x, double direction_y,
                          double& position_x, double& position_y );

};

#endif
