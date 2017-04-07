// Class declarations
#include "CsPixelGEMCluster.h"

// GEM headers
#include "CsPixelGEMPlane.h"

// C++ headers
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

//-----------------------------------------------------------------------------
// CsGEMCluster constructor
//-----------------------------------------------------------------------------
CsPixelGEMCluster::CsPixelGEMCluster(const std::list<CsGEMHit*>& _hits, bool _flags,
                                     const CsPixelGEMPlanePar* _par) : CsGEMCluster(_par)
{
  // Save reference to hits in cluster
  fHits = _hits;

  // does the cluster contain flagged channels
  fContainsFlaggedChannels = _flags;

  // calculate cluster amplitude from hits
  CalcAmps();

  // calculate cluster position from hits
  CalcCoG();

  // calculate cluster time
  CalcTime();
}

//-------------------------------------------------------------------------------------
// Calculate cluster pos from a list of hits
//-------------------------------------------------------------------------------------

void CsPixelGEMCluster::CalcCoG() {
  std::vector<float> sum(3,0.);
  double sumt = 0.;
  double wsumx = 0.;
  double wsumy = 0.;
  int maxX(0), minX(31), maxY(0), minY(31);

  // calculate position by centre-of-gravity method
  // Loop over hits
  std::list<CsGEMHit*>::iterator ithit;
  for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {

    // Detector channel, hemisphere, sigma of hit
    double posX = (*ithit)->GetChan()->GetId()->GetPosX();
    double posY = (*ithit)->GetChan()->GetId()->GetPosY();
    double sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();

    // determine maximal and minimal hit positions within one cluster
    if ((*ithit)->GetChan()->GetId()->GetPosX() > maxX) maxX = (*ithit)->GetChan()->GetId()->GetPosX();
    if ((*ithit)->GetChan()->GetId()->GetPosX() < minX) minX = (*ithit)->GetChan()->GetId()->GetPosX();
    if ((*ithit)->GetChan()->GetId()->GetPosY() > maxY) maxY = (*ithit)->GetChan()->GetId()->GetPosY();
    if ((*ithit)->GetChan()->GetId()->GetPosY() < minY) minY = (*ithit)->GetChan()->GetId()->GetPosY();

    // get the fraction of the hit amplitude that should be assigne to this cluster
    double frac = 1.;
    if (fClusterParams->GetShareHits() == 2)
      frac = 1. / (*ithit)->GetNrClusters();

    sumt   += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma);
    wsumx  += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma) * posX;
    wsumy  += frac*((*ithit)->GetAmp()[fClusterParams->GetSample()] - fClusterParams->GetThrHit()*sigma) * posY;
  }

  // Center of gravity
  double centerx = wsumx / sumt;
  double centery = wsumy / sumt;

  // Error on center of gravity
  unsigned int sizeX = maxX-minX+1;
  unsigned int sizeY = maxY-minY+1;
  double dcenterx(1./sqrt(12.)), dcentery(1./sqrt(12.));
  if (sizeX>fClusterParams->GetClusSizeRes().size()) {
    if (fClusterParams->GetClusSizeRes().size()>0)
      dcenterx = fClusterParams->GetClusSizeRes()[fClusterParams->GetClusSizeRes().size()-1];
  } else
    dcenterx = fClusterParams->GetClusSizeRes()[sizeX-1];
  if (sizeY>fClusterParams->GetClusSizeRes().size()) {
    if (fClusterParams->GetClusSizeRes().size()>0)
      dcentery = fClusterParams->GetClusSizeRes()[fClusterParams->GetClusSizeRes().size()-1];
  } else
    dcentery = fClusterParams->GetClusSizeRes()[sizeY-1];

  // Update cluster properties
  fPosX    = centerx;
  fPosY    = centery;
  fPosXErr = dcenterx;
  fPosYErr = dcentery;

  return;
}

//------------------------------------------------------------------------------------
// Correst position of clusters arcording to cluster category
//------------------------------------------------------------------------------------
void CsPixelGEMCluster::CorrectPos() {
  // get the position corrections from the clusterization parameters
  const CsPixelGEMPlanePar* params = dynamic_cast<const CsPixelGEMPlanePar*>(fClusterParams);
  const std::vector<float>& _corrs = params->GetPosCorrs();

  // in total we need 48 parameters, if we do not have all, we do not correct
  if ( _corrs.size() != 48 )
    return;

  // setting up parameters
  const float c010p0  = _corrs[ 0]; // category  10 - constant
  const float c010p1  = _corrs[ 1]; //                linear
  const float c010p2  = _corrs[ 2]; //                quadratic
  const float c010p3  = _corrs[ 3]; //                cubic
  const float c011p0  = _corrs[ 4]; // category  11 - constant
  const float c011p1  = _corrs[ 5]; //                linear
  const float c011p2  = _corrs[ 6]; //                quadratic
  const float c011p3  = _corrs[ 7]; //                cubic
  const float c104Xp0 = _corrs[ 8]; // category 104 - shift in X - constant
  const float c104Xp1 = _corrs[ 9]; //                             linear
  const float c104Xp2 = _corrs[10]; //                             quadratic
  const float c104Xp3 = _corrs[11]; //                             cubic
  const float c104Yp0 = _corrs[12]; // category 104 - shift in Y - constant
  const float c104Yp1 = _corrs[13]; //                             linear
  const float c104Yp2 = _corrs[14]; //                             quadratic
  const float c104Yp3 = _corrs[15]; //                             cubic
  const float c105Xp0 = _corrs[16]; // category 105 - shift in X - constant
  const float c105Xp1 = _corrs[17]; //                             linear
  const float c105Xp2 = _corrs[18]; //                             quadratic
  const float c105Xp3 = _corrs[19]; //                             cubic
  const float c105Yp0 = _corrs[20]; // category 105 - shift in Y - constant
  const float c105Yp1 = _corrs[21]; //                             linear
  const float c105Yp2 = _corrs[22]; //                             quadratic
  const float c105Yp3 = _corrs[23]; //                             cubic
  const float c106Xp0 = _corrs[24]; // category 106 - shift in X - constant
  const float c106Xp1 = _corrs[25]; //                             linear
  const float c106Xp2 = _corrs[26]; //                             quadratic
  const float c106Xp3 = _corrs[27]; //                             cubic
  const float c106Yp0 = _corrs[28]; // category 106 - shift in Y - constant
  const float c106Yp1 = _corrs[29]; //                             linear
  const float c106Yp2 = _corrs[30]; //                             quadratic
  const float c106Yp3 = _corrs[31]; //                             cubic
  const float c107Xp0 = _corrs[32]; // category 107 - shift in X - constant
  const float c107Xp1 = _corrs[33]; //                             linear
  const float c107Xp2 = _corrs[34]; //                             quadratic
  const float c107Xp3 = _corrs[35]; //                             cubic
  const float c107Yp0 = _corrs[36]; // category 107 - shift in Y - constant
  const float c107Yp1 = _corrs[37]; //                             linear
  const float c107Yp2 = _corrs[38]; //                             quadratic
  const float c107Yp3 = _corrs[39]; //                             cubic
  const float c200Xp0 = _corrs[40]; // category 200 - shift in X - constant
  const float c200Xp1 = _corrs[41]; //                             linear
  const float c200Xp2 = _corrs[42]; //                             quadratic
  const float c200Xp3 = _corrs[43]; //                             cubic
  const float c200Yp0 = _corrs[44]; // category 200 - shift in Y - constant
  const float c200Yp1 = _corrs[45]; //                             linear
  const float c200Yp2 = _corrs[46]; //                             quadratic
  const float c200Yp3 = _corrs[47]; //                             cubic

  // check if cluster is already categorized
  int fCategory = Categorize();

  int minX(fHits.front()->GetChan()->GetId()->GetPosX());
  int minY(fHits.front()->GetChan()->GetId()->GetPosY());
  std::list<CsGEMHit*>::iterator ithit;
  for (ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
    minX = std::min(minX, (*ithit)->GetChan()->GetId()->GetPosX());
    minY = std::min(minY, (*ithit)->GetChan()->GetId()->GetPosY());
  }

  // calculate eta
  // eta is the distance from the most bottom left hit to the cluster position
  const float etaX(fPosX - minX);
  const float etaY(fPosY - minY);

  // check for category and apply correction
  // positions are corrected with an res vs eta function calibrated pol3
  if (        fCategory ==  10 ) {
    fPosX -= (c010p0 + c010p1 * etaX + c010p2 * etaX*etaX + c010p3 * etaX*etaX*etaX);
  } else if ( fCategory ==  11 ) {
    fPosY -= (c011p0 + c011p1 * etaY + c011p2 * etaY*etaY + c011p3 * etaY*etaY*etaY);
  } else if ( fCategory == 104 ) {
    fPosX -= (c104Xp0 + c104Xp1 * etaX + c104Xp2 * etaX*etaX + c104Xp3 * etaX*etaX*etaX);
    fPosY -= (c104Yp0 + c104Yp1 * etaY + c104Yp2 * etaY*etaY + c104Yp3 * etaY*etaY*etaY);
  } else if ( fCategory == 105 ) {
    fPosX -= (c105Xp0 + c105Xp1 * etaX + c105Xp2 * etaX*etaX + c105Xp3 * etaX*etaX*etaX);
    fPosY -= (c105Yp0 + c105Yp1 * etaY + c105Yp2 * etaY*etaY + c105Yp3 * etaY*etaY*etaY);
  } else if ( fCategory == 106 ) {
    fPosX -= (c106Xp0 + c106Xp1 * etaX + c106Xp2 * etaX*etaX + c106Xp3 * etaX*etaX*etaX);
    fPosY -= (c106Yp0 + c106Yp1 * etaY + c106Yp2 * etaY*etaY + c106Yp3 * etaY*etaY*etaY);
  } else if ( fCategory == 107 ) {
    fPosX -= (c107Xp0 + c107Xp1 * etaX + c107Xp2 * etaX*etaX + c107Xp3 * etaX*etaX*etaX);
    fPosY -= (c107Yp0 + c107Yp1 * etaY + c107Yp2 * etaY*etaY + c107Yp3 * etaY*etaY*etaY);
  } else if ( fCategory == 200 ) {
    fPosX -= (c200Xp0 + c200Xp1 * etaX + c200Xp2 * etaX*etaX + c200Xp3 * etaX*etaX*etaX);
    fPosY -= (c200Yp0 + c200Yp1 * etaY + c200Yp2 * etaY*etaY + c200Yp3 * etaY*etaY*etaY);
  }
}

//-------------------------------------------------------------------------------------
// Determine cluster category
//-------------------------------------------------------------------------------------
int CsPixelGEMCluster::Categorize() {

  int fCategory(0);

  // restrain cluster size
  if ( fHits.size() <= 3 ) {

    // Generate 3x3 Matrix and
    // Vectors to hold X and Y coordinates and as iterator for Vectors
    int M[3][3];
    // initialise matrix
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
          M[x][y] = 0;

    // Save maximas
    int Xmin (33), Ymin(33);
    for (std::list<CsGEMHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ++ithit) {
      Xmin = std::min((*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
    }

    // loop over hits to fill matrix
    for ( std::list<CsGEMHit*>::iterator ithit=fHits.begin(); ithit != fHits.end(); ++ithit) {
      // consistency checks
      assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
      assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
      assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
      assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

      // fill matrix
      M[(*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
    }

    // binary coding of matrix
    int n = 1;
    for (int y = 0; y < 3; y++) {
      for (int x = 0; x < 3; x++) {
        fCategory += M[x][y] * n;
        n = n*2;
      }
    }


    // re-code binary information into teh final classification

    // clustersize 1
    if ( fCategory == 1 )       fCategory =  1;  // size 1

    // clustersize 2
    else if ( fCategory ==   03 ) fCategory = 10;  // horizontal
    else if ( fCategory ==  011 ) fCategory = 11;  // vertical
    else if ( fCategory ==  012 ) fCategory = 12;  // diagonal descending
    else if ( fCategory ==  021 ) fCategory = 13;  // diagonal ascending

    // clustersize 3
    else if ( fCategory ==   07 ) fCategory = 100; // horizontal
    else if ( fCategory == 0111 ) fCategory = 101; // vertical
    else if ( fCategory == 0124 ) fCategory = 102; // diagonal descending
    else if ( fCategory == 0421 ) fCategory = 103; // diagonal ascending
    else if ( fCategory ==  032 ) fCategory = 104; // -|
    else if ( fCategory ==  023 ) fCategory = 105; // _|
    else if ( fCategory ==  031 ) fCategory = 106; // |-
    else if ( fCategory ==  013 ) fCategory = 107; // |_
    else if ( fCategory == 0112 ) fCategory = 108; // |\ **
    else if ( fCategory == 0221 ) fCategory = 109; // |/
    else if ( fCategory == 0211 ) fCategory = 110; // /|
    else if ( fCategory == 0122 ) fCategory = 111; // \|
    else if ( fCategory ==  034 ) fCategory = 112; // -\ **
    else if ( fCategory ==  043 ) fCategory = 113; // -/
    else if ( fCategory ==  061 ) fCategory = 114; // /-
    else if ( fCategory ==  016 ) fCategory = 115; // \-
    else if ( fCategory == 0121 ) fCategory = 116; // >
    else if ( fCategory == 0212 ) fCategory = 117; // <
    else if ( fCategory ==  052 ) fCategory = 118; // \/
    else if ( fCategory ==  025 ) fCategory = 119; // /\ **
    // tell that an unknown form appeard
    //and push it out of normal coding range
    else {
      fCategory += 200;
      std::cout << "notcoded " << fCategory <<std::endl;
    }
  }

  if ( fHits.size() == 4 ) {
    int Xmin (33), Ymin(33), Xmax(0), Ymax(0);
    for ( std::list<CsGEMHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ++ithit) {
      Xmin = std::min((*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
      Xmax = std::max((*ithit)->GetChan()->GetId()->GetPosX(), Xmax);
      Ymax = std::max((*ithit)->GetChan()->GetId()->GetPosY(), Ymax);
    }
    if ( Xmax - Xmin == 1 && Ymax - Ymin == 1)
      fCategory = 200;
    else if ( Xmax - Xmin <= 2 && Ymax - Ymin <= 2) {
      // Generate 3x3 Matrix and
      // Vectors to hold X and Y coordinates and as iterator for Vectors
      int M[3][3];
      // initialise matrix
      for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++)
          M[x][y] = 0;

      // loop over hits to fill matrix
      for ( std::list<CsGEMHit*>::iterator ithit=fHits.begin(); ithit != fHits.end(); ++ithit) {
        // consistency checks
        assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
        assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
        assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
        assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

        // fill matrix
        M[(*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }

      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 3; y++) {
        for (int x = 0; x < 3; x++) {
          fCategory += M[x][y] * n;
          n = n*2;
        }
      }
      if ( fCategory ==  036 ) fCategory = 201; // -_
      else if ( fCategory ==  063 ) fCategory = 202; // _-
      else if ( fCategory == 0132 ) fCategory = 203; // |\|
      else if ( fCategory == 0231 ) fCategory = 204; // |/|
      else if ( fCategory ==  072 ) fCategory = 205; // T
      else if ( fCategory ==  027 ) fCategory = 206; // inv T
      else if ( fCategory == 0131 ) fCategory = 207; // rot left T or |-
      else if ( fCategory == 0232 ) fCategory = 208; // rot right T or -|
      else if ( fCategory == 0113 ) fCategory = 209; // L
      else if ( fCategory == 0223 ) fCategory = 210; // mir(vert) L
      else if ( fCategory ==  047 ) fCategory = 211; // rot left L
      else if ( fCategory ==  074 ) fCategory = 212; // rot left mir(vert) L
      else if ( fCategory ==  071 ) fCategory = 213; // rot right L
      else if ( fCategory ==  017 ) fCategory = 214; // rot right mir(vert) L
      else if ( fCategory == 0322 ) fCategory = 215; // mir(ho) L
      else if ( fCategory == 0311 ) fCategory = 216; // mir(ho+ver) L
      else fCategory = 217;

    }
    else fCategory = 218;
  }

  if ( fHits.size() == 5 ) {
    int Xmin (33), Ymin(33), Xmax(0), Ymax(0);
    for ( std::list<CsGEMHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ++ithit) {
      Xmin = std::min((*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
      Xmax = std::max((*ithit)->GetChan()->GetId()->GetPosX(), Xmax);
      Ymax = std::max((*ithit)->GetChan()->GetId()->GetPosY(), Ymax);
    }
    if ( Xmax - Xmin <= 2 && Ymax - Ymin <= 2) {
      // Generate 3x3 Matrix and
      // Vectors to hold X and Y coordinates and as iterator for Vectors
      int M[3][3];
      // initialise matrix
      for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++)
          M[x][y] = 0;

      // loop over hits to fill matrix
      for ( std::list<CsGEMHit*>::iterator ithit=fHits.begin(); ithit != fHits.end(); ++ithit) {
        // consistency checks
          assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
          assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
          assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
          assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

        // fill matrix
        M[(*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }

      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 3; y++) {
        for (int x = 0; x < 3; x++) {
          fCategory += M[x][y] * n;
          n = n*2;
        }
      }
      if ( fCategory == 0272 ) fCategory = 300; // +
      else if ( fCategory == 0133 ) fCategory = 301; //
      else if ( fCategory == 0233 ) fCategory = 302; //
      else if ( fCategory ==  073 ) fCategory = 303;
      else if ( fCategory ==  037 ) fCategory = 304;
      else if ( fCategory ==  067 ) fCategory = 305;
      else if ( fCategory ==  076 ) fCategory = 306;
      else if ( fCategory == 0331 ) fCategory = 307;
      else if ( fCategory == 0332 ) fCategory = 308;
      else if ( fCategory == 0334 ) fCategory = 309;
      else if ( fCategory == 0433 ) fCategory = 310;
      else if ( fCategory == 0166 ) fCategory = 311;
      else if ( fCategory == 0661 ) fCategory = 312;
      else fCategory = 313;
    } else
      fCategory = 314;
  }

  if ( fHits.size() == 6 ) {
    int Xmin (33), Ymin(33), Xmax(0), Ymax(0);
    for ( std::list<CsGEMHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ++ithit) {
      Xmin = std::min((*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
      Xmax = std::max((*ithit)->GetChan()->GetId()->GetPosX(), Xmax);
      Ymax = std::max((*ithit)->GetChan()->GetId()->GetPosY(), Ymax);
    }
    if ( Xmax - Xmin <= 2 && Ymax - Ymin <= 2) {
      // Generate 3x3 Matrix and
      // Vectors to hold X and Y coordinates and as iterator for Vectors
      int M[3][3];
      // initialise matrix
      for (int y = 0; y < 3; y++)
        for (int x = 0; x < 3; x++)
          M[x][y] = 0;

      // loop over hits to fill matrix
      for ( std::list<CsGEMHit*>::iterator ithit=fHits.begin(); ithit != fHits.end(); ++ithit) {
        // consistency checks
        assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
        assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
        assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
        assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

        // fill matrix
        M[(*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }

      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 3; y++) {
        for (int x = 0; x < 3; x++) {
          fCategory += M[x][y] * n;
          n = n*2;
        }
      }
      if ( fCategory ==  077 ) fCategory = 400;
      else if ( fCategory == 0333 ) fCategory = 401;
      else if ( fCategory == 0633 ) fCategory = 402;
      else if ( fCategory == 0363 ) fCategory = 403;
      else if ( fCategory == 0336 ) fCategory = 404;
      else if ( fCategory == 0636 ) fCategory = 405;
      else if ( fCategory == 0176 ) fCategory = 406;
      else if ( fCategory == 0671 ) fCategory = 407;
      else if ( fCategory == 0275 ) fCategory = 408;
      else if ( fCategory == 0572 ) fCategory = 409;
      else if ( fCategory == 0276 ) fCategory = 410;
      else if ( fCategory == 0672 ) fCategory = 411;
      else if ( fCategory == 0372 ) fCategory = 412;
      else if ( fCategory == 0273 ) fCategory = 413;
      else fCategory = 414;
    } else {
      // Generate 3x3 Matrix and
      // Vectors to hold X and Y coordinates and as iterator for Vectors
      int M[9][9];
      // initialise matrix
      for (int y = 0; y < 9; y++)
        for (int x = 0; x < 9; x++)
          M[x][y] = 0;

      // loop over hits to fill matrix
      for ( std::list<CsGEMHit*>::iterator ithit=fHits.begin(); ithit != fHits.end(); ++ithit) {
        // consistency checks
        assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 8);
        assert((*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
        assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 8);
        assert((*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

        // fill matrix
        M[(*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }

      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 9; y++) {
        if (false) std::cout << std::endl;
        for (int x = 0; x < 9; x++) {
          fCategory += M[x][y] * n;
          n = n*2;
          // output of matrix
          if ( false ){
            if ( M[x][y] )
              std::cout << M[x][y];
            else
              std::cout << ".\t";
          }
        }
      }
      if (false)
        std::cout << "\n\tBin 415:" << fCategory << std::endl;
      fCategory = 415;
    }
  }

  if ( fHits.size() > 6)
    fCategory = 1000;

  return fCategory;
}

bool CsPixelGEMCluster::AddHit(CsGEMHit* _hit) {
    bool insert = false; // should the current hit be added to the cluster
    bool changeflagged = false; // contains the expanded cluster flagged channels

    // read the clusterization configuration mask
    bool  clpossizesplit  = (fClusterParams->GetClusConfigMask() >> 4) & 1; // cut on the distance of hit to current CoG
    bool  clpossizemethod = (fClusterParams->GetClusConfigMask() >> 5) & 1; // cut radial (true) or else (false)

    // Add if list of hits in cluster empty
    if (fHits.empty())
        insert=true;
    else { // Check if new hit should be added to existing list of hits
        // Get list of active neighbours for new hit
        const std::list<CsGEMChan*>& neighbours = _hit->GetChan()->GetActiveNeighbours();
        for (std::list<CsGEMChan*>::const_iterator itchan=neighbours.begin(); itchan!=neighbours.end(); ++itchan) {
            // Check if any hit in cluster is a neighbour of the new hit
            for (std::list<CsGEMHit*>::const_iterator ithit=fHits.begin(); ithit!=fHits.end(); ++ithit) {
                if ((*itchan) == (*ithit)->GetChan()) {
                    insert = true;

//                    std::cout<<"CsGEMCluster::AddHit() : hit "
//                             <<_hit->GetChan()->GetId()->GetDetectorChannel()
//                             <<" neighbour of "
//                             <<(*ithit)->GetChan()->GetId()->GetDetectorChannel()
//                             <<std::endl;

                    // Check if there is a gap between neighbours
                    int distx = abs( (*itchan)->GetId()->GetPosX() - (*ithit)->GetChan()->GetId()->GetPosX());
                    int disty = abs( (*itchan)->GetId()->GetPosY() - (*ithit)->GetChan()->GetId()->GetPosY());
                    if ( distx>1 || disty>1 )
                        changeflagged = true;
                }
            }
        }

        // hit and cluster are neighbours, but does the hit fulfill all conditions
        if (insert) {
            // calculate current cluster center and distance to current hit
            CalcAmps();
            CalcCoG();
            float distx = fabs( _hit->GetChan()->GetId()->GetPosX() - fPosX );
            float disty = fabs( _hit->GetChan()->GetId()->GetPosY() - fPosY );
            float dists = distx*distx + disty*disty;

            // cluster position size splitting
            if (clpossizesplit) {
                if ( (clpossizemethod && dists>2.) || (!clpossizemethod && (distx>1. || disty>1.)) )
                    insert = false;
            }
        }
    }

    // Add hit to cluster
    if (insert) {
        fHits.push_back(_hit);
        if (changeflagged)
            fContainsFlaggedChannels = true;
    }

    return insert;
}

//-----------------------------------------------------------------------------
// Comparison of 2 clusters
//-----------------------------------------------------------------------------
bool CsPixelGEMCluster::operator< (const CsPixelGEMCluster &cluster) const {
    if (fPosX == cluster.GetPositionX())
        return (fPosY < cluster.GetPositionY());

    return (fPosX < cluster.GetPositionX());
}
