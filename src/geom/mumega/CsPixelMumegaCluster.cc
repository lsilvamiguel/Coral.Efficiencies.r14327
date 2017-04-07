// include C++ headers
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdlib>   // for abs()

#include "CsPixelMumegaCluster.h"
#include "CsPixelMumegaPlane.h"


CsPixelMumegaCluster::CsPixelMumegaCluster() {
  fPosX = 0.;
  fPosY = 0.;
  fPosXErr = 0.;
  fPosYErr = 0.;
  fNoise = 0.;
  fXTalk = 0.;
  fHasTime = false;
  fContainsFlaggedChannels = kFALSE;
}


CsPixelMumegaCluster::CsPixelMumegaCluster(const std::list<CsMumegaHit*> &_hits, Bool_t _flags, const CsPixelMumegaPlanePar &_par) {
  // Save reference to hits in cluster
  fHits = _hits;
  // does the cluster contain flagged channels ?
  fContainsFlaggedChannels = _flags;

  // get clusterization parameters
  Int_t _sample = _par.GetSample();
  Float_t _thrhit = _par.GetThrHit();
  Int_t _share = (_par.GetClusConfigMask() >> 6) & 3;
  std::vector<Float_t> _time_cal = _par.GetTimeCal();

  // Calculate cluster amplitude from hits
  CalcAmps(_sample, _thrhit, _share);

  // Calculate cluster position from hits
  CalcCoG(_sample, _thrhit, _share);

  // Calculate cluster time
  //CalcTime(_time_cal);
}


void CsPixelMumegaCluster::CalcCoG(int _sample, float _thr, int _share) {
  std::vector<float> sum(3,0.);
  double sumt = 0.;
  double wsumx = 0.;
  double wsumy = 0.;
  double dcenter = 0.;

  // Calculate position by center-of-gravity method
  // Loop over hits
  std::list<CsMumegaHit*>::iterator ithit;
  for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {

    // Detector channel, hemisphere, sigma of hit
    float posX = (*ithit)->GetChan()->GetId()->GetPosX();
    float posY = (*ithit)->GetChan()->GetId()->GetPosY();
    float sigma = (*ithit)->GetChan()->GetCal()->GetPedSigma();

    // Get the fraction of the hit amplitude that should be assigned to this cluster
    float frac = 1.;
    if (_share == 2)
      frac = 1. / (*ithit)->GetNrClusters();

    sumt += frac * ((*ithit)->GetAmp()[_sample] - _thr*sigma);
    wsumx += frac * ((*ithit)->GetAmp()[_sample] - _thr*sigma) * posX;
    wsumy += frac * ((*ithit)->GetAmp()[_sample] - _thr*sigma) * posY;
  }

  // Center of gravity
  double centerx = wsumx / sumt;
  double centery = wsumy / sumt;

  // Error on center of gravity
  dcenter = 0.0833; // s / sqrt(12) ^ 2
  dcenter = sqrt(dcenter);

  // Update cluster properties
  fPosX = centerx;
  fPosY = centery;
  fPosXErr = dcenter;
  fPosYErr = dcenter;

  return;
}


void CsPixelMumegaCluster::CorrectPos(const std::vector<float> &_corrs) {
  // in total, we need 16 parameters, if we do not have all, we do not correct
  if ( _corrs.size() != 16 )
    return;

  // setting up parameters
  const float c010p0 = _corrs[ 0]; // category  10 - constant
  const float c010p1 = _corrs[ 1]; //                linear
  const float c010p2 = _corrs[ 2]; //                quadratic
  const float c010p3 = _corrs[ 3]; //                cubic
  const float c011p0 = _corrs[ 4]; // category  11 - constant
  const float c011p1 = _corrs[ 5]; //                linear
  const float c011p2 = _corrs[ 6]; //                quadratic
  const float c011p3 = _corrs[ 7]; //                cubic
  const float c104pX = _corrs[ 8]; // category 104 - shift in X
  const float c104pY = _corrs[ 9]; //                      in Y
  const float c105pX = _corrs[10]; // category 105 - shift in X
  const float c105pY = _corrs[11]; //                      in Y
  const float c106pX = _corrs[12]; // category 106 - shift in X
  const float c106pY = _corrs[13]; //                      in Y
  const float c107pX = _corrs[14]; // category 107 - shift in X
  const float c107pY = _corrs[15]; //                      in Y

  // Check if cluster is already categorized
  int fCategory = Categorize();

  // clusters of size 2 and 3 are corrected
  if (fHits.size() == 2) {
    // clusters of size 2 are corrected with an res vs eta function calibrated po13
    // Check for category and apply correction
    if ( fCategory == 10 ) {
      // calculate eta
      // eta is the distance from the most bottom left hit to the cluster position
      float eta;
      if (fHits.back()->GetChan()->GetId()->GetPosX() > fHits.front()->GetChan()->GetId()->GetPosX())
	eta = fPosX - fHits.front()->GetChan()->GetId()->GetPosX();
      else
	eta = fPosX - fHits.back()->GetChan()->GetId()->GetPosX();

      fPosX -= (c010p0 + c010p1 * eta + c010p2 * eta*eta + c010p3 * eta*eta*eta);
    } else if ( fCategory == 11 ) {
      float eta;
      if (fHits.back()->GetChan()->GetId()->GetPosY() > fHits.front()->GetChan()->GetId()->GetPosY())
	eta = fPosY - fHits.front()->GetChan()->GetId()->GetPosY();
      else
	eta = fPosY - fHits.back()->GetChan()->GetId()->GetPosY();

      fPosY -= (c011p0 + c011p1 * eta + c011p2 * eta*eta + c011p3 * eta*eta*eta);
    }
  } else if (fHits.size() == 3) {
    // clusters of size 3 in a square of 2x2 pads are corrected by shifting the center towards the missing pad
    // check for category and apply correction
    if ( fCategory == 104 ) {
      fPosX -= c104pX;
      fPosY -= c104pY;
    } else if ( fCategory == 105 ) {
      fPosX -= c105pX;
      fPosY -= c105pY;
    } else if ( fCategory == 106 ) {
      fPosX -= c106pX;
      fPosY -= c106pY;
    } else if ( fCategory == 107 ) {
      fPosX -= c107pX;
      fPosY -= c107pY;
    }
  }
}


int CsPixelMumegaCluster::Categorize() {
  int fCategory = 0;

  // restrain cluster size;
  if ( fHits.size() <= 3 ) {

    fCategory = 0;

    // Generate 3x3 Matrix and
    // Vectors to hold X and Y coordinates and as iterator for Vectors
    int M[3][3], X[3], Y[3], i = 0;
    // intitialize matrix
    for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
	M[x][y] = 0;

    // loop over hits in cluster to set X, Y vectors
    for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
      X[i] = (int) (*ithit)->GetChan()->GetId()->GetPosX();
      Y[i] = (int) (*ithit)->GetChan()->GetId()->GetPosY();
      i++;
    }

    for (int k = i; k < 3; k++) {
      X[k] = 0;
      Y[k] = 0;
    }

    // Save maximas

    int Xmin (33), Ymin (33);
    for (unsigned int k = 0; k < fHits.size(); k++) {
      Xmin = std::min(X[k], Xmin);
      Ymin = std::min(Y[k], Ymin);
    }

    // loop over hits to fill matrix
    for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
      // consistency checks
      assert ((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
      assert ((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
      assert ((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
      assert ((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

      // fill matrix
      M[(int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
    }

    // binary coding matrix
    int n = 1;
    for (int y = 0; y < 3; y++) {
      for (int x = 0; x < 3; x++) {
	fCategory += M[x][y] * n;
	n = n * 2;
      }
    }

    // re-code binary information into the final classification

    // clustersize 1
    if ( fCategory == 1 ) fCategory = 1; // size 1

    // clustersize 2
    else if ( fCategory ==   03 ) fCategory =  10; // horizontal
    else if ( fCategory ==  011 ) fCategory =  11; // vertical
    else if ( fCategory ==  012 ) fCategory =  12; // diagonal descending
    else if ( fCategory ==  021 ) fCategory =  13; // diagonal ascending

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
    // tell that an unknown form appeared
    // and push it out of normal coding range
    else {
      fCategory += 200;
      std::cout << "notcoded " << fCategory << std::endl;
    }
  }

  if ( fHits.size() == 4 ) {
    int Xmin (33), Ymin (33), Xmax(0), Ymax(0);
    for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
      Xmin = std::min((int)(*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((int)(*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
      Xmax = std::max((int)(*ithit)->GetChan()->GetId()->GetPosX(), Xmax);
      Ymax = std::max((int)(*ithit)->GetChan()->GetId()->GetPosY(), Ymax);
    }
    if ( Xmax - Xmin == 1 && Ymax - Ymin == 1 )
      fCategory = 200;
    else if ( Xmax - Xmin <= 2 && Ymax - Ymin <= 2 ) {
      fCategory = 0;
      assert (fCategory == 0);
      assert (fHits.size() == 4);
      // Generate 3x3 matrix and
      // Vectors to hold X and Y coordinates and as iterator for Vectors
      int M[3][3], X[4], Y[4], i = 0;
      // initialize matrix
      for (int y = 0; y < 3; y++)
	for (int x = 0; x < 3; x++)
	  M[x][y] = 0;

      // loop over hits in cluster to fill matrix
      for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
	// consistency checks
	assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

	// fill matrix
	M[(int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }

      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 3; y++) {
	for (int x = 0; x < 3; x++) {
	  fCategory += M[x][y] * n;
	  n = n*2;
	}
      }
      assert (fHits.size() == 4);
      if ( fCategory == 036 ) fCategory = 201; // -_
      else if ( fCategory ==  063 ) fCategory = 202; // _-
      else if ( fCategory == 0132 ) fCategory = 203; // |\|
      else if ( fCategory == 0231 ) fCategory = 204; // |/|
      else if ( fCategory ==  072 ) fCategory = 205; // T
      else if ( fCategory ==  027 ) fCategory = 206; // inv T
      else if ( fCategory == 0131 ) fCategory = 207; // rot left T or |-
      else if ( fCategory == 0232 ) fCategory = 208; // rot rught T or -|
      else if ( fCategory == 0113 ) fCategory = 209; // L
      else if ( fCategory == 0223 ) fCategory = 210; // mir (vert) L
      else if ( fCategory ==  047 ) fCategory = 211; // rot left L
      else if ( fCategory ==  074 ) fCategory = 212; // rot left mir (vert) L
      else if ( fCategory ==  071 ) fCategory = 213; // rot right L
      else if ( fCategory ==  017 ) fCategory = 214; // rot right mir (vert) L
      else if ( fCategory == 0322 ) fCategory = 215; // mir (ho) L
      else if ( fCategory == 0311 ) fCategory = 216; // mir (ho + vert) L
      else fCategory = 217;

    }
    else fCategory = 218;
  }

  if ( fHits.size() == 5 ) {
    int Xmin (33), Ymin (33), Xmax (0), Ymax (0);
    for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
      Xmin = std::min((int)(*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((int)(*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
      Xmax = std::max((int)(*ithit)->GetChan()->GetId()->GetPosX(), Xmax);
      Ymax = std::max((int)(*ithit)->GetChan()->GetId()->GetPosY(), Ymax);
    }
    if ( Xmax - Xmin <= 2 && Ymax - Ymin <= 2 ) {
      fCategory =0;
      assert (fCategory == 0);
      assert (fHits.size() == 5);
      // Generate 3x3 matrix and
      // vectors to hold X and Y coordinates ans as iterator for vectors
      int M[3][3], X[5], Y[5], i = 0;
      // initialize matrix
      for (int y = 0; y < 3; y++)
	for (int x = 0; x < 3; x++)
	  M[x][y] = 0;

      // loop over hits in cluster to set X, Y vectors
      for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
	X[i] = (int) (*ithit)->GetChan()->GetId()->GetPosX();
	Y[i] = (int) (*ithit)->GetChan()->GetId()->GetPosY();
	i++;
      }

      // loop over hits to fill matrix
      for ( std::list<CsMumegaHit*>::iterator ithit=fHits.begin(); ithit != fHits.end(); ithit++) {
        // consistency checks
          assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
          assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
          assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
          assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

        // fill matrix
        M[(int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
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
      else if ( fCategory == 0133 ) fCategory = 301;
      else if ( fCategory == 0233 ) fCategory = 302;
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
    }
    else fCategory = 314;
  }

  if ( fHits.size() == 6 ) {
    int Xmin (33), Ymin (33), Xmax (0), Ymax (0);
    for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit!= fHits.end(); ithit++) {
      Xmin = std::min((int)(*ithit)->GetChan()->GetId()->GetPosX(), Xmin);
      Ymin = std::min((int)(*ithit)->GetChan()->GetId()->GetPosY(), Ymin);
      Xmax = std::max((int)(*ithit)->GetChan()->GetId()->GetPosX(), Xmax);
      Ymax = std::max((int)(*ithit)->GetChan()->GetId()->GetPosY(), Ymax);
    }
    if ( Xmax - Xmin <= 2 && Ymax -Ymin <= 2 ) {
      fCategory = 0;
      assert (fCategory == 0);
      assert (fHits.size() == 6);
      // Generate 3x3 matrix and
      // vectors to hold X and Y coordinates and as iterator for vectors
      int M[3][3], X[6], Y[6], i = 0;
      // initialize matrix
      for (int y = 0; y < 3; y++)
	for (int x = 0; x < 3; x++)
	  M[x][y] = 0;

      // loop over hits in cluster to set X, Y vectors
      for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
	X[i] = (int) (*ithit)->GetChan()->GetId()->GetPosX();
	Y[i] = (int) (*ithit)->GetChan()->GetId()->GetPosY();
	i++;
      }

      // loop over hits to fill matrix
      for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
	// consistency checks
	assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 2);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 2);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

	// fill matrix
	M[(int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }
      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 3; y++) {
	for (int x = 0; x < 3; x++) {
	  fCategory += M[x][y] * n;
	  n = n*2;
	}
      }
      if ( fCategory == 077 ) fCategory = 400;
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
      fCategory = 0;
      assert (fCategory == 0);
      assert (fHits.size() == 6);
      // Generate 3x3 matrix and
      // vectors to hold X and Y coordinates and as iterator for vectors
      int M[9][9], X[6], Y[6], i = 0;
      // initialize matrix
      for (int y = 0; y < 9; y++)
	for (int x = 0; x < 9; x++)
	  M[x][y] = 0;

      // loop over hits in cluster to set X, Y vectors
      for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
	X[i] = (int) (*ithit)->GetChan()->GetId()->GetPosX();
	Y[i] = (int) (*ithit)->GetChan()->GetId()->GetPosY();
	i++;
      }

      // loop over hits to fill matrix
      for ( std::list<CsMumegaHit*>::iterator ithit = fHits.begin(); ithit != fHits.end(); ithit++) {
	// consistency checks
	assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin <= 8);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin >= 0);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin <= 8);
	assert((int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin >= 0);

	// fill matrix
	M[(int) (*ithit)->GetChan()->GetId()->GetPosX() - Xmin][(int) (*ithit)->GetChan()->GetId()->GetPosY() - Ymin] = 1;
      }

      // binary coding of matrix
      int n = 1;
      for (int y = 0; y < 9; y++) {
	if (false) std::cout << std::endl;
	for (int x = 0; x < 9; x++) {
	  fCategory += M[x][y] * n;
	  n = n*2;
	  // output of matrix
	  if (false) {
	    if ( M[x][y] )
	      std::cout << M[x][y];
	    else
	      std::cout << ".\t";
	  }
	}
      }
      if (false)
	std::cout << "\n\tBin 415 : " << fCategory << std::endl;
      fCategory = 415;
    }
  }

  if ( fHits.size() > 6 )
    fCategory = 1000;

  return fCategory;
}


bool CsPixelMumegaCluster::AddHit(CsMumegaHit* _hit, Int_t _clusconfmask, Float_t _thr, Int_t _sample) {
  std::list<CsMumegaChan*>::iterator itchan;
  std::list<CsMumegaHit*>::iterator ithit;
  bool insert = false; // should the current hit be added to the cluster
  bool changeflagged = false; // contains the expanded cluster flagged channels

  // read the clusterization configuration mask
  bool clpossizesplit = (_clusconfmask >> 4) & 1; // cut on the distance of hit to current CoG
  bool clpossizemethod = (_clusconfmask >> 5) & 1; // cut radial (true) or else (false)
  Int_t share = (_clusconfmask >> 6) & 3; // hit sharing between clusters

  // Add if list of hits in cluster empty
  if ( fHits.size() == 0 )
    insert = true;
  else { // Check if new hit should be added to existing list of hits
    // Get list of active neighbours for new hit
    std::list<CsMumegaChan*> neighbours = _hit->GetChan()->GetActiveNeighbours();
    for ( itchan = neighbours.begin(); itchan != neighbours.end(); itchan++ ) {
      // Check if any hit in cluster is a neighbour of the new hit
      for ( ithit = fHits.begin(); ithit != fHits.end(); ithit++ ) {
	if ((*itchan) == (*ithit)->GetChan()) {
	  insert = true;

	  // Check if there is a gap between neighbours
	  register Float_t distx = fabs( (*itchan)->GetId()->GetPosX() - (*ithit)->GetChan()->GetId()->GetPosX() );
	  register Float_t disty = fabs( (*itchan)->GetId()->GetPosY() - (*ithit)->GetChan()->GetId()->GetPosY() );
	  if ( distx > 1 || disty > 1 )
	    changeflagged = true;
	}
      }
    }

    // hit and cluster are neighbours, but does the hit fulfill all conditions
    if (insert) {
      // calculate current cluster center and dictancce to current hit
      CalcAmps(_sample, _thr, share);
      CalcCoG(_sample, _thr, share);
      float distx = fabs( _hit->GetChan()->GetId()->GetPosX() - fPosX );
      float disty = fabs( _hit->GetChan()->GetId()->GetPosY() - fPosY );
      float dists = distx*distx + disty*disty;

      // cluster position size splitting
      if (clpossizesplit) {
	if ( ( clpossizemethod && dists > 2. ) || ( !clpossizemethod && (distx > 1. || disty > 1.) ) )
	  insert = false;
      }
    }
  }

  // Add hit to cluster
  if ( insert) {

    fHits.push_back(_hit);
    if (changeflagged)
      fContainsFlaggedChannels = true;
  }

  return insert;
}


Bool_t CsPixelMumegaCluster::operator< (const CsPixelMumegaCluster &cluster) const {
  if (fPosX == cluster.GetPositionX())
    return (fPosY < cluster.GetPositionY());

  return (fPosX < cluster.GetPositionX());
}
