// To apply the whole procedure 
// (i.e. going from CsI coordinates to lens and PMT number)
// include something
// like the following in your code:

#include "RayTrace.h"

// these variables must be read from the zebra bank:
// naming as in http://valexakh.home.cern.ch/valexakh/wwwcomg/comgeant.doc.html#ntu
double xr,  yr,  zr,  xm,  ym,  zm,  xd,  yd, et; 

static OptSyst Telescope; // must be static to avoid memory leaks

LightRay photon;

int ipin=-1, lens=-1;

lens=Telescope.FindFirstLens( xr,  yr,  zr,  xm,  ym,  zm,  xd,  yd,  et, photon);

// lens numbering scheme
//
// y
// ^
// | 133 ... 144
// |  .
// |  .
// |  13 ...  24 
// |   1 ...  12
// -----------------> x
//
// where x and y are the CsI coordinates (xd,yd)

if(lens>=0) ipin=Telescope.RayTrace(photon);

// ipin (PMT pixel) numbering scheme:
//  
//  y
//  ^  12 13 14 15
//  |   8  9 10 11
//  |   4  4  6  7
//  |   0  1  2  3
//  ---------------> x
//
// where x and y are the coordinates in the 
// ZEMAX local frame of the PMT 
// (which is rotated wrt the CsI frame and upside down
// for the bottom panels!)
// That means: x in the PMT frame corresponds to -x in the CsI frame
//             for the top panels (y correponds to y) and 
//             y in the PMT frame corresponds to -y in the CsI frame
//             for the bottom panels (x corresponds to x)...
//
// ATTENTION: The corresponding pseudo-pads are "inverted" 
//            by the optical imaging i.e.
// pixel 15 = pseudo-pad 0
// pixel 14 = pseudo-pad 1
// etc.




