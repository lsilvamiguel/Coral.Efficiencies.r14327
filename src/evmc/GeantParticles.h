// $Id: GeantParticles.h,v 1.2 2003/03/11 20:59:21 ybedfer Exp $

/*!
   \file    GeantParticles.h
   \brief   Geant particles table.
   \author  Benigno Gobbo
   \version $Revision: 1.2 $
   \date    $Date: 2003/03/11 20:59:21 $
*/

#ifndef GeantParticles_h
#define GeantParticles_h

const int nGeantParticles = 91;   //!< Total number of particles

/*!
   \struct  GeantParticles
   \brief   Geant particles table.
*/

const struct { 
  int    number;        //!< Particle Number (PDG convention)
  char   name[20];      //!< Particle Name
  double mass;          //!< Particle Mass (GeV/c^2)
  double width;         //!< Particle Width (GeV/c^2)
  int    charge;        //!< Particle Charge
  bool   antiparticle;  //!< True if sntiparticle
} GeantPart[nGeantParticles] = {
  {  1, "GAMMA",              0.00000, 0.1000e+16,  0, false },
  {  2, "POSITRON",           0.00051, 0.1000e+16,  1, true  },
  {  3, "ELECTRON",           0.00051, 0.1000e+16, -1, false },
  {  4, "NEUTRINO",           0.00000, 0.1000e+16,  0, false },
  {  5, "MUON +",             0.10566, 0.2197e-05,  1, true  },
  {  6, "MUON -",             0.10566, 0.2197e-05, -1, false },
  {  7, "PION 0",             0.13498, 0.8400e-16,  0, false },
  {  8, "PION +",             0.13957, 0.2603e-07,  1, false },
  {  9, "PION -",             0.13957, 0.2603e-07, -1, true  },
  { 10, "KAON 0 LONG",        0.49767, 0.5170e-07,  0, false },
  { 11, "KAON +",             0.49368, 0.1237e-07,  1, false },
  { 12, "KAON -",             0.49368, 0.1237e-07, -1, true  },
  { 13, "NEUTRON",            0.93957, 0.8870e+03,  0, false },
  { 14, "PROTON",             0.93827, 0.1000e+16,  1, false },
  { 15, "ANTIPROTON",         0.93827, 0.1000e+16, -1, true  },
  { 16, "KAON 0 SHORT",       0.49767, 0.8926e-10,  0, false },
  { 17, "ETA",                0.54745, 0.5485e-18,  0, false },
  { 18, "LAMBDA",             1.11568, 0.2632e-09,  0, false },
  { 19, "SIGMA +",            1.18937, 0.7990e-10,  1, false },
  { 20, "SIGMA 0",            1.19255, 0.7400e-19,  0, false },
  { 21, "SIGMA -",            1.19744, 0.1479e-09, -1, false },
  { 22, "XI 0",               1.31490, 0.2900e-09,  0, false },
  { 23, "XI -",               1.32132, 0.1639e-09, -1, false },
  { 24, "OMEGA -",            1.67245, 0.8220e-10, -1, false },
  { 25, "ANTINEUTRON",        0.93957, 0.8870e+03,  0, true  }, 
  { 26, "ANTILAMBDA",         1.11568, 0.2632e-09,  0, true  }, 
  { 27, "ANTISIGMA -",        1.18937, 0.7990e-10, -1, true  }, 
  { 28, "ANTISIGMA 0",        1.19255, 0.7400e-19,  0, true  }, 
  { 29, "ANTISIGMA +",        1.19744, 0.1479e-09,  1, true  }, 
  { 30, "ANTIXI 0",           1.31490, 0.2900e-09,  0, true  }, 
  { 31, "ANTIXI +",           1.32132, 0.1639e-09,  1, true  }, 
  { 32, "ANTIOMEGA +",        1.67245, 0.8220e-10,  1, true  }, 
  { 33, "",                   0.00000, 0.00000000,  0, false }, 
  { 34, "",                   0.00000, 0.00000000,  0, false }, 
  { 35, "D+",                 1.8693,  0.1062e-11,  1, false }, 
  { 36, "D-",                 1.8693,  0.1062e-11,  0, true  }, 
  { 37, "D0",                 1.8645,  0.4280e-12, -1, false }, 
  { 38, "Anti D0 (anti-c)",   1.8645,  0.4280e-12,  0, true  }, 
  { 39, "Ds+",                1.9685,  0.4670e-12,  1, false }, 
  { 40, "Ds-",                1.9685,  0.4670e-12, -1, false }, 
  { 41, "K0shortdetectable",  0.49770, 0.8920e-10,  0, false }, 
  { 42, "Lam0detectable",     1.11560, 0.2630e-09,  0, false }, 
  { 43, "Sigma- nodec,hadr",  1.19740, 0.1479e-09, -1, false }, 
  { 44, "Proton no hadr int", 0.93830, 0.1000e+23,  1, false }, 
  { 45, "DEUTERON",           1.87561, 0.1000e+16,  1, false }, 
  { 46, "TRITON",             2.80925, 0.1000e+16,  1, false }, 
  { 47, "ALPHA",              3.72742, 0.1000e+16,  2, false }, 
  { 48, "GEANTINO",           0.00000, 0.1000e+16,  0, false }, 
  { 49, "HE3",                2.80923, 0.1000e+16,  2, false }, 
  { 50, "Cerenkov",           0.00000, 0.1000e+16, -1, false }, 
  { 51, "Xi- spec dec",       1.32132, 0.1639e-09, -1, false }, 
  { 52, "Omega - spec dec",   1.67245, 0.8220e-10, -1, false }, 
  { 53, "Sigma+2",            1.18937, 0.7990e-10,  1, false }, 
  { 54, "Sigma-2",            1.19744, 0.1480e-09, -1, false }, 
  { 55, "D+*",                2.0100,  0.1000e-18,  1, false }, 
  { 56, "D-*",                2.0100,  0.1000e-18, -1, true  }, 
  { 57, "D0*",                2.0067,  0.1000e-18,  0, true  }, 
  { 58, "Anti D0* (anti-c)",  2.0067,  0.1000e-18,  0, false }, 
  { 59, "Ds+*",               2.1124,  0.1000e-18,  1, false }, 
  { 60, "Ds-*",               2.1124,  0.1000e-18, -1, true  }, 
  { 61, "",                   0.00000, 0.00000000,  0, false }, 
  { 62, "",                   0.00000, 0.00000000,  0, false }, 
  { 63, "",                   0.00000, 0.00000000,  0, false }, 
  { 64, "",                   0.00000, 0.00000000,  0, false }, 
  { 65, "",                   0.00000, 0.00000000,  0, false }, 
  { 66, "",                   0.00000, 0.00000000,  0, false }, 
  { 67, "",                   0.00000, 0.00000000,  0, false }, 
  { 68, "",                   0.00000, 0.00000000,  0, false }, 
  { 69, "",                   0.00000, 0.00000000,  0, false }, 
  { 70, "",                   0.00000, 0.00000000,  0, false }, 
  { 71, "LamC+",              2.28500, 0.2000e-12,  1, false }, 
  { 72, "XiC+",               2.46500, 0.3500e-12,  1, false }, 
  { 73, "XiC0",               2.47000, 0.1000e-12,  0, false }, 
  { 74, "OmC0",               2.70500, 0.8000e-13,  0, false }, 
  { 75, "SigC++",             2.45300, 0.1000e-18,  2, false }, 
  { 76, "SigC+",              2.45400, 0.1000e-18,  1, false }, 
  { 77, "SigC0",              2.45200, 0.1000e-18,  0, false }, 
  { 78, "",                   0.00000, 0.00000000,  0, false }, 
  { 79, "",                   0.00000, 0.00000000,  0, false }, 
  { 80, "",                   0.00000, 0.00000000,  0, false }, 
  { 81, "anti-LamC+",         2.28500, 0.2000e-12, -1, true  }, 
  { 82, "",                   0.00000, 0.00000000,  0, false }, 
  { 83, "",                   0.00000, 0.00000000,  0, false }, 
  { 84, "",                   0.00000, 0.00000000,  0, false }, 
  { 85, "",                   0.00000, 0.00000000,  0, false }, 
  { 86, "",                   0.00000, 0.00000000,  0, false }, 
  { 87, "",                   0.00000, 0.00000000,  0, false }, 
  { 88, "",                   0.00000, 0.00000000,  0, false }, 
  { 89, "",                   0.00000, 0.00000000,  0, false }, 
  { 90, "",                   0.00000, 0.00000000,  0, false }, 
  { 91, "Phi       ",         1.01941, 0.1490e-21,  0, false } 
};
#endif //GeantParticles_h
