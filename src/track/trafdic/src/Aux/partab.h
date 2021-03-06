//
// Geant particle ID table
//
// NOTE: Particles
// 41,42, 35-40, 55-60, 71-77
// are redefined in Comgeant 


#ifndef partab_h
#define partab_h

const int Id2q[101]={  // particle ID to charge map 
  0,
  0, 1,-1, 0, 1,  //  1
  -1, 0, 1,-1, 0, //  6
  1,-1, 0, 1,-1,  // 11
  0, 0, 0, 1, 0,  // 16
  -1, 0,-1,-1, 0, // 21
  0,-1, 0, 1, 0,  // 26
  1, 1, 1,-1, 1,  // 31
  -1, 0, 0, 1,-1, // 36
  0, 0, -1, 0, 1, // 41
  1, 2, 0, 2, 0,  // 46
  0, 0, 0, 0, 1,  // 51
  -1, 0, 0, 1,-1, // 56
  3, 3, 4, 4, 5,  // 61
  5, 6, 7, 8, 9,  // 66
  1, 1, 0, 0, 2,  // 71
  1, 0, 0, 0, 0,  // 76
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0 };
  

const double Id2mass[101]={ // particle ID to mass map
  0, 
  0.0,       0.000511, 0.000511, 0.0,       0.1056583,   // 1
  0.1056583, 0.134937, 0.139567, 0.139567,  0.49767,     // 6
  0.493646,  0.493646, 0.939566, 0.9382723, 0.9382723,   // 11
  0.49767,   0.5488,   1.11563,  1.18937,   1.19255,     // 16
  1.19743,   1.3149,   1.32132,  1.67243,   0.939566,    // 21
  1.11563,   1.18937,  1.19255,  1.19743,   1.3149,      // 26
  1.32132,   1.67243,  1.7841,   1.7841,    1.8693,      // 31
  1.8693,    1.8645,   1.8645,   1.9693,    1.9693,      // 36
  0.4977,    1.1156,  81.0,     92.4,       1.875613,    // 41
  2.81448,   3.727417, 0.,       2.81448,   0.,          // 46
  0.,        0.,       0.,       0.,        2.010,       // 51
  2.010,     2.007,    2.007,    2.110,     2.110,       // 56
  5.60305,   6.53536,  6.53622,  8.39479,   9.32699,     // 61
 10.25510,  11.17793, 13.04378, 14.89917,  17.69690,     // 66
  2.285,     2.465,    2.470,    2.705,     2.453,       // 71
  2.454,     2.452,    0.,       0.,        0.,          // 76
  0.,        0.,       0.,       0.,        0.,   
  0.,        0.,       0.,       0.,        0. };   

const char Id2nam[101][9] = { // particle ID to name map
  " ",                                                    // 0
  "Gamma   ","Positron","Electron","Neutrino","Muon +  ", // 1	
  "Muon -  ","Pion 0  ","Pion +  ","Pion-   ","Kaon 0 l", // 6
  "Kaon +  ","Kaon -  ","Neutron ","Proton  ","A-proton", // 11
  "Kaon 0 s","Eta     ","Lambda  ","Sigma + ","Sigma 0 ", // 16
  "Sigma - ","Xi 0    ","Xi -    ","Omega   ","A-nutron", // 21
  "A-lambda","A-sigma-","A-sigma0","A-sigma+","A-xi 0  ", // 26
  "A-xi +  ","A-omega+","Tau +   ","Tau -   ","D +     ", // 31
  "D -     ","D 0     ","Anti D 0","Ds +    ","Ds -    ", // 36
  "K0 spec ","Lam0spec","W -     ","Z 0     ","Deutron ", // 41
  "Tritium ","Alpha   ","Geantino","He3     ","  ---   ", // 46
  "  ---   ","  ---   ","  ---   ","  ---   ","D + *   ", // 51   
  "D - *   ","D *   0 ","Anti-D0*","Ds + *  ","Ds - *  ", // 56  
  "Li6     ","Li7     ","Be7     ","Be9     ","B10     ", // 61
  "B11     ","C12     ","N14     ","O16     ","F19     ", // 66
  "Lam C + ","Xi C +  ","Xi C 0  ","Om C 0  ","Sig C ++", // 71
  "Sig C + ","Sig C 0 ","        ","        ","        ", // 76
  "        ","        ","        ","        ","        ",
  "        ","        ","        ","        ","        ",
  "        ","        ","        ","        ","        ",
  "        ","        ","        ","        ","        " };

#endif








