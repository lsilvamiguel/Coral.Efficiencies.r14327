#ifndef  RCHISTOS_H
#define  RCHISTOS_H

/*!
   \file    CsRCHistos.h
   \brief   Histogram booking for CsRichOne.
   \author  Paolo Schiavon
   \version $Revision: 1.59 $
   \date    $Date: 2009/02/13 15:08:17 $
*/


  #include "CsHistograms.h"


  class CsRCHistos {


    public:


  static CsRCHistos& Ref();


  CsHist1D *hRC1000, *hRC1001;
  CsHist2D *hRC1011, *hRC1021, *hRC1022;

  CsHist2D *hRC1031;

  CsHist2D *hRC1061, *hRC1062, *hRC1063, *hRC1066, *hRC1068;

  CsHist1D *hRC1121, *hRC1131;

  CsHist2D *hRC1200, *hRC1201, *hRC1202, *hRC1203;
  CsHist2D *hRC1204, *hRC1205, *hRC1206, *hRC1207, *hRC1208, *hRC1209;

  CsHist2D *hRC1220;

  std::vector<CsHist2D*> vRC1350;
  std::vector<CsHist2D*> vRC1370;

  std::vector<CsHist2D*> vRC1500;
  CsHist2D *hRC1505, *hRC1506;
  CsHist1D *hRC1507, *hRC1508;
  CsHist2D *hRC1510, *hRC1511, *hRC1512, *hRC1513, *hRC1514, *hRC1515;
  CsHist2D *hRC1516, *hRC1517, *hRC1518, *hRC1519;
  CsHist2D *hRC1520, *hRC1522, *hRC1523, *hRC1524, *hRC1525, *hRC1527;
  CsHist1D *hRC1521, *hRC1526, *hRC1528;
  CsHist2D *hRC1529, *hRC1530;
  CsHist1D *hRC1531, *hRC1533, *hRC1535, *hRC1536, *hRC1537, *hRC1539;
  CsHist2D *hRC1540;

  CsHist1D *hRC1545, *hRC1547, *hRC1548;

  CsHist2D *hRC1556, *hRC1557;
  CsHist1D *hRC1558, *hRC1559;

  CsHist2D *hRC1551;
  CsHist1D *hRC1571, *hRC1572, *hRC1573;
  CsHist1D *hRC1574, *hRC1575, *hRC1576, *hRC1577, *hRC1578, *hRC1579;
  CsHist1D *hRC1581, *hRC1582, *hRC1583, *hRC1584, *hRC1585;

  CsHist1D *hRC1568;
  CsHist2D *hRC1569;
  std::vector<CsHist2D*> vRC1590;

  std::vector<CsHist1D*> vRC1600;
  CsHist1D *hRC1621, *hRC1630;
  CsHist2D *hRC1623, *hRC1624, *hRC1625, *hRC1626, *hRC1627;
  CsHist1D *hRC1631, *hRC1632, *hRC1633, *hRC1636, *hRC1639;
  CsHist2D *hRC1634;
  CsHist2D *hRC1635, *hRC1637, *hRC1638;
  CsHist2D *hRC1640, *hRC1641, *hRC1642, *hRC1643;

  CsHist2D *hRC1661, *hRC1664; 
  CsHist1D *hRC1662, *hRC1666, *hRC1667, *hRC1668;
  CsHist2D *hRC1692;
  CsHist2D *hRC1693, *hRC1694,*hRC1695, *hRC1696;
  std::vector<CsHist1D*> vRC1650;
  CsHist2D *hRC1655, *hRC1663, *hRC1665, *hRC1669, *hRC1691;
  std::vector<CsHist2D*> vRC1655;
  std::vector<CsHist2D*> vRC1670;
  std::vector<CsHist2D*> vRC1680;

  CsHist2D *hRC1701, *hRC1702, *hRC1703, *hRC1705, *hRC1708;

  //std::vector<CsHist2D*> vRC1800;
  //CsHist2D *hRC1980;

  std::vector<CsHist2D*> vRC2100;
  std::vector<CsHist2D*> vRC2105;
  std::vector<CsHist2D*> vRC2110;
  std::vector<CsHist1D*> vRC2120;

  CsHist2D *hRC3008, *hRC3009, *hRC3010, *hRC3011, *hRC3012, *hRC3013,
           *hRC3014, *hRC3028, *hRC3029, *hRC3034;
  CsHist1D *hRC3001, *hRC3002, *hRC3015, *hRC3016, *hRC3017, *hRC3019,
           *hRC3020, *hRC3022, *hRC3023, *hRC3025, *hRC3026, *hRC3031;
  CsHist2D *hRC3101, *hRC3111;

  CsHist1D *hRC3250, *hRC3252;
  CsHist2D *hRC3255, *hRC3256;

  CsHist2D *hRC3501;
  CsHist2D *hRC3502, *hRC3503, *hRC3504, *hRC3505, *hRC3507;
  CsHist2D *hRC3506, *hRC3508, *hRC3509, *hRC3510;

  CsHist1D *hRC3525, *hRC3529, *hRC3550, *hRC3555, *hRC3557;
  CsHist1D *hRC3559, *hRC3560, *hRC3562;
  CsHist2D *hRC3526, *hRC3527, *hRC3528, *hRC3530, *hRC3531, *hRC3533;
  CsHist2D *hRC3551, *hRC3553, *hRC3554, *hRC3556, *hRC3558, *hRC3561;
  std::vector<CsHist2D*> vRC3520;
  std::vector<CsHist2D*> vRC3535, vRC3540, vRC3545;
  std::vector<CsHist2D*> vRC3570, vRC3580, vRC3590;
  CsHist2D *hRC3597;
  std::vector<CsHist2D*> vRC3565;
  CsHist2D *hRC3563;
  CsHist2D *hRC3569;

  CsHist2D *hRC3519, *hRC3520, *hRC3601, *hRC3602, *hRC3603, *hRC3604;
  CsHist1D *hRC3006, *hRC3036, *hRC3037, *hRC3042, *hRC3043,
           *hRC3045, *hRC3046, *hRC3047, *hRC3050, *hRC3051;
  CsHist1D *hRC3052, *hRC3057, *hRC3058, *hRC3059, *hRC3060;
  CsHist1D *hRC3061, *hRC3062, *hRC3063, *hRC3065, *hRC3067;
  CsHist2D *hRC3053, *hRC3054, *hRC3055, *hRC3056, *hRC3064;
  CsHist1D *hRC3072, *hRC3073, *hRC3074, *hRC3075, *hRC3076,
           *hRC3077, *hRC3078, *hRC3079;
  CsHist1D *hRC3610, *hRC3611;

  CsHist2D *hRC3605, *hRC3606, *hRC3607, *hRC3608, *hRC3609;
  CsHist2D *hRC3613, *hRC3614, *hRC3618;
  CsHist1D *hRC3617;

  CsHist1D *hRC3631;
  CsHist2D *hRC3615, *hRC3632;
  std::vector<CsHist1D*> vRC3620;

  CsHist2D *hRC3635, *hRC3636, *hRC3637;

  CsHist1D *hRC3640;
  std::vector<CsHist2D*> vRC3640;

  std::vector<CsHist2D*>  vRC3660;
  CsHist2D *hRC3677, *hRC3678;
  std::vector<CsHist2D*>  vRC3680, vRC3690;
  CsHist2D *hRC3687, *hRC3688, *hRC3697, *hRC3698;

  std::vector<CsHist2D*> vRC3700, vRC3710, vRC3720, vRC3730;

  CsHist2D *hRC3751, *hRC3765;
  CsHist1D *hRC3750, *hRC3760;

  //std::vector<CsHist2D*> vRC3800;

  std::vector<CsHist2D*> vRC6550;


  inline  bool bookHis() const { return bookHis_; };
  inline  int levelBk() const { return levelBk_; };

  void print() const;


    protected:


  CsRCHistos();

  ~CsRCHistos();


    private:


  static CsRCHistos* Ptr_;

  bool bookHis_;
  int levelBk_;

 
  };

#endif
