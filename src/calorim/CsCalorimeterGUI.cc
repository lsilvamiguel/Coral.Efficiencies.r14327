
#include "Reco/Calorimeter.h"
#include "coral_config.h"
#include "CsECAL1.h"

#include "CsCalorimeterGUI.h"
#include "CsECAL1.h"

using namespace std;

using namespace Reco;

////////////////////////////////////////////////////////////////////////////////

CsCalorimeterGUI::CsCalorimeterGUI(Reco::Calorimeter &c, QWidget* parent, const char* name, Qt::WFlags fl) :
  GUICalorimeter(c,parent,name,fl,c.XYRegularGrid())
{
  Update();
}

////////////////////////////////////////////////////////////////////////////////

CsCalorimeterGUI:: ~CsCalorimeterGUI (void)
{
}

////////////////////////////////////////////////////////////////////////////////

CsECAL1GUI::CsECAL1GUI(Reco::Calorimeter &c, QWidget* parent, const char* name, Qt::WFlags fl) :
    CsCalorimeterGUI(c,parent,name,fl)
{

  bool debug=true;
  if( debug ) std::cout << " CsECAL1GUI::Draw() " << c.GetName() << std::endl;
#if USE_Qt
  QGridLayout* gridLayout = dynamic_cast<QGridLayout*>(new_select->layout());

  selection_3 = new QComboBox(new_select);
  selection_3->setObjectName(QString::fromUtf8("selection_3"));
  selection_3->setEnabled (TRUE);

  selection_3->insertItem(0, tr( "GAMS" ) );
  selection_3->insertItem(1, tr( "MaintzTop" ) );
  selection_3->insertItem(2, tr( "MaintzBottom" ) );
  selection_3->insertItem(3, tr( "OlgaSaleve" ) );
  selection_3->insertItem(4, tr( "OlgaJura" ) );

  gridLayout->addWidget(selection_3, gridLayout->rowCount(), 0, 1, gridLayout->columnCount());
#endif
  Update();
  xy_selection=false;
}

////////////////////////////////////////////////////////////////////////////////

CsECAL1GUI::~CsECAL1GUI (void)
{
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1GUI::Draw (void)
{

//   bool debug=true;
  bool debug=false;
  if( debug ) std::cout << " CsECAL1GUI::Draw() " << cal.GetName() << endl;
#if USE_Qt

  int min_value = cut1_value_min->value ();
  int max_value = cut1_value_max->value ();
  int units = selection_units->currentIndex();
  double units_scale=1.;
  if(units == 1 ) units_scale=0.1;
  else if(units == 2 ) units_scale=0.01;
  else if(units == 3 ) units_scale=10.;
  else if(units == 4 ) units_scale=100.;
  double fmin = min_value*units_scale;
  double fmax = max_value*units_scale;

  int selection_type= selection->currentIndex();
//   cout << " CsECAL1GUI::Draw selection " << selection_type  << endl;
  if(selection_type < 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  else if(selection_type == 5) info_type = Reco::Calorimeter::NOISE;
  else if(selection_type == 0) info_type = Reco::Calorimeter::RAW;

  int xmin = cut1_value_from->value ();
  int xhow_much = cut1_value_to->value ();
  int xstep = cut1_value_step->value ();

  int ymin = y1_value_from->value ();
  int yhow_much = y1_value_to->value ();
  int ystep = y1_value_step->value ();

  int subset = selection_3->currentIndex();
  int subset_index = 0;

  CsECAL1 *ecal1  = reinterpret_cast <CsECAL1 *> (&cal);
  if( subset == 0 )
    subset_index = ecal1->index_gams_;
  else if( subset == 1 )
    subset_index = ecal1->index_maintz_top_;
  else if( subset == 2 )
    subset_index = ecal1->index_maintz_bottom_;
  else if( subset == 3 )
    subset_index = ecal1->index_olga_saleve_;
  else if( subset == 4 )
    subset_index = ecal1->index_olga_jura_;

  cal.DrawFitCellsHistoXYSubset(subset_index, info_type,"DRAW",xmin,xhow_much,xstep,ymin,yhow_much,ystep,fmin,fmax);

#endif
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1GUI::Fit (void)
{
#if USE_Qt

  int min_value = cut1_value_min->value ();
  int max_value = cut1_value_max->value ();
  int units = selection_units->currentIndex();
  double units_scale=1.;
  if(units == 1 ) units_scale=0.1;
  else if(units == 2 ) units_scale=0.01;
  else if(units == 3 ) units_scale=10.;
  else if(units == 4 ) units_scale=100.;
  double fmin = min_value*units_scale;
  double fmax = max_value*units_scale;

  int selection_type= selection->currentIndex();
//   cout << " CsECAL1GUI::Draw selection " << selection_type  << endl;
  if(selection_type < 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  else if(selection_type == 5) info_type = Reco::Calorimeter::NOISE;
  else if(selection_type == 0) info_type = Reco::Calorimeter::RAW;

  if( selection_2->currentIndex() == 6 )   // auto fit
  {
    cal.FitAllCellsHisto(info_type,"AUTO");
    return;
  }

  int xmin = cut1_value_from->value ();
  int xhow_much = cut1_value_to->value ();
  int xstep = cut1_value_step->value ();

  int ymin = y1_value_from->value ();
  int yhow_much = y1_value_to->value ();
  int ystep = y1_value_step->value ();

  int subset = 0; selection_3->currentIndex();
  int subset_index = 0;

  CsECAL1 *ecal1  = reinterpret_cast <CsECAL1 *> (&cal);
  if( subset == 0 )
    subset_index = ecal1->index_gams_;
  else if( subset == 1 )
    subset_index = ecal1->index_maintz_top_;
  else if( subset == 2 )
    subset_index = ecal1->index_maintz_bottom_;
  else if( subset == 3 )
    subset_index = ecal1->index_olga_saleve_;
  else if( subset == 4 )
    subset_index = ecal1->index_olga_jura_;

  cal.DrawFitCellsHistoXYSubset(subset_index, info_type,"FIT",xmin,xhow_much,xstep,ymin,yhow_much,ystep,fmin,fmax);
#endif
}

////////////////////////////////////////////////////////////////////////////////

void CsECAL1GUI::FitOK (void)
{
#if USE_Qt
  int selection_type= selection->currentIndex();
//   cout << " CsECAL1GUI::Draw selection " << selection_type  << endl;
  if(selection_type < 0) return;
  Reco::Calorimeter::CellInfoType info_type=Reco::Calorimeter::CALIB;
  if(selection_type == 1) info_type = Reco::Calorimeter::CALIB;
  else if(selection_type == 2) info_type = Reco::Calorimeter::LED;
  else if(selection_type == 3) info_type = Reco::Calorimeter::PED;
  else if(selection_type == 4) info_type = Reco::Calorimeter::TIME;
  else if(selection_type == 5) info_type = Reco::Calorimeter::NOISE;
  else if(selection_type == 0) info_type = Reco::Calorimeter::RAW;


  int xmin = cut1_value_from->value ();
  int xhow_much = cut1_value_to->value ();
  int xstep = cut1_value_step->value ();

  int ymin = y1_value_from->value ();
  int yhow_much = y1_value_to->value ();
  int ystep = y1_value_step->value ();

  int subset = selection_3->currentIndex();
  int subset_index = 0;

  CsECAL1 *ecal1  = reinterpret_cast <CsECAL1 *> (&cal);
  if( subset == 0 )
    subset_index = ecal1->index_gams_;
  else if( subset == 1 )
    subset_index = ecal1->index_maintz_top_;
  else if( subset == 2 )
    subset_index = ecal1->index_maintz_bottom_;
  else if( subset == 3 )
    subset_index = ecal1->index_olga_saleve_;
  else if( subset == 4 )
    subset_index = ecal1->index_olga_jura_;

  cal.DrawFitCellsHistoXYSubset(subset_index, info_type,"FIT_ON_OFF",xmin,xhow_much,xstep,ymin,yhow_much,ystep);
#endif
}
