#ifndef CsCalorimeterGUI___include
#define CsCalorimeterGUI___include

#include "coral_config.h"


#include "Reco/src/GUICalorimeter.h"

#if USE_Qt
#  include <Qt/qobject.h>
#else
typedef int QWidget, WFlags;
class QComboBox;
#endif

class Calorimeter;

class CsCalorimeterGUI : public Reco::GUICalorimeter
{
#if USE_Qt
  Q_OBJECT
#endif

  // =========================================================================
  // Constructors and destructor
  // =========================================================================

 public:

  CsCalorimeterGUI       ( Reco::Calorimeter &c, QWidget* parent = 0, const char* name = 0, Qt::WFlags fl = 0 );

  // By making a base class destructor virtual, one ensures that the destructor of any derived class is executed (in addition to, and prior to, that of the base class).
  virtual ~CsCalorimeterGUI();

 private:
//  CsCalorimeterGUI       (const CsCalorimeterGUI &g);

  // =========================================================================
  // Operators
  // =========================================================================

  CsCalorimeterGUI   &operator =              (const CsCalorimeterGUI &g);

  // =========================================================================
  // Methods
  // =========================================================================
protected:
#if USE_Qt
    protected slots:
#endif

};

class CsECAL1GUI : public CsCalorimeterGUI
{
#if USE_Qt
  Q_OBJECT
#endif

  // =========================================================================
  // Constructors and destructor
  // =========================================================================

 public:

  CsECAL1GUI    ( Reco::Calorimeter &c,
			      QWidget* parent = 0, const char* name = 0, Qt::WFlags fl = 0 );

   virtual ~CsECAL1GUI();

 private:
  CsECAL1GUI                (const CsECAL1GUI &g);

  // =========================================================================
  // Operators
  // =========================================================================

  CsECAL1GUI           &operator =              (const CsECAL1GUI &g);

  // =========================================================================
  // Methods
  // =========================================================================
 protected:
#if USE_Qt
 protected slots:
#endif

  void        Draw                    (void);

  void        Fit                     (void);

  void        FitOK                   (void);

  // =========================================================================
  // Members
  // =========================================================================
 private:
  QComboBox* selection_3;
};

#endif // CsCalorimeterGUI___include
