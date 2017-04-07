#include "Reco_config.h"
#include "GUIBaseCalorimeter.h"

Reco::GUIBaseCalorimeter::GUIBaseCalorimeter(QWidget* parent, const char* name, Qt::WFlags fl, bool xy_selection)
    : QWidget(parent, fl) {
    setObjectName(QString::fromUtf8(name));
    setupUi(this);

    if (!xy_selection) {
        y1_value_from->setDisabled(true);
        y2_value_from->setDisabled(true);
        y1_value_to->setDisabled(true);
        y2_value_to->setDisabled(true);
        y1_value_step->setDisabled(true);
        y2_value_step->setDisabled(true);
    }

    connect(cut1_value_min,  SIGNAL(valueChanged(int)), this, SLOT(selMinChanged()));
    connect(cut1_value_max,  SIGNAL(valueChanged(int)), this, SLOT(selMaxChanged()));
    connect(cut1_value_from, SIGNAL(valueChanged(int)), this, SLOT(fitXMinChanged()));
    connect(cut1_value_to,   SIGNAL(valueChanged(int)), this, SLOT(fitXMaxChanged()));
    connect(y1_value_from,   SIGNAL(valueChanged(int)), this, SLOT(fitYMinChanged()));
    connect(y1_value_to,     SIGNAL(valueChanged(int)), this, SLOT(fitYMaxChanged()));
}

Reco::GUIBaseCalorimeter::~GUIBaseCalorimeter() {
}

void Reco::GUIBaseCalorimeter::selMinChanged() {
    const int min = cut1_value_min->value();
    const int max = cut1_value_max->value();
    if (min > max)
        cut1_value_max->setValue(min);
}

void Reco::GUIBaseCalorimeter::selMaxChanged() {
    const int min = cut1_value_min->value();
    const int max = cut1_value_max->value();
    if (max < min)
        cut1_value_min->setValue(max);
}

void Reco::GUIBaseCalorimeter::fitXMinChanged() {
    const int min = cut1_value_from->value();
    const int max = cut1_value_to->value();
    if (min > max)
        cut1_value_to->setValue(min);
}

void Reco::GUIBaseCalorimeter::fitXMaxChanged() {
    const int min = cut1_value_from->value();
    const int max = cut1_value_to->value();
    if (max < min)
        cut1_value_from->setValue(max);
}

void Reco::GUIBaseCalorimeter::fitYMinChanged() {
    const int min = y1_value_from->value();
    const int max = y1_value_to->value();
    if (min > max)
        y1_value_to->setValue(min);
}

void Reco::GUIBaseCalorimeter::fitYMaxChanged() {
    const int min = y1_value_from->value();
    const int max = y1_value_to->value();
    if (max < min)
        y1_value_from->setValue(max);
}
