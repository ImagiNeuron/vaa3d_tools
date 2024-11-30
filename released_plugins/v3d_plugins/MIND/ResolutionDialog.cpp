#include "ResolutionDialog.h"

ResolutionDialog::ResolutionDialog(float x, float y, float z, float radius,
                                   QWidget *parent)
    : QDialog(parent) {
  QFormLayout *formLayout = new QFormLayout;

  QLabel *landmarkLabel = new QLabel(
      QString("Landmark: x=%1, y=%2, z=%3, radius=%4\n\n"
              "Enter resolutions of the image pixel along the 3 axes:")
          .arg(x)
          .arg(y)
          .arg(z)
          .arg(radius),
      this);
  formLayout->addRow(landmarkLabel);

  xSpinBox = new QDoubleSpinBox(this);
  ySpinBox = new QDoubleSpinBox(this);
  zSpinBox = new QDoubleSpinBox(this);

  // Set default values
  xSpinBox->setValue(1.0);
  ySpinBox->setValue(1.0);
  zSpinBox->setValue(1.0);

  // Set range for spin boxes
  xSpinBox->setRange(0.01, 100.0);
  ySpinBox->setRange(0.01, 100.0);
  zSpinBox->setRange(0.01, 100.0);

  // Set step size
  xSpinBox->setSingleStep(1.0);
  ySpinBox->setSingleStep(1.0);
  zSpinBox->setSingleStep(1.0);

  formLayout->addRow("Resolution X:", xSpinBox);
  formLayout->addRow("Resolution Y:", ySpinBox);
  formLayout->addRow("Resolution Z:", zSpinBox);

  QDialogButtonBox *buttonBox = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
  connect(buttonBox, &QDialogButtonBox::accepted, this, &QDialog::accept);
  connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(formLayout);
  mainLayout->addWidget(buttonBox);

  setLayout(mainLayout);
  setWindowTitle("Enter Resolutions");
}

double ResolutionDialog::getXResolution() const { return xSpinBox->value(); }

double ResolutionDialog::getYResolution() const { return ySpinBox->value(); }

double ResolutionDialog::getZResolution() const { return zSpinBox->value(); }