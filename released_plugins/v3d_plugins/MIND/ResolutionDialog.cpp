#include "ResolutionDialog.h"

ResolutionDialog::ResolutionDialog(QWidget *parent)
    : QDialog(parent), parent(parent) {
  QFormLayout *formLayout = new QFormLayout;

  QLabel *helpDescriptionLabel = new QLabel(
      QString("Enter resolutions of the image voxel in µm along the 3 axes:"),
      this);
  formLayout->addRow(helpDescriptionLabel);

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

void ResolutionDialog::setResolutionOfImage(Image4DSimple *p4DImage) {
  // Landmark info and ask for desired resolution of a image pixel along the 3
  // axes for isotropic correction
  double rez[3];
  if (exec() == QDialog::Accepted) {
    rez[0] = getXResolution();
    rez[1] = getYResolution();
    rez[2] = getZResolution();
  } else {
    qDebug() << "error: failed to get resolution";
    return;
  }

  // Set resolution of the image
  p4DImage->setRezX(rez[0]);
  p4DImage->setRezY(rez[1]);
  p4DImage->setRezZ(rez[2]);
}