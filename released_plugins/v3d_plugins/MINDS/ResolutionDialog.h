 /**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#ifndef __RESOLUTIONDIALOG_H__
#define __RESOLUTIONDIALOG_H__

#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QVBoxLayout>

#include "basic_4dimage.h"

class ResolutionDialog : public QDialog {
  Q_OBJECT

 public:
  ResolutionDialog(QWidget *parent = nullptr);

  double getXResolution() const;
  double getYResolution() const;
  double getZResolution() const;

  void setResolutionOfImage(Image4DSimple *p4DImage);

 private:
  double xResolution;
  double yResolution;
  double zResolution;
  QWidget *parent;
  QDoubleSpinBox *xSpinBox;
  QDoubleSpinBox *ySpinBox;
  QDoubleSpinBox *zSpinBox;
};

#endif  // __RESOLUTIONDIALOG_H__