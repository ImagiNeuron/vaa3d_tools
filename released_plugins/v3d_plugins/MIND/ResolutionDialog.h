/* ResolutionDialog.h
 * This dialog allows the user to enter resolutions for a landmark and its image
 * to be used in the soma segmentation plugin to correct for isotropic and
 * z-thickness
 *
 * 2024-11-30 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette
 */

#ifndef __RESOLUTIONDIALOG_H__
#define __RESOLUTIONDIALOG_H__

#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QVBoxLayout>

class ResolutionDialog : public QDialog {
  Q_OBJECT

 public:
  ResolutionDialog(float x, float y, float z, float radius,
                   QWidget *parent = nullptr);

  double getXResolution() const;
  double getYResolution() const;
  double getZResolution() const;

 private:
  QDoubleSpinBox *xSpinBox;
  QDoubleSpinBox *ySpinBox;
  QDoubleSpinBox *zSpinBox;
};

#endif  // __RESOLUTIONDIALOG_H__