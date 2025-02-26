/* soma_segmentation_plugin.h
 * A plugin for analysis of neuron somas in the brain.
 * 2024-11-16: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette
 */

#ifndef __SOMA_SEGMENTATION_PLUGIN_H__
#define __SOMA_SEGMENTATION_PLUGIN_H__

#include <v3d_interface.h>

#include <QFileDialog>  // Add this include
#include <QtGui>
#include <fstream>  // Add this for file output

#include "cellSegmentation_plugin.h"

// A basic structure for a marker
struct MyMarker {
  float x, y, z;
  float radius;
};

// A basic 4D image container
struct MIND_4DImage {
  unsigned char *data;
  V3DLONG xdim, ydim, zdim, cdim;

  MIND_4DImage();
  MIND_4DImage(const MIND_4DImage &other);
  MIND_4DImage &operator=(const MIND_4DImage &other);
  ~MIND_4DImage();
};

// Plugin class
class SomaSegmentation : public QObject, public V3DPluginInterface2_1 {
  Q_OBJECT
  Q_INTERFACES(V3DPluginInterface2_1);
  Q_PLUGIN_METADATA(IID "com.janelia.v3d.V3DPluginInterface/2.1")

 public:
  float getPluginVersion() const { return 1.2f; }

  QStringList menulist() const;
  void domenu(const QString &menu_name, V3DPluginCallback2 &callback,
              QWidget *parent);

  QStringList funclist() const;
  bool dofunc(const QString &func_name, const V3DPluginArgList &input,
              V3DPluginArgList &output, V3DPluginCallback2 &callback,
              QWidget *parent);
};

// Main reconstruction function (invoked from menu or command-line)
struct input_PARA {
  QString inimg_file;
  V3DLONG channel;
};
MIND_4DImage *reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                                  input_PARA &PARA, bool bmenu);

// Isotropic correction function

void isotropic_correction_func(V3DPluginCallback2 &callback, QWidget *parent,
                               input_PARA &PARA, bool bmenu);

// PCA Analysis

void pca_func(V3DPluginCallback2 &callback, QWidget *parent, input_PARA &PARA,
              bool bmenu);

void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex);

void savePCAResultsToCSV(
    const QString &filename, int somaIndex, const LocationSimple &lm,
    double pc1, double pc2, double pc3, const double *vec1, const double *vec2,
    const double *vec3, double x_center, double y_center, double z_center,
    QWidget *parent = nullptr,
    bool *saveEnabled = nullptr);  // Add save flag parameter

void drawLine(Image4DSimple *image, double *from, double *to);

void visualizePCA_func(V3DPluginCallback2 &callback, QWidget *parent);

#endif  // __SOMA_SEGMENTATION_PLUGIN_H__
