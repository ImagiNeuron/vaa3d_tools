/* soma_segmentation_plugin.h
 * A plugin for soma segmentation using 3D watershed (Vincent & Soille style).
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette
 */

#ifndef __SOMA_SEGMENTATION_PLUGIN_H__
#define __SOMA_SEGMENTATION_PLUGIN_H__

#include <v3d_interface.h>

#include <QtGui>

#include "cellSegmentation_plugin.h"

// A basic structure for a marker
struct MyMarker {
  float x, y, z;
  float radius;
};

// A basic 4D image container
struct My4DImage {
  unsigned char *data;
  V3DLONG xdim, ydim, zdim, cdim;

  My4DImage();
  My4DImage(const My4DImage &other);
  My4DImage &operator=(const My4DImage &other);
  ~My4DImage();
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

/***********************************
 * Declaration of Helper Functions
 **********************************/

// Median filter
void applyMedianFilter(const unsigned char *inputData,
                       unsigned char *outputData, V3DLONG N, V3DLONG M,
                       V3DLONG P, int windowSize);

// Gaussian filter
void applyGaussianFilter(const unsigned char *inputData,
                         unsigned char *outputData, V3DLONG N, V3DLONG M,
                         V3DLONG P, float sigma);

// Subvolume extraction
void extractSubvolume(const unsigned char *inData, unsigned char *outData,
                      V3DLONG N, V3DLONG M, V3DLONG P, int x1, int x2, int y1,
                      int y2, int z1, int z2);

// Watershed (Vincent & Soille)
void applyWatershedVS(const unsigned char *subvol, unsigned short *labelOut,
                      int sx, int sy, int sz);

// (Minimal) BFS-based 3D distance transform for removing watershed lines
void dt3d_binary(const float *inData, V3DLONG *pix_index, const V3DLONG *sz,
                 float threshVal = 0);

// Simple 3D morphological opening (erosion + dilation) with a 3x3x3
// neighborhood
static void morphologicalOpen3D(unsigned char *vol, int sx, int sy, int sz);

// Compute Otsu threshold
static int computeOtsuThreshold(const unsigned char *data, int length);

// Main reconstruction function (invoked from menu or command-line)
struct input_PARA {
  QString inimg_file;
  V3DLONG channel;
};
My4DImage *reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                               input_PARA &PARA, bool bmenu);

#endif  // __SOMA_SEGMENTATION_PLUGIN_H__
