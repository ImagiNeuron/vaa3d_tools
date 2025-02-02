/* soma_segmentation_plugin.h
 * A plugin for soma segmentation using adaptive region growing with custom
 * smoothing.
 */

#ifndef __SOMA_SEGMENTATION_PLUGIN_H__
#define __SOMA_SEGMENTATION_PLUGIN_H__

#include <v3d_interface.h>

#include <QtGui>

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

// Plugin class declaration
class SomaSegmentation : public QObject, public V3DPluginInterface2_1 {
  Q_OBJECT
  Q_INTERFACES(V3DPluginInterface2_1);
  Q_PLUGIN_METADATA(IID "com.janelia.v3d.V3DPluginInterface/2.1")

 public:
  // Version number for this plugin
  float getPluginVersion() const { return 1.3f; }

  // Standard Vaa3D plugin interface methods
  QStringList menulist() const;
  void domenu(const QString &menu_name, V3DPluginCallback2 &callback,
              QWidget *parent);

  QStringList funclist() const;
  bool dofunc(const QString &func_name, const V3DPluginArgList &input,
              V3DPluginArgList &output, V3DPluginCallback2 &callback,
              QWidget *parent);
};

/***********************************
 * Declarations of helper functions
 ***********************************/

// (Internal functions for custom smoothing and region growing are defined in
// the cpp file) For example:
//   My4DImage *smoothMy4DImage(const My4DImage *input, int ITER);
//   My4DImage *performFloodFill(const My4DImage *image, const MyMarker &marker,
//   int unused,
//                               int &min_x, int &max_x, int &min_y, int &max_y,
//                               int &min_z, int &max_z);
//   My4DImage *mergeFloodedImages(const QList<My4DImage *> &floodedImages);

// Legacy function declarations (if needed, can be removed or commented out):
/*
void applyMedianFilter(const unsigned char *inputData,
                       unsigned char *outputData,
                       V3DLONG N, V3DLONG M, V3DLONG P,
                       int windowSize);

void applyGaussianFilter(const unsigned char *inputData,
                         unsigned char *outputData,
                         V3DLONG N, V3DLONG M, V3DLONG P,
                         float sigma);

void extractSubvolume(const unsigned char *inData,
                      unsigned char *outData,
                      V3DLONG N, V3DLONG M, V3DLONG P,
                      int x1, int x2,
                      int y1, int y2,
                      int z1, int z2);

void applyWatershedVS(const unsigned char *subvol,
                      unsigned short *labelOut,
                      int sx, int sy, int sz);
*/

// You may also declare any other helper functions if they need to be visible
// outside the cpp file.

#endif  // __SOMA_SEGMENTATION_PLUGIN_H__
