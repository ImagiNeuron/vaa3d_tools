/* soma_segmentation_plugin.h
 * This plugin allows the segmentation of individual somas with 3D flooding algorithms.
 * It generates ground truth data to train machine learning algorithms for automatic soma segmnentation in the brain
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and Thibaut Baguette
 */

#ifndef __SOMA_SEGMENTATION_PLUGIN_H__
#define __SOMA_SEGMENTATION_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>

// Define MyMarker structure
struct MyMarker
{
	float x, y, z;
	float radius;
};

// Modify My4DImage to prevent shallow copies
struct My4DImage
{
	unsigned char *data;
	V3DLONG xdim, ydim, zdim, cdim;

	My4DImage();
	My4DImage(const My4DImage &other);
	My4DImage &operator=(const My4DImage &other);
	~My4DImage();
};

// Corrected Flood Fill function declaration
My4DImage *performFloodFill(const My4DImage *image, const MyMarker &marker, int threshold, int &min_x, int &max_x, int &min_y, int &max_y, int &min_z, int &max_z);
My4DImage *mergeFloodedImages(const QList<My4DImage *> &floodedImages);
bool saveFloodedImage(const My4DImage *floodedImage, const QString &filename);

class SomaSegmentation : public QObject, public V3DPluginInterface2_1
{
	Q_OBJECT
	Q_INTERFACES(V3DPluginInterface2_1);
	Q_PLUGIN_METADATA(IID "com.janelia.v3d.V3DPluginInterface/2.1")

public:
	float getPluginVersion() const { return 1.1f; }

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);
};

#endif
