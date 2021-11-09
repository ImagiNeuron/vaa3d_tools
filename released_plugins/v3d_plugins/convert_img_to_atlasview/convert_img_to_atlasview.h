/* convert_img_to_atlasview.h
 * extract the channels of an image and convert it as an atlas-viewer vieweable format
 * 2011-02-28: by Hanchuan Peng
 */


#ifndef __CONVERT_IMG_TO_ATLAS_H__
#define __CONVERT_IMG_TO_ATLAS_H__

#include <QtGui>

#include "v3d_interface.h"

class ConvertImg2AtlasPlugin : public QObject, public V3DPluginInterface2_1
{
    Q_OBJECT
    Q_INTERFACES(V3DPluginInterface2_1);
    Q_PLUGIN_METADATA(IID"com.janelia.v3d.V3DPluginInterface/2.1")


public:
    float getPluginVersion() const {return 1.3f;}

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);
};



#endif

