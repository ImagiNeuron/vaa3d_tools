/* shape_extractor_plugin.h
 * A tool to extract cell shapes
 * 2015-2-17 : by Yujie Li
 */
 
#ifndef __SHAPE_EXTRACTOR_PLUGIN_H__
#define __SHAPE_EXTRACTOR_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>

class shape_extr_plugin : public QObject, public V3DPluginInterface2_1
{
	Q_OBJECT
	Q_INTERFACES(V3DPluginInterface2_1);
    Q_PLUGIN_METADATA(IID"com.janelia.v3d.V3DPluginInterface/2.1")


public:
	float getPluginVersion() const {return 1.1f;}

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const ;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);
    void shape_extract(V3DPluginCallback2 & callback, const V3DPluginArgList & input, V3DPluginArgList & output);


private:
    V3DLONG sz_img[4];
    unsigned char *image_data;
};

#endif

