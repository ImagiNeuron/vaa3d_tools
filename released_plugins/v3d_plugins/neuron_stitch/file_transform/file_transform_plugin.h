/* file_transform_plugin.h
 * This plugin will transform and combine image by given affine matrix. Some functions are based on littleQuickWarper.
 * 2015-3-18 : by Hanbo Chen
 */
 
#ifndef __FILE_TRANSFORM_PLUGIN_H__
#define __FILE_TRANSFORM_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>

class file_transform : public QObject, public V3DPluginInterface2_1
{
	Q_OBJECT
	Q_INTERFACES(V3DPluginInterface2_1);
    Q_PLUGIN_METADATA(IID"com.janelia.v3d.V3DPluginInterface/2.1")

public:
    float getPluginVersion() const {return 2.0f;}

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const ;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);

private:
    void doAffineImage(V3DPluginCallback2 & callback, QString fname_sub, QString fname_tar, QString fname_amat,
                       QString fname_output, int interpMethod, bool b_negativeShift, bool b_channelSeperate);
    int dotransform_swc(V3DPluginCallback2 &callback, QWidget *parent);
    int dotransform_marker(V3DPluginCallback2 &callback, QWidget *parent);
    void printHelp();
};


#endif

