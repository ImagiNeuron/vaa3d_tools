/* demo_opensurfacefile_plugin.h
 * demo for open a surface file
 * 2015-2-10 : by Hanchuan Peng
 */
 
#ifndef __DEMO_OPENSURFACEFILE_PLUGIN_H__
#define __DEMO_OPENSURFACEFILE_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>

class DemoOpenSurfaceFilePlugin : public QObject, public V3DPluginInterface2_1
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
};

#endif

