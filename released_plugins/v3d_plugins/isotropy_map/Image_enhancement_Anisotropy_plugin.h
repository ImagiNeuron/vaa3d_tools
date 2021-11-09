/* Image_enhancement_Anisotropy_plugin.h
 * This is a test plugin, you can use it as a demo.
 * 2013-01-28 : by Zhi
 */
 
#ifndef __IMAGE_ENHANCEMENT_ANISOTROPY_PLUGIN_H__
#define __IMAGE_ENHANCEMENT_ANISOTROPY_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>

class anisotropy_enhancement : public QObject, public V3DPluginInterface2_1
{
	Q_OBJECT
	Q_INTERFACES(V3DPluginInterface2_1);
    Q_PLUGIN_METADATA(IID"com.janelia.v3d.V3DPluginInterface/2.1")

public:
    float getPluginVersion() const {return 1.2f;}

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const ;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);
};

#endif

