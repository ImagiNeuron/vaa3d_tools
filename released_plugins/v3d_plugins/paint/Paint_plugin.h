/* Paint_plugin.h
 * This is a paint toolbox
 * 2015-02-04 : by Yujie Li
 */
 
#ifndef __PAINT_PLUGIN_H__
#define __PAINT_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>
#include <QImage>
#include <QPoint>
#include <QWidget>

class paint : public QObject, public V3DPluginInterface2_1
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


private:
    QString fileName;
    V3DPluginCallback2 * callback;



};


#endif

