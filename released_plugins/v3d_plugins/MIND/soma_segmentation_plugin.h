#ifndef __SOMA_SEGMENTATION_PLUGIN_H__
#define __SOMA_SEGMENTATION_PLUGIN_H__

#include <v3d_interface.h>

#include <QtGui>

#include "cellSegmentation_plugin.h"  // use cellSegmentation algorithm

class soma_segmentation : public QObject, public V3DPluginInterface2_1 {
  Q_OBJECT
  Q_INTERFACES(V3DPluginInterface2_1)
  Q_PLUGIN_METADATA(IID "com.janelia.v3d.V3DPluginInterface/2.1")

 public:
  float getPluginVersion() const { return 1.0f; }
  QStringList menulist() const { return QStringList() << "Soma Segmentation"; }
  QStringList funclist() const { return QStringList() << "somasegmentation"; }

  void domenu(const QString &menu_name, V3DPluginCallback2 &callback,
              QWidget *parent);
  bool dofunc(const QString &func_name, const V3DPluginArgList &input,
              V3DPluginArgList &output, V3DPluginCallback2 &callback,
              QWidget *parent);
};

#endif
