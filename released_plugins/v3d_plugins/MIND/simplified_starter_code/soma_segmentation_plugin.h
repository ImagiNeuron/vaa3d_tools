/* soma_segmentation_plugin.h
 *
 * Simplified plugin header for obtaining a 3D image and its landmarks.
 */

#ifndef __SOMA_SEGMENTATION_PLUGIN_H__
#define __SOMA_SEGMENTATION_PLUGIN_H__

#include <v3d_interface.h>

#include <QtGui>

class SomaSegmentation : public QObject, public V3DPluginInterface2_1 {
  Q_OBJECT
  Q_INTERFACES(V3DPluginInterface2_1);
  Q_PLUGIN_METADATA(IID "com.janelia.v3d.V3DPluginInterface/2.1")

 public:
  // Returns the plugin version.
  float getPluginVersion() const { return 1.0f; }

  // Returns the list of menu items.
  QStringList menulist() const;

  // Called when a menu item is chosen.
  void domenu(const QString &menu_name, V3DPluginCallback2 &callback,
              QWidget *parent);

  // Returns the list of functions.
  QStringList funclist() const;

  // Called when a function is invoked.
  bool dofunc(const QString &func_name, const V3DPluginArgList &input,
              V3DPluginArgList &output, V3DPluginCallback2 &callback,
              QWidget *parent);
};

#endif
