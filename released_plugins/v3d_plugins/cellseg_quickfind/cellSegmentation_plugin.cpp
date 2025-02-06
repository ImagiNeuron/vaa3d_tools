/* cellSegmentation.cpp
 * It aims to automatically segment cells;
 * 2014-10-12 : by Xiang Li (lindbergh.li@gmail.com);
 */

#include "cellSegmentation_plugin.h"

#include <basic_landmark.h>
#include <math.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "convert_type2uint8.h"
#include "sstream"
#include "string"
#include "v3d_message.h"

#pragma region "UI-related functions and function defination"
// Q_EXPORT_PLUGIN2(cellSegmentation, cellSegmentation);

QStringList cellSegmentation::menulist() const {
  return QStringList()
         //<<tr("Exemplar defination")
         << tr("Exemplar defination and propagation")
         //<<tr("Further segmentation")
         << tr("About");
}

QStringList cellSegmentation::funclist() const {
  return QStringList() << tr("cellsegmentation") << tr("help");
}

void cellSegmentation::domenu(
    const QString &menu_name,
    V3DPluginCallback2 &V3DPluginCallback2_currentCallback,
    QWidget *QWidget_parent) {
  if (menu_name == tr("Exemplar defination and propagation")) {
    cellSegmentation::interface_run(V3DPluginCallback2_currentCallback,
                                    QWidget_parent);
  } else {
    v3d_msg(
        tr("Cell segmentation and counting by exemplar matching\nby Xiang Li "
           "(lindbergh.li@gmail.com);"));
  }
}

bool cellSegmentation::dofunc(
    const QString &func_name, const V3DPluginArgList &input,
    V3DPluginArgList &output,
    V3DPluginCallback2 &V3DPluginCallback2_currentCallback,
    QWidget *QWidget_parent) {
  vector<char *> infiles, inparas, outfiles;
  if (input.size() >= 1) infiles = *((vector<char *> *)input.at(0).p);
  if (input.size() >= 2) inparas = *((vector<char *> *)input.at(1).p);
  if (output.size() >= 1) outfiles = *((vector<char *> *)output.at(0).p);

  if (func_name == tr("cellsegmentation")) {
  } else if (func_name == tr("help")) {
  } else {
    return false;
  }
  return true;
}
#pragma endregion
