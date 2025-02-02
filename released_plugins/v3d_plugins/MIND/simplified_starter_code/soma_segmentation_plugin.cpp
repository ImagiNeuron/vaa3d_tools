/* soma_segmentation_plugin.cpp
 *
 * This simplified plugin demonstrates how to obtain a 3D image and its
 * associated landmarks from the current image window. It then prints debugging
 * information (image dimensions and landmark details) to the console.
 *
 * 2025-02-01: Edited to only keep the image/landmark retrieval and debug code.
 */

#include "soma_segmentation_plugin.h"

#include <QDebug>  // For debug output
#include <QInputDialog>
#include <QMessageBox>

#include "basic_surf_objs.h"  // Contains Image4DSimple and LandmarkList definitions
#include "v3d_message.h"

using namespace std;

struct input_PARA {
  QString inimg_file;
  V3DLONG channel;
};

// This function obtains the current 3D image and the list of landmarks,
// then prints debugging information.
void demo_obtainImageAndLandmarks(V3DPluginCallback2 &callback, QWidget *parent,
                                  input_PARA &PARA, bool bmenu) {
  // Obtain the current image window
  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image window is currently open.", bmenu);
    return;
  }

  // Retrieve the image pointer
  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("Invalid image pointer.", bmenu);
    return;
  }

  // Debug: Print image dimensions and channel count
  V3DLONG xDim = p4DImage->getXDim();
  V3DLONG yDim = p4DImage->getYDim();
  V3DLONG zDim = p4DImage->getZDim();
  V3DLONG cDim = p4DImage->getCDim();
  qDebug() << "Image dimensions:" << xDim << "x" << yDim << "x" << zDim
           << "with" << cDim << "channel(s).";

  // Retrieve the landmarks
  LandmarkList landmarkList = callback.getLandmark(curwin);
  if (landmarkList.isEmpty()) {
    v3d_msg("No landmarks defined. Please define at least one landmark.",
            bmenu);
    return;
  }

  qDebug() << "Number of landmarks:" << landmarkList.size();

  // Debug: Print details of each landmark
  for (int i = 0; i < landmarkList.size(); ++i) {
    LocationSimple lm = landmarkList.at(i);
    qDebug() << "Landmark" << i << ": x =" << lm.x << ", y =" << lm.y
             << ", z =" << lm.z << ", radius =" << lm.radius;
  }

  v3d_msg(
      "Debug: 3D image and landmarks obtained. Check console output for "
      "details.",
      bmenu);
}

/////////////////////////////////////////////////////////////////////////
// Plugin interface implementations

QStringList SomaSegmentation::menulist() const {
  return QStringList() << tr("obtainImageAndLandmarks") << tr("about");
}

void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  if (menu_name == tr("obtainImageAndLandmarks")) {
    bool bmenu = true;
    input_PARA PARA;
    demo_obtainImageAndLandmarks(callback, parent, PARA, bmenu);
  } else {
    v3d_msg(
        "This plugin obtains a 3D image and its landmarks, printing debug "
        "info.",
        0);
  }
}

QStringList SomaSegmentation::funclist() const {
  // No command-line functions implemented in this simplified version.
  return QStringList();
}

bool SomaSegmentation::dofunc(const QString &func_name,
                              const V3DPluginArgList &input,
                              V3DPluginArgList &output,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  // For simplicity, this example only supports the menu interface.
  return false;
}
