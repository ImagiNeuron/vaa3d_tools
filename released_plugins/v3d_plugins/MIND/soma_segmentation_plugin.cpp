/* soma_segmentation_plugin.cpp
 * This plugin allows the segmentation of individual somas with 3D flooding
 * algorithms. It generates ground truth data to train machine learning
 * algorithms for automatic soma segmnentation in the brain
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"

#include <QInputDialog>
#include <QMessageBox>
#include <vector>

#include "ResolutionDialog.h"
#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

struct input_PARA {
  QString inimg_file;
  V3DLONG channel;
};

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu);

/**
 * @brief Menu option under the MIND plugins
 */
QStringList SomaSegmentation::menulist() const {
  return QStringList() << tr("soma_segmentation") << tr("about");
}

/**
 * @brief Function list for the soma segmentation plugin
 */
QStringList SomaSegmentation::funclist() const {
  return QStringList() << tr("segment_somas") << tr("help");
}

/**
 * @brief Display a description about the soma segmentation plugin and what it
 * does
 *
 * @param menu_name - the name of the menu being described
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 */
void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  if (menu_name == tr("soma_segmentation")) {
    bool bmenu = true;
    input_PARA PARA;
    reconstruction_func(callback, parent, PARA, bmenu);

  } else {
    v3d_msg(tr(
        "This plugin allows the segmentation of individual somas with 3D "
        "flooding algorithms. It generates ground truth data to train machine "
        "learning algorithms for automatic soma segmnentation in the brain. "
        "Developed by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous "
        "and Thibaut Baguette, 2024-11-16"));
  }
}

/**
 * @brief Function to segment somas
 *
 * Describe how function works here
 *
 * @param func_name - the name of the function
 * @param input - the input arguments
 * @param output - the output arguments
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @return true if the function was successful, false otherwise
 */
bool SomaSegmentation::dofunc(const QString &func_name,
                              const V3DPluginArgList &input,
                              V3DPluginArgList &output,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  if (func_name == tr("segment_somas")) {
    bool bmenu = false;
    input_PARA PARA;

    vector<char *> *pinfiles =
        (input.size() >= 1) ? (vector<char *> *)input[0].p : 0;
    vector<char *> *pparas =
        (input.size() >= 2) ? (vector<char *> *)input[1].p : 0;
    vector<char *> infiles = (pinfiles != 0) ? *pinfiles : vector<char *>();
    vector<char *> paras = (pparas != 0) ? *pparas : vector<char *>();

    if (infiles.empty()) {
      fprintf(stderr, "Need input image. \n");
      return false;
    } else
      PARA.inimg_file = infiles[0];
    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;
    reconstruction_func(callback, parent, PARA, bmenu);
  } else if (func_name == tr("help")) {
    ////HERE IS WHERE THE DEVELOPERS SHOULD UPDATE THE USAGE OF THE PLUGIN

    printf("**** Usage of soma_segmentation tracing **** \n");
    printf(
        "vaa3d -x soma_segmentation -f segment_somas -i <inimg_file> -p "
        "<channel> <other parameters>\n");
    printf("inimg_file       The input image\n");
    printf(
        "channel          Data channel for tracing. Start from 1 (default "
        "1).\n");

    printf(
        "outswc_file      Will be named automatically based on the input image "
        "file name, so you don't have to specify it.\n\n");

  } else
    return false;

  return true;
}

/**
 * @brief Function to reconstruct somas
 *
 * Describe how function works here
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param PARA - the input parameters
 * @param bmenu - whether the function is being called from the menu
 */
void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu) {
  unsigned char *data1d = 0;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];
  if (bmenu) {
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin) {
      QMessageBox::information(
          0, "", "You don't have any image open in the main window.");
      return;
    }

    Image4DSimple *p4DImage = callback.getImage(curwin);

    if (!p4DImage) {
      QMessageBox::information(0, "",
                               "The image pointer is invalid. Ensure your data "
                               "is valid and try again!");
      return;
    }

    data1d = p4DImage->getRawData();
    N = p4DImage->getXDim();
    M = p4DImage->getYDim();
    P = p4DImage->getZDim();
    sc = p4DImage->getCDim();

    bool ok1;

    if (sc == 1) {
      c = 1;
      ok1 = true;
    } else {
      c = QInputDialog::getInt(parent, "Channel", "Enter channel NO:", 1, 1, sc,
                               1, &ok1);
    }

    if (!ok1) return;

    in_sz[0] = N;
    in_sz[1] = M;
    in_sz[2] = P;
    in_sz[3] = sc;

    PARA.inimg_file = p4DImage->getFileName();
  } else {
    int datatype = 0;
    if (!simple_loadimage_wrapper(callback,
                                  PARA.inimg_file.toStdString().c_str(), data1d,
                                  in_sz, datatype)) {
      fprintf(stderr,
              "Error happens in reading the subject file [%s]. Exit. \n",
              PARA.inimg_file.toStdString().c_str());
      return;
    }
    if (PARA.channel < 1 || PARA.channel > in_sz[3]) {
      fprintf(stderr, "Invalid channel number. \n");
      return;
    }
    N = in_sz[0];
    M = in_sz[1];
    P = in_sz[2];
    sc = in_sz[3];
    c = PARA.channel;
  }

  //// THIS IS WHERE THE DEVELOPERS SHOULD ADD THEIR OWN NEURON TRACING CODE

  // get current window, image, and landmarks
  v3dhandle curwin = callback.currentImageWindow();
  Image4DSimple *p4DImage = callback.getImage(curwin);
  LandmarkList landmarkList = callback.getLandmark(curwin);

  // check if image is valid
  if (!p4DImage) {
    qDebug() << "error: invalid image";
    return;
  }

  // reset ROI
  ROIList roiList = callback.getROI(curwin);
  for (int j = 0; j < 3; j++) {
    roiList[j].clear();
  }

  // get 1st landmark
  if (landmarkList.size() < 1) {
    qDebug() << "error: no landmark found";
    return;
  }
  LocationSimple lm = landmarkList[0];
  float x, y, z;
  float radius;
  x = lm.x;
  y = lm.y;
  z = lm.z;
  radius = lm.radius;

  // Landmark info and ask for desired resolution of a image pixel along the 3
  // axes for isotropic correction and set resolution of the image
  ResolutionDialog dialog(x, y, z, radius, parent);
  dialog.setResolutionOfImage(p4DImage);

  // set ROI
  // ROIList being a QList<QPolygon>, and QPolygon being a QVector<QPoint>,
  // we represent the x, y, and z planes as polygons with 4 points each to
  // form a cube with the landmark at the center
  float x_min = x - radius * 2.0f;
  float x_max = x + radius * 2.0f;
  float y_min = y - radius * 2.0f;
  float y_max = y + radius * 2.0f;
  float z_min = z - radius * 2.0f;
  float z_max = z + radius * 2.0f;
  // x-y plane
  roiList[0] << QPoint(x_min, y_min);
  roiList[0] << QPoint(x_max, y_min);
  roiList[0] << QPoint(x_max, y_max);
  roiList[0] << QPoint(x_min, y_max);
  // z-y plane
  roiList[1] << QPoint(z_min, y_min);
  roiList[1] << QPoint(z_max, y_min);
  roiList[1] << QPoint(z_max, y_max);
  roiList[1] << QPoint(z_min, y_max);
  // x-z plane
  roiList[2] << QPoint(x_min, z_min);
  roiList[2] << QPoint(x_max, z_min);
  roiList[2] << QPoint(x_max, z_max);
  roiList[2] << QPoint(x_min, z_max);

  if (callback.setROI(curwin, roiList)) {
    callback.updateImageWindow(curwin);
  } else {
    qDebug() << "error: failed to set ROI";
    return;
  }

  callback.openROI3DWindow(curwin);

  // Update landmark
  View3DControl *v3dlocalcontrol = callback.getLocalView3DControl(curwin);
  if (v3dlocalcontrol) {
    v3dlocalcontrol->updateLandmark();
  } else {
    qDebug() << "error: failed to update 3D viewer";
    return;
  }

  return;
}
