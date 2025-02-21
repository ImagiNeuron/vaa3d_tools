/** soma_segmentation_plugin.cpp
 * This plugin supports the analaysis and segmentaiton of neuron somas in the
 * brain.
 *
 * It includes several functions:
 * - isotropic correction - to view 3D imagery of the brain isotropically
 * - soma segmentation - to segment individual somas using a user defined 3D
 * region-growing algorithm, and output a binary segmented image
 *
 * 2024-11-16: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"

#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <queue>
#include <vector>

#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

// Main reconstruction function (invoked from menu or command-line)
My4DImage *reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                               input_PARA &PARA, bool bmenu) {
  unsigned char *data1d = 0;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];
  if (bmenu) {
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin) {
      QMessageBox::information(
          0, "", "You don't have any image open in the main window.");
      return nullptr;
    }
    Image4DSimple *p4DImage = callback.getImage(curwin);
    if (!p4DImage) {
      QMessageBox::information(0, "",
                               "The image pointer is invalid. Ensure your data "
                               "is valid and try again!");
      return nullptr;
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
    } else
      c = QInputDialog::getInt(parent, "Channel", "Enter channel NO:", 1, 1, sc,
                               1, &ok1);
    if (!ok1) return nullptr;
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
      fprintf(stderr, "Error reading file [%s].\n",
              PARA.inimg_file.toStdString().c_str());
      return nullptr;
    }
    if (PARA.channel < 1 || PARA.channel > in_sz[3]) {
      fprintf(stderr, "Invalid channel number.\n");
      return nullptr;
    }
    N = in_sz[0];
    M = in_sz[1];
    P = in_sz[2];
    sc = in_sz[3];
    c = PARA.channel;
  }

  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image window is currently open.", bmenu);
    return nullptr;
  }
  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("Invalid image pointer.", bmenu);
    return nullptr;
  }
  LandmarkList landmarkList = callback.getLandmark(curwin);
  if (landmarkList.isEmpty()) {
    v3d_msg("No landmarks defined. Please define at least one landmark.",
            bmenu);
    return nullptr;
  }

  // Convert the LandmarkList to a QList of MyMarker objects.
  QList<MyMarker> myMarkerList;
  for (const auto &lm : landmarkList) {
    MyMarker marker;
    marker.x = lm.x;
    marker.y = lm.y;
    marker.z = lm.z;
    marker.radius = lm.radius;
    myMarkerList.append(marker);
  }
}

////////////////////////////////////////////////////////////////////////
// Standard Vaa3D plugin interface implementations

QStringList SomaSegmentation::menulist() const {
  return QStringList() << tr("soma_segmentation") << tr("about");
}

QStringList SomaSegmentation::funclist() const {
  return QStringList() << tr("segment_somas") << tr("help");
}

void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  if (menu_name == tr("soma_segmentation")) {
    cellSegmentation cellseg;
    cellseg.interface_run(callback, parent);
  } else {
    v3d_msg(tr("This plugin segments individual somas using a 3D "
               "region-growing algorithm "
               "with median filtering, and one of three possible tresholding "
               "methods: local or global otsu thresholding, as well as "
               "iterative tresholding."
               "All methods use landmarks as seeds. Developed by ImagiNeuron "
               "(2025-02-20)"),
            0);
  }
}

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
      fprintf(stderr, "Need input image.\n");
      return false;
    } else
      PARA.inimg_file = infiles[0];
    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;
    reconstruction_func(callback, parent, PARA, bmenu);
  } else if (func_name == tr("help")) {
    printf("**** Usage of soma_segmentation ****\n");
    printf(
        "vaa3d -x soma_segmentation -f segment_somas -i <inimg_file> -p "
        "<channel>\n");
    printf("inimg_file       The input image\n");
    printf(
        "channel          Data channel for processing (starting from 1, "
        "default 1)\n");
    printf(
        "The output segmented image is binary (white soma, black "
        "background).\n");
  } else
    return false;
  return true;
}

////////////////////////////////////////////////////////////////////////
// Implementation of My4DImage methods

My4DImage::My4DImage() : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0) {}

My4DImage::My4DImage(const My4DImage &other) {
  xdim = other.xdim;
  ydim = other.ydim;
  zdim = other.zdim;
  cdim = other.cdim;
  if (other.data) {
    data = new unsigned char[xdim * ydim * zdim * cdim];
    std::copy(other.data, other.data + xdim * ydim * zdim * cdim, data);
  } else {
    data = nullptr;
  }
}

My4DImage &My4DImage::operator=(const My4DImage &other) {
  if (this == &other) return *this;
  delete[] data;
  xdim = other.xdim;
  ydim = other.ydim;
  zdim = other.zdim;
  cdim = other.cdim;
  if (other.data) {
    data = new unsigned char[xdim * ydim * zdim * cdim];
    std::copy(other.data, other.data + xdim * ydim * zdim * cdim, data);
  } else {
    data = nullptr;
  }
  return *this;
}

My4DImage::~My4DImage() { delete[] data; }
