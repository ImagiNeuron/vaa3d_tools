/**
 * This is the main file of the MIND plugin that supports the analaysis and
 * segmentaiton of neuron somas in the brain. It contains the logic that
 * controls which functions are called when the user selects a plugin option.
 *
 * This file also contains the logic for isotropic correction, and intensity
 * Weighted PCA on an image
 *
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"

#include <QApplication>
#include <QFileDialog>
#include <QFont>
#include <QInputDialog>
#include <QMessageBox>
#include <QObject>
#include <QPainter>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <queue>
#include <vector>

#include "ResolutionDialog.h"
#include "basic_4dimage.h"
#include "basic_surf_objs.h"
#include "mainwindow.h"
#include "v3d_message.h"

using namespace std;

/**
 * @brief Function to reconstruct somas
 *
 * Default function from plugin generator - copy format and add functionality
 * afterwards
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param PARA - the input parameters
 * @param bmenu - whether the function is being called from the menu
 */
MIND_4DImage *reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
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
  // Put code here
}

/**************************************
 * Plugin Interface Methods
 **************************************/
QStringList SomaSegmentation::menulist() const {
  return QStringList() << tr("Isotropic Correction") << tr("Soma Segmentation")
                       << tr("PC Analysis") << tr("Visualize PCA")
                       << tr("Visualize Probability Model")
                       << tr("Probability Model Legend")
                       << tr("Create Background") << tr("Simulate Somas")
                       << tr("about");
}

QStringList SomaSegmentation::funclist() const {
  return QStringList() << tr("isotropic_correction") << tr("segment_somas")
                       << tr("simulate_soma_data") << tr("pc_analysis")
                       << tr("help");
}

/**
 * @brief Call the appropriate methods when the menue items are clicked
 *
 * @param menu_name - the name of the menu being described
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 */
void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  if (menu_name == tr("Soma Segmentation")) {
    bool bmenu = true;
    input_PARA PARA;
    cellSegmentation cellseg;
    cellseg.interface_run(callback, parent);
  } else if (menu_name == tr("Isotropic Correction")) {
    bool bmenu = true;
    input_PARA PARA;
    isotropic_correction_func(callback, parent, PARA, bmenu);
  } else if (menu_name == tr("PC Analysis")) {
    bool bmenu = true;
    input_PARA PARA;
    pca_func(callback, parent, PARA, bmenu);
  } else if (menu_name == tr("Visualize PCA")) {
    visualizePCA_func(callback, parent);
  } else if (menu_name == tr("Visualize Probability Model")) {
    visualizeProbabilityModel_func(callback, parent);
  } else if (menu_name == tr("Probability Model Legend")) {
    probabilityModelLegend_func(callback, parent);
  } else if (menu_name == tr("Create Background")) {
    create_background(callback, parent);
  } else if (menu_name == tr("Simulate Somas")) {
    bool bmenu = true;
    input_PARA PARA;
    simulate_soma_data(callback, parent, PARA, bmenu);
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

/**
 * @brief Call the appropriate methods when the function items are called from
 * the command line
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
      fprintf(stderr, "Need input image.\n");
      return false;
    } else {
      PARA.inimg_file = infiles[0];
    }

    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;
    cellSegmentation cellseg;
    cellseg.interface_run(callback, parent);
  } else if (func_name == tr("pc_analysis")) {
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
    } else {
      PARA.inimg_file = infiles[0];
    }

    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;

    pca_func(callback, parent, PARA, bmenu);
  } else if (func_name == tr("isotropic_correction")) {
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
    } else {
      PARA.inimg_file = infiles[0];
    }

    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;

    isotropic_correction_func(callback, parent, PARA, bmenu);
  } else if (func_name == tr("simulate_soma_data")) {
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
    } else {
      PARA.inimg_file = infiles[0];
    }

    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;

    simulate_soma_data(callback, parent, PARA, bmenu);
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
  } else {
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////
// Implementation of MIND_4DImage methods

MIND_4DImage::MIND_4DImage()
    : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0) {}

MIND_4DImage::MIND_4DImage(const MIND_4DImage &other) {
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

MIND_4DImage &MIND_4DImage::operator=(const MIND_4DImage &other) {
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

MIND_4DImage::~MIND_4DImage() { delete[] data; }

/**
 * @brief Function to correct isotropic resolution of an image
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param PARA - the input parameters
 * @param bmenu - whether the function is being called from the menu
 */
void isotropic_correction_func(V3DPluginCallback2 &callback, QWidget *parent,
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

  //// ISOTROPIC CORRECTION CODE GOES HERE

  // get current window, image
  v3dhandle curwin = callback.currentImageWindow();
  Image4DSimple *p4DImage = callback.getImage(curwin);

  // Ask for desired resolution of a image pixel along the 3
  // axes for isotropic correction and set resolution of the image
  ResolutionDialog dialog(parent);
  dialog.setResolutionOfImage(p4DImage);

  return;
}

/**************************************
 * PC Analysis for each soma
 **************************************/
void pca_func(V3DPluginCallback2 &callback, QWidget *parent, input_PARA &PARA,
              bool bmenu) {
  /*************************************
   * Load Image
   *************************************/
  unsigned char *data1d = 0;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];

  if (bmenu) {
    // From the current Vaa3D window
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin) {
      QMessageBox::information(0, "", "No image open in the main window.");
      return;
    }

    Image4DSimple *p4DImage = callback.getImage(curwin);
    if (!p4DImage) {
      QMessageBox::information(0, "", "Invalid image pointer.");
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
    // Command-line input
    int datatype = 0;
    if (!simple_loadimage_wrapper(callback,
                                  PARA.inimg_file.toStdString().c_str(), data1d,
                                  in_sz, datatype)) {
      fprintf(stderr, "Error loading file [%s].\n",
              PARA.inimg_file.toStdString().c_str());
      return;
    }
    if (PARA.channel < 1 || PARA.channel > in_sz[3]) {
      fprintf(stderr, "Invalid channel number.\n");
      return;
    }
    N = in_sz[0];
    M = in_sz[1];
    P = in_sz[2];
    sc = in_sz[3];
    c = PARA.channel;
  }

  if (!data1d) {
    v3d_msg("No valid image data!", bmenu);
    return;
  }

  if (sc > 1) {
    v3d_msg("For grayscale, ignoring additional channels");
    sc = 1;
    c = 1;
  }

  /*************************************
   * For each landmark, do PCA
   *************************************/
  v3dhandle curwin = callback.currentImageWindow();
  QString imageName = callback.getImageName(curwin);

  LandmarkList landmarkList;
  if (bmenu && curwin) {
    landmarkList = callback.getLandmark(curwin);
  } else {
    // In command-line or no open window, user must supply landmarks some other
    // way For now, if none, we just exit.
  }

  if (landmarkList.isEmpty()) {
    v3d_msg("No landmarks found. Please specify at least one landmark.", bmenu);
    delete[] data1d;
    return;
  }
  QString savePath =
      modifyFilePathForTeraFly(imageName) + "_pca_intensity_weighted.csv";
  for (int i = 0; i < landmarkList.size(); i++) {
    analyzeSomaPCA(data1d, N, M, P, landmarkList[i], i + 1, savePath);
  }
}
