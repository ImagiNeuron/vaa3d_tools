/** soma_segmentation_plugin.cpp
 * This plugin supports the analaysis and segmentaiton of neuron somas in the
 * brain.
 *
 * It includes several functions:
 * - isotropic correction - to view 3D imagery of the brain isotropically
 * - soma segmentation - to segment individual somas using a user defined 3D
 * region-growing algorithm, and output a binary segmented image
 * - pc_analysis - to perform Principal Component Analysis on the input image
 * and output the results
 *
 * 2024-11-16: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"

#include <QApplication>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
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
#include "cellSegmentation_plugin.h"
#include "v3d_message.h"

using namespace std;

/**
 * @brief Function to reconstruct somas
 *
 * Default function - copy format and add functionality afterwards
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

std::tuple<char, char, char> colormap(double value) {
  // Ensure value is between 0 and 1
  value = std::clamp(value, 0.0, 1.0);

  // Approximate viridis through polynomial fits
  // These are simplified approximations of the actual colormap
  double r = 0.267004 + value * (0.004974 + value * (0.9981 + value * -0.9988));
  double g = 0.004974 + value * (0.9186 + value * (0.1718 + value * -0.4439));
  double b = 0.329415 + value * (0.0875 + value * (-0.1234 + value * 0.1794));

  r = std::clamp(r, 0.0, 1.0);
  g = std::clamp(g, 0.0, 1.0);
  b = std::clamp(b, 0.0, 1.0);

  return {static_cast<char>(r * 255), static_cast<char>(g * 255),
          static_cast<char>(b * 255)};
}

void visualizeProbabilityModel_func(V3DPluginCallback2 &callback,
                                    QWidget *parent) {
  // Load data
  QString filename = QFileDialog::getOpenFileName(
      parent, "Open Probability Model", "", "Binary Files (*.bin)");

  if (filename.isEmpty()) {
    printf("No file selected.\n");
    return;
  }

  std::vector<double> data;
  V3DLONG dim_X, dim_Y, dim_Z;
  cellSegmentation::class_segmentationMain::loadProbabilityModel(
      filename.toStdString().c_str(), data, dim_X, dim_Y, dim_Z);

  Image4DSimple *p4DImage = new Image4DSimple();
  p4DImage->createBlankImage(dim_X, dim_Y, dim_Z, 3, V3D_UINT8);

  double min = data[0];
  double max = data[0];

  for (V3DLONG i = 0; i < dim_X * dim_Y * dim_Z; i++) {
    if (data[i] < min) min = data[i];
    if (data[i] > max) max = data[i];
  }

  unsigned char *pixels = p4DImage->getRawData();
  int channelSize = dim_X * dim_Y * dim_Z;
  for (V3DLONG i = 0; i < dim_X * dim_Y * dim_Z; i++) {
    auto [r, g, b] = colormap((data[i] - min) / (max - min));
    pixels[i] = r;
    pixels[i + channelSize] = g;
    pixels[i + 2 * channelSize] = b;
  }

  v3dhandle newwin = callback.newImageWindow("Probability Model");
  callback.setImage(newwin, p4DImage);
}

/**
 * @brief overlay the ground truth data on the new simulated image. Dimesnions
 * calculated based on current image
 */
void overlaySimulation(V3DPluginCallback2 &callback, QWidget *parent,
                       unsigned char *binarySegImage,
                       unsigned char *gradientImage,
                       unsigned char *simulatedImage) {
  v3dhandle curwin = callback.currentImageWindow();
  v3dhandle newwin = callback.newImageWindow();

  Image4DSimple *p4DImage = callback.getImage(curwin);
  unsigned char *newData = new unsigned char[p4DImage->getTotalBytes() * 3];

  memcpy(newData, simulatedImage, p4DImage->getTotalBytes());
  memcpy(newData + p4DImage->getTotalBytes(), binarySegImage,
         p4DImage->getTotalBytes());
  memcpy(newData + p4DImage->getTotalBytes() * 2, gradientImage,
         p4DImage->getTotalBytes());

  Image4DSimple *newImage = new Image4DSimple;
  newImage->setData(newData, p4DImage->getXDim(), p4DImage->getYDim(),
                    p4DImage->getZDim(), p4DImage->getCDim() * 3, V3D_UINT8);

  callback.setImage(newwin, newImage);
}

/**
 * @brief Function to simulate synthetic soma data
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param PARA - the input parameters
 * @param bmenu - whether the function is being called from the menu
 */

void simulate_soma_data(V3DPluginCallback2 &callback, QWidget *parent,
                        input_PARA &PARA, bool bmenu) {
  // Get current window and validate
  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image window open!");
    return;
  }

  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("Invalid image pointer!");
    return;
  }

  // Get dimensions of current image
  V3DLONG xDim = p4DImage->getXDim();
  V3DLONG yDim = p4DImage->getYDim();
  V3DLONG zDim = p4DImage->getZDim();

  printf("\nStarting soma simulation...\n");
  printf("Image dimensions: X=%ld, Y=%ld, Z=%ld\n", xDim, yDim, zDim);

  // Get the current image name and path
  QString imageName = callback.getImageName(curwin);
  QString currentImagePath = QFileInfo(imageName).absolutePath();
  QString baseImageName = QFileInfo(imageName).baseName();

  /*
   * Load original image data
   */

  unsigned char *originalData = p4DImage->getRawData();
  int channel = 0;  // Default to first channel (index 0)

  /*
   * Load soma segmentation data
   */

  unsigned char *segData = nullptr;
  V3DLONG sz[4];
  int datatype = 0;
  loadSegmentationFile(imageName, segData, sz, datatype, callback, parent);

  /*
   * Load segmentation image PCA and get distribution of soma properties
   */

  // Find the segmentation image PCA file
  QString segPcaFileName =
      modifyFilePathForTeraFly(imageName) + "_pca_binary_segmentation.csv";

  // Check if the PCA file exists
  if (!QFile::exists(segPcaFileName)) {
    v3d_msg(
        QString("segmentation image PCA file not found: %1\nPlease run soma "
                "segmentation first.")
            .arg(segPcaFileName),
        parent);
    delete[] segData;
    return;
  }

  printf("\nLoading segmentation image PCA data from: %s\n",
         segPcaFileName.toStdString().c_str());

  // Load PCA data from CSV file
  std::vector<double> eigenVectors;  // 9 eigenvector components
  std::vector<double> centerCoords;  // CenterMassX, CenterMassY, CenterMassZ
  std::vector<double> markerCoords;  // X, Y, Z
  std::vector<double> somaRadii;     // Soma radii

  std::ifstream segPcaFile(segPcaFileName.toStdString().c_str());
  if (!segPcaFile.is_open()) {
    v3d_msg("Could not open segmentation image PCA file!");
    delete[] segData;
    return;
  }

  std::string line;
  // Skip header line
  std::getline(segPcaFile, line);

  // Read PCA data
  int pcaRowCount = 0;
  while (std::getline(segPcaFile, line)) {
    std::stringstream ss(line);
    std::string value;
    std::vector<double> row;

    while (std::getline(ss, value, ',')) {
      row.push_back(std::stod(value));
    }

    // Expected columns in each CSV row:
    //  0: SomaID
    //  1: X
    //  2: Y
    //  3: Z
    //  4: Radius
    //  5: CenterMassX
    //  6: CenterMassY
    //  7: CenterMassZ
    //  8: eigenvalue1
    //  9: eigenvalue2
    // 10: eigenvalue3
    // 11: eigenvector1_x
    // 12: eigenvector1_y
    // 13: eigenvector1_z
    // 14: eigenvector2_x
    // 15: eigenvector2_y
    // 16: eigenvector2_z
    // 17: eigenvector3_x
    // 18: eigenvector3_y
    // 19: eigenvector3_z

    markerCoords.push_back(row[1]);  // X
    markerCoords.push_back(row[2]);  // Y
    markerCoords.push_back(row[3]);  // Z

    somaRadii.push_back(row[4]);  // Radius

    centerCoords.push_back(row[5]);  // CenterMassX
    centerCoords.push_back(row[6]);  // CenterMassY
    centerCoords.push_back(row[7]);  // CenterMassZ

    for (int iVec = 11; iVec < 20; iVec++) {
      eigenVectors.push_back(row[iVec]);
    }

    pcaRowCount++;
  }

  printf("Loaded %d soma segmentation image PCA records\n", pcaRowCount);

  // Calculate mean and standard deviation of center of mass coordinates,
  // and eigenvectors
  std::vector<double> meanCenter(3, 0.0);
  std::vector<double> stdCenter(3, 0.0);
  std::vector<double> meanEigenvectors(9, 0.0);
  std::vector<double> stdEigenvectors(9, 0.0);

  for (size_t i = 0; i < centerCoords.size(); i += 3) {
    for (int j = 0; j < 3; j++) {
      meanCenter[j] += centerCoords[i + j];
    }
  }

  int numSomas = centerCoords.size() / 3;
  for (int j = 0; j < 3; j++) {
    meanCenter[j] /= numSomas;
  }

  for (size_t i = 0; i < centerCoords.size(); i += 3) {
    for (int j = 0; j < 3; j++) {
      stdCenter[j] += pow(centerCoords[i + j] - meanCenter[j], 2);
    }
  }

  for (int j = 0; j < 3; j++) {
    stdCenter[j] = sqrt(stdCenter[j] / numSomas);
  }

  for (size_t i = 0; i < eigenVectors.size(); i += 9) {
    for (int j = 0; j < 9; j++) {
      meanEigenvectors[j] += eigenVectors[i + j];
    }
  }

  for (int j = 0; j < 9; j++) {
    meanEigenvectors[j] /= numSomas;
  }

  // Calculate standard deviations
  for (size_t i = 0; i < eigenVectors.size(); i += 9) {
    for (int j = 0; j < 9; j++) {
      stdEigenvectors[j] += pow(eigenVectors[i + j] - meanEigenvectors[j], 2);
    }
  }

  for (int j = 0; j < 9; j++) {
    stdEigenvectors[j] = sqrt(stdEigenvectors[j] / numSomas);
  }

  printf("\nPCA Statistics:\n");
  printf("Mean center: (%.2f, %.2f, %.2f)\n", meanCenter[0], meanCenter[1],
         meanCenter[2]);
  printf("Std dev: (%.2f, %.2f, %.2f)\n", stdCenter[0], stdCenter[1],
         stdCenter[2]);

  /*
   * Create synthetic somas and image
   */

  // Create output image
  V3DLONG totalSize = xDim * yDim * zDim;
  unsigned char *outSegData = new unsigned char[totalSize];
  memset(outSegData, 0, totalSize);

  // Create output image with original intensity values
  unsigned char *outIntensityData =
      create_background(callback, parent, xDim, yDim, zDim);

  // Generate random positions and place synthetic somas
  std::random_device rd;
  std::mt19937 gen(rd());

  // Ask user for the number of synthetic somas to generate
  bool ok;
  int numSynthetic =
      QInputDialog::getInt(parent, "Synthetic Soma Generation",
                           "Enter the number of synthetic somas to generate:",
                           20,    // Default value
                           1,     // Minimum value
                           1000,  // Maximum value
                           1,     // Step
                           &ok);

  if (!ok) {
    // User canceled the dialog
    v3d_msg("Synthetic soma generation canceled.", parent);
    delete[] segData;
    return;
  }

  v3d_msg(QString("Generating %1 synthetic somas...").arg(numSynthetic));

  int successfulPlacements = 0;

  // Structure to keep track of placed somas for overlap detection
  struct PlacedSoma {
    double x, y, z;  // Center coordinates
    double radius;   // Soma radius
  };
  std::vector<PlacedSoma> placedSomas;

  // Helper function to check if two somas overlap
  auto somasOverlap = [](const PlacedSoma &s1, const PlacedSoma &s2) -> bool {
    // Calculate squared distance between centers
    double dx = s1.x - s2.x;
    double dy = s1.y - s2.y;
    double dz = s1.z - s2.z;
    double distSq = dx * dx + dy * dy + dz * dz;

    // If distance is less than sum of radii, they overlap
    double minDist = s1.radius + s2.radius;
    return distSq < (minDist * minDist);
  };

  for (int i = 0; i < numSynthetic; i++) {
    // Choose a random soma from the available ones for this synthetic soma
    int randomSomaIndex = gen() % numSomas;

    // Get the radius and calculate appropriate cube size
    double radius = somaRadii[randomSomaIndex];
    V3DLONG cubeSize = static_cast<V3DLONG>(
        2.5 * radius);  // Use 2.5x radius to ensure we capture the whole soma

    // Make sure cubeSize is odd for centering purposes
    if (cubeSize % 2 == 0) cubeSize += 1;

    // Set boundary margin based on the cube size
    int boundaryMargin = cubeSize / 2;

    // Generate random position using normal distribution
    std::vector<double> newCenter(3);
    bool validPosition = false;
    int maxAttempts = 100;
    int attempts = 0;

    while (!validPosition && attempts < maxAttempts) {
      attempts++;
      validPosition = true;

      // Generate new potential position
      for (int j = 0; j < 3; j++) {
        std::normal_distribution<> d(meanCenter[j], stdCenter[j]);
        newCenter[j] = d(gen);

        // Check if position is within image boundaries
        if (j == 0 && (newCenter[j] < boundaryMargin ||
                       newCenter[j] > xDim - boundaryMargin)) {
          validPosition = false;
          break;
        } else if (j == 1 && (newCenter[j] < boundaryMargin ||
                              newCenter[j] > yDim - boundaryMargin)) {
          validPosition = false;
          break;
        } else if (j == 2 && (newCenter[j] < boundaryMargin ||
                              newCenter[j] > zDim - boundaryMargin)) {
          validPosition = false;
          break;
        }
      }

      // If position is within boundaries, check for overlap with existing somas
      if (validPosition) {
        PlacedSoma newSoma = {newCenter[0], newCenter[1], newCenter[2], radius};

        // Check against all previously placed somas
        for (const auto &existingSoma : placedSomas) {
          if (somasOverlap(newSoma, existingSoma)) {
            validPosition = false;
            printf(
                "Soma %d position attempt %d: Overlap detected with existing "
                "soma\n",
                i + 1, attempts);
            break;
          }
        }
      }
    }

    if (!validPosition) {
      printf(
          "Failed to find valid non-overlapping position for soma %d after %d "
          "attempts\n",
          i + 1, maxAttempts);
      continue;
    }

    // Add this soma to our placed somas list for future overlap checking
    placedSomas.push_back({newCenter[0], newCenter[1], newCenter[2], radius});

    printf(
        "Placed soma %d/%d at (%.1f, %.1f, %.1f) with radius %.2f and cube "
        "size %ld\n",
        i + 1, numSynthetic, newCenter[0], newCenter[1], newCenter[2], radius,
        cubeSize);
    successfulPlacements++;

    // Generate random PCA values based on the distribution
    double randomVec1[3], randomVec2[3], randomVec3[3];

    // Generate eigenvectors with normal distributions
    for (int j = 0; j < 3; j++) {
      std::normal_distribution<> dv1(meanEigenvectors[j], stdEigenvectors[j]);
      randomVec1[j] = dv1(gen);
    }

    // Normalize first vector
    double norm1 =
        sqrt(randomVec1[0] * randomVec1[0] + randomVec1[1] * randomVec1[1] +
             randomVec1[2] * randomVec1[2]);
    for (int j = 0; j < 3; j++) {
      randomVec1[j] /= norm1;
    }

    // Make second vector orthogonal to the first using Gram-Schmidt process
    bool validSecondVector = false;
    int maxRetries = 10;  // Prevent infinite loops
    int retryCount = 0;
    double norm2 = 0.0;

    while (!validSecondVector && retryCount < maxRetries) {
      // Generate eigenvectors with normal distributions
      for (int j = 0; j < 3; j++) {
        std::normal_distribution<> dv2(meanEigenvectors[j + 3],
                                       stdEigenvectors[j + 3]);
        randomVec2[j] = dv2(gen);
      }

      // Project randomVec2 onto randomVec1
      double dot_product1 = randomVec2[0] * randomVec1[0] +
                            randomVec2[1] * randomVec1[1] +
                            randomVec2[2] * randomVec1[2];

      // Subtract the projection from randomVec2
      for (int j = 0; j < 3; j++) {
        randomVec2[j] -= dot_product1 * randomVec1[j];
      }

      // Compute the norm of the resulting vector
      norm2 =
          sqrt(randomVec2[0] * randomVec2[0] + randomVec2[1] * randomVec2[1] +
               randomVec2[2] * randomVec2[2]);

      // Check if the resulting vector is not too small
      if (norm2 >= 1e-6) {
        validSecondVector = true;
      }

      retryCount++;
    }

    // If still no valid second vector after retries, generate an arbitrary
    // vector perpendicular to randomVec1
    if (!validSecondVector) {
      if (fabs(randomVec1[0]) < fabs(randomVec1[1]) &&
          fabs(randomVec1[0]) < fabs(randomVec1[2])) {
        randomVec2[0] = 1.0;
        randomVec2[1] = 0.0;
        randomVec2[2] = 0.0;
      } else if (fabs(randomVec1[1]) < fabs(randomVec1[2])) {
        randomVec2[0] = 0.0;
        randomVec2[1] = 1.0;
        randomVec2[2] = 0.0;
      } else {
        randomVec2[0] = 0.0;
        randomVec2[1] = 0.0;
        randomVec2[2] = 1.0;
      }

      // Make it orthogonal to randomVec1
      double dot_product1 = randomVec2[0] * randomVec1[0] +
                            randomVec2[1] * randomVec1[1] +
                            randomVec2[2] * randomVec1[2];

      for (int j = 0; j < 3; j++) {
        randomVec2[j] -= dot_product1 * randomVec1[j];
      }

      // Compute the norm of the resulting vector
      norm2 =
          sqrt(randomVec2[0] * randomVec2[0] + randomVec2[1] * randomVec2[1] +
               randomVec2[2] * randomVec2[2]);
    }

    // Normalize second vector
    for (int j = 0; j < 3; j++) {
      randomVec2[j] /= norm2;
    }

    // Make third vector orthogonal to first two using Gram-Schmidt
    bool validThirdVector = false;
    retryCount = 0;
    double norm3 = 0.0;

    while (!validThirdVector && retryCount < maxRetries) {
      // Generate eigenvectors with normal distributions
      for (int j = 0; j < 3; j++) {
        std::normal_distribution<> dv3(meanEigenvectors[j + 6],
                                       stdEigenvectors[j + 6]);
        randomVec3[j] = dv3(gen);
      }

      // Project randomVec3 onto randomVec1
      double dot_product3_1 = randomVec3[0] * randomVec1[0] +
                              randomVec3[1] * randomVec1[1] +
                              randomVec3[2] * randomVec1[2];

      // Project randomVec3 onto randomVec2
      double dot_product3_2 = randomVec3[0] * randomVec2[0] +
                              randomVec3[1] * randomVec2[1] +
                              randomVec3[2] * randomVec2[2];

      // Subtract both projections to make it orthogonal to both vectors
      for (int j = 0; j < 3; j++) {
        randomVec3[j] = randomVec3[j] - dot_product3_1 * randomVec1[j] -
                        dot_product3_2 * randomVec2[j];
      }

      // Normalize third vector
      norm3 =
          sqrt(randomVec3[0] * randomVec3[0] + randomVec3[1] * randomVec3[1] +
               randomVec3[2] * randomVec3[2]);

      // Check if the resulting vector is not too small
      if (norm3 >= 1e-6) {
        validThirdVector = true;
      }

      retryCount++;
    }

    // If still no valid third vector after retries, generate with cross product
    // method
    if (!validThirdVector) {
      randomVec3[0] =
          randomVec1[1] * randomVec2[2] - randomVec1[2] * randomVec2[1];
      randomVec3[1] =
          randomVec1[2] * randomVec2[0] - randomVec1[0] * randomVec2[2];
      randomVec3[2] =
          randomVec1[0] * randomVec2[1] - randomVec1[1] * randomVec2[0];
      norm3 =
          sqrt(randomVec3[0] * randomVec3[0] + randomVec3[1] * randomVec3[1] +
               randomVec3[2] * randomVec3[2]);
    }

    for (int j = 0; j < 3; j++) {
      randomVec3[j] /= norm3;
    }

    // Extract a soma from segmentation data
    // Get center of mass for the selected soma
    V3DLONG sourceCenterX =
        static_cast<V3DLONG>(centerCoords[randomSomaIndex * 3]);
    V3DLONG sourceCenterY =
        static_cast<V3DLONG>(centerCoords[randomSomaIndex * 3 + 1]);
    V3DLONG sourceCenterZ =
        static_cast<V3DLONG>(centerCoords[randomSomaIndex * 3 + 2]);

    V3DLONG totalVoxels = cubeSize * cubeSize * cubeSize;
    double *tempSegmentation = new double[totalVoxels];
    double *tempIntensity = new double[totalVoxels];
    memset(tempSegmentation, 0, totalVoxels * sizeof(double));
    memset(tempIntensity, 0, totalVoxels * sizeof(double));

    // Extract both segmentation and intensity data for the soma
    for (int z = 0; z < cubeSize; z++) {
      for (int y = 0; y < cubeSize; y++) {
        for (int x = 0; x < cubeSize; x++) {
          // Calculate positions relative to the soma's center
          V3DLONG sourceX = sourceCenterX + x - cubeSize / 2;
          V3DLONG sourceY = sourceCenterY + y - cubeSize / 2;
          V3DLONG sourceZ = sourceCenterZ + z - cubeSize / 2;

          // Target index in temporary buffer
          int targetIdx = z * cubeSize * cubeSize + y * cubeSize + x;

          // Check if coordinates are within the image bounds
          if (sourceX >= 0 && sourceX < xDim && sourceY >= 0 &&
              sourceY < yDim && sourceZ >= 0 && sourceZ < zDim) {
            // Calculate index in the source images
            V3DLONG sourceIdx =
                sourceZ * xDim * yDim + sourceY * xDim + sourceX;

            // Copy the segmentation value (0 or 255 for binary image)
            tempSegmentation[targetIdx] = segData[sourceIdx] > 0 ? 1 : 0;

            // Copy the original intensity value if this voxel is part of the
            // soma
            if (segData[sourceIdx] > 0) {
              tempIntensity[targetIdx] =
                  static_cast<double>(originalData[sourceIdx]);
            }
          }
        }
      }
    }

    // Apply random rotation to the synthetic soma
    cellSegmentation::class_segmentationMain segMain;
    segMain.rotateSegmentation(tempSegmentation, cubeSize, randomVec1,
                               randomVec2, randomVec3);
    segMain.rotateSegmentation(tempIntensity, cubeSize, randomVec1, randomVec2,
                               randomVec3);

    // Place rotated synthetic soma at generated position
    int centerX = static_cast<int>(newCenter[0]);
    int centerY = static_cast<int>(newCenter[1]);
    int centerZ = static_cast<int>(newCenter[2]);

    // Copy rotated soma to both output images
    for (int z = 0; z < cubeSize; z++) {
      for (int y = 0; y < cubeSize; y++) {
        for (int x = 0; x < cubeSize; x++) {
          int sourceIdx = z * cubeSize * cubeSize + y * cubeSize + x;

          // Convert to output image coordinates
          int targetX = centerX + x - cubeSize / 2;
          int targetY = centerY + y - cubeSize / 2;
          int targetZ = centerZ + z - cubeSize / 2;

          if (targetX >= 0 && targetX < xDim && targetY >= 0 &&
              targetY < yDim && targetZ >= 0 && targetZ < zDim) {
            V3DLONG targetIdx =
                targetZ * xDim * yDim + targetY * xDim + targetX;

            // Only set voxel if rotated model indicates soma presence
            if (tempSegmentation[sourceIdx] > 0) {
              outSegData[targetIdx] = 255;  // Binary segmentation
              outIntensityData[targetIdx] =
                  static_cast<unsigned char>(std::round(std::min(
                      255.0,
                      std::max(
                          0.0,
                          tempIntensity[sourceIdx]))));  // Original intensity
            }
          }
        }
      }
    }

    delete[] tempSegmentation;
    delete[] tempIntensity;
  }

  /*
   * Save the synthetic soma segmentation images
   */
  QString outSegFileName =
      modifyFilePathForTeraFly(imageName) + "_simulated_segmentation.tif";
  QString outIntensityFileName =
      modifyFilePathForTeraFly(imageName) + "_simulated_intensity.tif";

  // Create dimension array for saving images
  V3DLONG out_sz[4];
  out_sz[0] = xDim;
  out_sz[1] = yDim;
  out_sz[2] = zDim;
  out_sz[3] = 1;  // Single channel

  // Save the binary segmentation image
  simple_saveimage_wrapper(callback, outSegFileName.toStdString().c_str(),
                           outSegData, out_sz, V3D_UINT8);

  // Save the intensity image
  simple_saveimage_wrapper(callback, outIntensityFileName.toStdString().c_str(),
                           outIntensityData, out_sz, V3D_UINT8);

  printf("\nSimulation complete:\n");
  printf("Successfully placed %d/%d somas\n", successfulPlacements,
         numSynthetic);
  v3d_msg(QString("Simulation complete. Generated %1/%2 synthetic "
                  "somas. Saved images as %3 and %4.")
              .arg(successfulPlacements)
              .arg(numSynthetic)
              .arg(outSegFileName)
              .arg(outIntensityFileName),
          parent);

  delete[] segData;

  /*
   * Open new windows and display the synthetic soma data
   */

  unsigned char *gradientImage = new unsigned char[totalSize];
  cellSegmentation cellSeg;
  cellSeg.sobel3D(outSegData, gradientImage, xDim, yDim,
                                          zDim);

  // overlay
  overlaySimulation(callback, parent, outSegData, gradientImage,
                    outIntensityData);
  // // Create and show new window with binary simulated data. Now Obsolote
  // Image4DSimple outSegImage;
  // outSegImage.setData(outSegData, out_sz[0], out_sz[1], out_sz[2], out_sz[3],
  //                     V3D_UINT8);
  // v3dhandle segWin = callback.newImageWindow();
  // callback.setImage(segWin, &outSegImage);
  // callback.setImageName(segWin, outSegFileName);
  // callback.updateImageWindow(segWin);

  // // Create and show new window with intensity simulated data
  // Image4DSimple outIntensityImage;
  // outIntensityImage.setData(outIntensityData, out_sz[0], out_sz[1],
  // out_sz[2],
  //                           out_sz[3], V3D_UINT8);
  // v3dhandle intensityWin = callback.newImageWindow();
  // callback.setImage(intensityWin, &outIntensityImage);
  // callback.setImageName(intensityWin, outIntensityFileName);
  // callback.updateImageWindow(intensityWin);
}
