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
  return QStringList() << tr("isotropic_correction") << tr("soma_segmentation")
                       << tr("pc_analysis") << tr("Visualize PCA")
                       << tr("Simulate Somas") << tr("Create Background")
                       << tr("simulate_soma_data") << tr("about");
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
  if (menu_name == tr("soma_segmentation")) {
    bool bmenu = true;
    input_PARA PARA;
    cellSegmentation cellseg;
    cellseg.interface_run(callback, parent);
  } else if (menu_name == tr("isotropic_correction")) {
    bool bmenu = true;
    input_PARA PARA;
    isotropic_correction_func(callback, parent, PARA, bmenu);
  } else if (menu_name == tr("pc_analysis")) {
    bool bmenu = true;
    input_PARA PARA;
    pca_func(callback, parent, PARA, bmenu);
  } else if (menu_name == tr("Visualize PCA")) {
    visualizePCA_func(callback, parent);
  } else if (menu_name == tr("Simulate Somas")) {
    simulate_somas(callback, parent);
  } else if (menu_name == tr("Create Background")) {
    create_background(callback, parent);
  } else if (menu_name == tr("simulate_soma_data")) {
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
  QString savePath = imageName + "_pca.csv";
  for (int i = 0; i < landmarkList.size(); i++) {
    analyzeSomaPCA(data1d, N, M, P, landmarkList[i], i + 1, savePath);
  }
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

  // Construct segmentation filename (try different options)
  QStringList possibleSegFiles;
  possibleSegFiles << imageName + "_seg.tif"  // Original approach
                   << currentImagePath + "/" + baseImageName +
                          "_seg.tif"               // Full path + basename
                   << baseImageName + "_seg.tif";  // Just basename

  QString segFileName;
  bool foundSegFile = false;

  for (int i = 0; i < possibleSegFiles.size(); i++) {
    if (QFile::exists(possibleSegFiles[i])) {
      segFileName = possibleSegFiles[i];
      foundSegFile = true;
      printf("Found segmentation file: %s\n",
             segFileName.toStdString().c_str());
      break;
    }
  }

  // If segmentation file still not found, ask the user to select it
  if (!foundSegFile) {
    v3d_msg("No segmentation file found. Please segment the image first.",
            parent);
    return;
  }

  // Load the binary segmentation file
  unsigned char *segData = nullptr;
  V3DLONG sz[4];
  int datatype = 0;
  if (!simple_loadimage_wrapper(callback, segFileName.toStdString().c_str(),
                                segData, sz, datatype)) {
    v3d_msg("Failed to load segmentation file.", parent);
    return;
  }

  /*
   * Load segmentation image PCA and get distribution of soma properties
   */

  // Find the segmentation image PCA file
  QString segPcaFileName = imageName + "_seg_pca.csv";

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
  std::vector<double> pcValues;      // eigenvalues
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

    pcValues.push_back(row[8]);   // eigenvalue1
    pcValues.push_back(row[9]);   // eigenvalue2
    pcValues.push_back(row[10]);  // eigenvalue3

    for (int iVec = 11; iVec < 20; iVec++) {
      eigenVectors.push_back(row[iVec]);
    }

    pcaRowCount++;
  }

  printf("Loaded %d soma segmentation image PCA records\n", pcaRowCount);

  // Calculate mean and standard deviation of center of mass coordinates,
  // eigenvalues and eigenvectors
  std::vector<double> meanCenter(3, 0.0);
  std::vector<double> stdCenter(3, 0.0);
  std::vector<double> meanEigenvalues(3, 0.0);
  std::vector<double> stdEigenvalues(3, 0.0);
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

  for (size_t i = 0; i < pcValues.size(); i += 3) {
    for (int j = 0; j < 3; j++) {
      meanEigenvalues[j] += pcValues[i + j];
    }
  }

  for (size_t i = 0; i < eigenVectors.size(); i += 9) {
    for (int j = 0; j < 9; j++) {
      meanEigenvectors[j] += eigenVectors[i + j];
    }
  }

  for (int j = 0; j < 3; j++) {
    meanEigenvalues[j] /= numSomas;
  }

  for (int j = 0; j < 9; j++) {
    meanEigenvectors[j] /= numSomas;
  }

  // Calculate standard deviations
  for (size_t i = 0; i < pcValues.size(); i += 3) {
    for (int j = 0; j < 3; j++) {
      stdEigenvalues[j] += pow(pcValues[i + j] - meanEigenvalues[j], 2);
    }
  }

  for (size_t i = 0; i < eigenVectors.size(); i += 9) {
    for (int j = 0; j < 9; j++) {
      stdEigenvectors[j] += pow(eigenVectors[i + j] - meanEigenvectors[j], 2);
    }
  }

  for (int j = 0; j < 3; j++) {
    stdEigenvalues[j] = sqrt(stdEigenvalues[j] / numSomas);
  }

  for (int j = 0; j < 9; j++) {
    stdEigenvectors[j] = sqrt(stdEigenvectors[j] / numSomas);
  }

  printf("\nPCA Statistics:\n");
  printf("Mean center: (%.2f, %.2f, %.2f)\n", meanCenter[0], meanCenter[1],
         meanCenter[2]);
  printf("Std dev: (%.2f, %.2f, %.2f)\n", stdCenter[0], stdCenter[1],
         stdCenter[2]);
  printf("Mean eigenvalues: (%.2f, %.2f, %.2f)\n", meanEigenvalues[0],
         meanEigenvalues[1], meanEigenvalues[2]);
  printf("Std dev eigenvalues: (%.2f, %.2f, %.2f)\n\n", stdEigenvalues[0],
         stdEigenvalues[1], stdEigenvalues[2]);

  /*
   * Create synthetic somas and image
   */

  // Create output image
  V3DLONG totalSize = xDim * yDim * zDim;
  unsigned char *outSegData = new unsigned char[totalSize];
  memset(outSegData, 0, totalSize);

  // Create output image with original intensity values
  unsigned char *outIntensityData = new unsigned char[totalSize];
  memset(outIntensityData, 0, totalSize);

  // Generate random positions and place synthetic somas
  std::random_device rd;
  std::mt19937 gen(rd());

  // Number of synthetic somas to generate
  int numSynthetic = 20;

  v3d_msg(QString("Generating %1 synthetic somas...").arg(numSynthetic));

  int successfulPlacements = 0;

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

      for (int j = 0; j < 3; j++) {
        std::normal_distribution<> d(meanCenter[j], stdCenter[j]);
        newCenter[j] = d(gen);

        if (j == 0 && (newCenter[j] < boundaryMargin ||
                       newCenter[j] > xDim - boundaryMargin)) {
          validPosition = false;
        } else if (j == 1 && (newCenter[j] < boundaryMargin ||
                              newCenter[j] > yDim - boundaryMargin)) {
          validPosition = false;
        } else if (j == 2 && (newCenter[j] < boundaryMargin ||
                              newCenter[j] > zDim - boundaryMargin)) {
          validPosition = false;
        }
      }
    }

    if (!validPosition) {
      printf("Failed to find valid position for soma %d after %d attempts\n",
             i + 1, maxAttempts);
      continue;
    }

    printf(
        "Placed soma %d/%d at (%.1f, %.1f, %.1f) with radius %.2f and cube "
        "size %ld\n",
        i + 1, numSynthetic, newCenter[0], newCenter[1], newCenter[2], radius,
        cubeSize);
    successfulPlacements++;

    // Generate random PCA values based on the distribution
    double randomPC1, randomPC2, randomPC3;
    double randomVec1[3], randomVec2[3], randomVec3[3];

    // Generate eigenvalues with normal distribution
    std::normal_distribution<> d1(meanEigenvalues[0], stdEigenvalues[0]);
    std::normal_distribution<> d2(meanEigenvalues[1], stdEigenvalues[1]);
    std::normal_distribution<> d3(meanEigenvalues[2], stdEigenvalues[2]);
    randomPC1 = d1(gen);
    randomPC2 = d2(gen);
    randomPC3 = d3(gen);

    // Ensure eigenvalues are positive and in descending order
    randomPC1 = std::max(randomPC1, 0.1);
    randomPC2 = std::max(std::min(randomPC2, randomPC1 - 0.1), 0.1);
    randomPC3 = std::max(std::min(randomPC3, randomPC2 - 0.1), 0.1);

    // Generate eigenvectors with normal distributions
    for (int j = 0; j < 3; j++) {
      std::normal_distribution<> dv1(meanEigenvectors[j], stdEigenvectors[j]);
      std::normal_distribution<> dv2(meanEigenvectors[j + 3],
                                     stdEigenvectors[j + 3]);
      std::normal_distribution<> dv3(meanEigenvectors[j + 6],
                                     stdEigenvectors[j + 6]);
      randomVec1[j] = dv1(gen);
      randomVec2[j] = dv2(gen);
      randomVec3[j] = dv3(gen);
    }

    // Normalize eigenvectors
    double norm1 =
        sqrt(randomVec1[0] * randomVec1[0] + randomVec1[1] * randomVec1[1] +
             randomVec1[2] * randomVec1[2]);
    double norm2 =
        sqrt(randomVec2[0] * randomVec2[0] + randomVec2[1] * randomVec2[1] +
             randomVec2[2] * randomVec2[2]);
    double norm3 =
        sqrt(randomVec3[0] * randomVec3[0] + randomVec3[1] * randomVec3[1] +
             randomVec3[2] * randomVec3[2]);

    for (int j = 0; j < 3; j++) {
      randomVec1[j] /= norm1;
      randomVec2[j] /= norm2;
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
    int *tempSegmentation = new int[totalVoxels];
    unsigned char *tempIntensity = new unsigned char[totalVoxels];
    memset(tempSegmentation, 0, totalVoxels * sizeof(int));
    memset(tempIntensity, 0, totalVoxels * sizeof(unsigned char));

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
              tempIntensity[targetIdx] = originalData[sourceIdx];
            }
          }
        }
      }
    }

    // Apply random rotation based on PCA values
    // cellSegmentation::class_segmentationMain segMain;
    // segMain.rotateSegmentation(tempSegmentation, cubeSize, randomPC1,
    // randomPC2,
    //                            randomPC3, randomVec1, randomVec2,
    //                            randomVec3);

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
                  tempIntensity[sourceIdx];  // Original intensity
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
  QString outSegFileName = imageName + "_simulated_seg.tif";
  QString outIntensityFileName = imageName + "_simulated_intensity.tif";

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
                  "somas.\nSaved images as %3 and %4.")
              .arg(successfulPlacements)
              .arg(numSynthetic)
              .arg(outSegFileName)
              .arg(outIntensityFileName),
          parent);

  delete[] segData;

  /*
   * Open new windows and display the synthetic soma data
   */

  // Create and show new window with binary simulated data
  Image4DSimple outSegImage;
  outSegImage.setData(outSegData, out_sz[0], out_sz[1], out_sz[2], out_sz[3],
                      V3D_UINT8);
  v3dhandle segWin = callback.newImageWindow();
  callback.setImage(segWin, &outSegImage);
  callback.setImageName(segWin, outSegFileName);
  callback.updateImageWindow(segWin);

  // Create and show new window with intensity simulated data
  Image4DSimple outIntensityImage;
  outIntensityImage.setData(outIntensityData, out_sz[0], out_sz[1], out_sz[2],
                            out_sz[3], V3D_UINT8);
  v3dhandle intensityWin = callback.newImageWindow();
  callback.setImage(intensityWin, &outIntensityImage);
  callback.setImageName(intensityWin, outIntensityFileName);
  callback.updateImageWindow(intensityWin);
}
