/** soma_segmentation_plugin.cpp
 * This plugin supports the analaysis and segmentaiton of neuron somas in the
 * brain.
 *
 * It includes several functions:
 * - isotropic correction - to view 3D imagery of the brain isotropically
 * - soma segmentation - to segment individual somas using a user defined 3D
 * region-growing algorithm, and output a binary segmented image
 * - PCA analysis - to perform PCA analysis on the input image and output the
 * results
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
#include "compute_win_pca_wp.h"
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
                       << tr("PCA analysis") << tr("Visualize PCA")
                       << tr("about");
}

QStringList SomaSegmentation::funclist() const {
  return QStringList() << tr("isotropic_correction") << tr("segment_somas")
                       << tr("PCA analysis") << tr("help");
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
  } else if (menu_name == tr("PCA analysis")) {
    bool bmenu = true;
    input_PARA PARA;
    pca_func(callback, parent, PARA, bmenu);
  } else if (menu_name == tr("Visualize PCA")) {
    visualizePCA_func(callback, parent);
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
  } else if (func_name == tr("PCA analysis")) {
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
 * PCA analysis for each soma
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

  for (int i = 0; i < landmarkList.size(); i++) {
    analyzeSomaPCA(data1d, N, M, P, landmarkList[i], i + 1);
  }
}

void savePCAResultsToCSV(const QString &filename, int somaIndex,
                         const LocationSimple &lm, double pc1, double pc2,
                         double pc3, const double *vec1, const double *vec2,
                         const double *vec3, double x_center, double y_center,
                         double z_center, QWidget *parent, bool *saveEnabled) {
  static bool shouldSave = true;  // Default to true
  static QString actualFilename = filename;

  // Ask user if they want to save only for the first soma
  if (somaIndex == 1) {
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(
        parent, "Save Results",
        "Would you like to save the PCA results to a CSV file?",
        QMessageBox::Yes | QMessageBox::No);

    shouldSave = (reply == QMessageBox::Yes);
    if (saveEnabled) *saveEnabled = shouldSave;

    if (shouldSave) {
      QString suggestedName = QFileInfo(filename).fileName();
      actualFilename = QFileDialog::getSaveFileName(
          parent, "Save PCA Results", suggestedName, "CSV Files (*.csv)");
      if (actualFilename.isEmpty()) {
        printf("Save cancelled by user.\n");
        shouldSave = false;
        if (saveEnabled) *saveEnabled = false;
        return;
      }

      // Ensure it has .csv extension
      if (!actualFilename.endsWith(".csv", Qt::CaseInsensitive)) {
        actualFilename += ".csv";
      }

      bool fileExists = QFile::exists(actualFilename);
      if (fileExists) {
        // Remove and replace existing file
        if (QFile::remove(actualFilename)) {
          printf("Existing file removed: %s\n",
                 actualFilename.toStdString().c_str());
        } else {
          printf("Failed to remove existing file: %s\n",
                 actualFilename.toStdString().c_str());
          shouldSave = false;
          if (saveEnabled) *saveEnabled = false;
          return;
        }
      }

      // Write header if new file
      std::ofstream outFile(actualFilename.toStdString().c_str(),
                            std::ios::app);
      outFile << "SomaID,X,Y,Z,Radius,CenterMassX,CenterMassY,CenterMassZ,"
              << "eigenvalue1,eigenvalue2,eigenvalue3,"
              << "eigenvector1_x,eigenvector1_y,eigenvector1_z,"
              << "eigenvector2_x,eigenvector2_y,eigenvector2_z,"
              << "eigenvector3_x,eigenvector3_y,eigenvector3_z\n";
      outFile.close();
    } else {
      return;  // User chose not to save
    }
  }

  if (!shouldSave) return;  // Skip if user chose not to save

  std::ofstream outFile(actualFilename.toStdString().c_str(), std::ios::app);

  // Write data
  outFile << somaIndex << "," << lm.x << "," << lm.y << "," << lm.z << ","
          << lm.radius << "," << x_center << "," << y_center << "," << z_center
          << "," << pc1 << "," << pc2 << "," << pc3 << "," << vec1[0] << ","
          << vec1[1] << "," << vec1[2] << "," << vec2[0] << "," << vec2[1]
          << "," << vec2[2] << "," << vec3[0] << "," << vec3[1] << ","
          << vec3[2] << "\n";

  outFile.close();
}

void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex) {
  // Extract soma info
  float x = lm.x;
  float y = lm.y;
  float z = lm.z;
  float r = lm.radius > 0 ? lm.radius : 5.0f;

  // Create 3D array wrapper for the data
  unsigned char ***img3d = new unsigned char **[P];
  for (V3DLONG k = 0; k < P; k++) {
    img3d[k] = new unsigned char *[M];
    for (V3DLONG j = 0; j < M; j++) {
      img3d[k][j] = labeledData + (k * M * N + j * N);
    }
  }

  // Analyze each soma with appropriate window size
  double pc1, pc2, pc3;
  double vec1[3], vec2[3], vec3[3];
  double x_center, y_center, z_center;

  // Use sphere window type (1) and window size based on soma radius
  if (compute_sphere_win3d_pca_eigVec(
          img3d, N, M, P, x, y, z,  // Center on soma
          2 * r, 2 * r, 2 * r,      // Window size based on radius
          pc1, pc2, pc3, vec1, vec2, vec3, x_center, y_center, z_center)) {
    // Print results for this soma
    printf("\nSoma #%d PCA Results:\n", somaIndex);
    printf("  Center: (%.1f, %.1f, %.1f)\n", x, y, z);
    printf("  Radius: %.1f\n", r);
    printf("  Center of mass: (%f, %f, %f)\n", x_center, y_center, z_center);
    printf("  Eigenvalues:\n");
    printf("    pc1: %f\n", pc1);
    printf("    pc2: %f\n", pc2);
    printf("    pc3: %f\n", pc3);
    printf("  Principal axes:\n");
    printf("    pc1: [%f, %f, %f]\n", vec1[0], vec1[1], vec1[2]);
    printf("    pc2: [%f, %f, %f]\n", vec2[0], vec2[1], vec2[2]);
    printf("    pc3: [%f, %f, %f]\n\n\n", vec3[0], vec3[1], vec3[2]);

    // Save to CSV with save flag
    QString defaultFileName = "soma_pca_results.csv";
    QWidget *mainWin = QApplication::activeWindow();
    bool saveEnabled = false;
    savePCAResultsToCSV(defaultFileName, somaIndex, lm, pc1, pc2, pc3, vec1,
                        vec2, vec3, x_center, y_center, z_center, mainWin,
                        &saveEnabled);

    if (somaIndex == 1) {
      if (saveEnabled) {
        printf("PCA results will be saved to the selected file.\n");
      } else {
        printf("PCA results will not be saved to file.\n");
      }
    }
  } else {
    printf("\nSoma #%d PCA failed.\n", somaIndex);
  }

  // Cleanup 3D array wrapper
  for (V3DLONG k = 0; k < P; k++) {
    delete[] img3d[k];
  }
  delete[] img3d;
}

void visualizePCA_func(V3DPluginCallback2 &callback, QWidget *parent) {
  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image opened.", parent);
    return;
  }
  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("No image opened.", parent);
    return;
  }
  Image4DSimple *pcaVisualization = new Image4DSimple();
  pcaVisualization->createBlankImage(p4DImage->getXDim(), p4DImage->getYDim(),
                                     p4DImage->getZDim(), p4DImage->getCDim(),
                                     V3D_UINT8);
  pcaVisualization->setOriginX(p4DImage->getOriginX());
  pcaVisualization->setOriginY(p4DImage->getOriginY());
  pcaVisualization->setOriginZ(p4DImage->getOriginZ());
  pcaVisualization->setRezX(p4DImage->getRezX());
  pcaVisualization->setRezY(p4DImage->getRezY());
  pcaVisualization->setRezZ(p4DImage->getRezZ());

  // load pca data from csv
  QString filename = QFileDialog::getOpenFileName(parent, "Open PCA Results",
                                                  "", "CSV Files (*.csv)");
  if (filename.isEmpty()) {
    printf("No file selected.\n");
    return;
  }

  std::ifstream inFile(filename.toStdString().c_str());
  if (!inFile.is_open()) {
    printf("Failed to open file: %s\n", filename.toStdString().c_str());
    return;
  }

  // Skip header
  std::string line;
  std::getline(inFile, line);

  // Read data
  while (std::getline(inFile, line)) {
    std::istringstream ss(line);
    std::string token;

    int somaID;
    double center[3];
    double pc1, pc2, pc3;
    double vec1Pos[3], vec2Pos[3], vec3Pos[3];

    int col = 0;
    while (std::getline(ss, token, ',')) {
      switch (col) {
        case 0:
          somaID = std::stoi(token);
          break;
        case 1:
          center[0] = std::stod(token);
          break;
        case 2:
          center[1] = std::stod(token);
          break;
        case 3:
          center[2] = std::stod(token);
          break;
        // case 4: // radius
        // case 5,6,7: // center of mass
        case 8:
          pc1 = std::stod(token);
          break;
        case 9:
          pc2 = std::stod(token);
          break;
        case 10:
          pc3 = std::stod(token);
          break;
        case 11:
          vec1Pos[0] = pc1 * std::stod(token) + center[0];
          break;
        case 12:
          vec1Pos[1] = pc1 * std::stod(token) + center[1];
          break;
        case 13:
          vec1Pos[2] = pc1 * std::stod(token) + center[2];
          break;
        case 14:
          vec2Pos[0] = pc2 * std::stod(token) + center[0];
          break;
        case 15:
          vec2Pos[1] = pc2 * std::stod(token) + center[1];
          break;
        case 16:
          vec2Pos[2] = pc2 * std::stod(token) + center[2];
          break;
        case 17:
          vec3Pos[0] = pc3 * std::stod(token) + center[0];
          break;
        case 18:
          vec3Pos[1] = pc3 * std::stod(token) + center[1];
          break;
        case 19:
          vec3Pos[2] = pc3 * std::stod(token) + center[2];
          break;
      }
      col++;
    }

    printf(
        "Soma #%d PCA Results. Center (%f, %f, %f) | eigenvals (%f, %f, %f) | "
        "vec1Pos (%f, %f, %f) | vec2Pos (%f, %f, %f) | vec3Pos (%f, %f, "
        "%f)\n\n",
        somaID, center[0], center[1], center[2], pc1, pc2, pc3, vec1Pos[0],
        vec1Pos[1], vec1Pos[2], vec2Pos[0], vec2Pos[1], vec2Pos[2], vec3Pos[0],
        vec3Pos[1], vec3Pos[2]);

    // Visualize
    drawLine(pcaVisualization, center, vec1Pos);
    drawLine(pcaVisualization, center, vec2Pos);
    drawLine(pcaVisualization, center, vec3Pos);
  }

  inFile.close();

  // Show the PCA visualization
  v3dhandle newwin = callback.newImageWindow();
  callback.setImage(newwin, pcaVisualization);
}

void drawLine(Image4DSimple *image, double *from, double *to) {
  // convert points from world space to image space
  int x1 = (int)((from[0] - image->getOriginX()) / image->getRezX());
  int y1 = (int)((from[1] - image->getOriginY()) / image->getRezY());
  int z1 = (int)((from[2] - image->getOriginZ()) / image->getRezZ());
  int x2 = (int)((to[0] - image->getOriginX()) / image->getRezX());
  int y2 = (int)((to[1] - image->getOriginY()) / image->getRezY());
  int z2 = (int)((to[2] - image->getOriginZ()) / image->getRezZ());

  unsigned char *imgData = image->getRawData();

  auto fillPixel = [&](int x, int y, int z) {
    if (x >= 0 && x < image->getXDim() && y >= 0 && y < image->getYDim() &&
        z >= 0 && z < image->getZDim()) {
      imgData[z * image->getYDim() * image->getXDim() + y * image->getXDim() +
              x] = 255;
    }
  };

  // Bresenham's line algorithm
  int dx = abs(x2 - x1), xs = x2 > x1 ? 1 : -1;
  int dy = abs(y2 - y1), ys = y2 > y1 ? 1 : -1;
  int dz = abs(z2 - z1), zs = z2 > z1 ? 1 : -1;
  int drivingAxis;
  // Determine which difference is largest
  if (dx >= dy && dx >= dz)
    drivingAxis = 0;  // X-axis
  else if (dy >= dx && dy >= dz)
    drivingAxis = 1;  // Y-axis
  else
    drivingAxis = 2;  // Z-axis

  // Plot the initial point
  fillPixel(x1, y1, z1);

  if (drivingAxis == 0) {
    int p1 = 2 * dy - dx;
    int p2 = 2 * dz - dx;
    while (x1 != x2) {
      x1 += xs;
      if (p1 >= 0) {
        y1 += ys;
        p1 -= 2 * dx;
      }
      if (p2 >= 0) {
        z1 += zs;
        p2 -= 2 * dx;
      }
      p1 += 2 * dy;
      p2 += 2 * dz;
      fillPixel(x1, y1, z1);
    }
  } else if (drivingAxis == 1) {
    int p1 = 2 * dx - dy;
    int p2 = 2 * dz - dy;
    while (y1 != y2) {
      y1 += ys;
      if (p1 >= 0) {
        x1 += xs;
        p1 -= 2 * dy;
      }
      if (p2 >= 0) {
        z1 += zs;
        p2 -= 2 * dy;
      }
      p1 += 2 * dx;
      p2 += 2 * dz;
      fillPixel(x1, y1, z1);
    }
  } else {
    int p1 = 2 * dy - dz;
    int p2 = 2 * dx - dz;
    while (z1 != z2) {
      z1 += zs;
      if (p1 >= 0) {
        y1 += ys;
        p1 -= 2 * dz;
      }
      if (p2 >= 0) {
        x1 += xs;
        p2 -= 2 * dz;
      }
      p1 += 2 * dy;
      p2 += 2 * dx;
      fillPixel(x1, y1, z1);
    }
  }
}
