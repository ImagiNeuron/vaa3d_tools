 /**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */
#include "pcaAnalysis.h"

/**
 * @brief Function to save PCA results to a CSV file
 *
 * @param filename The name of the file to save the results to
 * @param somaIndex The index of the soma being analyzed
 * @param lm The location of the soma
 * @param pc1 The first principal component
 * @param pc2 The second principal component
 * @param pc3 The third principal component
 * @param vec1 The first eigenvector
 * @param vec2 The second eigenvector
 * @param vec3 The third eigenvector
 * @param x_center The x-coordinate of the center of mass
 * @param y_center The y-coordinate of the center of mass
 * @param z_center The z-coordinate of the center of mass
 * @param parent The parent widget for the file dialog
 */
void savePCAResultsToCSV(const QString &filename, int somaIndex,
                         const LocationSimple &lm, double pc1, double pc2,
                         double pc3, const double *vec1, const double *vec2,
                         const double *vec3, double x_center, double y_center,
                         double z_center, QWidget *parent) {
  static QString actualFilename = filename;

  // Ask user if they want to save only for the first soma
  if (somaIndex == 1) {
    actualFilename = filename;

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
        return;
      }
    }

    // Write header if new file
    std::ofstream outFile(actualFilename.toStdString().c_str(), std::ios::app);
    outFile << "SomaID,X,Y,Z,Radius,CenterMassX,CenterMassY,CenterMassZ,"
            << "eigenvalue1,eigenvalue2,eigenvalue3,"
            << "eigenvector1_x,eigenvector1_y,eigenvector1_z,"
            << "eigenvector2_x,eigenvector2_y,eigenvector2_z,"
            << "eigenvector3_x,eigenvector3_y,eigenvector3_z\n";
    outFile.close();
  }

  std::ofstream outFile(actualFilename.toStdString().c_str(), std::ios::app);

  // Write data
  outFile << somaIndex << "," << lm.x << "," << lm.y << "," << lm.z << ","
          << lm.radius << "," << x_center << "," << y_center << "," << z_center
          << "," << pc1 << "," << pc2 << "," << pc3 << "," << vec1[0] << ","
          << vec1[1] << "," << vec1[2] << "," << vec2[0] << "," << vec2[1]
          << "," << vec2[2] << "," << vec3[0] << "," << vec3[1] << ","
          << vec3[2] << "\n";

  outFile.close();
  printf("PCA results saved to file: %s\n",
         actualFilename.toStdString().c_str());
}

/**
 * * @brief Function to analyze a soma using PCA and save results to a CSV file
 *
 * @param labeledData The labeled data of the image
 * @param N The width of the image
 * @param M The height of the image
 * @param P The depth of the image
 * @param lm The location of the soma
 * @param somaIndex The index of the soma being analyzed
 * @param savePath The path to save the results
 * @param pc1 The pointer to the first principal component
 * @param pc2 The pointer to the second principal component
 * @param pc3 The pointer to the third principal component
 * @param vec1 The pointer to the first eigenvector
 * @param vec2 The pointer to the second eigenvector
 * @param vec3 The pointer to the third eigenvector
 * @param x_center The pointer to the x-coordinate of the center of mass
 * @param y_center The pointer to the y-coordinate of the center of mass
 * @param z_center The pointer to the z-coordinate of the center of mass
 */
void analyzeSomaPCAReturnResults(unsigned char *labeledData, V3DLONG N,
                                 V3DLONG M, V3DLONG P, const LocationSimple &lm,
                                 int somaIndex, QString savePath, double &pc1,
                                 double &pc2, double &pc3, double vec1[3],
                                 double vec2[3], double vec3[3],
                                 double &x_center, double &y_center,
                                 double &z_center) {
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

  // Use sphere window type (1) and window size based on soma radius
  if (compute_sphere_win3d_pca_eigVec(
          img3d, N, M, P, x, y, z,  // Center on soma
          2 * r, 2 * r, 2 * r,      // Window size based on radius
          pc1, pc2, pc3, vec1, vec2, vec3, x_center, y_center, z_center)) {
    // Save to CSV with save flag
    QWidget *mainWin = QApplication::activeWindow();
    savePCAResultsToCSV(savePath, somaIndex, lm, pc1, pc2, pc3, vec1, vec2,
                        vec3, x_center, y_center, z_center, mainWin);
  } else {
    printf("\nSoma #%d PCA failed.\n", somaIndex);
  }

  // Cleanup 3D array wrapper
  for (V3DLONG k = 0; k < P; k++) {
    delete[] img3d[k];
  }
  delete[] img3d;
}

/**
 * @brief Function to modify file path for TeraFly data structure
 *
 * @param fileName The original file name
 * @return The modified file path
 */
QString modifyFilePathForTeraFly(const QString &fileName) {
  QString modifiedFilePath = fileName;

  // we are using terafly if the fileName starts with ID
  if (modifiedFilePath.startsWith("ID")) {
    // check if the MIND folder already exists, if not, create it
    QDir dir("MIND");
    if (!dir.exists()) {
      dir.mkpath(".");
    }

    // replace ID(%), with ""
    modifiedFilePath.replace(QRegularExpression("ID\\(.*\\), "), "");
    modifiedFilePath.replace("1 channels_processed", "");

    // check if folder for this image exists
    QDir dir2("MIND/" + modifiedFilePath);
    if (!dir2.exists()) {
      dir2.mkpath(".");
    }

    modifiedFilePath = QDir::currentPath() + "/MIND/" + modifiedFilePath + "/";
  } else {
    // if not using terafly, just remove the file extension
    modifiedFilePath = modifiedFilePath.split(".")[0];
  }

  return modifiedFilePath;
}

/**
 * * @brief Function to load segmentation file
 *
 * @param imageName The name of the image
 * @param segData The segmentation data
 * @param sz The size of the image
 * @param datatype The data type of the image
 * @param callback The V3D plugin callback interface
 * @param parent The parent interface
 */
void loadSegmentationFile(const QString &imageName, unsigned char *&segData,
                          V3DLONG sz[4], int &datatype,
                          V3DPluginCallback2 &callback, QWidget *parent) {
  // Construct segmentation filename (try different options)
  QString segFileName =
      modifyFilePathForTeraFly(imageName) + "_binary_segmentation.tif";

  printf("Looking for segmentation file: %s\n",
         segFileName.toStdString().c_str());

  bool foundSegFile = false;
  if (QFile::exists(segFileName)) {
    foundSegFile = true;
    printf("Found segmentation file: %s\n", segFileName.toStdString().c_str());
  }

  // If segmentation file still not found, ask the user to select it
  if (!foundSegFile) {
    v3d_msg("No segmentation file found. Please segment the image first.",
            parent);
    return;
  }

  // Load the binary segmentation file
  if (!simple_loadimage_wrapper(callback, segFileName.toStdString().c_str(),
                                segData, sz, datatype)) {
    v3d_msg("Failed to load segmentation file.", parent);
    return;
  }
}

/**
 * * @brief Function to analyze a soma using PCA and save results to a CSV file,
 * wihtout returning the results
 *
 * @param labeledData The labeled data of the image
 * @param N The width of the image
 * @param M The height of the image
 * @param P The depth of the image
 * @param lm The location of the soma
 * @param somaIndex The index of the soma being analyzed
 * @param savePath The path to save the results
 */
void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex, QString savePath) {
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
    // Save to CSV with save flag
    QWidget *mainWin = QApplication::activeWindow();
    savePCAResultsToCSV(savePath, somaIndex, lm, pc1, pc2, pc3, vec1, vec2,
                        vec3, x_center, y_center, z_center, mainWin);
  } else {
    printf("\nSoma #%d PCA failed.\n", somaIndex);
  }

  // Cleanup 3D array wrapper
  for (V3DLONG k = 0; k < P; k++) {
    delete[] img3d[k];
  }
  delete[] img3d;
}

/**
 * @brief Function to visuzalize PCA results
 *
 * @param callback The V3D plugin callback interface
 * @param parent The parent interface
 */
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
                                     p4DImage->getZDim(), 3, V3D_UINT8);
  pcaVisualization->setOriginX(p4DImage->getOriginX());
  pcaVisualization->setOriginY(p4DImage->getOriginY());
  pcaVisualization->setOriginZ(p4DImage->getOriginZ());
  pcaVisualization->setRezX(p4DImage->getRezX());
  pcaVisualization->setRezY(p4DImage->getRezY());
  pcaVisualization->setRezZ(p4DImage->getRezZ());

  // copy original image to channel 1
  const int size = p4DImage->getTotalUnitNumber() * p4DImage->getUnitBytes();
  memcpy(pcaVisualization->getRawData(), p4DImage->getRawData(), size);
  memcpy(pcaVisualization->getRawData() + size, p4DImage->getRawData(), size);
  memcpy(pcaVisualization->getRawData() + 2 * size, p4DImage->getRawData(),
         size);

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
    double pc1Length, pc2Length, pc3Length;
    double vec1Pos[3], vec2Pos[3], vec3Pos[3];
    // double radius;

    int col = 0;
    while (std::getline(ss, token, ',')) {
      switch (col) {
        case 0:
          somaID = std::stoi(token);
          break;
        // case 1,2,3: // marker position
        // case 4:
        //   radius = std::stod(token);
        //   break;
        case 5:
          center[0] = std::stod(token);
          break;
        case 6:
          center[1] = std::stod(token);
          break;
        case 7:
          center[2] = std::stod(token);
          break;
        case 8:
          pc1Length = 2.35 * sqrt(std::stod(token));
          break;
        case 9:
          pc2Length = 2.35 * sqrt(std::stod(token));
          break;
        case 10:
          pc3Length = 2.35 * sqrt(std::stod(token));
          break;
        case 11:
          vec1Pos[0] = pc1Length * std::stod(token) + center[0];
          break;
        case 12:
          vec1Pos[1] = pc1Length * std::stod(token) + center[1];
          break;
        case 13:
          vec1Pos[2] = pc1Length * std::stod(token) + center[2];
          break;
        case 14:
          vec2Pos[0] = pc2Length * std::stod(token) + center[0];
          break;
        case 15:
          vec2Pos[1] = pc2Length * std::stod(token) + center[1];
          break;
        case 16:
          vec2Pos[2] = pc2Length * std::stod(token) + center[2];
          break;
        case 17:
          vec3Pos[0] = pc3Length * std::stod(token) + center[0];
          break;
        case 18:
          vec3Pos[1] = pc3Length * std::stod(token) + center[1];
          break;
        case 19:
          vec3Pos[2] = pc3Length * std::stod(token) + center[2];
          break;
      }
      col++;
    }

    // Visualize: green longest, red second, blue (cyan for better visibility)
    // third
    drawLine(pcaVisualization, 0, 255, 0, center, vec1Pos);
    drawLine(pcaVisualization, 255, 0, 0, center, vec2Pos);
    drawLine(pcaVisualization, 0, 255, 255, center, vec3Pos);
  }

  inFile.close();

  // Show the PCA visualization
  v3dhandle newwin = callback.newImageWindow();
  callback.setImage(newwin, pcaVisualization);
}

/**
 * @brief Function to draw a line in the image
 *
 * @param image - the image to draw on
 * @param r - red color value
 * @param g - green color value
 * @param b - blue color value
 * @param from - starting point of the line
 * @param to - ending point of the line
 */
void drawLine(Image4DSimple *image, unsigned char r, unsigned char g,
              unsigned char b, double *from, double *to) {
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
      const int channel_size =
          image->getXDim() * image->getYDim() * image->getZDim();
      const int subindex =
          z * image->getXDim() * image->getYDim() + y * image->getXDim() + x;
      imgData[subindex] = r;
      imgData[channel_size + subindex] = g;
      imgData[2 * channel_size + subindex] = b;
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
