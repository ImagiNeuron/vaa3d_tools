#include "pcaAnalysis.h"

// Function to perform PCA analysis

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
      // QString suggestedName = QFileInfo(filename).fileName();
      // actualFilename = QFileDialog::getSaveFileName(
      //     parent, "Save PCA Results", suggestedName, "CSV Files (*.csv)");
      // if (actualFilename.isEmpty()) {
      //   printf("Save cancelled by user.\n");
      //   shouldSave = false;
      //   if (saveEnabled) *saveEnabled = false;
      //   return;
      // }
      // actualFilename = suggestedName;

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
  printf("PCA results saved to file: %s\n",
         actualFilename.toStdString().c_str());
}

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
    QWidget *mainWin = QApplication::activeWindow();
    bool saveEnabled = false;
    savePCAResultsToCSV(savePath, somaIndex, lm, pc1, pc2, pc3, vec1, vec2,
                        vec3, x_center, y_center, z_center, mainWin,
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

QString modifyFileNameForTeraFly(const QString &fileName) {
  QString modifiedFileName = fileName;

  // we are using terafly if the fileName starts with ID
  if (modifiedFileName.startsWith("ID")) {
    // check if the MIND folder already exists, if not, create it
    QDir dir("MIND");
    if (!dir.exists()) {
      dir.mkpath(".");
    }

    // replace ID(%), with ""
    modifiedFileName.replace(QRegularExpression("ID\\(.*\\), "), "");
    modifiedFileName.replace("1 channels_processed", "");

    // check if folder for this image exists
    QDir dir2("MIND/" + modifiedFileName);
    if (!dir2.exists()) {
      dir2.mkpath(".");
    }

    modifiedFileName = "MIND/" + modifiedFileName + "/";
  }

  return modifiedFileName;
}

void loadSegmentationFile(const QString &imageName, unsigned char *&segData,
                          V3DLONG sz[4], int &datatype,
                          V3DPluginCallback2 &callback, QWidget *parent) {
  // Construct segmentation filename (try different options)
  QString segFileName =
      modifyFileNameForTeraFly(imageName + "_binary_segmentation.tif");

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
    QWidget *mainWin = QApplication::activeWindow();
    bool saveEnabled = false;
    savePCAResultsToCSV(savePath, somaIndex, lm, pc1, pc2, pc3, vec1, vec2,
                        vec3, x_center, y_center, z_center, mainWin,
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

unsigned char *create_background(V3DPluginCallback2 &callback, QWidget *parent,
                                 V3DLONG &dim_X, V3DLONG &dim_Y,
                                 V3DLONG &dim_Z) {
  // Get current image
  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image opened.", parent);
    return nullptr;
  }

  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("No image opened.", parent);
    return nullptr;
  }

  // Get image dimensions and data
  dim_X = p4DImage->getXDim();
  dim_Y = p4DImage->getYDim();
  dim_Z = p4DImage->getZDim();
  V3DLONG totalSize = dim_X * dim_Y * dim_Z;
  unsigned char *originalData = p4DImage->getRawData();

  // Calculate chunk dimensions (ensure at least 1)
  V3DLONG factor = 4;
  V3DLONG chunk_X = std::max(1L, dim_X / factor);
  V3DLONG chunk_Y = std::max(1L, dim_Y / factor);
  V3DLONG chunk_Z = std::max(1L, dim_Z / factor);

  // Calculate number of chunks in each dimension
  V3DLONG num_chunks_X = (dim_X + chunk_X - 1) / chunk_X;
  V3DLONG num_chunks_Y = (dim_Y + chunk_Y - 1) / chunk_Y;
  V3DLONG num_chunks_Z = (dim_Z + chunk_Z - 1) / chunk_Z;

  printf("Image dimensions: %ld x %ld x %ld\n", dim_X, dim_Y, dim_Z);
  printf("Chunk dimensions: %ld x %ld x %ld\n", chunk_X, chunk_Y, chunk_Z);
  printf("Number of chunks: %ld x %ld x %ld = %ld chunks total\n", num_chunks_X,
         num_chunks_Y, num_chunks_Z,
         num_chunks_X * num_chunks_Y * num_chunks_Z);

  // Allocate 3D array for background
  unsigned char *backgroundArray = new unsigned char[totalSize];

  // Seed random generator
  std::srand(std::time(nullptr));

  // Store statistics for each chunk for blending
  std::vector<std::vector<std::vector<std::pair<double, double>>>> chunkStats(
      num_chunks_Z, std::vector<std::vector<std::pair<double, double>>>(
                        num_chunks_Y, std::vector<std::pair<double, double>>(
                                          num_chunks_X, {0.0, 0.0})));

  // Calculate statistics for each chunk first
  for (V3DLONG cz = 0; cz < num_chunks_Z; cz++) {
    V3DLONG z_start = cz * chunk_Z;
    V3DLONG z_end = std::min(z_start + chunk_Z, dim_Z);

    for (V3DLONG cy = 0; cy < num_chunks_Y; cy++) {
      V3DLONG y_start = cy * chunk_Y;
      V3DLONG y_end = std::min(y_start + chunk_Y, dim_Y);

      for (V3DLONG cx = 0; cx < num_chunks_X; cx++) {
        V3DLONG x_start = cx * chunk_X;
        V3DLONG x_end = std::min(x_start + chunk_X, dim_X);

        // Create histogram for this chunk
        int hist[256] = {0};
        int chunkSize = 0;

        for (V3DLONG z = z_start; z < z_end; z++) {
          for (V3DLONG y = y_start; y < y_end; y++) {
            for (V3DLONG x = x_start; x < x_end; x++) {
              V3DLONG idx = z * dim_X * dim_Y + y * dim_X + x;
              hist[originalData[idx]]++;
              chunkSize++;
            }
          }
        }

        // Calculate Otsu threshold for this chunk using the new function
        int threshold = calculateOtsuThreshold(hist, chunkSize);

        // Calculate statistics for background voxels (below threshold)
        double mean = 0.0;
        double stdDev = 0.0;
        V3DLONG backgroundCount = 0;
        double sum_bg = 0.0;
        double sumSq_bg = 0.0;

        for (V3DLONG z = z_start; z < z_end; z++) {
          for (V3DLONG y = y_start; y < y_end; y++) {
            for (V3DLONG x = x_start; x < x_end; x++) {
              V3DLONG idx = z * dim_X * dim_Y + y * dim_X + x;
              if (originalData[idx] < threshold) {
                sum_bg += originalData[idx];
                sumSq_bg += originalData[idx] * originalData[idx];
                backgroundCount++;
              }
            }
          }
        }

        // If we found background voxels, calculate statistics
        if (backgroundCount > 0) {
          mean = sum_bg / backgroundCount;
          double variance = (sumSq_bg / backgroundCount) - (mean * mean);
          stdDev = sqrt(variance);

          // Ensure minimum stdDev to avoid issues with uniform regions
          if (stdDev < 1.0) stdDev = 1.0;
        } else {
          // If no background voxels found, use default values
          mean = 10.0;
          stdDev = 3.0;
        }

        // Store the mean and stdDev for this chunk
        chunkStats[cz][cy][cx] = {mean, stdDev};
      }
    }
  }

  // Define blending parameters - use a smaller radius for more localized
  // blending
  double blendRadius = 1.0;

  // Generate background values with blending between chunks
  for (V3DLONG z = 0; z < dim_Z; z++) {
    for (V3DLONG y = 0; y < dim_Y; y++) {
      for (V3DLONG x = 0; x < dim_X; x++) {
        // Get blended distribution parameters for this voxel
        double blendedMean, blendedStdDev;
        getBlendedDistributionParams(x, y, z, chunkStats, chunk_X, chunk_Y,
                                     chunk_Z, num_chunks_X, num_chunks_Y,
                                     num_chunks_Z, blendRadius, blendedMean,
                                     blendedStdDev);

        // Generate random value from the blended distribution
        double u1 = std::rand() / (RAND_MAX + 1.0);
        double u2 = std::rand() / (RAND_MAX + 1.0);

        // Avoid log(0)
        if (u1 < 1e-10) u1 = 1e-10;

        double randStdNormal =
            std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        double randNormal = blendedMean + blendedStdDev * randStdNormal;

        // Clamp to valid unsigned char range [0,255]
        int value = std::round(randNormal);
        value = std::max(0, std::min(255, value));

        backgroundArray[z * dim_X * dim_Y + y * dim_X + x] =
            (unsigned char)value;
      }
    }
  }

  printf(
      "Background generation complete with smooth blending (radius = %.1f).\n",
      blendRadius);
  return backgroundArray;
}

void free_mapped_arrays(unsigned char ***intensities,
                        unsigned char ***segmentation, V3DLONG dim_Z,
                        V3DLONG dim_Y) {
  // Free intensities array
  if (intensities) {
    for (V3DLONG z = 0; z < dim_Z; z++) {
      if (intensities[z]) {
        for (V3DLONG y = 0; y < dim_Y; y++) {
          if (intensities[z][y]) delete[] intensities[z][y];
        }
        delete[] intensities[z];
      }
    }
    delete[] intensities;
  }

  // Free segmentation array
  if (segmentation) {
    for (V3DLONG z = 0; z < dim_Z; z++) {
      if (segmentation[z]) {
        for (V3DLONG y = 0; y < dim_Y; y++) {
          if (segmentation[z][y]) delete[] segmentation[z][y];
        }
        delete[] segmentation[z];
      }
    }
    delete[] segmentation;
  }
}

void create_background(V3DPluginCallback2 &callback, QWidget *parent) {
  // Old version that creates and displays an image directly for debugging

  // Call the new version to get the background intensities
  V3DLONG dimX, dimY, dimZ;
  unsigned char *backgroundIntensities =
      create_background(callback, parent, dimX, dimY, dimZ);

  if (!backgroundIntensities) {
    v3d_msg("Failed to generate background intensities.", parent);
    return;
  }

  // Create new image for display
  Image4DSimple *backgroundImage = new Image4DSimple();
  backgroundImage->createBlankImage(dimX, dimY, dimZ, 1, V3D_UINT8);

  // Get image properties from current image for consistency
  v3dhandle curwin = callback.currentImageWindow();
  if (curwin) {
    Image4DSimple *p4DImage = callback.getImage(curwin);
    if (p4DImage) {
      backgroundImage->setOriginX(p4DImage->getOriginX());
      backgroundImage->setOriginY(p4DImage->getOriginY());
      backgroundImage->setOriginZ(p4DImage->getOriginZ());
      backgroundImage->setRezX(p4DImage->getRezX());
      backgroundImage->setRezY(p4DImage->getRezY());
      backgroundImage->setRezZ(p4DImage->getRezZ());
    }
  }

  // Copy the background intensities to the new image
  unsigned char *backgroundData = backgroundImage->getRawData();
  std::memcpy(backgroundData, backgroundIntensities, dimX * dimY * dimZ);

  // Free the background intensities array
  delete[] backgroundIntensities;

  // Show the generated background
  QString imageName = callback.getImageName(curwin);
  v3dhandle newwin = callback.newImageWindow();
  callback.setImage(newwin, backgroundImage);
  callback.setImageName(newwin, imageName + "_background");
  callback.updateImageWindow(newwin);

  v3d_msg("Background image generated successfully.");
}

/**
 * @brief Calculate Otsu threshold for a histogram
 * @param hist Array of 256 histogram values
 * @param totalPixels Total number of pixels in the region
 * @return The calculated Otsu threshold value
 */
int calculateOtsuThreshold(const int hist[256], int totalPixels) {
  double sum = 0;
  for (int t = 0; t < 256; t++) {
    sum += t * hist[t];
  }

  double sumB = 0;
  int wB = 0;
  double varMax = 0;
  int threshold = 0;

  for (int t = 0; t < 256; t++) {
    wB += hist[t];
    if (wB == 0) continue;

    int wF = totalPixels - wB;
    if (wF == 0) break;

    sumB += t * hist[t];
    double mB = sumB / wB;
    double mF = (sum - sumB) / wF;
    double varBetween = (double)wB * wF * (mB - mF) * (mB - mF);

    if (varBetween > varMax) {
      varMax = varBetween;
      threshold = t;
    }
  }

  return threshold;
}

void getBlendedDistributionParams(
    V3DLONG x, V3DLONG y, V3DLONG z,
    const std::vector<std::vector<std::vector<std::pair<double, double>>>>
        &chunkStats,
    V3DLONG chunk_X, V3DLONG chunk_Y, V3DLONG chunk_Z, V3DLONG num_chunks_X,
    V3DLONG num_chunks_Y, V3DLONG num_chunks_Z, double blendRadius,
    double &blendedMean, double &blendedStdDev) {
  // Calculate chunk coordinates and relative position
  double cx_pos = (double)x / chunk_X;
  V3DLONG cx = std::min(num_chunks_X - 1, (V3DLONG)cx_pos);

  double cy_pos = (double)y / chunk_Y;
  V3DLONG cy = std::min(num_chunks_Y - 1, (V3DLONG)cy_pos);

  double cz_pos = (double)z / chunk_Z;
  V3DLONG cz = std::min(num_chunks_Z - 1, (V3DLONG)cz_pos);

  // Variables for weighted blending
  double totalWeight = 0.0;
  double weightedMean = 0.0;
  double weightedVar = 0.0;

  // Iterate over neighboring chunks within blend radius
  for (int nz = std::max(0L, cz - (V3DLONG)blendRadius);
       nz <= std::min(num_chunks_Z - 1, cz + (V3DLONG)blendRadius); nz++) {
    for (int ny = std::max(0L, cy - (V3DLONG)blendRadius);
         ny <= std::min(num_chunks_Y - 1, cy + (V3DLONG)blendRadius); ny++) {
      for (int nx = std::max(0L, cx - (V3DLONG)blendRadius);
           nx <= std::min(num_chunks_X - 1, cx + (V3DLONG)blendRadius); nx++) {
        // Calculate distance to chunk center in chunk-space
        double dx = (cx_pos - nx - 0.5);
        double dy = (cy_pos - ny - 0.5);
        double dz = (cz_pos - nz - 0.5);
        double distSq = dx * dx + dy * dy + dz * dz;

        // Skip chunks that are too far
        if (distSq > blendRadius * blendRadius) continue;

        // Calculate weight based on distance (quadratic falloff)
        double weight = std::max(0.0, 1.0 - std::sqrt(distSq) / blendRadius);
        weight = weight * weight;  // Square the weight for smoother falloff

        // Get the statistics for this chunk
        double mean = chunkStats[nz][ny][nx].first;
        double stdDev = chunkStats[nz][ny][nx].second;

        // Accumulate weighted statistics
        weightedMean += weight * mean;
        weightedVar += weight * stdDev * stdDev;  // Weighted variance
        totalWeight += weight;
      }
    }
  }

  // Normalize by total weight
  if (totalWeight > 0) {
    weightedMean /= totalWeight;
    weightedVar /= totalWeight;
  } else {
    // Fallback to central chunk if no weights (should not happen)
    weightedMean = chunkStats[cz][cy][cx].first;
    weightedVar = chunkStats[cz][cy][cx].second * chunkStats[cz][cy][cx].second;
  }

  // Set output parameters
  blendedMean = weightedMean;
  blendedStdDev = std::sqrt(weightedVar);
}
