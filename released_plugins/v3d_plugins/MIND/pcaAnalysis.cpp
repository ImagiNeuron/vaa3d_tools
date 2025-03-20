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

void simulate_somas(V3DPluginCallback2 &callback, QWidget *parent) {
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

  // Get the current image name and path to construct the segmentation filename
  QString imageName = callback.getImageName(curwin);
  QString currentImagePath = QFileInfo(imageName).absolutePath();
  QString baseImageName = QFileInfo(imageName).baseName();

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

  // Get current landmarks (soma markers)
  LandmarkList markers = callback.getLandmark(curwin);
  if (markers.isEmpty()) {
    v3d_msg("No markers found. Please add markers to identify somas.", parent);
    if (segData) {
      delete[] segData;
    }
    return;
  }

  // Get the original image data for intensity extraction
  unsigned char *originalData = p4DImage->getRawData();
  int channel = 0;  // Default to first channel (index 0)

  // Use the PCA file created by analyzeSomaPCA (same naming convention)
  QString pcaFileName = imageName + "_pca.csv";

  // Check if the PCA file exists
  if (!QFile::exists(pcaFileName)) {
    v3d_msg(QString("PCA file not found: %1\nPlease run PC Analysis first.")
                .arg(pcaFileName),
            parent);
    if (segData) {
      delete[] segData;
    }
    return;
  }

  printf("Using PCA file: %s\n", pcaFileName.toStdString().c_str());

  // Read the PCA CSV file to get the center of mass
  std::ifstream pcaFile(pcaFileName.toStdString().c_str());
  if (!pcaFile.is_open()) {
    v3d_msg("Failed to open PCA file.", parent);
    if (segData) {
      delete[] segData;
    }
    return;
  }

  // Skip header line
  std::string line;
  std::getline(pcaFile, line);

  // Variables to store center of mass and selected marker info
  double x_center = 0, y_center = 0, z_center = 0;
  int somaIndex = -1;
  LocationSimple selectedMarker;
  bool foundMarker = false;

  // Read data lines
  while (std::getline(pcaFile, line)) {
    std::istringstream ss(line);
    std::string token;
    std::vector<double> values;

    // Parse each column from CSV
    while (std::getline(ss, token, ',')) {
      values.push_back(std::stod(token));
    }

    // Check if we have enough columns
    if (values.size() >= 8) {
      somaIndex = static_cast<int>(values[0]);

      // Find the corresponding marker
      for (int i = 0; i < markers.size(); i++) {
        double dx = markers[i].x - values[1];  // marker x - csv x
        double dy = markers[i].y - values[2];  // marker y - csv y
        double dz = markers[i].z - values[3];  // marker z - csv z

        // If this marker is close to the position in CSV
        if (dx * dx + dy * dy + dz * dz < 25) {  // Within 5 pixel distance
          selectedMarker = markers[i];
          foundMarker = true;

          // Get center of mass from columns 5, 6, 7
          x_center = values[5];  // CenterMassX
          y_center = values[6];  // CenterMassY
          z_center = values[7];  // CenterMassZ

          printf("Found soma #%d with center of mass (%f, %f, %f)\n", somaIndex,
                 x_center, y_center, z_center);
          break;
        }
      }

      if (foundMarker) break;  // Exit after finding the first matching marker
    }
  }
  pcaFile.close();

  // If no center of mass was found, select a random marker
  if (!foundMarker) {
    v3d_msg(
        "Could not find matching marker in PCA file. Selecting a random "
        "marker.",
        parent);
    int randIndex = rand() % markers.size();
    selectedMarker = markers[randIndex];
    somaIndex = randIndex + 1;

    // Use marker position as center of mass
    x_center = selectedMarker.x;
    y_center = selectedMarker.y;
    z_center = selectedMarker.z;

    printf(
        "Selected random marker #%d at position (%f, %f, %f) with radius %f\n",
        somaIndex, selectedMarker.x, selectedMarker.y, selectedMarker.z,
        selectedMarker.radius);
  }

  // Calculate cube size based on the soma radius (ensure adequate space)
  double radius = selectedMarker.radius > 0 ? selectedMarker.radius : 10.0;
  V3DLONG cubeSize = ((V3DLONG)ceil(radius) + 3) *
                     2;  // Same sizing as used in main segmentation
  V3DLONG totalVoxels = cubeSize * cubeSize * cubeSize;

  // Create the somaSegmentation array and somaIntensity array
  int *somaSegmentation = new int[totalVoxels];
  unsigned char *somaIntensities = new unsigned char[totalVoxels];
  memset(somaSegmentation, 0, totalVoxels * sizeof(int));
  memset(somaIntensities, 0, totalVoxels * sizeof(unsigned char));

  // Get the segmentation dimension information
  V3DLONG dimX = sz[0];
  V3DLONG dimY = sz[1];
  V3DLONG dimZ = sz[2];

  // Make sure the original image dimensions match segmentation
  if (dimX != p4DImage->getXDim() || dimY != p4DImage->getYDim() ||
      dimZ != p4DImage->getZDim()) {
    v3d_msg("Original image and segmentation dimensions don't match.", parent);
    if (segData) {
      delete[] segData;
    }
    delete[] somaSegmentation;
    delete[] somaIntensities;
    return;
  }

  printf("Using center of mass: (%f, %f, %f)\n", x_center, y_center, z_center);

  // Copy both the binary segmentation and intensities into arrays
  // centering the soma at the center of the cube
  int center = cubeSize / 2;
  V3DLONG sliceSize = dimX * dimY;
  V3DLONG channelOffset = channel * dimX * dimY * dimZ;

  for (V3DLONG z = 0; z < dimZ; z++) {
    for (V3DLONG y = 0; y < dimY; y++) {
      for (V3DLONG x = 0; x < dimX; x++) {
        V3DLONG idx = z * dimX * dimY + y * dimX + x;

        // Check if the voxel is part of the soma
        if (segData[idx] > 0) {
          // Calculate squared distance to the marker
          double dx = x - selectedMarker.x;
          double dy = y - selectedMarker.y;
          double dz = z - selectedMarker.z;
          double distSq = dx * dx + dy * dy + dz * dz;

          // Include this voxel if it's close to the selected marker
          if (distSq <= 4 * radius * radius) {
            // Compute coordinates relative to the center of mass
            int relX = x - (int)round(x_center) + center;
            int relY = y - (int)round(y_center) + center;
            int relZ = z - (int)round(z_center) + center;

            // Check if coordinates are within the cube bounds
            if (relX >= 0 && relX < cubeSize && relY >= 0 && relY < cubeSize &&
                relZ >= 0 && relZ < cubeSize) {
              // Get index in cube
              V3DLONG cubeIdx =
                  relZ * cubeSize * cubeSize + relY * cubeSize + relX;

              // Set binary segmentation
              somaSegmentation[cubeIdx] = 1;

              // Copy intensity value from original image
              V3DLONG origIdx = channelOffset + z * sliceSize + y * dimX + x;
              somaIntensities[cubeIdx] = originalData[origIdx];
            }
          }
        }
      }
    }
  }

  // Print the central slices of both extracted binary segmentation and
  // intensities
  int centralSlice = cubeSize / 2;
  printf("\nExtracted soma binary segmentation (central slice):\n");
  for (int y = 0; y < cubeSize; y++) {
    for (int x = 0; x < cubeSize; x++) {
      V3DLONG idx = centralSlice * cubeSize * cubeSize + y * cubeSize + x;
      printf("%d ", somaSegmentation[idx]);
    }
    printf("\n");
  }

  printf("\nExtracted soma intensities (central slice):\n");
  for (int y = 0; y < cubeSize; y++) {
    for (int x = 0; x < cubeSize; x++) {
      V3DLONG idx = centralSlice * cubeSize * cubeSize + y * cubeSize + x;
      printf("%3d ", somaIntensities[idx]);
    }
    printf("\n");
  }

  // Generate automatic filename from the base image name
  QString saveSomaPath =
      currentImagePath + "/" + baseImageName + "_soma_data.csv";

  // Offer to save the extracted soma information to a CSV file
  int saveResponse = QMessageBox::question(
      parent, "Save Extracted Soma Data",
      QString("Save extracted soma data to %1?").arg(saveSomaPath),
      QMessageBox::Yes | QMessageBox::No);

  if (saveResponse == QMessageBox::Yes) {
    std::ofstream csvFile(saveSomaPath.toStdString().c_str());
    if (csvFile.is_open()) {
      // Write header with metadata
      csvFile << "# Soma Data Extraction\n";
      csvFile << "# Image: " << baseImageName.toStdString() << "\n";
      csvFile << "# Marker position: " << selectedMarker.x << ","
              << selectedMarker.y << "," << selectedMarker.z << "\n";
      csvFile << "# Computed center of mass: " << x_center << "," << y_center
              << "," << z_center << "\n";
      csvFile << "# Cube size: " << cubeSize << "\n";
      csvFile << "# Soma radius: " << radius << "\n";
      csvFile << "# Center index in cube: " << center << "\n\n";

      // Write data headers
      csvFile << "x,y,z,segmentation,intensity\n";

      // Write all voxel data
      for (V3DLONG z = 0; z < cubeSize; z++) {
        for (V3DLONG y = 0; y < cubeSize; y++) {
          for (V3DLONG x = 0; x < cubeSize; x++) {
            V3DLONG idx = z * cubeSize * cubeSize + y * cubeSize + x;
            // Only save non-zero segmentation or intensity values to save space
            if (somaSegmentation[idx] > 0 || somaIntensities[idx] > 0) {
              csvFile << x << "," << y << "," << z << ","
                      << somaSegmentation[idx] << ","
                      << (int)somaIntensities[idx] << "\n";
            }
          }
        }
      }

      csvFile.close();
      v3d_msg(QString("Soma data saved to %1").arg(saveSomaPath));
    } else {
      v3d_msg("Failed to save soma data to CSV file.");
    }
  }

  // Clean up
  delete[] somaSegmentation;
  delete[] somaIntensities;
  if (segData) {
    delete[] segData;
  }
  v3d_msg("Soma extraction complete.");
}

unsigned char ***create_background(V3DPluginCallback2 &callback,
                                   QWidget *parent, V3DLONG &dim_X,
                                   V3DLONG &dim_Y, V3DLONG &dim_Z) {
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
  unsigned char ***backgroundArray = new unsigned char **[dim_Z];
  for (V3DLONG z = 0; z < dim_Z; z++) {
    backgroundArray[z] = new unsigned char *[dim_Y];
    for (V3DLONG y = 0; y < dim_Y; y++) {
      backgroundArray[z][y] = new unsigned char[dim_X];
      // Initialize to zero
      memset(backgroundArray[z][y], 0, dim_X * sizeof(unsigned char));
    }
  }

  // Seed random generator
  std::srand(std::time(nullptr));

  // Process each chunk
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

        // Generate background values for this chunk using normal distribution
        for (V3DLONG z = z_start; z < z_end; z++) {
          for (V3DLONG y = y_start; y < y_end; y++) {
            for (V3DLONG x = x_start; x < x_end; x++) {
              // Use Box-Muller transform to generate Gaussian distributed
              // random numbers
              double u1 = std::rand() / (RAND_MAX + 1.0);
              double u2 = std::rand() / (RAND_MAX + 1.0);

              // Avoid log(0)
              if (u1 < 1e-10) u1 = 1e-10;

              double randStdNormal =
                  std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
              double randNormal = mean + stdDev * randStdNormal;

              // Clamp to valid unsigned char range [0,255]
              int value = std::round(randNormal);
              value = std::max(0, std::min(255, value));

              backgroundArray[z][y][x] = (unsigned char)value;
            }
          }
        }

        // Print chunk info
        // printf(
        //     "Chunk [%ld,%ld,%ld]: Threshold=%d, Background mean=%.2f, "
        //     "stddev=%.2f\n",
        //     cx, cy, cz, threshold, mean, stdDev);
      }
    }
  }

  printf("Background generation complete.\n");
  return backgroundArray;
}

void free_3d_array(unsigned char ***array, V3DLONG dim_Z, V3DLONG dim_Y) {
  if (array) {
    for (V3DLONG z = 0; z < dim_Z; z++) {
      if (array[z]) {
        for (V3DLONG y = 0; y < dim_Y; y++) {
          if (array[z][y]) delete[] array[z][y];
        }
        delete[] array[z];
      }
    }
    delete[] array;
  }
}

bool get_PCA_info(int somaID, const QString &imageName, double &pc1,
                  double &pc2, double &pc3, double vec1[3], double vec2[3],
                  double vec3[3], double &x_center, double &y_center,
                  double &z_center) {
  // Construct the PCA CSV filename
  QString pcaFileName = imageName + "_pca.csv";

  // Check if the PCA file exists
  if (!QFile::exists(pcaFileName)) {
    printf("PCA file not found: %s\n", pcaFileName.toStdString().c_str());
    return false;
  }

  // Open the CSV file
  std::ifstream inFile(pcaFileName.toStdString().c_str());
  if (!inFile.is_open()) {
    printf("Failed to open PCA file: %s\n", pcaFileName.toStdString().c_str());
    return false;
  }

  // Skip the header row
  std::string line;
  std::getline(inFile, line);

  // Read each row until we find the matching soma ID
  bool found = false;
  while (std::getline(inFile, line)) {
    std::istringstream ss(line);
    std::string token;
    std::vector<double> values;

    // Parse the comma-separated values
    while (std::getline(ss, token, ',')) {
      values.push_back(std::stod(token));
    }

    // Check if we have enough columns and if this is the soma we want
    if (values.size() >= 20 && static_cast<int>(values[0]) == somaID) {
      // Extract values from the CSV columns
      // Format: SomaID,X,Y,Z,Radius,CenterMassX,CenterMassY,CenterMassZ,
      //         eigenvalue1,eigenvalue2,eigenvalue3,
      //         eigenvector1_x,eigenvector1_y,eigenvector1_z,
      //         eigenvector2_x,eigenvector2_y,eigenvector2_z,
      //         eigenvector3_x,eigenvector3_y,eigenvector3_z

      // Get center of mass (columns 5,6,7)
      x_center = values[5];
      y_center = values[6];
      z_center = values[7];

      // Get eigenvalues (columns 8,9,10)
      pc1 = values[8];
      pc2 = values[9];
      pc3 = values[10];

      // Get eigenvectors (columns 11-19)
      vec1[0] = values[11];
      vec1[1] = values[12];
      vec1[2] = values[13];

      vec2[0] = values[14];
      vec2[1] = values[15];
      vec2[2] = values[16];

      vec3[0] = values[17];
      vec3[1] = values[18];
      vec3[2] = values[19];

      found = true;
      break;
    }
  }

  // Close the file
  inFile.close();

  if (!found) {
    printf("Soma ID %d not found in PCA file: %s\n", somaID,
           pcaFileName.toStdString().c_str());
    return false;
  }

  printf("Successfully loaded PCA info for soma #%d:\n", somaID);
  printf("  Center of mass: (%f, %f, %f)\n", x_center, y_center, z_center);
  printf("  Eigenvalues: %f, %f, %f\n", pc1, pc2, pc3);
  printf("  Eigenvector 1: [%f, %f, %f]\n", vec1[0], vec1[1], vec1[2]);
  printf("  Eigenvector 2: [%f, %f, %f]\n", vec2[0], vec2[1], vec2[2]);
  printf("  Eigenvector 3: [%f, %f, %f]\n", vec3[0], vec3[1], vec3[2]);

  return true;
}

bool map_intensities(V3DPluginCallback2 &callback,
                     unsigned char ***&intensities,
                     unsigned char ***&segmentation, V3DLONG &dim_X,
                     V3DLONG &dim_Y, V3DLONG &dim_Z, int channel) {
  // Get the current image window
  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image window is open!");
    return false;
  }

  // Get the image name and path to find the segmentation file
  QString imageName = callback.getImageName(curwin);
  QString currentImagePath = QFileInfo(imageName).absolutePath();
  QString baseImageName = QFileInfo(imageName).baseName();

  // Get the image data
  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("Failed to get the image data!");
    return false;
  }

  // Get image dimensions
  dim_X = p4DImage->getXDim();
  dim_Y = p4DImage->getYDim();
  dim_Z = p4DImage->getZDim();
  V3DLONG dim_C = p4DImage->getCDim();

  // Check if the channel is valid
  if (channel < 0 || channel >= dim_C) {
    v3d_msg("Invalid channel index!");
    return false;
  }

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

  if (!foundSegFile) {
    v3d_msg("No segmentation file found. Please segment the image first.");
    return false;
  }

  // Load the binary segmentation file
  unsigned char *segData = nullptr;
  V3DLONG sz[4];
  int datatype = 0;
  if (!simple_loadimage_wrapper(callback, segFileName.toStdString().c_str(),
                                segData, sz, datatype)) {
    v3d_msg("Failed to load segmentation file.");
    return false;
  }

  // Check if segmentation dimensions match original image dimensions
  if (sz[0] != dim_X || sz[1] != dim_Y || sz[2] != dim_Z) {
    v3d_msg("Segmentation dimensions don't match the original image.");
    if (segData) delete[] segData;
    return false;
  }

  // Allocate memory for 3D arrays
  intensities = new unsigned char **[dim_Z];
  segmentation = new unsigned char **[dim_Z];
  for (V3DLONG z = 0; z < dim_Z; z++) {
    intensities[z] = new unsigned char *[dim_Y];
    segmentation[z] = new unsigned char *[dim_Y];
    for (V3DLONG y = 0; y < dim_Y; y++) {
      intensities[z][y] = new unsigned char[dim_X];
      segmentation[z][y] = new unsigned char[dim_X];
    }
  }

  // Get raw data from original image
  unsigned char *data1d = p4DImage->getRawData();

  // Copy data to 3D arrays
  V3DLONG offset_channel = dim_X * dim_Y * dim_Z * channel;
  for (V3DLONG z = 0; z < dim_Z; z++) {
    for (V3DLONG y = 0; y < dim_Y; y++) {
      for (V3DLONG x = 0; x < dim_X; x++) {
        V3DLONG idx = offset_channel + z * dim_X * dim_Y + y * dim_X + x;
        V3DLONG segIdx = z * dim_X * dim_Y + y * dim_X + x;

        intensities[z][y][x] = data1d[idx];
        segmentation[z][y][x] = segData[segIdx];
      }
    }
  }

  // Free the segmentation data as it's been copied to the 3D array
  if (segData) delete[] segData;

  return true;
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
  unsigned char ***backgroundIntensities =
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
  for (V3DLONG z = 0; z < dimZ; z++) {
    for (V3DLONG y = 0; y < dimY; y++) {
      for (V3DLONG x = 0; x < dimX; x++) {
        V3DLONG idx = z * dimX * dimY + y * dimX + x;
        backgroundData[idx] = backgroundIntensities[z][y][x];
      }
    }
  }

  // Free the 3D array now that we've copied the data
  free_3d_array(backgroundIntensities, dimZ, dimY);

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
