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
  int channel = 0;  // Default to first channel

  // If image has multiple channels, ask user which channel to use
  if (p4DImage->getCDim() > 1) {
    QDialog dialog(parent);
    QComboBox channelComboBox;
    QDialogButtonBox buttonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    QVBoxLayout layout;

    // Build channel selection
    for (int c = 0; c < p4DImage->getCDim(); c++) {
      channelComboBox.addItem(QString("Channel %1").arg(c + 1));
    }

    layout.addWidget(new QLabel("Select channel for intensity extraction:"));
    layout.addWidget(&channelComboBox);
    layout.addWidget(&buttonBox);
    dialog.setLayout(&layout);

    // Connect buttons
    QObject::connect(&buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
    QObject::connect(&buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

    if (dialog.exec() == QDialog::Accepted) {
      channel = channelComboBox.currentIndex();
    }
  }

  // Randomly select a soma
  int markerCount = markers.size();
  int randIndex = rand() % markerCount;
  LocationSimple selectedMarker = markers[randIndex];
  printf("Selected marker #%d at position (%f, %f, %f) with radius %f\n",
         randIndex + 1, selectedMarker.x, selectedMarker.y, selectedMarker.z,
         selectedMarker.radius);

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

  // Calculate center of mass of the soma
  double x_center = 0, y_center = 0, z_center = 0;
  double totalMass = 0;

  // First pass: find center of mass of the segmented soma
  for (V3DLONG z = 0; z < dimZ; z++) {
    for (V3DLONG y = 0; y < dimY; y++) {
      for (V3DLONG x = 0; x < dimX; x++) {
        V3DLONG idx = z * dimX * dimY + y * dimX + x;
        // Check if the voxel is part of the soma (value is 255 in the binary
        // segmentation)
        if (segData[idx] > 0) {
          // Calculate squared distance to the marker
          double dx = x - selectedMarker.x;
          double dy = y - selectedMarker.y;
          double dz = z - selectedMarker.z;
          double distSq = dx * dx + dy * dy + dz * dz;

          // Include this voxel if it's close to the selected marker (within
          // twice the radius)
          if (distSq <= 4 * radius * radius) {
            x_center += x;
            y_center += y;
            z_center += z;
            totalMass += 1;
          }
        }
      }
    }
  }

  // Compute center of mass
  if (totalMass > 0) {
    x_center /= totalMass;
    y_center /= totalMass;
    z_center /= totalMass;
  } else {
    // If no soma voxels were found, use the marker position
    x_center = selectedMarker.x;
    y_center = selectedMarker.y;
    z_center = selectedMarker.z;
    v3d_msg(
        "Warning: No soma voxels found near the marker. Using marker position "
        "as center.");
  }

  printf("Computed center of mass: (%f, %f, %f)\n", x_center, y_center,
         z_center);

  // Second pass: copy both the binary segmentation and intensities into arrays
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
