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
                                     p4DImage->getZDim(), 3,
                                     V3D_UINT8);
  pcaVisualization->setOriginX(p4DImage->getOriginX());
  pcaVisualization->setOriginY(p4DImage->getOriginY());
  pcaVisualization->setOriginZ(p4DImage->getOriginZ());
  pcaVisualization->setRezX(p4DImage->getRezX());
  pcaVisualization->setRezY(p4DImage->getRezY());
  pcaVisualization->setRezZ(p4DImage->getRezZ());

  // copy original image to channel 1
  memcpy(pcaVisualization->getRawData(), p4DImage->getRawData(),
         p4DImage->getTotalUnitNumber() * p4DImage->getUnitBytes());

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
        // case 1,2,3: // marker position
        // case 4: // radius
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
    if (pc1 > pc2 && pc1 > pc3) {
      drawLine(pcaVisualization, 1, center, vec1Pos);
      drawLine(pcaVisualization, 2, center, vec2Pos);
      drawLine(pcaVisualization, 2, center, vec3Pos);
    } else if (pc2 > pc1 && pc2 > pc3) {
      drawLine(pcaVisualization, 1, center, vec2Pos);
      drawLine(pcaVisualization, 2, center, vec1Pos);
      drawLine(pcaVisualization, 2, center, vec3Pos);
    } else {
      drawLine(pcaVisualization, 1, center, vec3Pos);
      drawLine(pcaVisualization, 2, center, vec1Pos);
      drawLine(pcaVisualization, 2, center, vec2Pos);
    }

  }

  inFile.close();

  // Show the PCA visualization
  v3dhandle newwin = callback.newImageWindow();
  callback.setImage(newwin, pcaVisualization);
}

void drawLine(Image4DSimple *image, int channel, double *from, double *to) {
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
      imgData[channel * image->getXDim() * image->getYDim() * image->getZDim() +
              z * image->getXDim() * image->getYDim() + y * image->getXDim() + x] = 255;
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
