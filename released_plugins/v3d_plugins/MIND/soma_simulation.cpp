/**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */
#include "soma_simulation.h"

/**
 * @brief Function to create a background image for the current image
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param dim_X - the X dimension of the image
 * @param dim_Y - the Y dimension of the image
 * @param dim_Z - the Z dimension of the image
 * @return - a pointer to the background image data
 */
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

/**
 * @brief Create a background image based on segmentation threshold
 * @param callback V3DPluginCallback2 reference
 * @param parent Parent widget
 */
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

/**
 * @brief Get blended distribution parameters for a voxel based on its
 * location
 *
 * @param x X coordinate of the voxel
 * @param y Y coordinate of the voxel
 * @param z Z coordinate of the voxel
 * @param chunkStats 3D vector containing mean and stdDev pairs for each chunk
 * @param chunk_X Chunk size in X dimension
 * @param chunk_Y Chunk size in Y dimension
 * @param chunk_Z Chunk size in Z dimension
 * @param num_chunks_X Number of chunks in X dimension
 * @param num_chunks_Y Number of chunks in Y dimension
 * @param num_chunks_Z Number of chunks in Z dimension
 * @param blendRadius Radius for blending (in chunk units)
 * @param blendedMean Output parameter for the blended mean value
 * @param blendedStdDev Output parameter for the blended standard deviation
 */
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

/**
 * @brief overlay the ground truth data on the new simulated image. Dimensions
 * calculated based on current image
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param binarySegImage - the binary segmentation image
 * @param gradientImage - the gradient image
 * @param simulatedImage - the simulated image
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
 */
void simulate_soma_data(V3DPluginCallback2 &callback, QWidget *parent) {
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
            .arg(segPcaFileName));
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
    printf("Synthetic soma generation canceled.");
    delete[] segData;
    return;
  }

  printf("\nGenerating %d synthetic somas...\n", numSynthetic);

  int successfulPlacements = 0;

  // Initialize simulatedLandmarks list
  LandmarkList simulatedLandmarks;

  // Helper function to check if two somas overlap
  auto somasOverlap = [](const LocationSimple &s1,
                         const LocationSimple &s2) -> bool {
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
    LocationSimple newSoma;
    newSoma.radius = round(radius);
    newSoma.name = qPrintable(QString("Sim_%1").arg(i + 1));
    newSoma.comments = "";
    newSoma.shape = pxSphere;

    bool validPosition = false;
    int maxAttempts = 100;
    int attempts = 0;

    while (!validPosition && attempts < maxAttempts) {
      attempts++;
      validPosition = true;

      // Generate new potential position
      for (int j = 0; j < 3; j++) {
        std::normal_distribution<> d(meanCenter[j], stdCenter[j]);
        double coord = round(d(gen));

        // Store coordinate in new soma
        if (j == 0) {
          newSoma.x = coord;
          // Check if position is within image boundaries
          if (coord < boundaryMargin || coord > xDim - boundaryMargin) {
            validPosition = false;
            break;
          }
        } else if (j == 1) {
          newSoma.y = coord;
          // Check if position is within image boundaries
          if (coord < boundaryMargin || coord > yDim - boundaryMargin) {
            validPosition = false;
            break;
          }
        } else if (j == 2) {
          newSoma.z = coord;
          // Check if position is within image boundaries
          if (coord < boundaryMargin || coord > zDim - boundaryMargin) {
            validPosition = false;
            break;
          }
        }
      }

      // If position is within boundaries, check for overlap with existing somas
      if (validPosition) {
        // Check against all previously placed somas
        for (int s = 0; s < simulatedLandmarks.size(); s++) {
          if (somasOverlap(newSoma, simulatedLandmarks[s])) {
            validPosition = false;
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

    // Add this soma to our landmarks list
    simulatedLandmarks.append(newSoma);
    successfulPlacements++;

    // Generate random PCA values based on the distribution
    double randomVec1[3], randomVec2[3], randomVec3[3];

    Eigen::Matrix3d A;
    for (int col = 0; col < 3; ++col) {
      for (int row = 0; row < 3; ++row) {
        std::normal_distribution<> dist(meanEigenvectors[col * 3 + row],
                                        stdEigenvectors[col * 3 + row]);
        A(row, col) = dist(gen);
      }
    }

    // Perform QR decomposition to obtain orthogonal vectors
    Eigen::HouseholderQR<Eigen::Matrix3d> qr(A);
    Eigen::Matrix3d Q = qr.householderQ();

    for (int i = 0; i < 3; ++i) {
      randomVec1[i] = Q(i, 0);
      randomVec2[i] = Q(i, 1);
      randomVec3[i] = Q(i, 2);
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
    int centerX = V3DLONG(newSoma.x);
    int centerY = V3DLONG(newSoma.y);
    int centerZ = V3DLONG(newSoma.z);

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

  // Generate timestamp for folder name
  QDateTime currentTime = QDateTime::currentDateTime();
  QString timestamp = currentTime.toString("yyyyMMdd_hhmmss");

  // Create timestamped directory to store all output files
  QString outputDirPath =
      modifyFilePathForTeraFly(imageName) + "_simulation_" + timestamp;
  QDir outputDir(outputDirPath);
  if (!outputDir.exists()) {
    outputDir.mkpath(".");
    printf("Created output directory: %s\n",
           outputDirPath.toStdString().c_str());
  }

  // Define output filenames in the new directory
  QString outSegFileName = outputDirPath + "/simulated_segmentation.tif";
  QString outIntensityFileName = outputDirPath + "/simulated_intensity.tif";

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

  /*
   * Save simulated soma landmarks as a marker file
   */

  QString markerFileName = outputDirPath + "/simulated_landmarks.marker";

  FILE *fp = fopen(markerFileName.toStdString().c_str(), "w");
  if (fp) {
    // Write header
    fprintf(fp, "#x, y, z, radius, shape, name, comment\n");

    // Write each landmark
    for (int i = 0; i < simulatedLandmarks.size(); i++) {
      // Format: x,y,z,radius,shape,name,comment
      fprintf(fp, "%ld,%ld,%ld,%ld,%ld,%s,%s\n",
              V3DLONG(simulatedLandmarks.at(i).x),
              V3DLONG(simulatedLandmarks.at(i).y),
              V3DLONG(simulatedLandmarks.at(i).z),
              V3DLONG(simulatedLandmarks.at(i).radius),
              V3DLONG(simulatedLandmarks.at(i).shape),
              simulatedLandmarks.at(i).name.c_str(),
              simulatedLandmarks.at(i).comments.c_str());
    }

    fclose(fp);
  } else {
    printf("Error: Could not save marker file: %s\n",
           markerFileName.toStdString().c_str());
  }

  /*
   * Perform PCA analysis on the simulated somas
   */

  QString pcaSimulatedFileName =
      outputDirPath + "/pca_simulated_segmentation.csv";

  // Perform PCA analysis on each simulated soma
  for (int i = 0; i < simulatedLandmarks.size(); i++) {
    analyzeSomaPCA(outSegData, xDim, yDim, zDim, simulatedLandmarks[i], i + 1,
                   pcaSimulatedFileName);
  }

  v3d_msg(QString("Simulation complete:\nGenerated %1/%2 synthetic "
                  "somas. Files saved to directory: %3\n")
              .arg(successfulPlacements)
              .arg(numSynthetic)
              .arg(outputDirPath));

  delete[] segData;

  /*
   * Open new windows and display the synthetic soma data
   */

  unsigned char *gradientImage = new unsigned char[totalSize];
  cellSegmentation cellSeg;
  cellSeg.sobel3D(outSegData, gradientImage, xDim, yDim, zDim);

  // overlay
  overlaySimulation(callback, parent, outSegData, gradientImage,
                    outIntensityData);
  // // Create and show new window with binary simulated data. Now Obsolete
  // Create an image for the binary and realistic simulation data
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