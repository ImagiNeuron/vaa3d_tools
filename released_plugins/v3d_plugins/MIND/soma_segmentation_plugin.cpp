/* soma_segmentation_plugin.cpp
 * This plugin segments individual somas with an improved 3D region growing
 * algorithm.
 *
 * The new implementation uses a manual intensity difference threshold (relative
 * to the seed) plus a gradient threshold to restrict growth. The image is first
 * smoothed and its gradient magnitude computed; during region growing (seeded
 * from a user-defined landmark) only voxels whose intensity lies within
 * [seedIntensity - diff, seedIntensity + diff] and whose gradient magnitude is
 * below a user-specified threshold are accepted. A spatial constraint (only
 * within a sphere of radius equal to the marker’s radius) is also enforced.
 *
 * The final segmented image is binary (white for soma, black for background).
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 * 2025-02-01 : Updated by Athmane for windowed segmentation view, saving the
 * segmentation as a TIFF file, prompting the user for parameters with input
 * boxes, and attempting to fix the segmentation algorithm (still problematic).
 */

#include "soma_segmentation_plugin.h"

#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <queue>
#include <vector>

#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// Helper function: 1D Gaussian kernel (centered, size = 2*radius+1)
static std::vector<double> gaussianKernel1D(int radius, double sigma) {
  int size = 2 * radius + 1;
  std::vector<double> kernel(size, 0.0);
  double sum = 0.0;
  double s2 = sigma * sigma;
  for (int i = -radius; i <= radius; ++i) {
    double val = exp(-(i * i) / (2 * s2));
    kernel[i + radius] = val;
    sum += val;
  }
  for (double &val : kernel) val /= sum;
  return kernel;
}

// 3D Gaussian smoothing (separable convolution)
// Assumes single-channel image (cdim == 1)
My4DImage *gaussianSmooth3D(const My4DImage *input, double sigma) {
  if (!input || !input->data || input->cdim != 1) return nullptr;

  V3DLONG xdim = input->xdim;
  V3DLONG ydim = input->ydim;
  V3DLONG zdim = input->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  int radius = static_cast<int>(ceil(3 * sigma));
  std::vector<double> kernel = gaussianKernel1D(radius, sigma);

  double *tempX = new double[totalSize]();
  double *tempY = new double[totalSize]();
  double *tempZ = new double[totalSize]();

  // Convolve along x
  for (V3DLONG z = 0; z < zdim; z++)
    for (V3DLONG y = 0; y < ydim; y++)
      for (V3DLONG x = 0; x < xdim; x++) {
        double sum = 0.0;
        for (int k = -radius; k <= radius; k++) {
          int xi = x + k;
          if (xi < 0)
            xi = 0;
          else if (xi >= xdim)
            xi = xdim - 1;
          sum +=
              kernel[k + radius] * input->data[z * ydim * xdim + y * xdim + xi];
        }
        tempX[z * ydim * xdim + y * xdim + x] = sum;
      }

  // Convolve along y
  for (V3DLONG z = 0; z < zdim; z++)
    for (V3DLONG x = 0; x < xdim; x++)
      for (V3DLONG y = 0; y < ydim; y++) {
        double sum = 0.0;
        for (int k = -radius; k <= radius; k++) {
          int yi = y + k;
          if (yi < 0)
            yi = 0;
          else if (yi >= ydim)
            yi = ydim - 1;
          sum += kernel[k + radius] * tempX[z * ydim * xdim + yi * xdim + x];
        }
        tempY[z * ydim * xdim + y * xdim + x] = sum;
      }

  // Convolve along z
  for (V3DLONG y = 0; y < ydim; y++)
    for (V3DLONG x = 0; x < xdim; x++)
      for (V3DLONG z = 0; z < zdim; z++) {
        double sum = 0.0;
        for (int k = -radius; k <= radius; k++) {
          int zi = z + k;
          if (zi < 0)
            zi = 0;
          else if (zi >= zdim)
            zi = zdim - 1;
          sum += kernel[k + radius] * tempY[zi * ydim * xdim + y * xdim + x];
        }
        tempZ[z * ydim * xdim + y * xdim + x] = sum;
      }

  My4DImage *smoothed = new My4DImage();
  smoothed->xdim = xdim;
  smoothed->ydim = ydim;
  smoothed->zdim = zdim;
  smoothed->cdim = 1;
  smoothed->data = new unsigned char[totalSize];
  for (V3DLONG i = 0; i < totalSize; i++) {
    int val = static_cast<int>(round(tempZ[i]));
    if (val < 0)
      val = 0;
    else if (val > 255)
      val = 255;
    smoothed->data[i] = static_cast<unsigned char>(val);
  }
  delete[] tempX;
  delete[] tempY;
  delete[] tempZ;
  return smoothed;
}

////////////////////////////////////////////////////////////////////////
// Compute gradient magnitude from a single-channel image using central
// differences.
My4DImage *computeGradientMagnitude(const My4DImage *input) {
  if (!input || !input->data || input->cdim != 1) return nullptr;
  V3DLONG xdim = input->xdim;
  V3DLONG ydim = input->ydim;
  V3DLONG zdim = input->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  My4DImage *gradMag = new My4DImage();
  gradMag->xdim = xdim;
  gradMag->ydim = ydim;
  gradMag->zdim = zdim;
  gradMag->cdim = 1;
  gradMag->data = new unsigned char[totalSize];

  for (V3DLONG z = 0; z < zdim; z++) {
    for (V3DLONG y = 0; y < ydim; y++) {
      for (V3DLONG x = 0; x < xdim; x++) {
        V3DLONG idx = z * ydim * xdim + y * xdim + x;
        double gx = 0, gy = 0, gz = 0;
        if (x > 0 && x < xdim - 1)
          gx = (input->data[idx + 1] - input->data[idx - 1]) / 2.0;
        if (y > 0 && y < ydim - 1)
          gy = (input->data[idx + xdim] - input->data[idx - xdim]) / 2.0;
        if (z > 0 && z < zdim - 1)
          gz = (input->data[idx + xdim * ydim] -
                input->data[idx - xdim * ydim]) /
               2.0;
        double mag = sqrt(gx * gx + gy * gy + gz * gz);
        if (mag > 255) mag = 255;
        gradMag->data[idx] = static_cast<unsigned char>(mag);
      }
    }
  }
  return gradMag;
}

////////////////////////////////////////////////////////////////////////
// Updated performFloodFill: now uses a relative (seed-based) intensity
// threshold. It computes lowerThreshold and upperThreshold from the seed
// intensity and a provided difference value, and also uses a gradient
// threshold.
My4DImage *performFloodFill(const My4DImage *image, const My4DImage *gradMag,
                            const MyMarker &marker, int intensityDiff,
                            int gradThreshold, int &min_x, int &max_x,
                            int &min_y, int &max_y, int &min_z, int &max_z) {
  if (!image || !image->data || !gradMag || !gradMag->data) return nullptr;

  V3DLONG xdim = image->xdim;
  V3DLONG ydim = image->ydim;
  V3DLONG zdim = image->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  My4DImage *segImage = new My4DImage();
  segImage->xdim = xdim;
  segImage->ydim = ydim;
  segImage->zdim = zdim;
  segImage->cdim = 1;
  segImage->data = new unsigned char[totalSize];
  memset(segImage->data, 0, totalSize * sizeof(unsigned char));

  int x0 = static_cast<int>(std::round(marker.x));
  int y0 = static_cast<int>(std::round(marker.y));
  int z0 = static_cast<int>(std::round(marker.z));
  if (x0 < 0 || x0 >= xdim || y0 < 0 || y0 >= ydim || z0 < 0 || z0 >= zdim) {
    delete segImage;
    return nullptr;
  }
  int seedIdx = z0 * ydim * xdim + y0 * xdim + x0;
  unsigned char seedIntensity = image->data[seedIdx];

  // Use the manual "intensityDiff" to set a relative threshold band.
  int lowerThreshold = std::max(0, (int)seedIntensity - intensityDiff);
  int upperThreshold = std::min(255, (int)seedIntensity + intensityDiff);

  size_t maxRegionSize =
      static_cast<size_t>((4.0 / 3.0) * M_PI * std::pow(marker.radius, 3));
  size_t minRegionSize = 5;

  struct Node {
    V3DLONG idx;
    double cost;
  };
  auto cmp = [](const Node &a, const Node &b) { return a.cost > b.cost; };
  std::priority_queue<Node, std::vector<Node>, decltype(cmp)> heap(cmp);
  heap.push({seedIdx, 0.0});

  std::vector<bool> visited(totalSize, false);
  min_x = x0;
  max_x = x0;
  min_y = y0;
  max_y = y0;
  min_z = z0;
  max_z = z0;
  size_t regionSize = 0;

  std::vector<int> neighborOffsets;
  for (int dz = -1; dz <= 1; dz++)
    for (int dy = -1; dy <= 1; dy++)
      for (int dx = -1; dx <= 1; dx++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        int offset = dz * ydim * xdim + dy * xdim + dx;
        neighborOffsets.push_back(offset);
      }

  while (!heap.empty()) {
    Node current = heap.top();
    heap.pop();
    if (visited[current.idx]) continue;
    visited[current.idx] = true;

    int cz = current.idx / (xdim * ydim);
    int rem = current.idx % (xdim * ydim);
    int cy = rem / xdim;
    int cx = rem % xdim;

    double distance = std::sqrt((cx - x0) * (cx - x0) + (cy - y0) * (cy - y0) +
                                (cz - z0) * (cz - z0));
    if (distance > marker.radius) continue;

    unsigned char currIntensity = image->data[current.idx];
    if (currIntensity < lowerThreshold || currIntensity > upperThreshold)
      continue;
    if (gradMag->data[current.idx] > gradThreshold) continue;

    segImage->data[current.idx] = 255;
    regionSize++;

    min_x = std::min(min_x, cx);
    max_x = std::max(max_x, cx);
    min_y = std::min(min_y, cy);
    max_y = std::max(max_y, cy);
    min_z = std::min(min_z, cz);
    max_z = std::max(max_z, cz);

    if (regionSize > maxRegionSize) break;

    for (int offset : neighborOffsets) {
      V3DLONG neighborIdx = current.idx + offset;
      int nz = neighborIdx / (xdim * ydim);
      int rem2 = neighborIdx % (xdim * ydim);
      int ny = rem2 / xdim;
      int nx = rem2 % xdim;
      if (nx < 0 || nx >= xdim || ny < 0 || ny >= ydim || nz < 0 || nz >= zdim)
        continue;
      if (visited[neighborIdx]) continue;
      unsigned char neighborVal = image->data[neighborIdx];
      if (neighborVal < lowerThreshold || neighborVal > upperThreshold)
        continue;
      if (gradMag->data[neighborIdx] > gradThreshold) continue;
      heap.push({neighborIdx, 0.0});
    }
  }

  if (regionSize < minRegionSize) {
    delete segImage;
    return nullptr;
  }
  return segImage;
}

////////////////////////////////////////////////////////////////////////
// Merge function: combine flood-filled images with a bitwise OR.
My4DImage *mergeFloodedImages(const QList<My4DImage *> &floodedImages) {
  if (floodedImages.isEmpty()) return nullptr;
  My4DImage *mergedImage = new My4DImage(*floodedImages[0]);
  V3DLONG totalSize = mergedImage->xdim * mergedImage->ydim * mergedImage->zdim;
  for (int i = 1; i < floodedImages.size(); ++i)
    for (V3DLONG idx = 0; idx < totalSize; ++idx)
      mergedImage->data[idx] |= floodedImages[i]->data[idx];
  for (V3DLONG idx = 0; idx < totalSize; ++idx)
    mergedImage->data[idx] = (mergedImage->data[idx] == 255 ? 255 : 0);
  return mergedImage;
}

////////////////////////////////////////////////////////////////////////
// Main reconstruction function (invoked from menu or command-line)
My4DImage *reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
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

  // Prompt user for intensity difference threshold (relative to the seed
  // intensity).
  bool okIntensity = false;
  int intensityDiff =
      QInputDialog::getInt(parent, "Intensity Difference Threshold",
                           "Enter allowed difference from seed intensity:", 20,
                           0, 255, 1, &okIntensity);
  if (!okIntensity) return nullptr;
  // Prompt user for gradient threshold.
  bool okGrad = false;
  int gradThreshold = QInputDialog::getInt(
      parent, "Gradient Threshold",
      "Enter maximum allowed gradient magnitude:", 20, 0, 255, 1, &okGrad);
  if (!okGrad) return nullptr;

  // Convert the LandmarkList to a QList of MyMarker objects.
  QList<MyMarker> myMarkerList;
  for (const auto &lm : landmarkList) {
    MyMarker marker;
    marker.x = lm.x;
    marker.y = lm.y;
    marker.z = lm.z;
    marker.radius = lm.radius;
    myMarkerList.append(marker);
  }

  // ============================
  // Processing steps:
  // 1. Create a global segmentation image (initialized to zero).
  // 2. Extract the desired channel image.
  // 3. Compute the smoothed image and gradient magnitude ONCE.
  // 4. Iterate over landmarks and add each flood-filled region onto the global
  // segmentation.
  // ============================

  V3DLONG totalSize =
      p4DImage->getXDim() * p4DImage->getYDim() * p4DImage->getZDim();

  // (1) Create a global segmentation image.
  My4DImage *globalSegImage = new My4DImage();
  globalSegImage->xdim = p4DImage->getXDim();
  globalSegImage->ydim = p4DImage->getYDim();
  globalSegImage->zdim = p4DImage->getZDim();
  globalSegImage->cdim = 1;
  globalSegImage->data = new unsigned char[totalSize];
  memset(globalSegImage->data, 0, totalSize * sizeof(unsigned char));

  // (2) Create a My4DImage holding the channel data for processing.
  My4DImage *channelImage = new My4DImage();
  channelImage->xdim = p4DImage->getXDim();
  channelImage->ydim = p4DImage->getYDim();
  channelImage->zdim = p4DImage->getZDim();
  channelImage->cdim = 1;
  channelImage->data = new unsigned char[totalSize];
  unsigned char *channelData = p4DImage->getRawDataAtChannel(PARA.channel - 1);
  if (!channelData) {
    v3d_msg("Failed to get image data for the specified channel.", bmenu);
    delete channelImage;
    return nullptr;
  }
  memcpy(channelImage->data, channelData, totalSize * sizeof(unsigned char));

  // (3) Compute a smoothed image and then its gradient magnitude.
  double sigma = 2.0;  // use a larger sigma for better smoothing
  My4DImage *smoothedImage = gaussianSmooth3D(channelImage, sigma);
  delete channelImage;
  if (!smoothedImage) {
    v3d_msg("Gaussian smoothing failed.", bmenu);
    delete globalSegImage;
    return nullptr;
  }

  My4DImage *gradMag = computeGradientMagnitude(smoothedImage);
  if (!gradMag) {
    v3d_msg("Gradient magnitude computation failed.", bmenu);
    delete smoothedImage;
    delete globalSegImage;
    return nullptr;
  }

  // (4) For each landmark, perform flood fill using the preprocessed images and
  // merge the result.
  for (const auto &marker : myMarkerList) {
    int min_x, max_x, min_y, max_y, min_z, max_z;
    My4DImage *localSeg = performFloodFill(smoothedImage, gradMag, marker,
                                           intensityDiff, gradThreshold, min_x,
                                           max_x, min_y, max_y, min_z, max_z);
    if (localSeg) {
      // Merge the current segmentation into the global segmentation.
      for (V3DLONG idx = 0; idx < totalSize; ++idx) {
        if (localSeg->data[idx] == 255) globalSegImage->data[idx] = 255;
      }
      delete localSeg;
    }
  }

  // Clean up the temporary images.
  delete smoothedImage;
  delete gradMag;

  // For diagnostic purposes, print the number of segmented (flooded) voxels.
  size_t flooded_voxels = 0;
  for (V3DLONG idx = 0; idx < totalSize; ++idx) {
    if (globalSegImage->data[idx] == 255) flooded_voxels++;
  }
  printf("Flooded voxels: %zu\n", flooded_voxels);

  // Open a new image window to display the global segmented image.
  if (globalSegImage) {
    Image4DSimple p4DImageSeg;
    p4DImageSeg.setData(globalSegImage->data, globalSegImage->xdim,
                        globalSegImage->ydim, globalSegImage->zdim, 1,
                        V3D_UINT8);
    // Detach the data pointer so that it isn’t deleted when globalSegImage is
    // destroyed.
    globalSegImage->data = nullptr;
    v3dhandle newwin = callback.newImageWindow();
    callback.setImage(newwin, &p4DImageSeg);
    callback.setImageName(newwin, QString("Segmented Image"));
    callback.updateImageWindow(newwin);
  }

  // (Commented out for now) Save the binary segmented image as a TIFF file.
  // QString savePath = QFileDialog::getSaveFileName(parent, "Save binary
  // segmented image", "",
  //                                                 "TIFF Files (*.tiff
  //                                                 *.tif)");
  // if (!savePath.isEmpty() && globalSegImage)
  // {
  //   V3DLONG outSZ[4] = {globalSegImage->xdim, globalSegImage->ydim,
  //   globalSegImage->zdim, 1}; simple_saveimage_wrapper(callback,
  //   savePath.toStdString().c_str(),
  //                            globalSegImage->data, outSZ, 1);
  //   v3d_msg("Binary segmented image saved.", bmenu);
  // }

  return globalSegImage;
}

////////////////////////////////////////////////////////////////////////
// Standard Vaa3D plugin interface implementations

QStringList SomaSegmentation::menulist() const {
  return QStringList() << tr("soma_segmentation") << tr("cellsegmentation") << tr("about");
}

QStringList SomaSegmentation::funclist() const {
  return QStringList() << tr("segment_somas") << tr("help");
}

void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback, QWidget *parent) {
  if (menu_name == tr("soma_segmentation")) {
    bool bmenu = true;
    input_PARA PARA;
    reconstruction_func(callback, parent, PARA, bmenu);
  } else if (menu_name == tr("cellsegmentation")) {
    cellSegmentation cellseg;
    cellseg.interface_run(callback, parent);
  } else {
    v3d_msg(tr("This plugin segments individual somas using a 3D "
               "region-growing algorithm "
               "with enhanced smoothing, relative intensity thresholding and "
               "gradient control, "
               "with landmarks as seeds. Developed by ImagiNeuron (2024-11-16) "
               "and updated for improved segmentation."),
            0);
  }
}

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
    } else
      PARA.inimg_file = infiles[0];
    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;
    reconstruction_func(callback, parent, PARA, bmenu);
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
  } else
    return false;
  return true;
}

////////////////////////////////////////////////////////////////////////
// Implementation of My4DImage methods

My4DImage::My4DImage() : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0) {}

My4DImage::My4DImage(const My4DImage &other) {
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

My4DImage &My4DImage::operator=(const My4DImage &other) {
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

My4DImage::~My4DImage() { delete[] data; }
