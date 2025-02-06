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
// Helper: Merge two binary segmentation images (logical OR)
void mergeSegmentations(My4DImage *globalSeg, const My4DImage *localSeg) {
  V3DLONG totalSize = globalSeg->xdim * globalSeg->ydim * globalSeg->zdim;
  for (V3DLONG i = 0; i < totalSize; i++) {
    if (localSeg->data[i] == 255) globalSeg->data[i] = 255;
  }
}

////////////////////////////////////////////////////////////////////////
// Segmentation routine using region growing.
// This function uses cell segmentation using the provided landmark as the seed
My4DImage *segmentSoma(const My4DImage *input, const MyMarker &marker) {
  // Create a binary output image of the same dimensions.
  My4DImage *seg = new My4DImage();
  seg->xdim = input->xdim;
  seg->ydim = input->ydim;
  seg->zdim = input->zdim;
  seg->cdim = 1;
  V3DLONG totalSize = input->xdim * input->ydim * input->zdim;
  seg->data = new unsigned char[totalSize];
  memset(seg->data, 0, totalSize);

  // Use the landmark as the seed.
  int seedX = static_cast<int>(round(marker.x));
  int seedY = static_cast<int>(round(marker.y));
  int seedZ = static_cast<int>(round(marker.z));
  if (seedX < 0 || seedX >= input->xdim || seedY < 0 || seedY >= input->ydim ||
      seedZ < 0 || seedZ >= input->zdim) {
    delete seg;
    return nullptr;
  }
  V3DLONG seedIdx =
      seedZ * input->ydim * input->xdim + seedY * input->xdim + seedX;
  unsigned char seedIntensity = input->data[seedIdx];

  // Set an intensity threshold.
  // (You may choose a more sophisticated parameter derivation here.)
  int intensityThreshold =
      seedIntensity - 20;  // example: allow intensities close to seed
  if (intensityThreshold < 0) intensityThreshold = 0;

  // Create a queue and a visited flag vector for region growing.
  std::queue<V3DLONG> q;
  vector<bool> visited(totalSize, false);

  // Begin at the seed.
  q.push(seedIdx);
  visited[seedIdx] = true;
  seg->data[seedIdx] = 255;

  // Define neighbor offsets for a 26-connected neighborhood.
  int dims[3] = {(int)input->xdim, (int)input->ydim, (int)input->zdim};
  vector<int> neighborOffsets;
  for (int dz = -1; dz <= 1; dz++) {
    for (int dy = -1; dy <= 1; dy++) {
      for (int dx = -1; dx <= 1; dx++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        int offset = dz * dims[0] * dims[1] + dy * dims[0] + dx;
        neighborOffsets.push_back(offset);
      }
    }
  }

  // Region growing loop.
  while (!q.empty()) {
    V3DLONG idx = q.front();
    q.pop();

    // Recover (x,y,z) coordinate from linear index.
    int z = idx / (dims[0] * dims[1]);
    int rem = idx % (dims[0] * dims[1]);
    int y = rem / dims[0];
    int x = rem % dims[0];

    // Enforce a spatial constraint: do not grow beyond the landmark radius.
    double distance =
        sqrt((x - marker.x) * (x - marker.x) + (y - marker.y) * (y - marker.y) +
             (z - marker.z) * (z - marker.z));
    if (distance > marker.radius) continue;

    // Examine all 26 neighbors.
    for (auto offset : neighborOffsets) {
      V3DLONG nIdx = idx + offset;
      // Compute neighbor coordinate
      int nz = nIdx / (dims[0] * dims[1]);
      int rem2 = nIdx % (dims[0] * dims[1]);
      int ny = rem2 / dims[0];
      int nx = rem2 % dims[0];

      // Check bounds.
      if (nx < 0 || nx >= dims[0] || ny < 0 || ny >= dims[1] || nz < 0 ||
          nz >= dims[2])
        continue;
      if (visited[nIdx]) continue;
      visited[nIdx] = true;

      // If the neighbor’s intensity meets the criteria, include it.
      if (input->data[nIdx] >= intensityThreshold) {
        seg->data[nIdx] = 255;
        q.push(nIdx);
      }
    }
  }
  return seg;
}

////////////////////////////////////////////////////////////////////////
// Main reconstruction function
// Iterate over the given landmark list and call our segmentation for each soma.
My4DImage *reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                               input_PARA &PARA, bool bmenu) {
  // (1) Get the current image window.
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

  // (2) Get image dimensions and raw data.
  V3DLONG N = p4DImage->getXDim();
  V3DLONG M = p4DImage->getYDim();
  V3DLONG P = p4DImage->getZDim();
  V3DLONG sc = p4DImage->getCDim();
  V3DLONG totalSize = N * M * P;
  unsigned char *data1d = p4DImage->getRawData();

  // (3) Get the list of landmarks (soma seed markers).
  LandmarkList landmarkList = callback.getLandmark(curwin);
  if (landmarkList.isEmpty()) {
    v3d_msg("No soma landmarks provided. Please mark at least one soma.",
            bmenu);
    return nullptr;
  }

  // (4) Choose the channel to process.
  int channel = 1;
  if (sc > 1) {
    bool ok;
    channel = QInputDialog::getInt(parent, "Channel",
                                   "Enter channel number:", 1, 1, sc, 1, &ok);
    if (!ok) return nullptr;
  }

  // (5) Extract the desired channel image into a My4DImage.
  My4DImage *channelImage = new My4DImage();
  channelImage->xdim = N;
  channelImage->ydim = M;
  channelImage->zdim = P;
  channelImage->cdim = 1;
  channelImage->data = new unsigned char[totalSize];
  unsigned char *channelData = p4DImage->getRawDataAtChannel(channel - 1);
  memcpy(channelImage->data, channelData, totalSize * sizeof(unsigned char));

  // (6) Create a global segmentation image (binary; initialized to zero).
  My4DImage *globalSeg = new My4DImage();
  globalSeg->xdim = N;
  globalSeg->ydim = M;
  globalSeg->zdim = P;
  globalSeg->cdim = 1;
  globalSeg->data = new unsigned char[totalSize];
  memset(globalSeg->data, 0, totalSize * sizeof(unsigned char));

  // (7) For each landmark, perform segmentation and merge.
  for (int i = 0; i < landmarkList.count(); i++) {
    MyMarker marker;
    marker.x = landmarkList.at(i).x;
    marker.y = landmarkList.at(i).y;
    marker.z = landmarkList.at(i).z;
    marker.radius = landmarkList.at(i).radius;  // if defined

    My4DImage *localSeg = segmentSoma(channelImage, marker);
    if (localSeg) {
      mergeSegmentations(globalSeg, localSeg);
      delete localSeg;
    }
  }

  // (8) Clean up temporary images.
  delete channelImage;

  // (9) Display the global segmentation result in a new window.
  Image4DSimple p4DImageSeg;
  p4DImageSeg.setData(globalSeg->data, N, M, P, 1, V3D_UINT8);
  v3dhandle newwin = callback.newImageWindow();
  callback.setImage(newwin, &p4DImageSeg);
  callback.setImageName(newwin, QString("Soma Segmentation"));
  callback.updateImageWindow(newwin);

  // Detach the pointer so that it isn’t deleted when globalSeg is destroyed.
  globalSeg->data = nullptr;
  return globalSeg;
}

////////////////////////////////////////////////////////////////////////
// Standard Vaa3D plugin interface implementations

QStringList SomaSegmentation::menulist() const {
  return QStringList() << tr("soma_segmentation") << tr("about");
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
