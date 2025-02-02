/* soma_segmentation_plugin.cpp
 *
 * This plugin segments individual somas using a 3D region growing algorithm.
 * It performs a custom smoothing (using a weighted-average filter), then
 * applies adaptive thresholding with a dynamic region–mean constraint to
 * prevent overflooding.
 *
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
// --- Smoothing code using a custom weighted-average filter ---
//
// We first provide helper functions to allocate/free a 3D float array.
static float ***Falloc3d(int x_size, int y_size, int z_size) {
  float ***array = new float **[z_size];
  for (int z = 0; z < z_size; z++) {
    array[z] = new float *[x_size];
    for (int x = 0; x < x_size; x++) {
      array[z][x] = new float[y_size];
      for (int y = 0; y < y_size; y++) array[z][x][y] = 0.0f;
    }
  }
  return array;
}

static void Ffree3d(float ***array, int z_size, int x_size) {
  for (int z = 0; z < z_size; z++) {
    for (int x = 0; x < x_size; x++) delete[] array[z][x];
    delete[] array[z];
  }
  delete[] array;
}

// The provided smoothing routine. The filter uses a weighted average
// with 0.4 weight for the center voxel and 0.1 for each of its 6 direct
// neighbors.
void smooth_image(float ***image, int ITER, int x_size, int y_size,
                  int z_size) {
  int x, y, z, i;
  float ***simage = Falloc3d(x_size, y_size, z_size);

  for (i = 0; i < ITER; i++) {
    // Zero the temporary image
    for (z = 0; z < z_size; z++)
      for (x = 0; x < x_size; x++)
        for (y = 0; y < y_size; y++) simage[z][x][y] = 0;

    // Apply smoothing (avoid borders)
    for (z = 1; z < z_size - 1; z++)
      for (x = 1; x < x_size - 1; x++)
        for (y = 1; y < y_size - 1; y++) {
          simage[z][x][y] = 0.4f * image[z][x][y] +
                            0.1f * (image[z - 1][x][y] + image[z + 1][x][y] +
                                    image[z][x - 1][y] + image[z][x + 1][y] +
                                    image[z][x][y - 1] + image[z][x][y + 1]);
        }
    // Copy back smoothed values.
    for (z = 0; z < z_size; z++)
      for (x = 0; x < x_size; x++)
        for (y = 0; y < y_size; y++) image[z][x][y] = simage[z][x][y];
  }

  Ffree3d(simage, z_size, x_size);
}

////////////////////////////////////////////////////////////////////////
// --- Helper: Smooth a My4DImage using the above routine ---
//
// This function converts the 1D unsigned char image into a 3D float array,
// applies the smoothing filter, and then converts back.
My4DImage *smoothMy4DImage(const My4DImage *input, int ITER) {
  if (!input || !input->data || input->cdim != 1) return nullptr;

  V3DLONG xdim = input->xdim, ydim = input->ydim, zdim = input->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  // Allocate a 3D float array and convert the data.
  float ***fimg = Falloc3d(xdim, ydim, zdim);
  for (int z = 0; z < zdim; z++)
    for (int x = 0; x < xdim; x++)
      for (int y = 0; y < ydim; y++) {
        int idx = z * ydim * xdim + x * ydim + y;
        fimg[z][x][y] = static_cast<float>(input->data[idx]);
      }

  // Smooth the image.
  smooth_image(fimg, ITER, xdim, ydim, zdim);

  // Create a new My4DImage and convert the float values back to unsigned char.
  My4DImage *smoothed = new My4DImage();
  smoothed->xdim = xdim;
  smoothed->ydim = ydim;
  smoothed->zdim = zdim;
  smoothed->cdim = 1;
  smoothed->data = new unsigned char[totalSize];
  for (V3DLONG z = 0; z < zdim; z++)
    for (V3DLONG x = 0; x < xdim; x++)
      for (V3DLONG y = 0; y < ydim; y++) {
        int idx = z * ydim * xdim + x * ydim + y;
        int val = static_cast<int>(round(fimg[z][x][y]));
        val = std::min(255, std::max(0, val));
        smoothed->data[idx] = static_cast<unsigned char>(val);
      }

  Ffree3d(fimg, zdim, xdim);
  return smoothed;
}

////////////////////////////////////////////////////////////////////////
// --- Helper: Compute local Otsu threshold over a cube ---
//
// This function computes a threshold over a cube of half–size “radius”
// centered at (x0,y0,z0). This adaptive threshold is used in region growing.
static int localOtsuThreshold(const My4DImage *image, int x0, int y0, int z0,
                              int radius) {
  V3DLONG xdim = image->xdim, ydim = image->ydim, zdim = image->zdim;
  int startX = std::max(0, x0 - radius);
  int endX = std::min(static_cast<int>(xdim) - 1, x0 + radius);
  int startY = std::max(0, y0 - radius);
  int endY = std::min(static_cast<int>(ydim) - 1, y0 + radius);
  int startZ = std::max(0, z0 - radius);
  int endZ = std::min(static_cast<int>(zdim) - 1, z0 + radius);

  vector<int> histogram(256, 0);
  int count = 0;
  for (int z = startZ; z <= endZ; z++)
    for (int x = startX; x <= endX; x++)
      for (int y = startY; y <= endY; y++) {
        int idx = z * ydim * xdim + x * ydim + y;
        histogram[image->data[idx]]++;
        count++;
      }

  double sum = 0;
  for (int t = 0; t < 256; t++) sum += t * histogram[t];

  double sumB = 0;
  int wB = 0, wF = 0;
  double maxVar = 0;
  int threshold = 0;
  for (int t = 0; t < 256; t++) {
    wB += histogram[t];
    if (wB == 0) continue;
    wF = count - wB;
    if (wF == 0) break;
    sumB += t * histogram[t];
    double mB = sumB / static_cast<double>(wB);
    double mF = (sum - sumB) / static_cast<double>(wF);
    double varBetween = static_cast<double>(wB) * wF * (mB - mF) * (mB - mF);
    if (varBetween > maxVar) {
      maxVar = varBetween;
      threshold = t;
    }
  }
  return threshold;
}

////////////////////////////////////////////////////////////////////////
// --- Improved Region Growing with Adaptive Threshold and Region Mean ---
//
// In this version we maintain a running average (regionMean) of the accepted
// voxels. A candidate voxel is accepted only if its intensity is at least
//     max( localOtsuThreshold, regionMean - margin ).
My4DImage *performFloodFill(const My4DImage *image, const MyMarker &marker,
                            int /*unused*/, int &min_x, int &max_x, int &min_y,
                            int &max_y, int &min_z, int &max_z) {
  if (!image || !image->data) return nullptr;

  V3DLONG xdim = image->xdim, ydim = image->ydim, zdim = image->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  // Create an empty binary image.
  My4DImage *segImage = new My4DImage();
  segImage->xdim = xdim;
  segImage->ydim = ydim;
  segImage->zdim = zdim;
  segImage->cdim = 1;
  segImage->data = new unsigned char[totalSize];
  memset(segImage->data, 0, totalSize);

  // Use the marker location (rounded) as the seed.
  int x0 = static_cast<int>(std::round(marker.x));
  int y0 = static_cast<int>(std::round(marker.y));
  int z0 = static_cast<int>(std::round(marker.z));
  if (x0 < 0 || x0 >= xdim || y0 < 0 || y0 >= ydim || z0 < 0 || z0 >= zdim) {
    delete segImage;
    return nullptr;
  }
  int seedIdx = z0 * ydim * xdim + x0 * ydim + y0;
  unsigned char seedValue = image->data[seedIdx];

  // Compute a local Otsu threshold over a cube (half–size = marker.radius).
  int localThresh =
      localOtsuThreshold(image, x0, y0, z0, static_cast<int>(marker.radius));

  // Set maximum allowed region size (roughly the sphere volume).
  size_t maxRegionSize =
      static_cast<size_t>((4.0 / 3.0) * M_PI * std::pow(marker.radius, 3));
  size_t minRegionSize = 5;  // Reject too small regions

  // For improved control, we maintain the running sum and mean of accepted
  // voxels.
  double regionSum = 0;
  int regionCount = 0;
  double regionMean = seedValue;  // initial guess
  const int margin = 10;          // margin subtracted from regionMean

  // Define the node for the priority queue.
  struct Node {
    V3DLONG idx;
    double cost;
  };
  auto cmp = [](const Node &a, const Node &b) { return a.cost > b.cost; };
  std::priority_queue<Node, vector<Node>, decltype(cmp)> heap(cmp);
  heap.push({seedIdx, 0.0});

  vector<bool> visited(totalSize, false);

  // Initialize bounding box.
  min_x = x0;
  max_x = x0;
  min_y = y0;
  max_y = y0;
  min_z = z0;
  max_z = z0;
  size_t regionSize = 0;

  // Precompute 26-neighborhood offsets.
  vector<int> neighborOffsets;
  for (int dz = -1; dz <= 1; dz++)
    for (int dx = -1; dx <= 1; dx++)
      for (int dy = -1; dy <= 1; dy++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        int offset = dz * ydim * xdim + dx * ydim + dy;
        neighborOffsets.push_back(offset);
      }

  // Region growing loop.
  while (!heap.empty()) {
    Node current = heap.top();
    heap.pop();
    if (visited[current.idx]) continue;
    visited[current.idx] = true;

    // Compute voxel coordinates.
    int z = current.idx / (xdim * ydim);
    int rem = current.idx % (xdim * ydim);
    int x = rem / ydim;
    int y = rem % ydim;

    // Enforce spatial constraint: only grow within a sphere of radius
    // marker.radius.
    double distance = std::sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) +
                                (z - z0) * (z - z0));
    if (distance > marker.radius) continue;

    unsigned char currIntensity = image->data[current.idx];
    // For all voxels except the very first, require that the intensity be at
    // least the higher of the local Otsu threshold and (regionMean - margin).
    if (regionCount > 0 &&
        currIntensity <
            std::max(localThresh, static_cast<int>(regionMean - margin)))
      continue;

    // Accept voxel.
    segImage->data[current.idx] = 255;
    regionSize++;
    // Update region running sum and mean.
    regionSum += currIntensity;
    regionCount++;
    regionMean = regionSum / regionCount;

    // Update bounding box.
    min_x = std::min(min_x, x);
    max_x = std::max(max_x, x);
    min_y = std::min(min_y, y);
    max_y = std::max(max_y, y);
    min_z = std::min(min_z, z);
    max_z = std::max(max_z, z);

    if (regionSize > maxRegionSize) break;

    // Process 26 neighbors.
    for (int offset : neighborOffsets) {
      V3DLONG nIdx = current.idx + offset;
      int nz = nIdx / (xdim * ydim);
      int rem2 = nIdx % (xdim * ydim);
      int nx = rem2 / ydim;
      int ny = rem2 % ydim;
      if (nx < 0 || nx >= xdim || ny < 0 || ny >= ydim || nz < 0 || nz >= zdim)
        continue;
      if (visited[nIdx]) continue;
      unsigned char neighborVal = image->data[nIdx];
      if (neighborVal <
          std::max(localThresh, static_cast<int>(regionMean - margin)))
        continue;
      double cost = std::abs(neighborVal - seedValue);
      heap.push({nIdx, cost});
    }
  }

  if (regionSize < minRegionSize) {
    delete segImage;
    return nullptr;
  }
  return segImage;
}

////////////////////////////////////////////////////////////////////////
// --- Merge function (bitwise OR) for multiple flood–filled images ---
My4DImage *mergeFloodedImages(const QList<My4DImage *> &floodedImages) {
  if (floodedImages.isEmpty()) return nullptr;
  My4DImage *merged = new My4DImage(*floodedImages[0]);
  V3DLONG totalSize = merged->xdim * merged->ydim * merged->zdim;
  for (int i = 1; i < floodedImages.size(); i++)
    for (V3DLONG idx = 0; idx < totalSize; idx++)
      merged->data[idx] |= floodedImages[i]->data[idx];
  for (V3DLONG idx = 0; idx < totalSize; idx++)
    merged->data[idx] = (merged->data[idx] == 255 ? 255 : 0);
  return merged;
}

////////////////////////////////////////////////////////////////////////
// --- Main Reconstruction Function ---
struct input_PARA {
  QString inimg_file;
  V3DLONG channel;
};

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu) {
  unsigned char *data1d = nullptr;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];

  if (bmenu) {
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin) {
      QMessageBox::information(0, "", "No image is open in the main window.");
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
    if (sc == 1) {
      c = 1;
    } else {
      bool ok;
      c = QInputDialog::getInt(parent, "Channel", "Enter channel number:", 1, 1,
                               sc, 1, &ok);
      if (!ok) return;
    }
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
      fprintf(stderr, "Error reading image file [%s].\n",
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

  v3dhandle curwin = callback.currentImageWindow();
  if (!curwin) {
    v3d_msg("No image window is open.", bmenu);
    return;
  }
  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("Invalid image pointer.", bmenu);
    return;
  }
  LandmarkList landmarkList = callback.getLandmark(curwin);
  if (landmarkList.isEmpty()) {
    v3d_msg("No landmarks defined. Please define at least one.", bmenu);
    return;
  }

  // Convert landmarks to our MyMarker list.
  QList<MyMarker> myMarkerList;
  for (const auto &lm : landmarkList) {
    MyMarker marker;
    marker.x = lm.x;
    marker.y = lm.y;
    marker.z = lm.z;
    marker.radius = lm.radius;
    myMarkerList.append(marker);
  }

  QList<My4DImage *> floodedImages;
  int dummy_threshold = 2;  // Unused now.

  // For each marker, perform segmentation.
  for (const auto &marker : myMarkerList) {
    My4DImage *tempImage = new My4DImage();
    tempImage->xdim = p4DImage->getXDim();
    tempImage->ydim = p4DImage->getYDim();
    tempImage->zdim = p4DImage->getZDim();
    tempImage->cdim = 1;
    V3DLONG totalSize = tempImage->xdim * tempImage->ydim * tempImage->zdim;
    tempImage->data = new unsigned char[totalSize];

    unsigned char *channelData =
        p4DImage->getRawDataAtChannel(PARA.channel - 1);
    if (!channelData) {
      v3d_msg("Failed to retrieve channel data.", bmenu);
      delete tempImage;
      continue;
    }
    memcpy(tempImage->data, channelData, totalSize);

    // Smooth the image using our custom filter (adjust ITER as needed).
    My4DImage *smoothedImage = smoothMy4DImage(tempImage, 2);
    delete tempImage;
    if (!smoothedImage) {
      v3d_msg("Smoothing failed.", bmenu);
      continue;
    }

    int min_x, max_x, min_y, max_y, min_z, max_z;
    My4DImage *flooded =
        performFloodFill(smoothedImage, marker, dummy_threshold, min_x, max_x,
                         min_y, max_y, min_z, max_z);
    delete smoothedImage;
    if (flooded) floodedImages.append(flooded);
  }

  // Merge all flood-filled regions.
  My4DImage *finalFloodedImage = mergeFloodedImages(floodedImages);
  if (finalFloodedImage) {
    printf("Final image: %ld x %ld x %ld, with %ld voxels segmented.\n",
           finalFloodedImage->xdim, finalFloodedImage->ydim,
           finalFloodedImage->zdim,
           std::count(finalFloodedImage->data,
                      finalFloodedImage->data + finalFloodedImage->xdim *
                                                    finalFloodedImage->ydim *
                                                    finalFloodedImage->zdim,
                      (unsigned char)255));
  } else {
    printf("No region was segmented.\n");
  }

  for (auto img : floodedImages) delete img;

  // Save final binary segmentation.
  QString savePath = QFileDialog::getSaveFileName(
      parent, "Save segmented image", "", "TIFF Files (*.tiff *.tif)");
  if (!savePath.isEmpty()) {
    V3DLONG outSZ[4] = {finalFloodedImage->xdim, finalFloodedImage->ydim,
                        finalFloodedImage->zdim, 1};
    simple_saveimage_wrapper(callback, savePath.toStdString().c_str(),
                             finalFloodedImage->data, outSZ, 1);
    v3d_msg("Segmented image saved.", bmenu);
  }

  delete finalFloodedImage;
  v3d_msg("Segmentation complete.", bmenu);
  return;
}

////////////////////////////////////////////////////////////////////////
// --- Plugin interface implementations ---
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
    v3d_msg(
        tr("This plugin segments somas using region growing with custom "
           "smoothing "
           "and adaptive thresholding (with a dynamic region–mean check) to "
           "prevent overflooding.\n"
           "Developed by ImagiNeuron and updated for improved segmentation."),
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
        (input.size() >= 1) ? (vector<char *> *)input[0].p : nullptr;
    vector<char *> *pparas =
        (input.size() >= 2) ? (vector<char *> *)input[1].p : nullptr;
    vector<char *> infiles = (pinfiles ? *pinfiles : vector<char *>());
    vector<char *> paras = (pparas ? *pparas : vector<char *>());

    if (infiles.empty()) {
      fprintf(stderr, "Need input image.\n");
      return false;
    }
    PARA.inimg_file = infiles[0];
    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;
    reconstruction_func(callback, parent, PARA, bmenu);
  } else if (func_name == tr("help")) {
    printf("Usage: \n");
    printf(
        "vaa3d -x soma_segmentation -f segment_somas -i <inimg_file> -p "
        "<channel>\n");
    printf(
        "The output is a binary segmented image (white somas, black "
        "background).\n");
  } else
    return false;

  return true;
}

////////////////////////////////////////////////////////////////////////
// --- Implementation of My4DImage methods ---
My4DImage::My4DImage() : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0) {}

My4DImage::My4DImage(const My4DImage &other) {
  xdim = other.xdim;
  ydim = other.ydim;
  zdim = other.zdim;
  cdim = other.cdim;
  if (other.data) {
    data = new unsigned char[xdim * ydim * zdim * cdim];
    std::copy(other.data, other.data + xdim * ydim * zdim * cdim, data);
  } else
    data = nullptr;
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
  } else
    data = nullptr;
  return *this;
}

My4DImage::~My4DImage() { delete[] data; }
