/* soma_segmentation_plugin.cpp
 * This plugin segments individual somas using a marker–controlled 3D watershed.
 *
 * The new implementation computes a smoothed image and its gradient magnitude.
 * Then, using the user–defined landmarks as markers (with a spatial constraint
 * given by each marker’s radius), a bucket–based watershed algorithm is applied
 * to “flood” the gradient image. When two marker regions meet the voxel is set
 * as a watershed (boundary). Finally the segmentation is binarized (white soma,
 * black background).
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 * 2025-02-01 : Updated by Athmane for marker–controlled watershed segmentation.
 */

#include "soma_segmentation_plugin.h"

#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <algorithm>
#include <cassert>
#include <cfloat>  // For FLT_MAX
#include <cmath>
#include <cstring>
#include <queue>
#include <tuple>
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
// Assumes single–channel image (cdim == 1)
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
// Compute gradient magnitude from a single–channel image using central
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
// Optimized: Marker–controlled 3D watershed segmentation using bucket–based
// flooding This version uses 256 buckets (one per gradient level) and squared
// distances to avoid expensive sqrt calls. It also stores labels in a vector of
// short values.
static My4DImage *markerControlledWatershed(const My4DImage *gradMag,
                                            const QList<MyMarker> &markers) {
  if (!gradMag || !gradMag->data) return nullptr;
  int xdim = gradMag->xdim;
  int ydim = gradMag->ydim;
  int zdim = gradMag->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  // Define label codes using a 16-bit type.
  const short UNVISITED = -1;
  const short WSHED = 0;
  vector<short> labels(totalSize, UNVISITED);
  // Distance (squared) from the marker seed for each voxel.
  vector<float> dist(totalSize, FLT_MAX);

  // Precompute marker information.
  // (Each marker gets a unique label: marker label = index+1)
  struct MarkerInfo {
    int mx, my, mz;
    float r2;  // squared radius
  };
  vector<MarkerInfo> markerInfos;
  for (int i = 0; i < markers.size(); i++) {
    int mx = static_cast<int>(std::round(markers[i].x));
    int my = static_cast<int>(std::round(markers[i].y));
    int mz = static_cast<int>(std::round(markers[i].z));
    if (mx < 0 || mx >= xdim || my < 0 || my >= ydim || mz < 0 || mz >= zdim)
      continue;
    MarkerInfo mi;
    mi.mx = mx;
    mi.my = my;
    mi.mz = mz;
    mi.r2 = markers[i].radius * markers[i].radius;
    markerInfos.push_back(mi);
  }
  if (markerInfos.empty()) return nullptr;

  // Bucket element for our bucket-based queue.
  struct BucketElement {
    V3DLONG idx;
    short markerLabel;  // (1-indexed, corresponding to markerInfos)
    float dist2;        // squared distance from marker seed
  };

  // Create 256 buckets corresponding to gradient intensity values (0-255).
  vector<vector<BucketElement> > buckets(256);

  // Precompute 26-neighborhood offsets.
  vector<int> neighborOffsets;
  for (int dz = -1; dz <= 1; dz++)
    for (int dy = -1; dy <= 1; dy++)
      for (int dx = -1; dx <= 1; dx++) {
        if (dx == 0 && dy == 0 && dz == 0) continue;
        int offset = dz * (xdim * ydim) + dy * xdim + dx;
        neighborOffsets.push_back(offset);
      }

  // Initialize: For each marker, label its seed voxel and add its neighbors.
  for (int i = 0; i < markerInfos.size(); i++) {
    short markerLabel = i + 1;
    MarkerInfo mi = markerInfos[i];
    V3DLONG seedIdx = mi.mz * (xdim * ydim) + mi.my * xdim + mi.mx;
    labels[seedIdx] = markerLabel;
    dist[seedIdx] = 0.0f;
    // For each neighbor of the seed:
    for (int offset : neighborOffsets) {
      V3DLONG nIdx = seedIdx + offset;
      int nz = nIdx / (xdim * ydim);
      int rem = nIdx % (xdim * ydim);
      int ny = rem / xdim;
      int nx = rem % xdim;
      if (nx < 0 || nx >= xdim || ny < 0 || ny >= ydim || nz < 0 || nz >= zdim)
        continue;
      int dx_ = nx - mi.mx;
      int dy_ = ny - mi.my;
      int dz_ = nz - mi.mz;
      float d2 = float(dx_ * dx_ + dy_ * dy_ + dz_ * dz_);
      if (d2 > mi.r2) continue;
      if (labels[nIdx] == UNVISITED) {
        labels[nIdx] = markerLabel;
        dist[nIdx] = d2;
        unsigned char intensity = gradMag->data[nIdx];
        BucketElement be;
        be.idx = nIdx;
        be.markerLabel = markerLabel;
        be.dist2 = d2;
        buckets[intensity].push_back(be);
      } else if (labels[nIdx] != markerLabel && labels[nIdx] != WSHED) {
        labels[nIdx] = WSHED;
      }
    }
  }

  // Process buckets in order of increasing gradient intensity.
  for (int level = 0; level < 256; level++) {
    size_t pos = 0;
    while (pos < buckets[level].size()) {
      BucketElement be = buckets[level][pos];
      pos++;
      V3DLONG idx = be.idx;
      int z = idx / (xdim * ydim);
      int rem = idx % (xdim * ydim);
      int y = rem / xdim;
      int x = rem % xdim;
      short currentLabel = be.markerLabel;
      MarkerInfo mi = markerInfos[currentLabel - 1];
      // Check all 26 neighbors.
      for (int offset : neighborOffsets) {
        V3DLONG nIdx = idx + offset;
        int nz = nIdx / (xdim * ydim);
        int rem2 = nIdx % (xdim * ydim);
        int ny = rem2 / xdim;
        int nx = rem2 % xdim;
        if (nx < 0 || nx >= xdim || ny < 0 || ny >= ydim || nz < 0 ||
            nz >= zdim)
          continue;
        int dx_ = nx - mi.mx;
        int dy_ = ny - mi.my;
        int dz_ = nz - mi.mz;
        float nd2 = float(dx_ * dx_ + dy_ * dy_ + dz_ * dz_);
        if (nd2 > mi.r2) continue;  // outside the allowed sphere.
        if (labels[nIdx] == UNVISITED) {
          labels[nIdx] = currentLabel;
          dist[nIdx] = nd2;
          unsigned char nIntensity = gradMag->data[nIdx];
          BucketElement newBe;
          newBe.idx = nIdx;
          newBe.markerLabel = currentLabel;
          newBe.dist2 = nd2;
          buckets[nIntensity].push_back(newBe);
        } else if (labels[nIdx] != currentLabel && labels[nIdx] != WSHED) {
          labels[nIdx] = WSHED;
        }
      }
    }
  }

  // Build the binary segmentation image: any voxel with label > 0 becomes 255.
  My4DImage *segImage = new My4DImage();
  segImage->xdim = xdim;
  segImage->ydim = ydim;
  segImage->zdim = zdim;
  segImage->cdim = 1;
  segImage->data = new unsigned char[totalSize];
  for (V3DLONG i = 0; i < totalSize; i++) {
    segImage->data[i] = (labels[i] > 0 ? 255 : 0);
  }
  return segImage;
}

////////////////////////////////////////////////////////////////////////
// Main reconstruction function (invoked from menu or command-line)
void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
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
    } else
      c = QInputDialog::getInt(parent, "Channel", "Enter channel NO:", 1, 1, sc,
                               1, &ok1);
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
      fprintf(stderr, "Error reading file [%s].\n",
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
    v3d_msg("No image window is currently open.", bmenu);
    return;
  }
  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage) {
    v3d_msg("Invalid image pointer.", bmenu);
    return;
  }
  LandmarkList landmarkList = callback.getLandmark(curwin);
  if (landmarkList.isEmpty()) {
    v3d_msg("No landmarks defined. Please define at least one landmark.",
            bmenu);
    return;
  }

  // ---- In the watershed approach we do not use intensity difference or
  // gradient threshold ---- Convert the LandmarkList to a QList of MyMarker
  // objects.
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
  // 1. Create a My4DImage holding the desired channel.
  // 2. Compute the smoothed image and gradient magnitude ONCE.
  // 3. Run the optimized marker–controlled watershed on the gradient image.
  // ============================

  V3DLONG totalSize =
      p4DImage->getXDim() * p4DImage->getYDim() * p4DImage->getZDim();

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
    return;
  }
  memcpy(channelImage->data, channelData, totalSize * sizeof(unsigned char));

  // (3) Compute a smoothed image and then its gradient magnitude.
  double sigma = 2.0;  // use a larger sigma for better smoothing
  My4DImage *smoothedImage = gaussianSmooth3D(channelImage, sigma);
  delete channelImage;
  if (!smoothedImage) {
    v3d_msg("Gaussian smoothing failed.", bmenu);
    return;
  }

  My4DImage *gradMag = computeGradientMagnitude(smoothedImage);
  if (!gradMag) {
    v3d_msg("Gradient magnitude computation failed.", bmenu);
    delete smoothedImage;
    return;
  }

  // (4) Run the optimized marker–controlled watershed segmentation.
  My4DImage *globalSegImage = markerControlledWatershed(gradMag, myMarkerList);

  // Clean up temporary images.
  delete smoothedImage;
  delete gradMag;

  // For diagnostic purposes, count the number of segmented (labeled) voxels.
  size_t segmented_voxels = 0;
  for (V3DLONG idx = 0; idx < totalSize; ++idx) {
    if (globalSegImage->data[idx] == 255) segmented_voxels++;
  }
  printf("Segmented voxels: %zu\n", segmented_voxels);

  // Open a new image window to display the segmentation.
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
    callback.setImageName(newwin, QString("Optimized Watershed Segmentation"));
    callback.updateImageWindow(newwin);
  }

  // (Optional) Save the binary segmented image as a TIFF file.
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

  delete globalSegImage;
  return;
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
    v3d_msg(
        tr("This plugin segments individual somas using a marker–controlled "
           "3D watershed segmentation algorithm. Developed by ImagiNeuron "
           "(2024-11-16) "
           "and updated by Athmane (2025-02-01), then optimized for speed and "
           "memory usage."),
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
