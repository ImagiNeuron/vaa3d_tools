/* soma_segmentation_plugin.cpp
 * This plugin segments individual somas using an improved 3D region-growing
 * algorithm. The algorithm uses a manual intensity difference threshold
 * (relative to the seed) plus additional constraints. In this merged version
 * the core region-growing routines are directly included.
 */

#include "soma_segmentation_plugin.h"

#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <vector>

#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

/***********************************
 * Implementation of My4DImage methods
 ***********************************/

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

/***********************************
 * Gaussian Smoothing and Gradient Computation
 ***********************************/

// Helper: 1D Gaussian kernel (centered).
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
// Assumes a single-channel image.
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

  // Convolve along x.
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

  // Convolve along y.
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

  // Convolve along z.
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

// Compute gradient magnitude using central differences.
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

/***********************************
 * Implementation of Merged Algorithm Class
 ***********************************/

SomaSegmentationAlgorithm::SomaSegmentationAlgorithm()
    : Image1D_page(nullptr),
      dim_X(0),
      dim_Y(0),
      dim_Z(0),
      size_page(0),
      offset_Y(0),
      offset_Z(0) {}

SomaSegmentationAlgorithm::~SomaSegmentationAlgorithm() {
  // Do not free Image1D_page here.
}

void SomaSegmentationAlgorithm::initialize(unsigned char *img, V3DLONG _dim_X,
                                           V3DLONG _dim_Y, V3DLONG _dim_Z) {
  Image1D_page = img;
  dim_X = _dim_X;
  dim_Y = _dim_Y;
  dim_Z = _dim_Z;
  size_page = dim_X * dim_Y * dim_Z;
  offset_Y = dim_X;
  offset_Z = dim_X * dim_Y;
  initializeConstants();
}

void SomaSegmentationAlgorithm::initializeConstants() {
  poss_neighborRelative.clear();
  point_neighborRelative.clear();
  for (V3DLONG z = -1; z <= 1; z++) {
    for (V3DLONG y = -1; y <= 1; y++) {
      for (V3DLONG x = -1; x <= 1; x++) {
        if (x == 0 && y == 0 && z == 0) continue;
        V3DLONG offset = z * offset_Z + y * offset_Y + x;
        poss_neighborRelative.push_back(offset);
        point_neighborRelative.push_back(double3D(x, y, z));
      }
    }
  }
}

vector<V3DLONG> SomaSegmentationAlgorithm::index2Coordinate(V3DLONG idx) {
  vector<V3DLONG> coord(3, 0);
  coord[2] = idx / offset_Z;
  V3DLONG rem = idx % offset_Z;
  coord[1] = rem / offset_Y;
  coord[0] = rem % offset_Y;
  return coord;
}

V3DLONG SomaSegmentationAlgorithm::coordinate2Index(V3DLONG x, V3DLONG y,
                                                    V3DLONG z) {
  return z * offset_Z + y * offset_Y + x;
}

bool SomaSegmentationAlgorithm::checkValidity(V3DLONG idx) {
  return (idx >= 0 && idx < size_page);
}

vector<V3DLONG> SomaSegmentationAlgorithm::regionGrowOnPos(
    V3DLONG pos_seed, V3DLONG threshold_voxelValue,
    double threshold_valueChangeRatio, V3DLONG uThreshold_regionSize,
    unsigned char *mask_input) {
  vector<V3DLONG> poss_result;
  vector<V3DLONG> poss_growing;
  poss_growing.push_back(pos_seed);
  poss_result.push_back(pos_seed);
  V3DLONG count_voxel = 1;
  mask_input[pos_seed] = 0;  // Mark as visited.
  V3DLONG min_voxelValue = Image1D_page[pos_seed];

  while (!poss_growing.empty()) {
    V3DLONG pos_current = poss_growing.back();
    poss_growing.pop_back();
    vector<V3DLONG> xyz_current = index2Coordinate(pos_current);
    for (size_t j = 0; j < poss_neighborRelative.size(); j++) {
      vector<V3DLONG> neighbor_coord = {
          xyz_current[0] + static_cast<V3DLONG>(point_neighborRelative[j].x),
          xyz_current[1] + static_cast<V3DLONG>(point_neighborRelative[j].y),
          xyz_current[2] + static_cast<V3DLONG>(point_neighborRelative[j].z)};
      if (neighbor_coord[0] < 0 || neighbor_coord[0] >= dim_X ||
          neighbor_coord[1] < 0 || neighbor_coord[1] >= dim_Y ||
          neighbor_coord[2] < 0 || neighbor_coord[2] >= dim_Z)
        continue;

      V3DLONG pos_neighbor = pos_current + poss_neighborRelative[j];
      if (!checkValidity(pos_neighbor)) continue;
      if (mask_input[pos_neighbor] > 0) {
        V3DLONG value_neighbor = Image1D_page[pos_neighbor];
        if (value_neighbor > threshold_voxelValue &&
            ((min_voxelValue - value_neighbor) <
             (min_voxelValue * threshold_valueChangeRatio))) {
          mask_input[pos_neighbor] = 0;
          poss_growing.push_back(pos_neighbor);
          poss_result.push_back(pos_neighbor);
          count_voxel++;
          if (value_neighbor < min_voxelValue) min_voxelValue = value_neighbor;
          if (count_voxel > (uThreshold_regionSize + 2)) return poss_result;
        }
      }
    }
  }
  return poss_result;
}

V3DLONG SomaSegmentationAlgorithm::getCenterByMass(
    const vector<V3DLONG> &voxels) {
  if (voxels.empty()) return -1;
  double sum_X = 0, sum_Y = 0, sum_Z = 0, sum_mass = 0;
  for (size_t i = 0; i < voxels.size(); i++) {
    vector<V3DLONG> coord = index2Coordinate(voxels[i]);
    double value_voxel = Image1D_page[voxels[i]];
    sum_X += coord[0] * value_voxel;
    sum_Y += coord[1] * value_voxel;
    sum_Z += coord[2] * value_voxel;
    sum_mass += value_voxel;
  }
  double cx = sum_X / sum_mass;
  double cy = sum_Y / sum_mass;
  double cz = sum_Z / sum_mass;
  return coordinate2Index(static_cast<V3DLONG>(round(cx)),
                          static_cast<V3DLONG>(round(cy)),
                          static_cast<V3DLONG>(round(cz)));
}

vector<V3DLONG> SomaSegmentationAlgorithm::getBoundBox(
    const vector<V3DLONG> &indices) {
  vector<V3DLONG> bbox(6, 0);  // {minX, maxX, minY, maxY, minZ, maxZ}
  if (indices.empty()) return bbox;
  V3DLONG minX = 1e9, minY = 1e9, minZ = 1e9;
  V3DLONG maxX = -1e9, maxY = -1e9, maxZ = -1e9;
  for (size_t i = 0; i < indices.size(); i++) {
    vector<V3DLONG> coord = index2Coordinate(indices[i]);
    minX = min(minX, coord[0]);
    maxX = max(maxX, coord[0]);
    minY = min(minY, coord[1]);
    maxY = max(maxY, coord[1]);
    minZ = min(minZ, coord[2]);
    maxZ = max(maxZ, coord[2]);
  }
  bbox[0] = minX;
  bbox[1] = maxX;
  bbox[2] = minY;
  bbox[3] = maxY;
  bbox[4] = minZ;
  bbox[5] = maxZ;
  return bbox;
}

void SomaSegmentationAlgorithm::poss2Image1D(const vector<V3DLONG> &poss,
                                             unsigned char *Image1D_output,
                                             V3DLONG value) {
  for (size_t i = 0; i < poss.size(); i++) {
    if (checkValidity(poss[i]))
      Image1D_output[poss[i]] = static_cast<unsigned char>(value);
  }
}

/***********************************
 * New Segmentation Routine: segmentSoma
 ***********************************/

My4DImage *segmentSoma(const My4DImage *input, const MyMarker &marker) {
  if (!input || !input->data) return nullptr;

  V3DLONG xdim = input->xdim;
  V3DLONG ydim = input->ydim;
  V3DLONG zdim = input->zdim;
  V3DLONG totalSize = xdim * ydim * zdim;

  // Create output binary image.
  My4DImage *segImage = new My4DImage();
  segImage->xdim = xdim;
  segImage->ydim = ydim;
  segImage->zdim = zdim;
  segImage->cdim = 1;
  segImage->data = new unsigned char[totalSize];
  memset(segImage->data, 0, totalSize * sizeof(unsigned char));

  // Determine seed from marker.
  int x0 = static_cast<int>(std::round(marker.x));
  int y0 = static_cast<int>(std::round(marker.y));
  int z0 = static_cast<int>(std::round(marker.z));
  if (x0 < 0 || x0 >= xdim || y0 < 0 || y0 >= ydim || z0 < 0 || z0 >= zdim)
    return nullptr;
  V3DLONG seedIdx = z0 * ydim * xdim + y0 * xdim + x0;

  // Prepare a mask (255 = unvisited).
  unsigned char *mask = new unsigned char[totalSize];
  memset(mask, 255, totalSize);

  // Use an intensity difference threshold relative to the seed.
  int intensityDiff = 20;  // This value may be tuned or prompted.
  int seedIntensity = static_cast<int>(input->data[seedIdx]);
  int lowerThreshold = std::max(0, seedIntensity - intensityDiff);
  int upperThreshold = std::min(255, seedIntensity + intensityDiff);
  // (In this simple merge we use lowerThreshold.)

  // Initialize the algorithm.
  SomaSegmentationAlgorithm algo;
  algo.initialize(input->data, xdim, ydim, zdim);

  // Perform region growing.
  vector<V3DLONG> region =
      algo.regionGrowOnPos(seedIdx, lowerThreshold, 0.1, totalSize, mask);

  // Mark segmented voxels.
  for (V3DLONG idx : region) segImage->data[idx] = 255;

  delete[] mask;
  return segImage;
}

/***********************************
 * Plugin Reconstruction Function
 ***********************************/

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

  // Prompt for intensity difference threshold.
  bool okIntensity = false;
  int intensityDiff =
      QInputDialog::getInt(parent, "Intensity Difference Threshold",
                           "Enter allowed difference from seed intensity:", 20,
                           0, 255, 1, &okIntensity);
  if (!okIntensity) return nullptr;

  // Convert LandmarkList into QList<MyMarker>.
  QList<MyMarker> myMarkerList;
  for (const auto &lm : landmarkList) {
    MyMarker marker;
    marker.x = lm.x;
    marker.y = lm.y;
    marker.z = lm.z;
    marker.radius = lm.radius;
    myMarkerList.append(marker);
  }

  V3DLONG totalSize =
      p4DImage->getXDim() * p4DImage->getYDim() * p4DImage->getZDim();
  // Create global segmentation image.
  My4DImage *globalSegImage = new My4DImage();
  globalSegImage->xdim = p4DImage->getXDim();
  globalSegImage->ydim = p4DImage->getYDim();
  globalSegImage->zdim = p4DImage->getZDim();
  globalSegImage->cdim = 1;
  globalSegImage->data = new unsigned char[totalSize];
  memset(globalSegImage->data, 0, totalSize * sizeof(unsigned char));

  // (2) Extract the desired channel.
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

  // (3) Preprocess: Gaussian smoothing.
  double sigma = 2.0;
  My4DImage *smoothedImage = gaussianSmooth3D(channelImage, sigma);
  delete channelImage;
  if (!smoothedImage) {
    v3d_msg("Gaussian smoothing failed.", bmenu);
    delete globalSegImage;
    return nullptr;
  }

  // (Optional) Compute gradient magnitude.
  My4DImage *gradMag = computeGradientMagnitude(smoothedImage);
  if (!gradMag) {
    v3d_msg("Gradient magnitude computation failed.", bmenu);
    delete smoothedImage;
    delete globalSegImage;
    return nullptr;
  }

  // (4) For each landmark, use the new segmentSoma routine and merge the
  // result.
  for (const auto &marker : myMarkerList) {
    My4DImage *localSeg = segmentSoma(smoothedImage, marker);
    if (localSeg) {
      for (V3DLONG idx = 0; idx < totalSize; ++idx) {
        if (localSeg->data[idx] == 255) globalSegImage->data[idx] = 255;
      }
      delete localSeg;
    }
  }

  // Clean up temporary images.
  delete smoothedImage;
  delete gradMag;

  // Diagnostic: print number of segmented voxels.
  size_t flooded_voxels = 0;
  for (V3DLONG idx = 0; idx < totalSize; ++idx) {
    if (globalSegImage->data[idx] == 255) flooded_voxels++;
  }
  printf("Flooded voxels: %zu\n", flooded_voxels);

  // Display the segmentation.
  if (globalSegImage) {
    Image4DSimple p4DImageSeg;
    p4DImageSeg.setData(globalSegImage->data, globalSegImage->xdim,
                        globalSegImage->ydim, globalSegImage->zdim, 1,
                        V3D_UINT8);
    // Detach the pointer.
    globalSegImage->data = nullptr;
    v3dhandle newwin = callback.newImageWindow();
    callback.setImage(newwin, &p4DImageSeg);
    callback.setImageName(newwin, QString("Segmented Image"));
    callback.updateImageWindow(newwin);
  }

  return globalSegImage;
}

/***********************************
 * Plugin Interface Implementations
 ***********************************/

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
               "with landmarks as seeds. Developed by ImagiNeuron."),
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
    } else {
      PARA.inimg_file = infiles[0];
    }
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
