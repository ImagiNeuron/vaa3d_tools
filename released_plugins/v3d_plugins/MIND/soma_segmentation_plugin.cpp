/* soma_segmentation_plugin.cpp
 * A plugin for seeded soma segmentation using 3D watershed.
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"

#include <QInputDialog>
#include <QMessageBox>
#include <algorithm>  // For std::sort and std::nth_element
#include <cmath>
#include <cstring>  // For memcpy
#include <queue>
#include <vector>

#include "../../v3d/compute_win_pca_wp.h"
#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

/**************************************
 * My4DImage implementations
 **************************************/
My4DImage::My4DImage() : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0) {}

My4DImage::My4DImage(const My4DImage &other)
    : xdim(other.xdim), ydim(other.ydim), zdim(other.zdim), cdim(other.cdim) {
  V3DLONG totalSize = xdim * ydim * zdim * cdim;
  if (other.data && totalSize > 0) {
    data = new unsigned char[totalSize];
    memcpy(data, other.data, totalSize);
  } else {
    data = nullptr;
  }
}

My4DImage &My4DImage::operator=(const My4DImage &other) {
  if (this != &other) {
    delete[] data;
    xdim = other.xdim;
    ydim = other.ydim;
    zdim = other.zdim;
    cdim = other.cdim;

    V3DLONG totalSize = xdim * ydim * zdim * cdim;
    if (other.data && totalSize > 0) {
      data = new unsigned char[totalSize];
      memcpy(data, other.data, totalSize);
    } else {
      data = nullptr;
    }
  }
  return *this;
}

My4DImage::~My4DImage() { delete[] data; }

/**************************************
 * Internal struct to store params
 **************************************/
struct input_PARA {
  QString inimg_file;
  V3DLONG channel;
};

/**************************************
 * Forward declaration of main
 **************************************/
void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu);

/**************************************
 * Plugin Interface Methods
 **************************************/
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
        tr("Soma segmentation plugin using 3D watershed.\n"
           "Developed by ImagiNeuron, 2024-11-16"));
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
    printf("**** Usage of soma_segmentation (3D watershed) ****\n");
    printf(
        "vaa3d -x soma_segmentation -f segment_somas -i <inimg_file> -p "
        "<channel>\n");
  } else {
    return false;
  }
  return true;
}

/**************************************
 * The main function
 **************************************/
void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu) {
  /*************************************
   * 1) Load Image
   *************************************/
  unsigned char *data1d = 0;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];

  if (bmenu) {
    // From the current Vaa3D window
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin) {
      QMessageBox::information(0, "", "No image open in the main window.");
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

    bool ok1;
    if (sc == 1) {
      c = 1;
      ok1 = true;
    } else {
      c = QInputDialog::getInt(parent, "Channel", "Enter channel NO:", 1, 1, sc,
                               1, &ok1);
    }
    if (!ok1) return;

    in_sz[0] = N;
    in_sz[1] = M;
    in_sz[2] = P;
    in_sz[3] = sc;

    PARA.inimg_file = p4DImage->getFileName();
  } else {
    // Command-line input
    int datatype = 0;
    if (!simple_loadimage_wrapper(callback,
                                  PARA.inimg_file.toStdString().c_str(), data1d,
                                  in_sz, datatype)) {
      fprintf(stderr, "Error loading file [%s].\n",
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

  if (!data1d) {
    v3d_msg("No valid image data!", bmenu);
    return;
  }

  if (sc > 1) {
    v3d_msg("For grayscale, ignoring additional channels");
    sc = 1;
    c = 1;
  }

  /*************************************
   * 2) Preprocessing (Median + Gaussian)
   *************************************/
  {
    // Extract the desired channel
    V3DLONG totalSize = N * M * P;
    unsigned char *channelData = new unsigned char[totalSize];
    {
      // Usually channels are sequential. Copy channel c-1.
      unsigned char *cptr = data1d + (c - 1) * (N * M * P);
      memcpy(channelData, cptr, totalSize);
    }

    // 2.1) Median Filter
    unsigned char *medianData = new unsigned char[totalSize];
    applyMedianFilter(channelData, medianData, N, M, P, /*windowSize=*/3);
    delete[] channelData;

    // 2.2) Gaussian Filter
    unsigned char *gaussData = new unsigned char[totalSize];
    applyGaussianFilter(medianData, gaussData, N, M, P, /*sigma=*/1.0f);
    delete[] medianData;

    // Now data1d points to the new, single-channel preprocessed data
    data1d = gaussData;
  }

  /*************************************
   * 3) For each landmark, do watershed
   *************************************/
  v3dhandle curwin = callback.currentImageWindow();
  LandmarkList landmarkList;
  if (bmenu && curwin) {
    landmarkList = callback.getLandmark(curwin);
  } else {
    // In command-line or no open window, user must supply landmarks some other
    // way For now, if none, we just exit.
  }

  if (landmarkList.isEmpty()) {
    v3d_msg("No landmarks found. Please specify at least one landmark.", bmenu);
    delete[] data1d;
    return;
  }

  // Prepare a global label volume (16-bit) to store watershed results
  V3DLONG totalSize = N * M * P;
  unsigned short *labelVolume = new unsigned short[totalSize];
  memset(labelVolume, 0, totalSize * sizeof(unsigned short));

  // For printing
  int somaIndex = 1;

  for (int i = 0; i < landmarkList.size(); i++) {
    LocationSimple lm = landmarkList[i];
    float x = lm.x;
    float y = lm.y;
    float z = lm.z;
    float r = lm.radius > 0 ? lm.radius : 5.0f;

    // bounding box
    int x1 = std::max(0, (int)floor(x - r));
    int x2 = std::min((int)N - 1, (int)ceil(x + r));
    int y1 = std::max(0, (int)floor(y - r));
    int y2 = std::min((int)M - 1, (int)ceil(y + r));
    int z1 = std::max(0, (int)floor(z - r));
    int z2 = std::min((int)P - 1, (int)ceil(z + r));

    int sx = x2 - x1 + 1;
    int sy = y2 - y1 + 1;
    int sz = z2 - z1 + 1;

    // Extract subvolume
    unsigned char *subvol = new unsigned char[sx * sy * sz];
    extractSubvolume(data1d, subvol, N, M, P, x1, x2, y1, y2, z1, z2);

    // Compute Otsu threshold
    int otsuValue = computeOtsuThreshold(subvol, sx * sy * sz);

    for (int idx = 0; idx < sx * sy * sz; idx++) {
      if (subvol[idx] < otsuValue) {
        subvol[idx] = 0;
      }
    }

    // 2) Morphological opening to remove small noise or bridging artifacts
    morphologicalOpen3D(subvol, sx, sy, sz);

    // Run 3D watershed
    unsigned short *wsLabels = new unsigned short[sx * sy * sz];
    memset(wsLabels, 0, sx * sy * sz * sizeof(unsigned short));
    applyWatershedVS(subvol, wsLabels, sx, sy, sz);

    // Copy watershed result to global label volume
    V3DLONG voxelCount = 0;
    int somaMinX = x2, somaMaxX = x1;
    int somaMinY = y2, somaMaxY = y1;
    int somaMinZ = z2, somaMaxZ = z1;

    for (int zz = 0; zz < sz; zz++) {
      for (int yy = 0; yy < sy; yy++) {
        for (int xx = 0; xx < sx; xx++) {
          V3DLONG idxSub = (V3DLONG)zz * (sx * sy) + yy * sx + xx;
          unsigned short labelVal = wsLabels[idxSub];
          if (labelVal > 0) {
            V3DLONG gx = x1 + xx;
            V3DLONG gy = y1 + yy;
            V3DLONG gz = z1 + zz;
            V3DLONG idxGlobal = gz * (N * M) + gy * N + gx;
            labelVolume[idxGlobal] = somaIndex;  // simple approach

            // track bounding box
            if (gx < somaMinX) somaMinX = gx;
            if (gx > somaMaxX) somaMaxX = gx;
            if (gy < somaMinY) somaMinY = gy;
            if (gy > somaMaxY) somaMaxY = gy;
            if (gz < somaMinZ) somaMinZ = gz;
            if (gz > somaMaxZ) somaMaxZ = gz;

            voxelCount++;
          }
        }
      }
    }

    // Print bounding box info
    printf("\nSoma #%d at landmark (%.1f, %.1f, %.1f, r=%.1f):\n", somaIndex, x,
           y, z, r);
    printf("   X range: [%d, %d]\n", somaMinX, somaMaxX);
    printf("   Y range: [%d, %d]\n", somaMinY, somaMaxY);
    printf("   Z range: [%d, %d]\n", somaMinZ, somaMaxZ);
    printf("   Volume (# of voxels): %lld\n", (long long)voxelCount);

    delete[] subvol;
    delete[] wsLabels;
    somaIndex++;
  }

  // Create a new image to show or save, I will use this later for visualization
  {
    My4DImage *labeledResult = new My4DImage();
    labeledResult->xdim = N;
    labeledResult->ydim = M;
    labeledResult->zdim = P;
    labeledResult->cdim = 1;
    labeledResult->data = new unsigned char[totalSize];

    // Convert 16-bit labels to 8-bit (clamp).
    for (V3DLONG i = 0; i < totalSize; i++) {
      unsigned short val = labelVolume[i];
      labeledResult->data[i] = (val > 255) ? 255 : (unsigned char)(val);
    }

    printf("\nFinal labeled image (clamped to 8-bit) is ready.\n");

    // Apply PCA analysis to each soma
    for (int i = 0; i < landmarkList.size(); i++) {
      analyzeSomaPCA(labeledResult->data, N, M, P, landmarkList[i], i + 1);
    }

    delete labeledResult;
  }

  delete[] labelVolume;
  delete[] data1d;

  v3d_msg("3D watershed segmentation complete.", bmenu);
}

/*************************************
 *        MEDIAN FILTER
 *************************************/
void applyMedianFilter(const unsigned char *inputData,
                       unsigned char *outputData, V3DLONG N, V3DLONG M,
                       V3DLONG P, int windowSize) {
  int halfWindow = windowSize / 2;
  vector<unsigned char> window;
  window.reserve(windowSize * windowSize * windowSize);

  for (V3DLONG z = 0; z < P; z++) {
    for (V3DLONG y = 0; y < M; y++) {
      for (V3DLONG x = 0; x < N; x++) {
        window.clear();
        for (int dz = -halfWindow; dz <= halfWindow; dz++) {
          for (int dy = -halfWindow; dy <= halfWindow; dy++) {
            for (int dx = -halfWindow; dx <= halfWindow; dx++) {
              V3DLONG nx = x + dx;
              V3DLONG ny = y + dy;
              V3DLONG nz = z + dz;
              if (nx >= 0 && nx < N && ny >= 0 && ny < M && nz >= 0 && nz < P) {
                V3DLONG idx = nz * (M * N) + ny * N + nx;
                window.push_back(inputData[idx]);
              }
            }
          }
        }
        // find median
        std::nth_element(window.begin(), window.begin() + window.size() / 2,
                         window.end());
        unsigned char medVal = window[window.size() / 2];
        V3DLONG outIdx = z * (M * N) + y * N + x;
        outputData[outIdx] = medVal;
      }
    }
  }
}

/*************************************
 *        GAUSSIAN FILTER
 *************************************/
void applyGaussianFilter(const unsigned char *inputData,
                         unsigned char *outputData, V3DLONG N, V3DLONG M,
                         V3DLONG P, float sigma) {
  // For simplicity, use a fixed kernel size 5x5x5
  int windowSize = 5;
  int halfWindow = windowSize / 2;

  // build kernel
  vector<float> kernel(windowSize * windowSize * windowSize);
  float sumKernel = 0.0f;
  int idx = 0;
  for (int dz = -halfWindow; dz <= halfWindow; dz++) {
    for (int dy = -halfWindow; dy <= halfWindow; dy++) {
      for (int dx = -halfWindow; dx <= halfWindow; dx++) {
        float val = expf(-(dx * dx + dy * dy + dz * dz) / (2 * sigma * sigma));
        kernel[idx++] = val;
        sumKernel += val;
      }
    }
  }
  // normalize
  for (size_t i = 0; i < kernel.size(); i++) kernel[i] /= sumKernel;

  // convolve
  for (V3DLONG z = 0; z < P; z++) {
    for (V3DLONG y = 0; y < M; y++) {
      for (V3DLONG x = 0; x < N; x++) {
        float accum = 0.0f;
        idx = 0;
        for (int dz = -halfWindow; dz <= halfWindow; dz++) {
          for (int dy = -halfWindow; dy <= halfWindow; dy++) {
            for (int dx = -halfWindow; dx <= halfWindow; dx++) {
              V3DLONG nx = x + dx;
              V3DLONG ny = y + dy;
              V3DLONG nz = z + dz;
              if (nx >= 0 && nx < N && ny >= 0 && ny < M && nz >= 0 && nz < P) {
                V3DLONG inIdx = nz * (M * N) + ny * N + nx;
                accum += inputData[inIdx] * kernel[idx];
              }
              idx++;
            }
          }
        }
        V3DLONG outIdx = z * (M * N) + y * N + x;
        outputData[outIdx] = (unsigned char)(accum + 0.5f);
      }
    }
  }
}

/*************************************
 *     EXTRACT SUBVOLUME
 *************************************/
void extractSubvolume(const unsigned char *inData, unsigned char *outData,
                      V3DLONG N, V3DLONG M, V3DLONG P, int x1, int x2, int y1,
                      int y2, int z1, int z2) {
  int sx = x2 - x1 + 1;
  int sy = y2 - y1 + 1;
  int sz = z2 - z1 + 1;

  for (int zz = z1; zz <= z2; zz++) {
    for (int yy = y1; yy <= y2; yy++) {
      for (int xx = x1; xx <= x2; xx++) {
        V3DLONG gx = zz * (N * M) + yy * N + xx;
        int subz = zz - z1;
        int suby = yy - y1;
        int subx = xx - x1;
        V3DLONG lx = subz * (sx * sy) + suby * sx + subx;
        outData[lx] = inData[gx];
      }
    }
  }
}

/*************************************************
 * 3D Watershed (Vincent & Soille) Implementation
 *   - "watershed_vs" & "remove_watershed_lines"
 *************************************************/
#define INIT -1
#define MASK -2
#define WSHED 0
#define FICTITIOUS -1

static void remove_watershed_lines(float *&label_data, const V3DLONG *sz,
                                   const V3DLONG ndims) {
  V3DLONG num_elements = 1;
  for (int i = 0; i < ndims; i++) num_elements *= sz[i];

  // Prepare arrays
  float *tmp_data = new float[num_elements];
  V3DLONG *pix_idx = new V3DLONG[num_elements];
  // Binarize label_data: 1 for region, 0 for lines
  for (V3DLONG i = 0; i < num_elements; i++) {
    if (label_data[i] == 1)
      tmp_data[i] = 0;  // watershed line
    else if (label_data[i] > 1)
      tmp_data[i] = 1;  // region
    else
      tmp_data[i] = 0;  // region=1, ws=1 => 0 means line
  }

  // run distance transform to find nearest region pixel
  dt3d_binary(tmp_data, pix_idx, sz, 0);

  // Re-assign watershed lines to one of the regions
  float *copy_data = new float[num_elements];
  memcpy(copy_data, label_data, num_elements * sizeof(float));

  for (V3DLONG i = 0; i < num_elements; i++) {
    if (copy_data[i] == 0)  // was watershed line
    {
      label_data[i] = copy_data[pix_idx[i]];  // adopt label
    }
  }

  // Decrement all labels by 1 (so background becomes 0)
  for (V3DLONG i = 0; i < num_elements; i++) label_data[i] = label_data[i] - 1;

  delete[] tmp_data;
  delete[] pix_idx;
  delete[] copy_data;
}

// Internal BFS queue
template <class T>
class SimpleQueue {
 public:
  queue<T> q;
  void put(const T &val) { q.push(val); }
  T get() {
    T frontval = q.front();
    q.pop();
    return frontval;
  }
  bool empty() const { return q.empty(); }
};

// We define a minimal "6-connected" neighbor approach
static const int NB_6[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                               {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};

// A helper to do sorting with a separate index array
static void sort2(V3DLONG n, float *data, float *idx) {
  struct Info {
    float val;
    float ind;
  };
  vector<Info> all(n);
  for (V3DLONG i = 0; i < n; i++) {
    all[i].val = data[i + 1];
    all[i].ind = idx[i + 1];
  }
  std::sort(all.begin(), all.end(),
            [](const Info &a, const Info &b) { return (a.val < b.val); });
  for (V3DLONG i = 0; i < n; i++) {
    data[i + 1] = all[i].val;
    idx[i + 1] = all[i].ind;
  }
}

// The main watershed
template <class T>
void compute_watershed(T *indata, float *sortidx, V3DLONG num_elements,
                       const V3DLONG *sz, float *&label_data) {
  // label_data must be uninitialized
  for (V3DLONG i = 0; i < num_elements; i++) label_data[i] = INIT;

  V3DLONG current_label = 0;
  V3DLONG *dist = new V3DLONG[num_elements];
  memset(dist, 0, num_elements * sizeof(V3DLONG));

  SimpleQueue<V3DLONG> pixelQueue;
  V3DLONG count = 0;

  while (count < num_elements) {
    // group pixels with same intensity
    V3DLONG k1 = count;
    float current_level = indata[(V3DLONG)sortidx[k1]];
    V3DLONG k2 = k1;
    while ((k2 + 1 < num_elements) &&
           (indata[(V3DLONG)sortidx[k2 + 1]] == current_level)) {
      k2++;
    }

    // mask all pixels with value == current_level
    for (V3DLONG k = k1; k <= k2; k++) {
      V3DLONG p = (V3DLONG)sortidx[k];
      label_data[p] = MASK;

      // check neighbors that are labeled or WSHED
      // if so, dist[p] = 1, queue put(p)
      // minimal 6-connect
      V3DLONG z = p / (sz[0] * sz[1]);
      V3DLONG r = p % (sz[0] * sz[1]);
      V3DLONG y = r / sz[0];
      V3DLONG x = r % sz[0];

      for (int ni = 0; ni < 6; ni++) {
        V3DLONG nx = x + NB_6[ni][0];
        V3DLONG ny = y + NB_6[ni][1];
        V3DLONG nz = z + NB_6[ni][2];
        if (nx < 0 || nx >= sz[0] || ny < 0 || ny >= sz[1] || nz < 0 ||
            nz >= sz[2])
          continue;
        V3DLONG q = nz * (sz[0] * sz[1]) + ny * sz[0] + nx;

        if ((label_data[q] > 0) || (label_data[q] == WSHED)) {
          dist[p] = 1;
          pixelQueue.put(p);
          break;
        }
      }
      count++;
    }

    V3DLONG current_distance = 1;
    pixelQueue.put(FICTITIOUS);

    while (true) {
      V3DLONG p = pixelQueue.get();
      if (p == FICTITIOUS) {
        if (pixelQueue.empty()) {
          break;
        } else {
          pixelQueue.put(FICTITIOUS);
          current_distance++;
          p = pixelQueue.get();
        }
      }

      // find neighbors with the closest distance
      V3DLONG z = p / (sz[0] * sz[1]);
      V3DLONG r = p % (sz[0] * sz[1]);
      V3DLONG y = r / sz[0];
      V3DLONG x = r % sz[0];

      V3DLONG closest_dist = 999999;
      V3DLONG closest_label_val = 0;
      bool unique_label = true;

      // check neighbors
      for (int ni = 0; ni < 6; ni++) {
        V3DLONG nx = x + NB_6[ni][0];
        V3DLONG ny = y + NB_6[ni][1];
        V3DLONG nz = z + NB_6[ni][2];
        if (nx < 0 || nx >= sz[0] || ny < 0 || ny >= sz[1] || nz < 0 ||
            nz >= sz[2])
          continue;
        V3DLONG q = nz * (sz[0] * sz[1]) + ny * sz[0] + nx;

        if ((label_data[q] > 0) || (label_data[q] == WSHED)) {
          if (dist[q] < closest_dist) {
            closest_dist = dist[q];
            if (label_data[q] > 0) closest_label_val = (V3DLONG)label_data[q];
          } else if (dist[q] == closest_dist) {
            if (label_data[q] > 0 && (closest_label_val > 0) &&
                (label_data[q] != closest_label_val)) {
              unique_label = false;
            }
            if (label_data[q] > 0) closest_label_val = (V3DLONG)label_data[q];
          }
        } else if ((label_data[q] == MASK) && (dist[q] == 0)) {
          dist[q] = current_distance + 1;
          pixelQueue.put(q);
        }
      }

      // label p
      if ((closest_dist < current_distance) && (closest_label_val > 0)) {
        if (unique_label && (label_data[p] == MASK || label_data[p] == WSHED)) {
          label_data[p] = (float)closest_label_val;
        } else if (!unique_label || (label_data[p] != closest_label_val)) {
          label_data[p] = WSHED;
        }
      } else if (label_data[p] == MASK) {
        label_data[p] = WSHED;
      }
    }

    // detect new minima at current_level
    for (V3DLONG k = k1; k <= k2; k++) {
      V3DLONG p = (V3DLONG)sortidx[k];
      dist[p] = 0;
      if (label_data[p] == MASK) {
        current_label++;
        SimpleQueue<V3DLONG> Q2;
        Q2.put(p);
        label_data[p] = (float)current_label;
        while (!Q2.empty()) {
          V3DLONG q = Q2.get();
          // minimal 6-connect
          V3DLONG zq = q / (sz[0] * sz[1]);
          V3DLONG rq = q % (sz[0] * sz[1]);
          V3DLONG yq = rq / sz[0];
          V3DLONG xq = rq % sz[0];

          for (int ni = 0; ni < 6; ni++) {
            V3DLONG nx = xq + NB_6[ni][0];
            V3DLONG ny = yq + NB_6[ni][1];
            V3DLONG nz = zq + NB_6[ni][2];
            if (nx < 0 || nx >= sz[0] || ny < 0 || ny >= sz[1] || nz < 0 ||
                nz >= sz[2])
              continue;
            V3DLONG rr = nz * (sz[0] * sz[1]) + ny * sz[0] + nx;

            if (label_data[rr] == MASK) {
              Q2.put(rr);
              label_data[rr] = (float)current_label;
            }
          }
        }
      }
    }
  }

  delete[] dist;
}

// The main driver for watershed + remove lines
template <class T>
void watershed_vs(T *indata, float *&label_data, const V3DLONG *sz,
                  V3DLONG ndims, V3DLONG conn_code) {
  // Prepare
  V3DLONG num_elements = 1;
  for (int i = 0; i < ndims; i++) num_elements *= sz[i];

  label_data = new float[num_elements];
  // We must sort intensities
  float *sortidx = new float[num_elements + 1];
  float *valdata = new float[num_elements + 1];

  // fill [1..num_elements] for sorting
  valdata[0] = -999;
  sortidx[0] = -999;
  for (V3DLONG i = 0; i < num_elements; i++) {
    valdata[i + 1] = (float)indata[i];
    sortidx[i + 1] = (float)i;
  }

  // sort
  sort2(num_elements, valdata, sortidx);

  // run watershed
  compute_watershed(indata, sortidx, num_elements, sz, label_data);

  // remove lines
  remove_watershed_lines(label_data, sz, ndims);

  delete[] sortidx;
  delete[] valdata;
}

/**************************************
 *  applyWatershedVS()
 *   - wrap the snippet for subvolume
 **************************************/
void applyWatershedVS(const unsigned char *subvol, unsigned short *labelOut,
                      int sx, int sy, int sz) {
  // The snippet operates on float output
  float *label_f = nullptr;
  V3DLONG dims[3];
  dims[0] = sx;
  dims[1] = sy;
  dims[2] = sz;

  // Run watershed
  watershed_vs(subvol, label_f, dims, /*ndims=*/3, /*conn_code=*/6);

  // We then convert label_f to 16-bit for returning
  V3DLONG totalSz = (V3DLONG)sx * sy * sz;
  for (V3DLONG i = 0; i < totalSz; i++) {
    float val = label_f[i];
    // we clamp to [0, 65535]
    if (val < 0) val = 0;
    if (val > 65535) val = 65535;
    labelOut[i] = (unsigned short)(val + 0.5f);
  }
  delete[] label_f;
}

/*********************************************************
 * Minimal BFS-based 3D distance transform for watershed line removal
 *   - pix_index[i] gets index of the nearest "object" voxel
 *********************************************************/
void dt3d_binary(const float *inData, V3DLONG *pix_index, const V3DLONG *sz,
                 float threshVal) {
  // multi-source BFS:
  //   inData[i] > threshVal => source (object)
  //   inData[i] <= threshVal => background
  // For each background, find nearest object via BFS
  V3DLONG Nx = sz[0], Ny = sz[1], Nz = sz[2];
  V3DLONG numel = Nx * Ny * Nz;
  // arrays
  vector<int> dist(numel, INT_MAX);
  // pix_index[i] will store the index of nearest object voxel
  // We'll push all "object" voxels in a queue with dist=0
  queue<V3DLONG> Q;

  // Initialize
  for (V3DLONG i = 0; i < numel; i++) {
    if (inData[i] > threshVal) {
      dist[i] = 0;
      pix_index[i] = i;  // nearest to itself
      Q.push(i);
    }
  }

  // BFS
  while (!Q.empty()) {
    V3DLONG curr = Q.front();
    Q.pop();
    V3DLONG z = curr / (Nx * Ny);
    V3DLONG r = curr % (Nx * Ny);
    V3DLONG y = r / Nx;
    V3DLONG x = r % Nx;

    // 6 neighbors
    for (int i = 0; i < 6; i++) {
      int nx = x + NB_6[i][0];
      int ny = y + NB_6[i][1];
      int nz = z + NB_6[i][2];
      if (nx < 0 || nx >= Nx || ny < 0 || ny >= Ny || nz < 0 || nz >= Nz)
        continue;
      V3DLONG idx = nz * (Nx * Ny) + ny * Nx + nx;
      if (dist[idx] > dist[curr] + 1) {
        dist[idx] = dist[curr] + 1;
        pix_index[idx] = pix_index[curr];  // inherits the same nearest object
        Q.push(idx);
      }
    }
  }
}

// Simple 3D morphological opening (erosion + dilation) with a 3x3x3
// neighborhood
static void morphologicalOpen3D(unsigned char *vol, int sx, int sy, int sz) {
  // Erode
  unsigned char *eroded = new unsigned char[sx * sy * sz];
  memset(eroded, 0, sx * sy * sz);
  for (int z = 1; z < sz - 1; z++) {
    for (int y = 1; y < sy - 1; y++) {
      for (int x = 1; x < sx - 1; x++) {
        int idx = z * sx * sy + y * sx + x;
        unsigned char minVal = 255;
        for (int dz = -1; dz <= 1; dz++) {
          for (int dy = -1; dy <= 1; dy++) {
            for (int dx_ = -1; dx_ <= 1; dx_++) {
              int nx = x + dx_;
              int ny = y + dy;
              int nz = z + dz;
              int nidx = nz * sx * sy + ny * sx + nx;
              if (vol[nidx] < minVal) minVal = vol[nidx];
            }
          }
        }
        eroded[idx] = minVal;
      }
    }
  }
  // Dilate
  for (int z = 1; z < sz - 1; z++) {
    for (int y = 1; y < sy - 1; y++) {
      for (int x = 1; x < sx - 1; x++) {
        int idx = z * sx * sy + y * sx + x;
        unsigned char maxVal = 0;
        for (int dz = -1; dz <= 1; dz++) {
          for (int dy = -1; dy <= 1; dy++) {
            for (int dx_ = -1; dx_ <= 1; dx_++) {
              int nx = x + dx_;
              int ny = y + dy;
              int nz = z + dz;
              int nidx = nz * sx * sy + ny * sx + nx;
              if (eroded[nidx] > maxVal) maxVal = eroded[nidx];
            }
          }
        }
        vol[idx] = maxVal;  // final opened volume
      }
    }
  }
  delete[] eroded;
}

static int computeOtsuThreshold(const unsigned char *data, int length) {
  // Compute histogram
  int hist[256];
  memset(hist, 0, sizeof(hist));
  for (int i = 0; i < length; i++) hist[data[i]]++;

  // Calculate total number of pixels
  int total = length;

  float sum = 0.0f;
  for (int t = 0; t < 256; t++) sum += t * hist[t];

  float sumB = 0.0f;
  int wB = 0;
  int wF = 0;
  float varMax = 0.0f;
  int threshold = 0;

  for (int t = 0; t < 256; t++) {
    wB += hist[t];
    if (wB == 0) continue;
    wF = total - wB;
    if (wF == 0) break;
    sumB += (float)(t * hist[t]);
    float mB = sumB / wB;
    float mF = (sum - sumB) / wF;
    float varBetween = (float)wB * (float)wF * (mB - mF) * (mB - mF);
    if (varBetween > varMax) {
      varMax = varBetween;
      threshold = t;
    }
  }
  return threshold;
}

void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex) {
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

  // Use sphere window type (1) and window size based on soma radius
  if (compute_sphere_win3d_pca_eigVec(
          img3d, N, M, P, x, y, z,  // Center on soma
          2 * r, 2 * r, 2 * r,      // Window size based on radius
          pc1, pc2, pc3, vec1, vec2, vec3)) {
    // Print results for this soma
    printf("\nSoma #%d PCA Results:\n", somaIndex);
    printf("  Center: (%.1f, %.1f, %.1f)\n", x, y, z);
    printf("  Eigenvalues:\n");
    printf("    pc1: %f\n", pc1);
    printf("    pc2: %f\n", pc2);
    printf("    pc3: %f\n", pc3);
    printf("  Principal axes:\n");
    printf("    v1: [%f, %f, %f]\n", vec1[0], vec1[1], vec1[2]);
    printf("    v2: [%f, %f, %f]\n", vec2[0], vec2[1], vec2[2]);
    printf("    v3: [%f, %f, %f]\n\n\n", vec3[0], vec3[1], vec3[2]);
  } else {
    printf("\nSoma #%d PCA failed.\n", somaIndex);
  }

  // Cleanup 3D array wrapper
  for (V3DLONG k = 0; k < P; k++) {
    delete[] img3d[k];
  }
  delete[] img3d;
}
