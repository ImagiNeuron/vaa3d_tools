/* soma_segmentation_plugin.cpp
 * A plugin for seeded soma segmentation using threshold + region growing.
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"
#include <QInputDialog>
#include <QMessageBox>
#include <vector>
#include <cmath>
#include <algorithm> // For std::sort and std::nth_element
#include <queue>
#include <cstring> // For memcpy

#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

/**************************************
 * My4DImage implementations
 **************************************/
My4DImage::My4DImage()
    : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0)
{
}

My4DImage::My4DImage(const My4DImage &other)
    : xdim(other.xdim), ydim(other.ydim), zdim(other.zdim), cdim(other.cdim)
{
  V3DLONG totalSize = xdim * ydim * zdim * cdim;
  if (other.data && totalSize > 0)
  {
    data = new unsigned char[totalSize];
    memcpy(data, other.data, totalSize);
  }
  else
  {
    data = nullptr;
  }
}

My4DImage &My4DImage::operator=(const My4DImage &other)
{
  if (this != &other)
  {
    delete[] data;
    xdim = other.xdim;
    ydim = other.ydim;
    zdim = other.zdim;
    cdim = other.cdim;

    V3DLONG totalSize = xdim * ydim * zdim * cdim;
    if (other.data && totalSize > 0)
    {
      data = new unsigned char[totalSize];
      memcpy(data, other.data, totalSize);
    }
    else
    {
      data = nullptr;
    }
  }
  return *this;
}

My4DImage::~My4DImage()
{
  delete[] data;
}

/**************************************
 * Internal struct to store params
 **************************************/
struct input_PARA
{
  QString inimg_file;
  V3DLONG channel;
};

/**************************************
 * Forward declaration of main
 **************************************/
void reconstruction_func(V3DPluginCallback2 &callback,
                         QWidget *parent,
                         input_PARA &PARA,
                         bool bmenu);

/**************************************
 * Plugin Interface Methods
 **************************************/
QStringList SomaSegmentation::menulist() const
{
  return QStringList()
         << tr("soma_segmentation")
         << tr("about");
}

QStringList SomaSegmentation::funclist() const
{
  return QStringList()
         << tr("segment_somas")
         << tr("help");
}

void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback,
                              QWidget *parent)
{
  if (menu_name == tr("soma_segmentation"))
  {
    bool bmenu = true;
    input_PARA PARA;
    reconstruction_func(callback, parent, PARA, bmenu);
  }
  else
  {
    v3d_msg(tr("Soma segmentation plugin using threshold + region growing.\n"
               "Originally by ImagiNeuron, 2024-11-16"));
  }
}

bool SomaSegmentation::dofunc(const QString &func_name,
                              const V3DPluginArgList &input,
                              V3DPluginArgList &output,
                              V3DPluginCallback2 &callback,
                              QWidget *parent)
{
  if (func_name == tr("segment_somas"))
  {
    bool bmenu = false;
    input_PARA PARA;

    vector<char *> *pinfiles = (input.size() >= 1) ? (vector<char *> *)input[0].p : 0;
    vector<char *> *pparas = (input.size() >= 2) ? (vector<char *> *)input[1].p : 0;
    vector<char *> infiles = (pinfiles != 0) ? *pinfiles : vector<char *>();
    vector<char *> paras = (pparas != 0) ? *pparas : vector<char *>();

    if (infiles.empty())
    {
      fprintf(stderr, "Need input image.\n");
      return false;
    }
    else
    {
      PARA.inimg_file = infiles[0];
    }

    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;

    reconstruction_func(callback, parent, PARA, bmenu);
  }
  else if (func_name == tr("help"))
  {
    printf("**** Usage of soma_segmentation ****\n");
    printf("vaa3d -x soma_segmentation -f segment_somas -i <inimg_file> -p <channel>\n");
  }
  else
  {
    return false;
  }
  return true;
}

/**************************************
 *         HELPER DECLARATIONS
 **************************************/
void applyMedianFilter(const unsigned char *inputData,
                       unsigned char *outputData,
                       V3DLONG N, V3DLONG M, V3DLONG P,
                       int windowSize);

void applyGaussianFilter(const unsigned char *inputData,
                         unsigned char *outputData,
                         V3DLONG N, V3DLONG M, V3DLONG P,
                         float sigma);

static int computeOtsuThreshold(const unsigned char *data, int length);

void morphologicalOpen3D(unsigned char *vol, int sx, int sy, int sz);

void extractSubvolume(const unsigned char *inData,
                      unsigned char *outData,
                      V3DLONG N, V3DLONG M, V3DLONG P,
                      int x1, int x2,
                      int y1, int y2,
                      int z1, int z2);

/// Region-growing from a seed (seedX, seedY, seedZ) in subvol
void regionGrow3D(const unsigned char *subvol,
                  unsigned char *label3D,
                  int sx, int sy, int sz,
                  int seedX, int seedY, int seedZ);

/**************************************
 * The main function
 **************************************/
void reconstruction_func(V3DPluginCallback2 &callback,
                         QWidget *parent,
                         input_PARA &PARA,
                         bool bmenu)
{
  /*************************************
   * 1) Load Image
   *************************************/
  unsigned char *data1d = 0;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];

  if (bmenu)
  {
    // From the current Vaa3D window
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin)
    {
      QMessageBox::information(0, "", "No image open in the main window.");
      return;
    }

    Image4DSimple *p4DImage = callback.getImage(curwin);
    if (!p4DImage)
    {
      QMessageBox::information(0, "", "Invalid image pointer.");
      return;
    }

    data1d = p4DImage->getRawData();
    N = p4DImage->getXDim();
    M = p4DImage->getYDim();
    P = p4DImage->getZDim();
    sc = p4DImage->getCDim();

    bool ok1;
    if (sc == 1)
    {
      c = 1;
      ok1 = true;
    }
    else
    {
      c = QInputDialog::getInt(parent, "Channel",
                               "Enter channel NO:",
                               1, 1, sc, 1, &ok1);
    }
    if (!ok1)
      return;

    in_sz[0] = N;
    in_sz[1] = M;
    in_sz[2] = P;
    in_sz[3] = sc;

    PARA.inimg_file = p4DImage->getFileName();
  }
  else
  {
    // Command-line input
    int datatype = 0;
    if (!simple_loadimage_wrapper(callback,
                                  PARA.inimg_file.toStdString().c_str(),
                                  data1d, in_sz, datatype))
    {
      fprintf(stderr, "Error loading file [%s].\n",
              PARA.inimg_file.toStdString().c_str());
      return;
    }
    if (PARA.channel < 1 || PARA.channel > in_sz[3])
    {
      fprintf(stderr, "Invalid channel number.\n");
      return;
    }
    N = in_sz[0];
    M = in_sz[1];
    P = in_sz[2];
    sc = in_sz[3];
    c = PARA.channel;
  }

  if (!data1d)
  {
    v3d_msg("No valid image data!", bmenu);
    return;
  }

  // If multiple channels, keep only the one of interest
  if (sc > 1)
  {
    v3d_msg("Grayscale only: ignoring extra channels except the chosen one.");
    sc = 1;
    c = 1;
  }

  /*************************************
   * 2) Preprocessing (Median + Gaussian)
   *************************************/
  {
    // Extract the desired channel if needed
    V3DLONG totalSize = N * M * P;
    unsigned char *channelData = new unsigned char[totalSize];
    {
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

    data1d = gaussData; // now we use this preprocessed version
  }

  /*************************************
   * 3) For each landmark, do local segmentation
   *************************************/
  v3dhandle curwin = callback.currentImageWindow();
  LandmarkList landmarkList;
  if (bmenu && curwin)
  {
    landmarkList = callback.getLandmark(curwin);
  }
  else
  {
    // In command-line or no open window, you would specify landmarks some other way
  }

  if (landmarkList.isEmpty())
  {
    v3d_msg("No landmarks found. Please provide at least one landmark.", bmenu);
    delete[] data1d;
    return;
  }

  // Prepare a global label volume (unsigned char) for the final binary result
  V3DLONG totalSize = N * M * P;
  unsigned char *labelVolume = new unsigned char[totalSize];
  memset(labelVolume, 0, totalSize * sizeof(unsigned char)); // 0 = background

  for (int i = 0; i < landmarkList.size(); i++)
  {
    LocationSimple lm = landmarkList[i];
    float x = lm.x;
    float y = lm.y;
    float z = lm.z;
    float r = lm.radius > 0 ? lm.radius : 5.0f; // fallback radius

    // bounding box
    int x1 = std::max(0, (int)floor(x - r));
    int x2 = std::min((int)N - 1, (int)ceil(x + r));
    int y1 = std::max(0, (int)floor(y - r));
    int y2 = std::min((int)M - 1, (int)ceil(y + r));
    int z1 = std::max(0, (int)floor(z - r));
    int z2 = std::min((int)P - 1, (int)ceil(z + r));

    int sx = x2 - x1 + 1;
    int sy = y2 - y1 + 1;
    int sz_ = z2 - z1 + 1;

    // Extract subvolume
    unsigned char *subvol = new unsigned char[sx * sy * sz_];
    extractSubvolume(data1d, subvol, N, M, P, x1, x2, y1, y2, z1, z2);

    // Zero out voxels that lie outside the sphere defined by landmark radius
    int seedX = int(x - x1 + 0.5f);
    int seedY = int(y - y1 + 0.5f);
    int seedZ = int(z - z1 + 0.5f);

    for (int zz = 0; zz < sz_; zz++)
    {
      for (int yy = 0; yy < sy; yy++)
      {
        for (int xx = 0; xx < sx; xx++)
        {
          double dx = xx - seedX;
          double dy = yy - seedY;
          double dz = zz - seedZ;
          if ((dx * dx + dy * dy + dz * dz) > (r * r))
          {
            subvol[zz * sx * sy + yy * sx + xx] = 0;
          }
        }
      }
    }

    // Compute Otsu threshold in the subvolume
    int otsuValue = computeOtsuThreshold(subvol, sx * sy * sz_);
    int seedIdx = (seedZ * sx * sy) + (seedY * sx) + seedX;
    unsigned char seedVal = subvol[seedIdx];

    // Binarize depending on whether the seed is below or above the threshold
    for (int idx = 0; idx < sx * sy * sz_; idx++)
    {
      bool isForeground;
      if (seedVal < otsuValue)
      {
        isForeground = (subvol[idx] <= otsuValue);
      }
      else
      {
        isForeground = (subvol[idx] >= otsuValue);
      }
      subvol[idx] = isForeground ? 255 : 0;
    }

    // If the seed is background after binarization, skip
    if (subvol[seedIdx] != 255)
    {
      continue;
    }

    // Morphological opening
    morphologicalOpen3D(subvol, sx, sy, sz_);

    // Prepare a local label array
    unsigned char *subLabel = new unsigned char[sx * sy * sz_];
    memset(subLabel, 0, sx * sy * sz_);

    // Region-growing from the landmark center
    if (seedX >= 0 && seedX < sx &&
        seedY >= 0 && seedY < sy &&
        seedZ >= 0 && seedZ < sz_)
    {
      regionGrow3D(subvol, subLabel, sx, sy, sz_, seedX, seedY, seedZ);
    }

    // Copy the region back to the global labelVolume as 255
    for (int zz = 0; zz < sz_; zz++)
    {
      for (int yy = 0; yy < sy; yy++)
      {
        for (int xx = 0; xx < sx; xx++)
        {
          int subIdx = zz * (sx * sy) + yy * sx + xx;
          if (subLabel[subIdx] > 0)
          {
            V3DLONG gx = x1 + xx;
            V3DLONG gy = y1 + yy;
            V3DLONG gz = z1 + zz;
            V3DLONG idxGlobal = gz * (N * M) + gy * N + gx;
            labelVolume[idxGlobal] = 255;
          }
        }
      }
    }

    delete[] subvol;
    delete[] subLabel;
  }

  /*************************************
   * 4) Save the final segmentation
   *************************************/
  QString savePath = QFileDialog::getSaveFileName(parent,
                                                  "Save binary segmented image",
                                                  "",
                                                  "TIFF Files (*.tiff *.tif)");
  if (!savePath.isEmpty())
  {
    V3DLONG outSZ[4] = {N, M, P, 1};
    simple_saveimage_wrapper(callback,
                             savePath.toStdString().c_str(),
                             labelVolume,
                             outSZ,
                             1);
    v3d_msg("Binary segmented image saved.", bmenu);
  }

  delete[] labelVolume;
  delete[] data1d;

  v3d_msg("Local threshold + region-growing segmentation complete.", bmenu);
}

/*************************************
 *        MEDIAN FILTER
 *************************************/
void applyMedianFilter(const unsigned char *inputData,
                       unsigned char *outputData,
                       V3DLONG N, V3DLONG M, V3DLONG P,
                       int windowSize)
{
  int halfWindow = windowSize / 2;
  vector<unsigned char> window;
  window.reserve(windowSize * windowSize * windowSize);

  for (V3DLONG z = 0; z < P; z++)
  {
    for (V3DLONG y = 0; y < M; y++)
    {
      for (V3DLONG x = 0; x < N; x++)
      {
        window.clear();
        for (int dz = -halfWindow; dz <= halfWindow; dz++)
        {
          for (int dy = -halfWindow; dy <= halfWindow; dy++)
          {
            for (int dx_ = -halfWindow; dx_ <= halfWindow; dx_++)
            {
              V3DLONG nx = x + dx_;
              V3DLONG ny = y + dy;
              V3DLONG nz = z + dz;
              if (nx >= 0 && nx < N && ny >= 0 && ny < M && nz >= 0 && nz < P)
              {
                V3DLONG idx = nz * (M * N) + ny * N + nx;
                window.push_back(inputData[idx]);
              }
            }
          }
        }
        // find median
        std::nth_element(window.begin(),
                         window.begin() + window.size() / 2,
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
                         unsigned char *outputData,
                         V3DLONG N, V3DLONG M, V3DLONG P,
                         float sigma)
{
  // For simplicity, use a 5x5x5 kernel
  int windowSize = 5;
  int halfWindow = windowSize / 2;

  vector<float> kernel(windowSize * windowSize * windowSize);
  float sumKernel = 0.0f;
  int idx = 0;
  for (int dz = -halfWindow; dz <= halfWindow; dz++)
  {
    for (int dy = -halfWindow; dy <= halfWindow; dy++)
    {
      for (int dx_ = -halfWindow; dx_ <= halfWindow; dx_++)
      {
        float val = expf(-(dx_ * dx_ + dy * dy + dz * dz) / (2 * sigma * sigma));
        kernel[idx++] = val;
        sumKernel += val;
      }
    }
  }
  // normalize
  for (size_t i = 0; i < kernel.size(); i++)
    kernel[i] /= sumKernel;

  // convolve
  for (V3DLONG z = 0; z < P; z++)
  {
    for (V3DLONG y = 0; y < M; y++)
    {
      for (V3DLONG x = 0; x < N; x++)
      {
        float accum = 0.0f;
        idx = 0;
        for (int dz = -halfWindow; dz <= halfWindow; dz++)
        {
          for (int dy = -halfWindow; dy <= halfWindow; dy++)
          {
            for (int dx_ = -halfWindow; dx_ <= halfWindow; dx_++)
            {
              V3DLONG nx = x + dx_;
              V3DLONG ny = y + dy;
              V3DLONG nz = z + dz;
              if (nx >= 0 && nx < N && ny >= 0 && ny < M && nz >= 0 && nz < P)
              {
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
 *    EXTRACT SUBVOLUME
 *************************************/
void extractSubvolume(const unsigned char *inData,
                      unsigned char *outData,
                      V3DLONG N, V3DLONG M, V3DLONG P,
                      int x1, int x2,
                      int y1, int y2,
                      int z1, int z2)
{
  int sx = x2 - x1 + 1;
  int sy = y2 - y1 + 1;
  int sz = z2 - z1 + 1;

  for (int zz = z1; zz <= z2; zz++)
  {
    for (int yy = y1; yy <= y2; yy++)
    {
      for (int xx = x1; xx <= x2; xx++)
      {
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

/*************************************
 * 3D MORPHOLOGICAL OPENING
 *************************************/
static void dilate3D(const unsigned char *src, unsigned char *dst,
                     int sx, int sy, int sz)
{
  memset(dst, 0, sx * sy * sz);
  for (int z = 1; z < sz - 1; z++)
  {
    for (int y = 1; y < sy - 1; y++)
    {
      for (int x = 1; x < sx - 1; x++)
      {
        int idx = z * sx * sy + y * sx + x;
        unsigned char maxVal = 0;
        for (int dz = -1; dz <= 1; dz++)
        {
          for (int dy = -1; dy <= 1; dy++)
          {
            for (int dx_ = -1; dx_ <= 1; dx_++)
            {
              int nx = x + dx_;
              int ny = y + dy;
              int nz = z + dz;
              int nidx = nz * sx * sy + ny * sx + nx;
              if (src[nidx] > maxVal)
                maxVal = src[nidx];
            }
          }
        }
        dst[idx] = maxVal;
      }
    }
  }
}

static void erode3D(const unsigned char *src, unsigned char *dst,
                    int sx, int sy, int sz)
{
  memset(dst, 255, sx * sy * sz);
  for (int z = 1; z < sz - 1; z++)
  {
    for (int y = 1; y < sy - 1; y++)
    {
      for (int x = 1; x < sx - 1; x++)
      {
        int idx = z * sx * sy + y * sx + x;
        unsigned char minVal = 255;
        for (int dz = -1; dz <= 1; dz++)
        {
          for (int dy = -1; dy <= 1; dy++)
          {
            for (int dx_ = -1; dx_ <= 1; dx_++)
            {
              int nx = x + dx_;
              int ny = y + dy;
              int nz = z + dz;
              int nidx = nz * sx * sy + ny * sx + nx;
              if (src[nidx] < minVal)
                minVal = src[nidx];
            }
          }
        }
        dst[idx] = minVal;
      }
    }
  }
}

void morphologicalOpen3D(unsigned char *vol, int sx, int sy, int sz)
{
  unsigned char *eroded = new unsigned char[sx * sy * sz];
  unsigned char *dilated = new unsigned char[sx * sy * sz];

  erode3D(vol, eroded, sx, sy, sz);
  dilate3D(eroded, dilated, sx, sy, sz);

  memcpy(vol, dilated, sx * sy * sz);

  delete[] eroded;
  delete[] dilated;
}

/********************************************
 *         OTSU THRESHOLD
 ********************************************/
static int computeOtsuThreshold(const unsigned char *data, int length)
{
  int hist[256];
  memset(hist, 0, sizeof(hist));
  for (int i = 0; i < length; i++)
    hist[data[i]]++;

  int total = length;

  float sum = 0.0f;
  for (int t = 0; t < 256; t++)
    sum += (float)t * (float)hist[t];

  float sumB = 0.0f;
  int wB = 0;
  int wF = 0;
  float varMax = 0.0f;
  int threshold = 0;

  for (int t = 0; t < 256; t++)
  {
    wB += hist[t];
    if (wB == 0)
      continue;
    wF = total - wB;
    if (wF == 0)
      break;

    sumB += (float)(t * hist[t]);
    float mB = sumB / wB;
    float mF = (sum - sumB) / wF;
    float varBetween = (float)wB * (float)wF * (mB - mF) * (mB - mF);
    if (varBetween > varMax)
    {
      varMax = varBetween;
      threshold = t;
    }
  }
  return threshold;
}

/**********************************************
 *  SIMPLE 6-CONNECTED REGION GROWING
 **********************************************/
static const int NB_6[6][3] = {
    {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};

void regionGrow3D(const unsigned char *subvol,
                  unsigned char *label3D,
                  int sx, int sy, int sz,
                  int seedX, int seedY, int seedZ)
{
  std::queue<long> Q;

  auto idx3D = [&](int xx, int yy, int zz)
  {
    return (long)zz * (sx * sy) + (long)yy * sx + (long)xx;
  };

  long seedIdx = idx3D(seedX, seedY, seedZ);
  if (subvol[seedIdx] != 255)
    return;

  label3D[seedIdx] = 1;
  Q.push(seedIdx);

  while (!Q.empty())
  {
    long cur = Q.front();
    Q.pop();

    int z = (int)(cur / (sx * sy));
    int r = (int)(cur % (sx * sy));
    int y = r / sx;
    int x = r % sx;

    for (int i = 0; i < 6; i++)
    {
      int nx = x + NB_6[i][0];
      int ny = y + NB_6[i][1];
      int nz = z + NB_6[i][2];
      if (nx < 0 || nx >= sx ||
          ny < 0 || ny >= sy ||
          nz < 0 || nz >= sz)
        continue;

      long nIdx = idx3D(nx, ny, nz);
      if (label3D[nIdx] == 0 && subvol[nIdx] == 255)
      {
        label3D[nIdx] = 1;
        Q.push(nIdx);
      }
    }
  }
}
