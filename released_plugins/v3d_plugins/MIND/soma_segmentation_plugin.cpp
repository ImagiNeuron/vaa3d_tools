/* soma_segmentation_plugin.cpp
 * This plugin allows the segmentation of individual somas with 3D flooding
 * algorithms. It generates ground truth data to train machine learning
 * algorithms for automatic soma segmnentation in the brain
 *
 * 2024-11-16 : by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */

#include "soma_segmentation_plugin.h"

#include <QInputDialog>
#include <QMessageBox>
#include <vector>
#include <queue>
#include <cmath>

#include "basic_surf_objs.h"
#include "v3d_message.h"

using namespace std;

struct input_PARA
{
  QString inimg_file;
  V3DLONG channel;
};

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu);

/**
 * @brief Menu option under the MIND plugins
 */
QStringList SomaSegmentation::menulist() const
{
  return QStringList() << tr("soma_segmentation") << tr("about");
}

/**
 * @brief Function list for the soma segmentation plugin
 */
QStringList SomaSegmentation::funclist() const
{
  return QStringList() << tr("segment_somas") << tr("help");
}

/**
 * @brief Display a description about the soma segmentation plugin and what it
 * does
 *
 * @param menu_name - the name of the menu being described
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 */
void SomaSegmentation::domenu(const QString &menu_name,
                              V3DPluginCallback2 &callback, QWidget *parent)
{
  if (menu_name == tr("soma_segmentation"))
  {
    bool bmenu = true;
    input_PARA PARA;
    reconstruction_func(callback, parent, PARA, bmenu);
  }
  else
  {
    v3d_msg(tr(
        "This plugin allows the segmentation of individual somas with 3D "
        "flooding algorithms. It generates ground truth data to train machine "
        "learning algorithms for automatic soma segmnentation in the brain. "
        "Developed by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous "
        "and Thibaut Baguette, 2024-11-16"));
  }
}

/**
 * @brief Function to segment somas
 *
 * Describe how function works here
 *
 * @param func_name - the name of the function
 * @param input - the input arguments
 * @param output - the output arguments
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @return true if the function was successful, false otherwise
 */
bool SomaSegmentation::dofunc(const QString &func_name,
                              const V3DPluginArgList &input,
                              V3DPluginArgList &output,
                              V3DPluginCallback2 &callback, QWidget *parent)
{
  if (func_name == tr("segment_somas"))
  {
    bool bmenu = false;
    input_PARA PARA;

    vector<char *> *pinfiles =
        (input.size() >= 1) ? (vector<char *> *)input[0].p : 0;
    vector<char *> *pparas =
        (input.size() >= 2) ? (vector<char *> *)input[1].p : 0;
    vector<char *> infiles = (pinfiles != 0) ? *pinfiles : vector<char *>();
    vector<char *> paras = (pparas != 0) ? *pparas : vector<char *>();

    if (infiles.empty())
    {
      fprintf(stderr, "Need input image. \n");
      return false;
    }
    else
      PARA.inimg_file = infiles[0];
    int k = 0;
    PARA.channel = (paras.size() >= k + 1) ? atoi(paras[k]) : 1;
    k++;
    reconstruction_func(callback, parent, PARA, bmenu);
  }
  else if (func_name == tr("help"))
  {
    ////HERE IS WHERE THE DEVELOPERS SHOULD UPDATE THE USAGE OF THE PLUGIN

    printf("**** Usage of soma_segmentation tracing **** \n");
    printf(
        "vaa3d -x soma_segmentation -f segment_somas -i <inimg_file> -p "
        "<channel> <other parameters>\n");
    printf("inimg_file       The input image\n");
    printf(
        "channel          Data channel for tracing. Start from 1 (default "
        "1).\n");

    printf(
        "outswc_file      Will be named automatically based on the input image "
        "file name, so you don't have to specify it.\n\n");
  }
  else
    return false;

  return true;
}

/**
 * @brief Function to reconstruct somas
 *
 * Describe how function works here
 *
 * @param callback - the V3D plugin callback interface
 * @param parent - the parent interface
 * @param PARA - the input parameters
 * @param bmenu - whether the function is being called from the menu
 */
void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent,
                         input_PARA &PARA, bool bmenu)
{
  unsigned char *data1d = 0;
  V3DLONG N, M, P, sc, c;
  V3DLONG in_sz[4];
  if (bmenu)
  {
    v3dhandle curwin = callback.currentImageWindow();
    if (!curwin)
    {
      QMessageBox::information(
          0, "", "You don't have any image open in the main window.");
      return;
    }

    Image4DSimple *p4DImage = callback.getImage(curwin);

    if (!p4DImage)
    {
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

    if (sc == 1)
    {
      c = 1;
      ok1 = true;
    }
    else
    {
      c = QInputDialog::getInt(parent, "Channel", "Enter channel NO:", 1, 1, sc,
                               1, &ok1);
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
    int datatype = 0;
    if (!simple_loadimage_wrapper(callback,
                                  PARA.inimg_file.toStdString().c_str(), data1d,
                                  in_sz, datatype))
    {
      fprintf(stderr,
              "Error happens in reading the subject file [%s]. Exit. \n",
              PARA.inimg_file.toStdString().c_str());
      return;
    }
    if (PARA.channel < 1 || PARA.channel > in_sz[3])
    {
      fprintf(stderr, "Invalid channel number. \n");
      return;
    }
    N = in_sz[0];
    M = in_sz[1];
    P = in_sz[2];
    sc = in_sz[3];
    c = PARA.channel;
  }

  //// THIS IS WHERE THE DEVELOPERS SHOULD ADD THEIR OWN NEURON TRACING CODE

  // get current window, image, and landmarks
  v3dhandle curwin = callback.currentImageWindow();

  // Add null checks for curwin
  if (!curwin)
  {
    v3d_msg("No image window is currently open.", bmenu);
    return;
  }

  Image4DSimple *p4DImage = callback.getImage(curwin);
  if (!p4DImage)
  {
    v3d_msg("Invalid image pointer.", bmenu);
    return;
  }

  LandmarkList landmarkList = callback.getLandmark(curwin);

  // Check if landmarkList is empty
  if (landmarkList.isEmpty())
  {
    v3d_msg("No landmarks defined. Please define at least one landmark.", bmenu);
    return;
  }

  // reset ROI
  ROIList roiList = callback.getROI(curwin);
  for (int j = 0; j < 3; j++)
  {
    roiList[j].clear();
  }

  // get 1st landmark
  LocationSimple lm = landmarkList[0];
  float x, y, z;
  float radius;
  x = lm.x;
  y = lm.y;
  z = lm.z;
  radius = lm.radius;

  v3d_msg(QString("Landmark: x=%1, y=%2, z=%3, radius=%4")
              .arg(x)
              .arg(y)
              .arg(z)
              .arg(radius),
          bmenu);

  // set ROI
  // ROIList being a QList<QPolygon>, and QPolygon being a QVector<QPoint>,
  // we represent the x, y, and z planes as polygons with 4 points each to
  // form a cube with the landmark at the center
  float x_min = x - radius * 2.0f;
  float x_max = x + radius * 2.0f;
  float y_min = y - radius * 2.0f;
  float y_max = y + radius * 2.0f;
  float z_min = z - radius * 2.0f;
  float z_max = z + radius * 2.0f;
  // x-y plane
  roiList[0] << QPoint(x_min, y_min);
  roiList[0] << QPoint(x_max, y_min);
  roiList[0] << QPoint(x_max, y_max);
  roiList[0] << QPoint(x_min, y_max);
  // z-y plane
  roiList[1] << QPoint(z_min, y_min);
  roiList[1] << QPoint(z_max, y_min);
  roiList[1] << QPoint(z_max, y_max);
  roiList[1] << QPoint(z_min, y_max);
  // x-z plane
  roiList[2] << QPoint(x_min, z_min);
  roiList[2] << QPoint(x_max, z_min);
  roiList[2] << QPoint(x_max, z_max);
  roiList[2] << QPoint(x_min, z_max);

  if (callback.setROI(curwin, roiList))
  {
    callback.updateImageWindow(curwin);
  }
  else
  {
    qDebug() << "error: failed to set ROI";
    return;
  }

  callback.openROI3DWindow(curwin);

  // Perform flood fill on all landmarks
  QList<MyMarker> myMarkerList;

  // Convert LandmarkList to QList<MyMarker>
  for (const auto &lm : landmarkList)
  {
    MyMarker marker;
    marker.x = lm.x;
    marker.y = lm.y;
    marker.z = lm.z;
    marker.radius = lm.radius;
    myMarkerList.append(marker);
  }

  QList<My4DImage *> floodedImages;
  int threshold_value = 25; // Increased base threshold

  for (const auto &marker : myMarkerList)
  {
    // Create a deep copy of the image data
    My4DImage *myImage = new My4DImage();
    myImage->xdim = p4DImage->getXDim();
    myImage->ydim = p4DImage->getYDim();
    myImage->zdim = p4DImage->getZDim();
    myImage->cdim = 1; // Assuming single channel for processing

    V3DLONG totalSize = myImage->xdim * myImage->ydim * myImage->zdim * myImage->cdim;
    myImage->data = new unsigned char[totalSize];

    // Copy the data of the selected channel
    unsigned char *channelData = p4DImage->getRawDataAtChannel(c - 1);
    if (!channelData)
    {
      v3d_msg("Failed to get image data for the specified channel.", bmenu);
      delete myImage;
      continue;
    }
    memcpy(myImage->data, channelData, totalSize);

    // Declare variables to hold edges
    int min_x, max_x, min_y, max_y, min_z, max_z;

    // Perform flood fill and get edges
    My4DImage *flooded = performFloodFill(myImage, marker, threshold_value, min_x, max_x, min_y, max_y, min_z, max_z);
    delete myImage; // Safe to delete since we have copied the data

    // Check if flood fill was successful
    if (flooded)
    {
      floodedImages.append(flooded);

      // Print edges of the detected soma
      printf("Detected Soma Edges:\n");
      printf("X: [%d, %d]\n", min_x, max_x);
      printf("Y: [%d, %d]\n", min_y, max_y);
      printf("Z: [%d, %d]\n", min_z, max_z);
    }
    else
    {
      printf("Flood fill failed for marker at (%f, %f, %f).\n", marker.x, marker.y, marker.z);
    }
  }

  My4DImage *finalFloodedImage = mergeFloodedImages(floodedImages);

  // Instead of saving the flooded image, print its details
  if (finalFloodedImage)
  {
    printf("Final Flooded Image Details:\n");
    printf("Dimensions: %ld x %ld x %ld x %ld\n",
           finalFloodedImage->xdim,
           finalFloodedImage->ydim,
           finalFloodedImage->zdim,
           finalFloodedImage->cdim);
    // Optionally, print a summary of the data
    size_t flooded_voxels = 0;
    for (V3DLONG idx = 0; idx < finalFloodedImage->xdim * finalFloodedImage->ydim * finalFloodedImage->zdim * finalFloodedImage->cdim; ++idx)
    {
      if (finalFloodedImage->data[idx] == 255)
        flooded_voxels++;
    }
    printf("Number of Flooded Voxels: %zu\n", flooded_voxels);
  }
  else
  {
    printf("No flooded image generated.\n");
  }

  // Clean up
  for (auto img : floodedImages)
    delete img;
  delete finalFloodedImage;

  v3d_msg("Flood fill completed. Flooded image details printed to console.", bmenu);

  return;
}

My4DImage *performFloodFill(const My4DImage *image, const MyMarker &marker, int threshold, int &min_x, int &max_x, int &min_y, int &max_y, int &min_z, int &max_z)
{
  if (!image || !image->data)
    return nullptr;

  V3DLONG xdim = image->xdim;
  V3DLONG ydim = image->ydim;
  V3DLONG zdim = image->zdim;
  V3DLONG cdim = image->cdim;

  if (cdim < 1)
    return nullptr;

  My4DImage *floodedImage = new My4DImage(*image);
  unsigned char *floodData = floodedImage->data;

  int x0 = static_cast<int>(std::round(marker.x));
  int y0 = static_cast<int>(std::round(marker.y));
  int z0 = static_cast<int>(std::round(marker.z));

  if (x0 < 0 || x0 >= xdim || y0 < 0 || y0 >= ydim || z0 < 0 || z0 >= zdim)
  {
    delete floodedImage;
    return nullptr;
  }

  int seedIdx = z0 * ydim * xdim + y0 * xdim + x0;
  unsigned char seedValue = floodData[seedIdx];

  int neighborhood = 5;
  double localMean = 0;
  double localMax = 0;
  double localMin = 255;
  int count = 0;

  for (int dz = -neighborhood; dz <= neighborhood; dz++)
  {
    for (int dy = -neighborhood; dy <= neighborhood; dy++)
    {
      for (int dx = -neighborhood; dx <= neighborhood; dx++)
      {
        int nx = x0 + dx;
        int ny = y0 + dy;
        int nz = z0 + dz;

        if (nx >= 0 && nx < xdim && ny >= 0 && ny < ydim && nz >= 0 && nz < zdim)
        {
          unsigned char val = floodData[nz * ydim * xdim + ny * xdim + nx];
          localMean += val;
          localMax = std::max(localMax, (double)val);
          localMin = std::min(localMin, (double)val);
          count++;
        }
      }
    }
  }
  localMean /= count;

  printf("Local statistics: mean=%.2f, min=%.2f, max=%.2f\n", localMean, localMin, localMax);

  double dynamicThreshold = std::max(threshold, static_cast<int>(marker.radius * 0.8));
  int lowerThreshold = std::max(0, static_cast<int>(localMin - dynamicThreshold));
  int upperThreshold = std::min(255, static_cast<int>(localMax + dynamicThreshold));

  printf("Thresholds: lower=%d, upper=%d, dynamic=%.2f\n", lowerThreshold, upperThreshold, dynamicThreshold);

  size_t maxRegionSize = static_cast<size_t>(4.0 / 3.0 * M_PI * pow(marker.radius * 3, 3));
  size_t minRegionSize = 5;

  std::vector<bool> visited(xdim * ydim * zdim, false);

  std::queue<V3DLONG> q;
  q.push(seedIdx);
  visited[seedIdx] = true;

  int directions[6][3] = {
      {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};

  min_x = x0;
  max_x = x0;
  min_y = y0;
  max_y = y0;
  min_z = z0;
  max_z = z0;

  size_t regionSize = 0;

  while (!q.empty())
  {
    V3DLONG current = q.front();
    q.pop();

    int cz = current / (xdim * ydim);
    int cy = (current % (xdim * ydim)) / xdim;
    int cx = current % xdim;

    double distance = std::sqrt(pow(cx - x0, 2) + pow(cy - y0, 2) + pow(cz - z0, 2));
    if (distance > marker.radius * 4)
      continue;

    for (int i = 0; i < 6; ++i)
    {
      int nx = cx + directions[i][0];
      int ny = cy + directions[i][1];
      int nz = cz + directions[i][2];

      if (nx < 0 || nx >= xdim || ny < 0 || ny >= ydim || nz < 0 || nz >= zdim)
        continue;

      V3DLONG neighborIdx = nz * ydim * xdim + ny * xdim + nx;

      if (!visited[neighborIdx])
      {
        unsigned char neighborValue = floodData[neighborIdx];

        double gradient = std::abs(static_cast<double>(neighborValue) - floodData[current]);

        if (neighborValue >= lowerThreshold &&
            neighborValue <= upperThreshold &&
            gradient < dynamicThreshold * 2)
        {

          q.push(neighborIdx);
          visited[neighborIdx] = true;
          floodData[neighborIdx] = 255;

          if (nx < min_x)
            min_x = nx;
          if (nx > max_x)
            max_x = nx;
          if (ny < min_y)
            min_y = ny;
          if (ny > max_y)
            max_y = ny;
          if (nz < min_z)
            min_z = nz;
          if (nz > max_z)
            max_z = nz;

          regionSize++;
          if (regionSize > maxRegionSize)
          {
            printf("Region size exceeded for marker at (%f, %f, %f). Size: %zu, Max: %zu\n",
                   marker.x, marker.y, marker.z, regionSize, maxRegionSize);
            continue;
          }
        }
      }
    }
  }

  if (regionSize < minRegionSize)
  {
    printf("Region too small: %zu voxels (minimum: %zu)\n", regionSize, minRegionSize);
    delete floodedImage;
    return nullptr;
  }

  printf("Successfully flooded region with %zu voxels\n", regionSize);
  return floodedImage;
}

My4DImage *mergeFloodedImages(const QList<My4DImage *> &floodedImages)
{
  if (floodedImages.isEmpty())
    return nullptr;

  My4DImage *mergedImage = new My4DImage(*floodedImages[0]);

  for (int i = 1; i < floodedImages.size(); ++i)
  {
    for (V3DLONG idx = 0; idx < mergedImage->xdim * mergedImage->ydim * mergedImage->zdim * mergedImage->cdim; ++idx)
    {
      mergedImage->data[idx] |= floodedImages[i]->data[idx];
    }
  }

  return mergedImage;
}

// Ensure libtiff is linked correctly

My4DImage::My4DImage() : data(nullptr), xdim(0), ydim(0), zdim(0), cdim(0) {}

My4DImage::My4DImage(const My4DImage &other)
    : xdim(other.xdim), ydim(other.ydim), zdim(other.zdim), cdim(other.cdim)
{
  if (other.data)
  {
    data = new unsigned char[xdim * ydim * zdim * cdim];
    std::copy(other.data, other.data + xdim * ydim * zdim * cdim, data);
  }
  else
  {
    data = nullptr;
  }
}

My4DImage::~My4DImage()
{
  delete[] data;
}

My4DImage &My4DImage::operator=(const My4DImage &other)
{
  if (this != &other)
  {
    delete[] data; // Free existing data

    xdim = other.xdim;
    ydim = other.ydim;
    zdim = other.zdim;
    cdim = other.cdim;

    if (other.data)
    {
      V3DLONG totalSize = xdim * ydim * zdim * cdim;
      data = new unsigned char[totalSize];
      std::copy(other.data, other.data + totalSize, data);
    }
    else
    {
      data = nullptr;
    }
  }
  return *this;
}
