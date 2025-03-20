#ifndef __MIND_PCA_ANALYSIS_H__
#define __MIND_PCA_ANALYSIS_H__

#include <v3d_interface.h>

#include <QApplication>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QLabel>
#include <QMessageBox>
#include <QVBoxLayout>
#include <QtGui>
#include <fstream>

#include "basic_4dimage.h"
#include "compute_win_pca_wp.h"

// Function prototypes

/**
 * @brief Analyze soma using PCA and save results to CSV file - return results
 */
void analyzeSomaPCAReturnResults(unsigned char *labeledData, V3DLONG N,
                                 V3DLONG M, V3DLONG P, const LocationSimple &lm,
                                 int somaIndex, QString savePath, double &pc1,
                                 double &pc2, double &pc3, double vec1[3],
                                 double vec2[3], double vec3[3],
                                 double &x_center, double &y_center,
                                 double &z_center);

/**
 * @brief Analyze soma using PCA and save results to CSV file - results not
 * returned
 */
void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex, QString savePath);

/**
 * @brief Save PCA results to CSV file
 */
void savePCAResultsToCSV(const QString &filename, int somaIndex,
                         const LocationSimple &lm, double pc1, double pc2,
                         double pc3, const double *vec1, const double *vec2,
                         const double *vec3, double x_center, double y_center,
                         double z_center, QWidget *parent, bool *saveEnabled);

/**
 * @brief draws a line to represent a PC
 */
void drawLine(Image4DSimple *image, double *from, double *to);

/**
 * @brief Visualize PCA results
 */
void visualizePCA_func(V3DPluginCallback2 &callback, QWidget *parent);

/**
 * @brief Simulate somas based on segmentation
 */
void simulate_somas(V3DPluginCallback2 &callback, QWidget *parent);

/**
 * @brief Create a background image based on segmentation threshold and return
 * as 3D array
 * @param callback V3DPluginCallback2 reference
 * @param parent Parent widget
 * @param dim_X Output parameter for X dimension
 * @param dim_Y Output parameter for Y dimension
 * @param dim_Z Output parameter for Z dimension
 * @return 3D array of background intensities, or nullptr on failure
 */
unsigned char ***create_background(V3DPluginCallback2 &callback,
                                   QWidget *parent, V3DLONG &dim_X,
                                   V3DLONG &dim_Y, V3DLONG &dim_Z);

/**
 * @brief Create background intensity values based on segmentation threshold
 * @param callback V3DPluginCallback2 reference
 * @param parent Parent widget
 * @param dimX Output parameter for X dimension
 * @param dimY Output parameter for Y dimension
 * @param dimZ Output parameter for Z dimension
 * @return 3D array of background intensity values (not binary segmentation).
 *         The values follow Gaussian distribution with mean and standard
 * deviation calculated from background voxels in the original image. The caller
 * is responsible for freeing this memory using free_3d_array(array, dimZ, dimY)
 */
unsigned char ***create_background(V3DPluginCallback2 &callback,
                                   QWidget *parent, V3DLONG &dimX,
                                   V3DLONG &dimY, V3DLONG &dimZ);

/**
 * @brief Create a background image based on segmentation threshold
 * @param callback V3DPluginCallback2 reference
 * @param parent Parent widget
 * @deprecated Use the version that returns unsigned char*** instead
 */
void create_background(V3DPluginCallback2 &callback, QWidget *parent);

/**
 * @brief Free memory allocated for a single 3D array
 * @param array 3D array to free
 * @param dim_Z Number of Z slices
 * @param dim_Y Number of rows
 */
void free_3d_array(unsigned char ***array, V3DLONG dim_Z, V3DLONG dim_Y);

/**
 * @brief Retrieve PCA information for a specific soma from a CSV file
 *
 * @param somaID The ID of the soma to retrieve information for
 * @param imageName The base image name to construct the CSV filename
 * @param pc1 Output parameter for first principal component eigenvalue
 * @param pc2 Output parameter for second principal component eigenvalue
 * @param pc3 Output parameter for third principal component eigenvalue
 * @param vec1 Output array for first principal component eigenvector (size 3)
 * @param vec2 Output array for second principal component eigenvector (size 3)
 * @param vec3 Output array for third principal component eigenvector (size 3)
 * @param x_center Output parameter for x-coordinate of center of mass
 * @param y_center Output parameter for y-coordinate of center of mass
 * @param z_center Output parameter for z-coordinate of center of mass
 * @return bool True if successful, false otherwise
 */
bool get_PCA_info(int somaID, const QString &imageName, double &pc1,
                  double &pc2, double &pc3, double vec1[3], double vec2[3],
                  double vec3[3], double &x_center, double &y_center,
                  double &z_center);

/**
 * @brief Map intensities from an image and its segmentation into 3D arrays
 * @param callback V3DPluginCallback2 reference
 * @param intensities Output parameter for 3D array of intensity values
 * @param segmentation Output parameter for 3D array of binary segmentation
 * values
 * @param dim_X Output parameter for X dimension
 * @param dim_Y Output parameter for Y dimension
 * @param dim_Z Output parameter for Z dimension
 * @param channel Channel to extract intensity values from (0-based index)
 * @return true if successful, false otherwise
 */
bool map_intensities(V3DPluginCallback2 &callback,
                     unsigned char ***&intensities,
                     unsigned char ***&segmentation, V3DLONG &dim_X,
                     V3DLONG &dim_Y, V3DLONG &dim_Z, int channel = 0);

/**
 * @brief Free memory allocated for two 3D arrays (intensities and segmentation)
 * @param intensities 3D array of intensity values
 * @param segmentation 3D array of segmentation values
 * @param dim_Z Number of Z slices
 * @param dim_Y Number of rows
 */
void free_mapped_arrays(unsigned char ***intensities,
                        unsigned char ***segmentation, V3DLONG dim_Z,
                        V3DLONG dim_Y);

/**
 * @brief Calculate Otsu threshold for a histogram
 * @param hist Array of 256 histogram values
 * @param totalPixels Total number of pixels in the region
 * @return The calculated Otsu threshold value
 */
int calculateOtsuThreshold(const int hist[256], int totalPixels);

/**
 * @brief Get blended Gaussian distribution parameters for a given position
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
    double &blendedMean, double &blendedStdDev);

#endif  // __MIND_PCA_ANALYSIS_H__