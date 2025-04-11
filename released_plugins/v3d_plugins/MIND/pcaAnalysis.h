#ifndef __MIND_PCA_ANALYSIS_H__
#define __MIND_PCA_ANALYSIS_H__

#include <v3d_interface.h>

#include <QApplication>
#include <QDir>
#include <QFileDialog>
#include <QLabel>
#include <QMessageBox>
#include <QRegularExpression>
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
                         double z_center, QWidget *parent);

/**
 * @brief draws a line to represent a PC
 */
void drawLine(Image4DSimple *image, unsigned char r, unsigned char g, unsigned char b, double *from, double *to);

/**
 * @brief Visualize PCA results
 */
void visualizePCA_func(V3DPluginCallback2 &callback, QWidget *parent);

/**
 * @brief change names according to folder structure for storing things when
 * TeraFly data is used
 */
QString modifyFilePathForTeraFly(const QString &filename);

/**
 * @brief Find and load the segmentation file for the current image
 */
void loadSegmentationFile(const QString &imageName, unsigned char *&segData,
                          V3DLONG sz[4], int &datatype,
                          V3DPluginCallback2 &callback, QWidget *parent);

/**
 * @brief Create a background image based on segmentation threshold and return
 * as 1D array
 * @param callback V3DPluginCallback2 reference
 * @param parent Parent widget
 * @param dim_X Output parameter for X dimension
 * @param dim_Y Output parameter for Y dimension
 * @param dim_Z Output parameter for Z dimension
 * @return 1D array of background image
 */
unsigned char *create_background(V3DPluginCallback2 &callback, QWidget *parent,
                                 V3DLONG &dim_X, V3DLONG &dim_Y,
                                 V3DLONG &dim_Z);

/**
 * @brief Create a background image based on segmentation threshold
 * @param callback V3DPluginCallback2 reference
 * @param parent Parent widget
 * @deprecated Use the version that returns unsigned char*** instead
 */
void create_background(V3DPluginCallback2 &callback, QWidget *parent);

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