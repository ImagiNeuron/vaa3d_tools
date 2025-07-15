/**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */
#ifndef __MIND_SOMA_SIMULATION_H__
#define __MIND_SOMA_SIMULATION_H__

#include <v3d_interface.h>

#include <QDateTime>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <QLabel>
#include <QMessageBox>
#include <QRegularExpression>
#include <QtGui>
#include <algorithm>
#include <cmath>
#include <eigen/Dense>
#include <fstream>
#include <memory>
#include <random>
#include <vector>

#include "basic_4dimage.h"
#include "cellSegmentation_plugin.h"
#include "compute_win_pca_wp.h"
#include "soma_segmentation_plugin.h"

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

// Simulation of synthetic soma data
void simulate_soma_data(V3DPluginCallback2 &callback, QWidget *parent);

// display simulation results
void overlaySimulation(V3DPluginCallback2 &callback, QWidget *parent,
                       unsigned char *binarySegImage,
                       unsigned char *gradientImage,
                       unsigned char *simulatedImage);

// Extract and deform an existing soma shape
void extractAndDeformSomaShape(
    unsigned char *segData, unsigned char *originalData, V3DLONG xDim,
    V3DLONG yDim, V3DLONG zDim, V3DLONG sourceCenterX, V3DLONG sourceCenterY,
    V3DLONG sourceCenterZ, V3DLONG cubeSize, const double somaEigenvector1[3],
    const double somaEigenvector2[3], const double somaEigenvector3[3],
    const std::vector<double> &probabilisticModel,
    V3DLONG probabilisticModelDim_X, V3DLONG probabilisticModelDim_Y,
    V3DLONG probabilisticModelDim_Z, double radius, std::mt19937 &gen,
    double *tempSegmentation, double *tempIntensity);

#endif  // __MIND_SOMA_SIMULATION_H__