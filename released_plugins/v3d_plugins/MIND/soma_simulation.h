/**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */
#ifndef __MIND_SOMA_SIMULATION_H__
#define __MIND_SOMA_SIMULATION_H__

#include <v3d_interface.h>

#include <QDir>
#include <QFileDialog>
#include <QLabel>
#include <QMessageBox>
#include <QRegularExpression>
#include <QtCore>
#include <QtGui>
#include <QInputDialog>
#include <eigen/Dense>
#include <fstream>

#include "basic_4dimage.h"
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

// Simulation of synthetic soma data
void simulate_soma_data(V3DPluginCallback2 &callback, QWidget *parent,
                        input_PARA &PARA, bool bmenu);

// display simulation results
void overlaySimulation(V3DPluginCallback2 &callback, QWidget *parent,
                       unsigned char *binarySegImage,
                       unsigned char *gradientImage,
                       unsigned char *simulatedImage);

#endif  // __MIND_SOMA_SIMULATION_H__