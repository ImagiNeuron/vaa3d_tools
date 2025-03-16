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
 * @brief Create a background image based on segmentation threshold
 */
void create_background(V3DPluginCallback2 &callback, QWidget *parent);

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

#endif  // __MIND_PCA_ANALYSIS_H__