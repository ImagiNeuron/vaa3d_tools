#ifndef __MIND_PCA_ANALYSIS_H__
#define __MIND_PCA_ANALYSIS_H__

#include <v3d_interface.h>

#include <QApplication>
#include <QDir>
#include <QFileDialog>
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
 * @brief change names according to folder structure for storing things when
 * TeraFly data is used
 */
QString modifyFileNameForTeraFly(const QString &filename);

#endif  // __MIND_PCA_ANALYSIS_H__