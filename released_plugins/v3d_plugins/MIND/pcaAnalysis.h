/**
 * 2025-04-18: by ImagiNeuron: Shidan Javaheri, Siger Ma, Athmane Benarous and
 * Thibaut Baguette (McGill University)
 */
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

// function to analyze a soma using PCA and return results
void analyzeSomaPCAReturnResults(unsigned char *labeledData, V3DLONG N,
                                 V3DLONG M, V3DLONG P, const LocationSimple &lm,
                                 int somaIndex, QString savePath,
                                 double vec1[3], double vec2[3], double vec3[3],
                                 double &pc1, double &pc2, double &pc3,
                                 double &x_center, double &y_center,
                                 double &z_center);

// function to analyze a soma using PCA and save results to a CSV file
void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex, QString savePath);

// function to save PCA results to a CSV file
void savePCAResultsToCSV(const QString &filename, int somaIndex,
                         const LocationSimple &lm, const double *vec1,
                         const double *vec2, const double *vec3, double pc1,
                         double pc2, double pc3, double x_center,
                         double y_center, double z_center, QWidget *parent);

// function to draw a line on the image
void drawLine(Image4DSimple *image, unsigned char r, unsigned char g,
              unsigned char b, double *from, double *to);

// function to visualize PCA results
void visualizePCA_func(V3DPluginCallback2 &callback, QWidget *parent);

// change names according to folder structure for storing things when TeraFly
// data is used
QString modifyFilePathForTeraFly(const QString &filename);

// find and load the segmentation file for the current image
void loadSegmentationFile(const QString &imageName, unsigned char *&segData,
                          V3DLONG sz[4], int &datatype,
                          V3DPluginCallback2 &callback, QWidget *parent);

#endif  // __MIND_PCA_ANALYSIS_H__