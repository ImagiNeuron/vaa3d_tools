#ifndef __MIND_PCA_ANALYSIS_H__
#define __MIND_PCA_ANALYSIS_H__

#include <v3d_interface.h>

#include <QApplication>
#include <QFileDialog>
#include <QMessageBox>
#include <QtGui>
#include <fstream>

#include "basic_4dimage.h"
#include "compute_win_pca_wp.h"

// Function prototypes

void analyzeSomaPCA(unsigned char *labeledData, V3DLONG N, V3DLONG M, V3DLONG P,
                    const LocationSimple &lm, int somaIndex);

void savePCAResultsToCSV(const QString &filename, int somaIndex,
                         const LocationSimple &lm, double pc1, double pc2,
                         double pc3, const double *vec1, const double *vec2,
                         const double *vec3, double x_center, double y_center,
                         double z_center, QWidget *parent, bool *saveEnabled);

void drawLine(Image4DSimple *image, double *from, double *to);

void visualizePCA_func(V3DPluginCallback2 &callback, QWidget *parent);

#endif  // __MIND_PCA_ANALYSIS_H__