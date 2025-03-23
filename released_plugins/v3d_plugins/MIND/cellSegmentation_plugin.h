/* adapted from cellseg_quickfind - cellSegmentation.cpp
 * 2014-10-12: by Xiang Li (lindbergh.li@gmail.com);
 * 2025-02-10: By ImagiNeuron - Thibaut Baguette, Shidan Javaheri, Siger Ma and
 * Athmane Benarous. Performs 3D flood filling on marked cells using either
 * otsu, local otsu or iterative thresholding. */

#ifndef __CELLSEGMENTATION_PLUGIN_H__
#define __CELLSEGMENTATION_PLUGIN_H__

#pragma region "includes and constants"
#include <basic_landmark.h>
#include <math.h>
#include <time.h>
#include <v3d_interface.h>

#include <QCheckBox>
#include <QComboBox>
#include <QCommonStyle>
#include <QGridLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QRadioButton>
#include <QRegularExpression>
#include <QtGui>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "cellSegmentation_plugin.h"
#include "compute_win_pca_wp.h"
#include "convert_type2uint8.h"
#include "pcaAnalysis.h"
#include "sstream"
#include "string"
#include "v3d_message.h"
using namespace std;
const int const_length_histogram = 256;
const double const_max_voxelValue = 255;
// 27 directions -1
const int const_count_neighbors = 26;
// small enough global value as a last resort
const double default_threshold_global = 15;
// cube of voxels of length 2
const int default_threshold_regionSize = 8;
const double const_infinitesimal = 0.000000001;
#define INF 1E9
#define NINF -1E9
#define PI 3.14159265
enum enum_shape_t { sphere, cube };
#pragma endregion

#pragma region "dialogInitialization"
/**
 * @class dialogRun
 * @brief Dialog box for the inputs into cell segmentation which include:
 * - Color channel
 * - Restrictions on the exemplar labels (max deviations from mass center and
 * marker position)
 * - The shape of the cells being considered
 * - The type of tresholding to be used (iterative, global Otsu, local Otsu)
 */
class dialogRun : public QDialog {
  Q_OBJECT
 public:
  QComboBox *QComboBox_mode_selection;  // combo box for mode selection
  int segmentationMode;  // 1: local otsu, 2: global Otsu, 3: iterative
  QCheckBox *QCheckBox_medianFiltering;     // checkbox for median filtering
  bool applyMedianFiltering;                // flag read from the checkbox
  QCheckBox *QCheckBox_markerConstraint;    // checkbox for marker constraint
  bool applyMarkerConstraint;               // flag read from the checkbox
  QCheckBox *QCheckBox_manualThresholding;  // checkbox for manual thresholding
  bool manualThresholding;                  //  flag for manual thresholding
  QLineEdit *QLineEdit_localOtsuRadius;     // line edit for local Otsu radius
  QLineEdit *
      QLineEdit_medianFilteringRadius;  // line edit for median filtering radius
  double medianFilteringRadius;         // radius for median filtering radius
  dialogRun(V3DPluginCallback2 &V3DPluginCallback2_currentCallback,
            QWidget *QWidget_parent, int int_channelDim) {
    // channel
    QStringList QStringList_channel_items;
    if (int_channelDim == 1) {
      QStringList_channel_items << "1";
    } else if (int_channelDim == 3) {
      QStringList_channel_items << "1 - red";
      QStringList_channel_items << "2 - green";
      QStringList_channel_items << "3 - blue";
    } else {
      for (int i = 1; i <= int_channelDim; i++) {
        QStringList_channel_items << QString().setNum(i);
      }
    }
    QComboBox_channel_selection = new QComboBox();
    QComboBox_channel_selection->addItems(QStringList_channel_items);
    if (QStringList_channel_items.size() > 1) {
      QComboBox_channel_selection->setCurrentIndex(0);
    }
    QGroupBox *QGroupBox_channel_main = new QGroupBox("Color channel");
    QGroupBox_channel_main->setStyle(new QCommonStyle());
    QGridLayout *QGridLayout_channel_main = new QGridLayout();
    QGroupBox_channel_main->setStyle(new QCommonStyle());
    QGridLayout_channel_main->addWidget(QComboBox_channel_selection, 1, 1, 1,
                                        1);
    QGroupBox_channel_main->setLayout(QGridLayout_channel_main);

    // exemplar
    QGroupBox *QGroupBox_exemplar_main = new QGroupBox("Exemplar definition");
    QGridLayout *QGridLayout_exemplar_main = new QGridLayout();
    QLabel *QLabel_exemplar_maxMovement1 =
        new QLabel(QObject::tr("Max movement from\nmass center"));
    QLineEdit_exemplar_maxMovement1 = new QLineEdit("2", QWidget_parent);
    QLabel *QLabel_exemplar_maxMovement2 =
        new QLabel(QObject::tr("Max movement from\nmarker position"));
    QLineEdit_exemplar_maxMovement2 = new QLineEdit("4", QWidget_parent);
    QGridLayout_exemplar_main->addWidget(QLabel_exemplar_maxMovement1, 1, 1, 1,
                                         1);
    QGridLayout_exemplar_main->addWidget(QLineEdit_exemplar_maxMovement1, 1, 2,
                                         1, 1);
    QGridLayout_exemplar_main->addWidget(QLabel_exemplar_maxMovement2, 1, 3, 1,
                                         1);
    QGridLayout_exemplar_main->addWidget(QLineEdit_exemplar_maxMovement2, 1, 4,
                                         1, 1);
    QGroupBox_exemplar_main->setLayout(QGridLayout_exemplar_main);

    // shape
    QGroupBox *QGroupBox_shape_main = new QGroupBox("Geometry stat");
    QGridLayout *QGridLayout_shape_main = new QGridLayout();
    QRadioButton_shape_sphere = new QRadioButton("sphere-like", QWidget_parent);
    QRadioButton_shape_sphere->setChecked(true);
    QRadioButton_shape_cube = new QRadioButton("cube-like", QWidget_parent);
    QRadioButton_shape_cube->setChecked(false);
    QGridLayout_shape_main->addWidget(QRadioButton_shape_sphere, 1, 1, 1, 1);
    QGridLayout_shape_main->addWidget(QRadioButton_shape_cube, 1, 2, 1, 1);
    QLabel *QLabel_shape_delta =
        new QLabel(QObject::tr("max anisotropic\ndeviation:"));
    QLineEdit_Shape_delta = new QLineEdit("1", QWidget_parent);
    QGridLayout_shape_main->addWidget(QLabel_shape_delta, 1, 3, 1, 1);
    QGridLayout_shape_main->addWidget(QLineEdit_Shape_delta, 1, 4, 1, 1);
    QLabel *QLabel_shape_thresholdRegionSize =
        new QLabel(QObject::tr("Min region size\nvs. exemplar ratio:"));
    QLineEdit_shape_thresholdRegionSize = new QLineEdit("0.1", QWidget_parent);
    QLabel *QLabel_shape_uThresholdRegionSize =
        new QLabel(QObject::tr("Max region size\nvs. exemplar ratio:"));
    QLineEdit_shape_uThresholdRegionSize = new QLineEdit("20", QWidget_parent);
    QGridLayout_shape_main->addWidget(QLabel_shape_thresholdRegionSize, 2, 1, 1,
                                      1);
    QGridLayout_shape_main->addWidget(QLineEdit_shape_thresholdRegionSize, 2, 2,
                                      1, 1);
    QGridLayout_shape_main->addWidget(QLabel_shape_uThresholdRegionSize, 2, 3,
                                      1, 1);
    QGridLayout_shape_main->addWidget(QLineEdit_shape_uThresholdRegionSize, 2,
                                      4, 1, 1);
    QGroupBox_shape_main->setLayout(QGridLayout_shape_main);

    // Create a new group box for segmentation preferences
    QGroupBox *QGroupBox_segmentationMode =
        new QGroupBox("Segmentation Preferences", this);
    QVBoxLayout *vLayout_segmentation = new QVBoxLayout();
    // Row 1: Segmentation mode, Local Otsu radius, Manual thresholding
    QHBoxLayout *hLayout_segmentationTop = new QHBoxLayout();
    // Segmentation Mode Label
    QLabel *label_segMode = new QLabel("Segmentation Mode:", this);
    hLayout_segmentationTop->addWidget(label_segMode, 0);
    // Segmentation Mode Dropdown
    QComboBox_mode_selection = new QComboBox(this);
    QComboBox_mode_selection->addItem("Local Otsu");
    QComboBox_mode_selection->addItem("Global Otsu");
    QComboBox_mode_selection->addItem("Iterative Threshold");
    hLayout_segmentationTop->addWidget(QComboBox_mode_selection, 2);
    // Manual Thresholding Checkbox
    QCheckBox_manualThresholding = new QCheckBox("Manual Thresholding", this);
    QCheckBox_manualThresholding->setChecked(false);  // default off
    hLayout_segmentationTop->addWidget(QCheckBox_manualThresholding, 0);
    vLayout_segmentation->addLayout(hLayout_segmentationTop);
    // Row 2: Median filtering, its radius, and marker constraint
    QHBoxLayout *hLayout_segmentationBottom = new QHBoxLayout();
    // Median Filtering Checkbox
    QCheckBox_medianFiltering = new QCheckBox("Median Filtering", this);
    QCheckBox_medianFiltering->setChecked(true);  // default is enabled
    hLayout_segmentationBottom->addWidget(QCheckBox_medianFiltering, 0);
    // Median Filtering Radius Label
    QLabel *label_medianRadius =
        new QLabel("Median Filtering Radius: (Voxels)", this);
    hLayout_segmentationBottom->addWidget(label_medianRadius, 0);
    // Median Filtering Radius Input
    QLineEdit_medianFilteringRadius =
        new QLineEdit("1", this);  // default value
    hLayout_segmentationBottom->addWidget(QLineEdit_medianFilteringRadius, 1);
    // Marker Constraint Checkbox
    QCheckBox_markerConstraint = new QCheckBox("Marker Constraint", this);
    QCheckBox_markerConstraint->setChecked(true);  // default on
    hLayout_segmentationBottom->addWidget(QCheckBox_markerConstraint, 0);
    vLayout_segmentation->addLayout(hLayout_segmentationBottom);
    // Set the layout for the group box
    QGroupBox_segmentationMode->setLayout(vLayout_segmentation);

    // control
    QPushButton *QPushButton_control_start =
        new QPushButton(QObject::tr("Run"));
    QPushButton *QPushButton_control_close =
        new QPushButton(QObject::tr("Close"));
    QWidget *QWidget_control_bar = new QWidget();
    QGridLayout *QGridLayout_control_bar = new QGridLayout();
    QGridLayout_control_bar->addWidget(QPushButton_control_start, 1, 1, 1, 1);
    QGridLayout_control_bar->addWidget(QPushButton_control_close, 1, 2, 1, 1);
    QWidget_control_bar->setLayout(QGridLayout_control_bar);

    // main panel
    QGridLayout *QGridLayout_main = new QGridLayout();
    QGridLayout_main->addWidget(QGroupBox_channel_main);
    QGridLayout_main->addWidget(QGroupBox_shape_main);
    QGridLayout_main->addWidget(QGroupBox_exemplar_main);
    QGridLayout_main->addWidget(QGroupBox_segmentationMode);
    QGridLayout_main->addWidget(QWidget_control_bar);
    setLayout(QGridLayout_main);
    setWindowTitle(QString("Soma Segmentation"));
    // event evoking - connects defined slots to signals that should run
    // when they are pressed
    connect(QPushButton_control_start, SIGNAL(clicked()), this,
            SLOT(_slot_start()));
    connect(QPushButton_control_close, SIGNAL(clicked()), this, SLOT(reject()));
    update();
  }

  // desctructor - doesn't do anything
  ~dialogRun() {}

  // member variables
  QComboBox *QComboBox_channel_selection;
  V3DLONG channel_idx_selection;
  QLineEdit *QLineEdit_Shape_delta;
  QLineEdit *QLineEdit_shape_thresholdRegionSize;
  QLineEdit *QLineEdit_shape_uThresholdRegionSize;
  QLineEdit *QLineEdit_exemplar_maxMovement1;
  QLineEdit *QLineEdit_exemplar_maxMovement2;
  QRadioButton *QRadioButton_shape_sphere;
  QRadioButton *QRadioButton_shape_cube;
  enum_shape_t shape_type_selection;
  double shape_para_delta;
  double shape_multiplier_thresholdRegionSize;
  double shape_multiplier_uThresholdRegionSize;
  V3DLONG exemplar_maxMovement1;
  V3DLONG exemplar_maxMovement2;

  // definition of slots to be connected to events
 public slots:

  // retrieve the values from the dialog box and store them in the member
  // variables
  void _slot_start() {
    channel_idx_selection = QComboBox_channel_selection->currentIndex() + 1;
    shape_para_delta = this->QLineEdit_Shape_delta->text().toDouble();
    shape_multiplier_thresholdRegionSize =
        this->QLineEdit_shape_thresholdRegionSize->text().toDouble();
    shape_multiplier_uThresholdRegionSize =
        this->QLineEdit_shape_uThresholdRegionSize->text().toDouble();
    exemplar_maxMovement1 =
        this->QLineEdit_exemplar_maxMovement1->text().toUInt();
    exemplar_maxMovement2 =
        this->QLineEdit_exemplar_maxMovement2->text().toUInt();
    // retrieve the segmentation mode:
    segmentationMode = QComboBox_mode_selection->currentIndex() + 1;

    // retrieve the median filtering flag
    applyMedianFiltering = QCheckBox_medianFiltering->isChecked();
    // Retrieve the median filtering radius only if median filtering is enabled.
    if (applyMedianFiltering) {
      medianFilteringRadius =
          QLineEdit_medianFilteringRadius->text().toDouble();
      // Enforce a valid range of 1 to 9.
      if (medianFilteringRadius < 1)
        medianFilteringRadius = 1;
      else if (medianFilteringRadius > 9)
        medianFilteringRadius = 9;
    }
    // retrieve the marker constraint flag
    applyMarkerConstraint = QCheckBox_markerConstraint->isChecked();
    // Retrieve manual threshold flag:
    manualThresholding = QCheckBox_manualThresholding->isChecked();

    if (this->QRadioButton_shape_sphere->isChecked()) {
      this->shape_type_selection = sphere;
    } else if (this->QRadioButton_shape_cube->isChecked()) {
      this->shape_type_selection = cube;
    }
    accept();
  }
};
#pragma endregion

class cellSegmentation : public QObject {
 public:
#pragma region "class: class_segmentationMain"
  class class_segmentationMain {
#pragma region "class member"
   public:
    struct double3D {
      double x;
      double y;
      double z;
      double3D(double _x = 0, double _y = 0, double _z = 0) {
        x = _x;
        y = _y;
        z = _z;
      }
    };

    struct long3D {
      V3DLONG x;
      V3DLONG y;
      V3DLONG z;
      long3D(V3DLONG _x = 0, V3DLONG _y = 0, V3DLONG _z = 0) {
        x = _x;
        y = _y;
        z = _z;
      }
    };

    // constant
    vector<V3DLONG> poss_neighborRelative;
    vector<double3D> point_neighborRelative;
    vector<vector<V3DLONG> > colors_simpleTable;

    // Input or directly derived;
    bool is_initialized;
    bool errorOccurred = false;
    unsigned char *Image1D_page;
    unsigned char *Image1D_mask;
    unsigned char ***Image3D_page;
    V3DLONG dim_X;
    V3DLONG dim_Y;
    V3DLONG dim_Z;
    V3DLONG size_page;
    V3DLONG size_page3;
    V3DLONG offset_channel;
    V3DLONG offset_Z;
    V3DLONG offset_Y;
    int idx_channel;
    int idx_shape;
    double threshold_deltaShapeStat;
    double multiplier_thresholdRegionSize;
    double multiplier_uThresholdRegionSize;
    V3DLONG max_movment1;
    V3DLONG max_movment2;
    QString name_currentWindow;

    // Exemplar (or learn from it);
    unsigned char *Image1D_exemplar;

    // segmentation - where the results are stored
    vector<vector<V3DLONG> > possVct_segmentationResult;
    unsigned char *binarySegImage;

    vector<vector<V3DLONG> > possVct_seed;
    unsigned char *Image1D_segmentationResult;

    LandmarkList LandmarkList_segmentationResult;
    vector<V3DLONG> poss_segmentationResultCenter;

    // segmentation mode
    int segmentationMode;  // 1: local otsu, 2: global Otsu, 3: iterative
    // median filtering check
    bool applyMedianFiltering;
    // median filtering radius input
    double medianFilteringRadius;
    // marker constraint check
    bool applyMarkerConstraint;
    // manual thresholding check
    bool manualThresholding;

// vector<V3DLONG> poss_segmentationResultCenterMerged;
#pragma endregion
    class_segmentationMain() {
      is_initialized = false;
      applyMedianFiltering = true;    // default to true
      applyMarkerConstraint = false;  // default: no marker constraint
    }
    ~class_segmentationMain() {}

    void printSomaSlice(double *data, int size, int padding = 1) {
      for (int y = 0; y < size; y++) {
        for (int x = 0; x < size; x++) {
          int idx = y * size + x;
          // Round the data value to the nearest integer
          int rounded = (int)round(data[idx]);

          switch (padding) {
            case 0:
              printf("%2d", rounded);
              break;
            case 1:
              printf("%4d", rounded);
              break;
            default:
              printf("%2d", rounded);
              break;
          }
        }
        printf("\n");
      }
      printf("\n");
    }

#pragma region "control-run"
    /**
     * @brief - Main function that goes over landmarks and floods them
     */
    bool control_run(unsigned char *_Image1D_original, V3DLONG _dim_X,
                     V3DLONG _dim_Y, V3DLONG _dim_Z, int _idx_channel,
                     LandmarkList &_LandmarkList_exemplar, int _idx_shape,
                     double _threshold_deltaShapeStat,
                     double _multiplier_thresholdRegionSize,
                     double _multiplier_uThresholdRegionSize,
                     QString _name_currentWindow, V3DLONG _maxMovement1,
                     V3DLONG _maxMovement2, QString fileName, int mode = 1) {
      // if (!this->is_initialized) // Temporally solution for the "parameter
      // window not popped up" problem;
      {
        // initialize parameters
        this->dim_X = _dim_X;
        this->dim_Y = _dim_Y;
        this->dim_Z = _dim_Z;
        this->idx_channel = _idx_channel;
        this->size_page = dim_X * dim_Y * dim_Z;
        this->size_page3 = this->size_page + this->size_page + this->size_page;
        this->offset_channel = (idx_channel - 1) * size_page;
        this->offset_Z = dim_X * dim_Y;
        this->offset_Y = dim_X;
        segmentationMode = mode;  // set the segmentation mode (1, 2, or 3)

        // allocate memory to store segmentation results
        this->Image1D_page = memory_allocate_uchar1D(this->size_page);
        this->Image1D_mask = memory_allocate_uchar1D(this->size_page);
        this->Image3D_page = memory_allocate_uchar3D(this->dim_Y, this->dim_X,
                                                     this->dim_Z);  // tricky!
        this->Image1D_segmentationResult =
            memory_allocate_uchar1D(this->size_page3);
        this->Image1D_exemplar = memory_allocate_uchar1D(this->size_page3);

        // iterate over each voxel in the image
        for (V3DLONG i = 0; i < this->size_page; i++) {
          // copy image data
          this->Image1D_page[i] = _Image1D_original[i + offset_channel];
          // converts an index in a 1D array to a 3D coordinate, and store
          // data in 3D array
          vector<V3DLONG> xyz_i = this->index2Coordinate(i);
          this->Image3D_page[xyz_i[2]][xyz_i[1]][xyz_i[0]] =
              this->Image1D_page[i];

          // storing the same data in 3 different positions to support colored
          // cell segmentation for visualization.
          this->Image1D_segmentationResult[i] = this->Image1D_page[i];
          this->Image1D_segmentationResult[i + size_page] =
              this->Image1D_page[i];
          this->Image1D_segmentationResult[i + size_page + size_page] =
              this->Image1D_page[i];
        }

        // this sets all voxels to the maximum value to leave them ready for
        // preprocessing
        memset(this->Image1D_mask, const_max_voxelValue, this->size_page);
        this->idx_shape = _idx_shape;
        // get potential seeds for other cells
        // not necessary for now - for now we are just doing flooding on
        // markers
        this->categorizeVoxelsByValue();

        // optionally apply the median filter as preprocessing
        if (this->applyMedianFiltering) {
          // Use the user-provided median filtering radius.
          this->filter_Median((V3DLONG)(this->medianFilteringRadius));
        }

        // set parameters that contorl the segmentation
        this->threshold_deltaShapeStat = _threshold_deltaShapeStat;
        this->multiplier_thresholdRegionSize = _multiplier_thresholdRegionSize;
        this->multiplier_uThresholdRegionSize =
            _multiplier_uThresholdRegionSize;
        this->name_currentWindow = _name_currentWindow;
        this->max_movment1 = _maxMovement1 * _maxMovement1;
        this->max_movment2 = _maxMovement2 * _maxMovement2;
        this->is_initialized = false;
      }

      // stores the tresholds that were used in segmentation
      vector<double> thresholds_valueChangeRatio;
      vector<V3DLONG> thresholds_voxelValue;
      vector<V3DLONG> thresholds_regionSize;
      vector<V3DLONG> uThresholds_regionSize;
      vector<V3DLONG> thresholds_radius;

      // constants for segmentation
      this->initializeConstants();
      vector<vector<V3DLONG> > possVct_exemplarRegion;
      vector<vector<vector<double> > > valueVctVct_exemplarShapeStat;
      vector<V3DLONG> poss_exemplar =
          landMarkList2IndexList(_LandmarkList_exemplar);
      V3DLONG count_exemplar = poss_exemplar.size();
      vector<V3DLONG> poss_exemplarNew;

      // create list of indexes of successfully segmented labels
      vector<int> segmentedLabels;

      // keep track of the largest labelled radius
      double largestRadius = 0;

      // main loop - iterating over each exemplar
      for (V3DLONG idx_exemplar = 0; idx_exemplar < count_exemplar;
           idx_exemplar++) {
        V3DLONG pos_exemplar = poss_exemplar[idx_exemplar];
        // make sure this voxel has not already been processed (255 is not
        // processed)
        if (this->Image1D_mask[pos_exemplar] < 1) {
          // soma skipped as flooding of other soma overlapped with it
          printf(
              "Skipping soma %d - due to overlap with other soma "
              "segmentations\n",
              idx_exemplar + 1);
          continue;
        }
        // variables to track region growing and see if the center of mass has
        // moved
        V3DLONG marker_intensity = this->Image1D_page[pos_exemplar];
        V3DLONG count_step = (marker_intensity - default_threshold_global);
        V3DLONG pos_massCenterOld = -1;
        V3DLONG pos_massCenterNew = 0;
        vector<V3DLONG> poss_exemplarRegionNew;
        vector<V3DLONG> poss_exemplarRegionOld;
        // the change in the threshold from the marker intensity at each
        // flooding interation
        V3DLONG idx_step = 0;
        double value_centerMovement2 = 0;

        double radius_marker = _LandmarkList_exemplar[idx_exemplar].radius;

        // update the largest radius
        if (radius_marker > largestRadius) {
          largestRadius = radius_marker;
        }

        // Retrieve the current landmark comment.
        std::string comment = _LandmarkList_exemplar[idx_exemplar].comments;

        int manualThresh = -1;
        if (this->manualThresholding) {
          size_t pos = comment.find("threshold:");
          if (pos != std::string::npos) {
            pos += 10;  // length of "threshold:"
            // Skip any whitespace.
            while (pos < comment.size() && isspace(comment[pos])) pos++;
            try {
              manualThresh = std::stoi(comment.substr(pos));
            } catch (...) {
              manualThresh = -1;  // parsing failed, ignore
            }
            // Validate the extracted value.
            if (manualThresh < default_threshold_global || manualThresh > 255)
              manualThresh = -1;
          }
        }
        // region growing on each examplar label, trying different tresholds
        // until conditions are broken. Conditions on the size of the region,
        // and the distance from the center of mass / marker. The final
        // flooded region is the one that is closest to the center of mass of
        // the original cell Depending on the segmentation mode, choose the
        // threshold:
        V3DLONG threshold_exemplarRegion;
        // If manualThresh is provided, skip any threshold calculations.
        if (manualThresh != -1) {
          threshold_exemplarRegion = manualThresh;
          // Use the manual threshold in region growing.
          poss_exemplarRegionOld = this->regionGrowOnPos(
              pos_exemplar, threshold_exemplarRegion, INF,
              this->size_page / 1000, this->Image1D_mask, radius_marker);
          pos_massCenterOld = getCenterByMass(poss_exemplarRegionOld);
          value_centerMovement2 =
              this->getEuclideanDistance2(pos_exemplar, pos_massCenterOld);
        } else {
          // No manual threshold provided: choose segmentation mode to compute
          // threshold.
          if (segmentationMode == 3) {
            for (idx_step = 0; idx_step < count_step; idx_step++) {
              threshold_exemplarRegion = marker_intensity - idx_step;
              poss_exemplarRegionNew = this->regionGrowOnPos(
                  pos_exemplar, threshold_exemplarRegion, INF,
                  this->size_page / 1000, this->Image1D_mask, radius_marker);
              this->poss2Image1D(poss_exemplarRegionNew, this->Image1D_mask,
                                 const_max_voxelValue);
              pos_massCenterNew = this->getCenterByMass(poss_exemplarRegionNew);
              double value_centerMovement1 = this->getEuclideanDistance2(
                  pos_massCenterOld, pos_massCenterNew);
              value_centerMovement2 =
                  this->getEuclideanDistance2(pos_exemplar, pos_massCenterNew);
              if (value_centerMovement1 > max_movment1) {
                printf(
                    "final threshold: %d, for marker number %d (moved too far "
                    "from center of mass)\n",
                    threshold_exemplarRegion, idx_exemplar);
                break;
              }
              if (value_centerMovement2 > max_movment2) {
                printf(
                    "final threshold: %d, for marker number %d (moved too far "
                    "from marker)\n",
                    threshold_exemplarRegion, idx_exemplar);
                break;
              }
              pos_massCenterOld = pos_massCenterNew;
              poss_exemplarRegionOld = poss_exemplarRegionNew;
            }
            // global otsu method
          } else if (segmentationMode == 2) {
            threshold_exemplarRegion = globalOtsuThreshold();
            printf("Global Otsu threshold computed: %d\n",
                   threshold_exemplarRegion);
            poss_exemplarRegionOld =
                regionGrowOnPos(pos_exemplar, threshold_exemplarRegion, INF,
                                size_page / 1000, Image1D_mask, radius_marker);
            pos_massCenterOld = getCenterByMass(poss_exemplarRegionOld);
            value_centerMovement2 =
                this->getEuclideanDistance2(pos_exemplar, pos_massCenterOld);
          } else if (segmentationMode == 1) {
            // Get the otsu threshold around the soma
            threshold_exemplarRegion =
                localOtsuThreshold(pos_exemplar, (V3DLONG)(radius_marker));
            vector<V3DLONG> xyz_exemplar = this->index2Coordinate(pos_exemplar);
            printf(
                "Local Otsu threshold computed at landmark (%ld, %ld, %ld): "
                "%d\n",
                xyz_exemplar[0], xyz_exemplar[1], xyz_exemplar[2],
                threshold_exemplarRegion);
            poss_exemplarRegionOld =
                regionGrowOnPos(pos_exemplar, threshold_exemplarRegion, INF,
                                size_page / 1000, Image1D_mask, radius_marker);
            pos_massCenterOld = getCenterByMass(poss_exemplarRegionOld);
            value_centerMovement2 =
                this->getEuclideanDistance2(pos_exemplar, pos_massCenterOld);
          }
        }
        // Update the landmark's comment with the threshold value
        _LandmarkList_exemplar[idx_exemplar].comments =
            "threshold:" + std::to_string(threshold_exemplarRegion);
        // no need to attempt a second time

        // heuristics to remove floodings that are bad, as well as bad markers

        // initial threshold didn't lead to a grown region. This check is not
        // necessary without too small region check
        if (segmentationMode == 3 && idx_step < 1 && manualThresh == -1) {
          printf("Marker number %d failed - index step did not change%d\n",
                 idx_exemplar);
          continue;
        }

        // // on final iteration, the center of mass moved too far
        if (segmentationMode == 3 &&
            value_centerMovement2 > (max_movment2 * 4)) {
          printf(
              "Marker number %d failed - value_centerMovement2 was %f (too "
              "high)\n",
              idx_exemplar, value_centerMovement2);
          continue;
        }

        // region is too small
        if (poss_exemplarRegionOld.size() < default_threshold_regionSize) {
          printf("Marker number %d failed - poss_exemplarRegionOld was %d\n",
                 idx_exemplar, poss_exemplarRegionOld.size());
          continue;
        }

        // region is too large
        if (poss_exemplarRegionOld.size() > (this->size_page / 1000)) {
          printf("Marker number %d failed - poss_exemplarRegionOld was %d\n",
                 idx_exemplar, poss_exemplarRegionOld.size());
          continue;
        }  // failed;

        // store index of successful label
        segmentedLabels.push_back(idx_exemplar);

        // analyse the properties of the shape of the cell
        vector<V3DLONG> boundBox_exemplarRegion =
            this->getBoundBox(poss_exemplarRegionOld);
        V3DLONG radius_exemplarRegion =
            getMinDimension(boundBox_exemplarRegion) / 2;
        vector<V3DLONG> xyz_exemplarRegionCenter =
            this->index2Coordinate(pos_massCenterOld);
        // the result is not saved if the shape stats are empty
        // shape stats are empty if for some reason the PCA analysis cannot be
        // completed
        vector<vector<double> > valuesVct_shapeStatExemplarRegion =
            this->getShapeStat(
                xyz_exemplarRegionCenter[0], xyz_exemplarRegionCenter[1],
                xyz_exemplarRegionCenter[2], radius_exemplarRegion);
        if (valuesVct_shapeStatExemplarRegion.empty()) {
          printf(
              "Marker number %d failed - valuesVct_shapeStatExemplarRegion "
              "was "
              "empty\n",
              idx_exemplar);
          continue;
        }  // failed;

        // mark processed voxels
        this->poss2Image1D(poss_exemplarRegionOld, this->Image1D_mask, 0);

        // store segmentation results
        possVct_exemplarRegion.push_back(poss_exemplarRegionOld);
        poss_exemplarNew.push_back(pos_massCenterOld);

        // calculate and store statistics, including threshold that was used
        V3DLONG min_exemplarRegionValue = this->getMin(poss_exemplarRegionOld);
        V3DLONG threshold_exemplarRegionValue = marker_intensity - idx_step;
        thresholds_valueChangeRatio.push_back(
            (double)(min_exemplarRegionValue - threshold_exemplarRegionValue) /
            (double)min_exemplarRegionValue);
        thresholds_voxelValue.push_back(threshold_exemplarRegionValue);
        V3DLONG count_voxel = poss_exemplarRegionOld.size();
        V3DLONG size_upper =
            count_voxel * this->multiplier_uThresholdRegionSize;
        V3DLONG size_lower = count_voxel * this->multiplier_thresholdRegionSize;
        if (size_lower < default_threshold_regionSize) {
          size_lower = default_threshold_regionSize;
        }
        thresholds_regionSize.push_back(size_lower);
        uThresholds_regionSize.push_back(size_upper);
        valueVctVct_exemplarShapeStat.push_back(
            valuesVct_shapeStatExemplarRegion);
        thresholds_radius.push_back(radius_exemplarRegion);
      }

      if (possVct_exemplarRegion.empty()) {
        return false;
      }

      poss_exemplar.clear();
      poss_exemplar = poss_exemplarNew;
      count_exemplar = poss_exemplar.size();
      memset(this->Image1D_exemplar, 0, this->size_page3);
      // store segmented regions in image format
      this->possVct2Image1DC(possVct_exemplarRegion, this->Image1D_exemplar);

      // stores the tresholds used for segmentation
      vector<V3DLONG> mapping_exemplar =
          this->sort(thresholds_voxelValue);  // in ascending order;
      V3DLONG count_seedCategory = this->possVct_seed.size();
      unsigned char **masks_page =
          this->memory_allocate_uchar2D(count_exemplar, this->size_page);

      // code that adds further markers and segmentations

      // for (V3DLONG
      // idx_exemplar=0;idx_exemplar<count_exemplar;idx_exemplar++)
      // {
      // 	memset(masks_page[idx_exemplar], const_max_voxelValue,
      // this->size_page); 	for (V3DLONG i=0;i<this->size_page;i++)
      // 	{
      // 		if (this->Image1D_mask[i]<1)
      // {masks_page[idx_exemplar][i]=0;}
      // 	}
      // }
      // for (V3DLONG
      // idx_seedCategoy=0;idx_seedCategoy<count_seedCategory;idx_seedCategoy++)
      // {
      // 	V3DLONG count_seed = this->possVct_seed[idx_seedCategoy].size();
      // 	cout<<"at value: "<<(const_max_voxelValue-idx_seedCategoy)<<",
      // totally: "<<count_seed<<" seeds;"<<endl; 	for (V3DLONG
      // idx_seed=0;idx_seed<count_seed;idx_seed++)
      // 	{
      // 		V3DLONG pos_seed =
      // this->possVct_seed[idx_seedCategoy][idx_seed]; 		for
      // (V3DLONG idx_exemplar=0;idx_exemplar<count_exemplar;idx_exemplar++)
      // 		{
      // 			if (masks_page[idx_exemplar][pos_seed]<1)
      // {continue;} 			V3DLONG value_seed =
      // this->Image1D_page[pos_seed]; 			V3DLONG
      // idx_exemplarMapped = mapping_exemplar[idx_exemplar];
      // V3DLONG threshold_backgroundValue =
      // thresholds_voxelValue[idx_exemplarMapped]; 			double
      // threshold_valueChangeRatio =
      // thresholds_valueChangeRatio[idx_exemplarMapped];
      // V3DLONG threshold_regionSize =
      // thresholds_regionSize[idx_exemplarMapped]; 			if
      // (value_seed<threshold_backgroundValue) { break; }
      // V3DLONG uThreshold_regionSize =
      // uThresholds_regionSize[idx_exemplarMapped]; 			V3DLONG
      // threshold_radius = thresholds_radius[idx_exemplarMapped];
      // vector<V3DLONG> poss_region = this->regionGrowOnPos(pos_seed,
      // threshold_backgroundValue, threshold_valueChangeRatio,
      // uThreshold_regionSize, masks_page[idx_exemplar]);
      // V3DLONG count_voxel = poss_region.size(); 			if
      // (count_voxel>uThreshold_regionSize) {continue; }
      // else if (count_voxel<threshold_regionSize)
      // 			{
      // 				for (V3DLONG
      // idx_exemplar=0;idx_exemplar<count_exemplar;idx_exemplar++)
      // 				{
      // 					this->poss2Image1D(poss_region,
      // masks_page[idx_exemplar], 0);
      // 				}
      // 				break;
      // 			}
      // 			vector<V3DLONG> boundBox_region =
      // this->getBoundBox(poss_region); 			V3DLONG
      // size_radius = this->getMinDimension(boundBox_region)/2;
      // if
      // (size_radius<(threshold_radius*this->multiplier_thresholdRegionSize))
      // 			{
      // 				for (V3DLONG
      // idx_exemplar=0;idx_exemplar<count_exemplar;idx_exemplar++)
      // 				{
      // 					this->poss2Image1D(poss_region,
      // masks_page[idx_exemplar], 0);
      // 				}
      // 				break;
      // 			}
      // 			else if
      // (size_radius>(threshold_radius*this->multiplier_uThresholdRegionSize))
      // {continue;} 			V3DLONG pos_center =
      // this->getCenterByMass(poss_region); 			vector<V3DLONG>
      // xyz_center = this->index2Coordinate(pos_center);
      // V3DLONG x = V3DLONG(xyz_center[0]); V3DLONG y =
      // V3DLONG(xyz_center[1]); V3DLONG z = V3DLONG(xyz_center[2]);
      // vector<vector<double> > valuesVct_regionShapeStat =
      // this->getShapeStat(x, y, z, size_radius);
      // //consisted of 3 vectors with length 4; 			if
      // (valuesVct_regionShapeStat.empty())
      // 			{
      // 				for (V3DLONG
      // idx_exemplar=0;idx_exemplar<count_exemplar;idx_exemplar++)
      // 				{
      // 					this->poss2Image1D(poss_region,
      // masks_page[idx_exemplar], 0);
      // 				}
      // 				break;
      // 			}
      // 			vector<double> values_PC1 =
      // valuesVct_regionShapeStat[0]; vector<double> values_PC2 =
      // valuesVct_regionShapeStat[1]; vector<double> values_PC3 =
      // valuesVct_regionShapeStat[2]; 			bool
      // is_passedShapeTest = true; 			for (int m=0; m<4; m++)
      // 			{
      // 				double value_anisotropy =
      // valueVctVct_exemplarShapeStat[idx_exemplarMapped][0][m];
      // if
      // (fabs(values_PC1[m]-value_anisotropy)>(this->threshold_deltaShapeStat*value_anisotropy))
      // 				{is_passedShapeTest = false; break;}
      // 				value_anisotropy =
      // valueVctVct_exemplarShapeStat[idx_exemplarMapped][1][m];
      // if
      // (fabs(values_PC2[m]-value_anisotropy)>(this->threshold_deltaShapeStat*value_anisotropy))
      // 				{is_passedShapeTest = false; break;}
      // 				value_anisotropy =
      // valueVctVct_exemplarShapeStat[idx_exemplarMapped][2][m];
      // if
      // (fabs(values_PC3[m]-value_anisotropy)>(this->threshold_deltaShapeStat*value_anisotropy))
      // 				{is_passedShapeTest = false; break;}
      // 			}
      // 			if (is_passedShapeTest)
      // 			{
      // 				this->possVct_segmentationResult.push_back(poss_region);
      // 				this->poss_segmentationResultCenter.push_back(pos_center);
      // 				for (V3DLONG
      // idx_exemplar=0;idx_exemplar<count_exemplar;idx_exemplar++)
      // 				{
      // 					this->poss2Image1D(poss_region,
      // masks_page[idx_exemplar], 0);
      // 				}
      // 				for (V3DLONG i=0;i<count_voxel;i++)
      // 				{
      // 					vector<V3DLONG> xyz_i =
      // this->index2Coordinate(poss_region[i]);
      // 					this->Image3D_page[xyz_i[2]][xyz_i[1]][xyz_i[0]]
      // = 0;
      // 				}
      // 				break;
      // 			}
      // 		}
      // 	}
      // }

      // leave uncommented
      // memset(this->Image1D_mask, const_max_voxelValue, this->size_page);

      // merge all segmentation results into one
      this->possVct_segmentationResult = this->mergePossVector(
          possVct_exemplarRegion, this->possVct_segmentationResult);
      // update image mask (necessary after other segmentations)
      this->possVct2Image1D(this->possVct_segmentationResult,
                            this->Image1D_mask, 0);

      // merge all segmentation centers into one
      this->poss_segmentationResultCenter =
          this->mergePoss(poss_exemplar, this->poss_segmentationResultCenter);

      // creates landmarks from the centers of the resulting segmentation
      this->LandmarkList_segmentationResult =
          this->poss2LandMarkList(this->poss_segmentationResultCenter);

      // this doesn't work
      // QString filename="quickfind_test.v3draw";
      // simple_saveimage_wrapper(this->_V3DPluginCallback2_currentCallback,filename.toAscii(),this->Image1D_segmentationResult,this->dim_X,this->dim_Y,this->dim_Z);

      // get the 3D image of the segmentation result and save it in
      // Image1D_segmentationResult
      this->possVct2Image1DC(this->possVct_segmentationResult,
                             this->Image1D_segmentationResult);

      // create the binary segmentation image
      V3DLONG size_page = this->dim_X * this->dim_Y * this->dim_Z;
      binarySegImage = new unsigned char[size_page];
      memset(binarySegImage, 0,
             size_page);  // Initialize to background (black)
      for (const auto &region : this->possVct_segmentationResult) {
        for (V3DLONG idx : region) {
          binarySegImage[idx] = 255;  // Set somas to white
        }
      }

      QString savePath = fileName + "_pca_binary_segmentation.csv";

      // make an array to store the counts of each voxel being part of a
      // soma size of the array is based on the largest radius bounding the
      // somas
      V3DLONG cubeSize = ((V3DLONG)ceil(largestRadius) + 3) * 2;
      V3DLONG centralSlice = (cubeSize / 2) - 1;
      V3DLONG totalVoxels = cubeSize * cubeSize * cubeSize;

      // store the binary segmentation of each soma
      double *somaSegmentation = new double[totalVoxels];
      memset(somaSegmentation, 0, totalVoxels * sizeof(double));

      // store the counts of each voxel being part of a soma
      double *probabilityModel = new double[totalVoxels];
      memset(probabilityModel, 0, totalVoxels * sizeof(double));

      // for each index inside the segmentedLabels vector
      int segmentationCount = 0;

      for (int idx_exemplar : segmentedLabels) {
        // perform PCA analysis on the binary segmentation

        // results of PCA for alignment of soma
        double pc1, pc2, pc3;
        double vec1[3], vec2[3], vec3[3];
        double x_center, y_center, z_center;

        analyzeSomaPCAReturnResults(
            this->binarySegImage, this->dim_X, this->dim_Y, this->dim_Z,
            _LandmarkList_exemplar[idx_exemplar], idx_exemplar + 1, savePath,
            pc1, pc2, pc3, vec1, vec2, vec3, x_center, y_center, z_center);

        // for the label at idx_exemplar, get the segmentation (square around
        // marker center)
        vector<V3DLONG> binarySomaIndicies =
            possVct_exemplarRegion[segmentationCount];
        segmentationCount++;

        // Adjust the segmentation indices so that the center-of-mass aligns
        // with the volume center.
        adjustSegmentationCenter(binarySomaIndicies, cubeSize, x_center,
                                 y_center, z_center, somaSegmentation);

        // Rotate the segmentation so that its principal axes align with the x,
        // y, and z axes.
        rotateSegmentation(somaSegmentation, cubeSize, pc1, pc2, pc3, vec1,
                           vec2, vec3);

        // Accumulate the binary segmentation into the probability model.
        for (V3DLONG i = 0; i < totalVoxels; i++) {
          probabilityModel[i] += somaSegmentation[i];
        }

        // Compute the central slice index
        if (centralSlice < 0 || centralSlice >= cubeSize) {
          printf("Central slice out of bounds\n");
        } else {
          // Print the central slice of the soma segmentation if there was an
          // error
          if (somaSegmentation[centralSlice * cubeSize * cubeSize +
                               centralSlice * cubeSize + centralSlice] == 0) {
            errorOccurred = true;
            printf(
                "Error: Soma %d is not centered after rotation. Check "
                "segmentation\n",
                idx_exemplar + 1);
            printf("Soma segmentation (central slice) after rotation: \n");
            printSomaSlice(
                somaSegmentation + (centralSlice * cubeSize * cubeSize),
                cubeSize, 0);
            printf("Probability model (central slice): \n");
            printSomaSlice(
                probabilityModel + (centralSlice * cubeSize * cubeSize),
                cubeSize);
          }
        }

        // Clear somaSegmentation for the next exemplar.
        memset(somaSegmentation, 0, totalVoxels * sizeof(double));
      }

      // print value at the center of the probability model to see if it is
      // working

      printf("Final probability model (central slice) \n");
      printSomaSlice(probabilityModel + (centralSlice * cubeSize * cubeSize),
                     cubeSize);

      QString saveModelPath = fileName + "_probability_model.bin";

      if (!saveProbabilityModel(saveModelPath.toStdString(), probabilityModel,
                                totalVoxels)) {
        printf("Failed to save probability model\n");
      }

      // Free the allocated memory for the probability model and segmentation.
      delete[] somaSegmentation;
      delete[] probabilityModel;
      this->memory_free_uchar2D(masks_page, count_exemplar);
      return true;
    }

    /**
     * @brief Helper function to adjust coordinates of a segmented soma to align
     * with the center of mass.
     */
    void adjustSegmentationCenter(const vector<V3DLONG> &indices,
                                  V3DLONG cubeSize, double x_center,
                                  double y_center, double z_center,
                                  double *segmentation) {
      // Calculate the target center index of the cube.
      int center = cubeSize / 2;
      // Compute integer shifts (rounding the center-of-mass coordinates).
      int x_shift = center - static_cast<int>(round(x_center));
      int y_shift = center - static_cast<int>(round(y_center));
      int z_shift = center - static_cast<int>(round(z_center));

      // Loop over each voxel index in the provided segmentation region.
      for (size_t i = 0; i < indices.size(); i++) {
        V3DLONG idx = indices[i];
        // Convert the flat index to (x, y, z)
        int x = idx % dim_X;
        int y = (idx / dim_X) % dim_Y;
        int z = idx / (dim_X * dim_Y);

        // Adjust coordinates based on the computed
        int newX = x + x_shift;
        int newY = y + y_shift;
        int newZ = z + z_shift;

        // Check that the new coordinates lie within bounds.
        if (newX >= 0 && newX < cubeSize && newY >= 0 && newY < cubeSize &&
            newZ >= 0 && newZ < cubeSize) {
          V3DLONG newIdx = newZ * cubeSize * cubeSize + newY * cubeSize + newX;
          segmentation[newIdx] = 1.0;
          // optional debugging
          // printf(
          //     "Index %zu: original (%d, %d, %d) adjusted to (%d, %d, %d) -> "
          //     "newIdx = %ld\n",
          //     i, x, y, z, newX, newY, newZ, newIdx);
        }
      }
    }

    /**
     * @brief Helper function to rotate a segmented soma using PCA results.
     */
    void rotateSegmentationToAxes(double *segmentation, V3DLONG cubeSize,
                                  double pc1, double pc2, double pc3,
                                  double ev1[3], double ev2[3], double ev3[3],
                                  double ax1[3], double ax2[3], double ax3[3]) {
      V3DLONG totalVoxels = cubeSize * cubeSize * cubeSize;
      // Use a std::vector for temporary storage instead of raw new[]:
      std::vector<double> rotated(totalVoxels, 0);
      int center = cubeSize / 2;

      // construct the rotation matrix
      // R = B A^T
      // where B is the rotation matrix for canonical-to-new
      // and A is the rotation matrix for canonical-to-eigenvector
      //     [ev1[0] ev2[0] ev3[0]]
      // B = [ev1[1] ev2[1] ev3[1]]
      //     [ev1[2] ev2[2] ev3[2]]
      //     [ax1[0] ax2[0] ax3[0]]       [ax1[0] ax1[1] ax1[2]]
      // A = [ax1[1] ax2[1] ax3[1]] A^T = [ax2[0] ax2[1] ax2[2]]
      //     [ax1[2] ax2[2] ax3[2]]       [ax3[0] ax3[1] ax3[2]]
      double R[3][3] = {
          {ev1[0] * ax1[0] + ev2[0] * ax2[0] + ev3[0] * ax3[0],
           ev1[0] * ax1[1] + ev2[0] * ax2[1] + ev3[0] * ax3[1],
           ev1[0] * ax1[2] + ev2[0] * ax2[2] + ev3[0] * ax3[2]},
          {ev1[1] * ax1[0] + ev2[1] * ax2[0] + ev3[1] * ax3[0],
           ev1[1] * ax1[1] + ev2[1] * ax2[1] + ev3[1] * ax3[1],
           ev1[1] * ax1[2] + ev2[1] * ax2[2] + ev3[1] * ax3[2]},
          {ev1[2] * ax1[0] + ev2[2] * ax2[0] + ev3[2] * ax3[0],
           ev1[2] * ax1[1] + ev2[2] * ax2[1] + ev3[2] * ax3[1],
           ev1[2] * ax1[2] + ev2[2] * ax2[2] + ev3[2] * ax3[2]},
      };

      // Define supersampling resolution per axis.
      const int samplesPerAxis = 2;  // 2x2x2 grid => 8 samples per voxel.
      const int numSamples = samplesPerAxis * samplesPerAxis * samplesPerAxis;
      // Precompute the sub-voxel offsets (center of each sub-cube)
      std::vector<double> offsets(samplesPerAxis);
      for (int i = 0; i < samplesPerAxis; i++) {
        offsets[i] = (i + 0.5) / samplesPerAxis;  // e.g. for 2: 0.25, 0.75
      }

      // Inverse mapping: iterate over every voxel in the output (rotated)
      // volume.
      for (int z = 0; z < cubeSize; z++) {
        for (int y = 0; y < cubeSize; y++) {
          for (int x = 0; x < cubeSize; x++) {
            V3DLONG outIdx = z * cubeSize * cubeSize + y * cubeSize + x;
            double sum = 0;
            // Loop over sub-voxel samples.
            for (int dz = 0; dz < samplesPerAxis; dz++) {
              for (int dy = 0; dy < samplesPerAxis; dy++) {
                for (int dx = 0; dx < samplesPerAxis; dx++) {
                  // Compute sub-voxel coordinate in output volume.
                  // Adding the sub-voxel offset to the integer coordinate.
                  double sampleX = x + offsets[dx];
                  double sampleY = y + offsets[dy];
                  double sampleZ = z + offsets[dz];

                  // the rotated (output) basis. This is the vector r =
                  // [rx,ry,rz]. We will consider the original vector in the
                  // unrotated basis as o = [ox, oy, oz]
                  double rx = sampleX - center;
                  double ry = sampleY - center;
                  double rz = sampleZ - center;

                  // The rotation matrix A has the eigenvectors as its columns
                  // We know that o = A * r and r = A^T * o
                  // Horizontal axis (new X axis) (first col of A): Second
                  // longest PC Vertical axis (new Y axis) (second col of A):
                  // Longest PC Depth axis (new Z axis) (last col of A):
                  // Shortest PC allows a veritcal slice to show the most
                  // information

                  // multiply by the inverse of the rotation matrix
                  double ox = R[0][0] * rx + R[0][1] * ry + R[0][2] * rz;
                  double oy = R[1][0] * rx + R[1][1] * ry + R[1][2] * rz;
                  double oz = R[2][0] * rx + R[2][1] * ry + R[2][2] * rz;

                  // Convert back to original volume coordinates.
                  int src_x = static_cast<int>(round(ox)) + center;
                  int src_y = static_cast<int>(round(oy)) + center;
                  int src_z = static_cast<int>(round(oz)) + center;

                  // If the computed source coordinates are valid, sample the
                  // input segmentation.
                  if (src_x >= 0 && src_x < cubeSize && src_y >= 0 &&
                      src_y < cubeSize && src_z >= 0 && src_z < cubeSize) {
                    V3DLONG srcIdx =
                        src_z * cubeSize * cubeSize + src_y * cubeSize + src_x;
                    sum += segmentation[srcIdx];
                  }
                  // If out-of-bounds, we treat the sample as 0.
                }
              }
            }
            // Set the output voxel to 1 if the majority of sub-samples are 1.
            rotated[outIdx] = sum / numSamples;
          }
        }
      }

      // Copy the rotated volume back to the original segmentation array.
      memcpy(segmentation, rotated.data(), totalVoxels * sizeof(double));
    }

    /**
     * @brief Helper function to rotate a segmented soma using PCA results.
     */
    void rotateSegmentation(double *segmentation, V3DLONG cubeSize, double pc1,
                            double pc2, double pc3, double vec1[3],
                            double vec2[3], double vec3[3]) {
      // align first (longest) principal component with the y-axis
      double ax1[] = {0.0, 1.0, 0.0};
      double ax2[] = {1.0, 0.0, 0.0};
      double ax3[] = {0.0, 0.0, 1.0};
      rotateSegmentationToAxes(segmentation, cubeSize, pc1, pc2, pc3, vec1,
                               vec2, vec3, ax1, ax2, ax3);
    }

    /**
     * @brief Helper function to save the probability model to a binary file.
     */
    static bool saveProbabilityModel(const std::string &filename, const double *probabilityModel, V3DLONG dimX, V3DLONG dimY, V3DLONG dimZ) {
      std::ofstream outFile(filename, std::ios::binary);

      if (!outFile) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
      }

      outFile.write(reinterpret_cast<const char *>(&dimX), sizeof(V3DLONG));
      outFile.write(reinterpret_cast<const char *>(&dimY), sizeof(V3DLONG));
      outFile.write(reinterpret_cast<const char *>(&dimZ), sizeof(V3DLONG));

      // Write the entire array as binary.
      outFile.write(reinterpret_cast<const char *>(probabilityModel), dimX * dimY * dimZ * sizeof(double));
      if (!outFile.good()) {
        std::cerr << "Error: Failed to write data to file " << filename << "." << std::endl;
        return false;
      }

      outFile.close();
      return true;
    }

    /**
     * @brief Helper function to load a probability model from a binary file.
     */
    static bool loadProbabilityModel(const std::string &filename, std::vector<double>& probabilityModelOut, V3DLONG& dimXOut, V3DLONG& dimYOut, V3DLONG& dimZOut) {
      std::ifstream inFile(filename, std::ios::binary);

      if (!inFile) {
        std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
        return false;
      }

      inFile.read(reinterpret_cast<char *>(&dimXOut), sizeof(V3DLONG));
      inFile.read(reinterpret_cast<char *>(&dimYOut), sizeof(V3DLONG));
      inFile.read(reinterpret_cast<char *>(&dimZOut), sizeof(V3DLONG));

      // Read the entire array from the binary file.
      probabilityModelOut.resize(dimXOut * dimYOut * dimZOut);
      inFile.read(reinterpret_cast<char *>(probabilityModelOut.data()), dimXOut * dimYOut * dimZOut * sizeof(double));

      if (!inFile.good() && !inFile.eof()) {
        std::cerr << "Error: Failed to read data from file " << filename << "." << std::endl;
        return false;
      }

      inFile.close();
      return true;
    }

#pragma endregion

    /**
     * @brief Function to indentify voxels that can be potential seeds for
     * cells based on their intensity
     */
#pragma region "regionGrow"
    void categorizeVoxelsByValue() {
      this->possVct_seed.clear();
      vector<V3DLONG> poss_empty(0, 0);
      for (V3DLONG i = default_threshold_global; i < const_length_histogram;
           i++) {
        this->possVct_seed.push_back(poss_empty);
      }
      for (V3DLONG i = 0; i < this->size_page; i++) {
        V3DLONG value_i = this->Image1D_page[i];
        if (value_i > default_threshold_global) {
          V3DLONG offset_i = const_max_voxelValue - value_i;
          this->possVct_seed[offset_i].push_back(i);
        }
      }
    }

    /**
     * @brief given a seed voxel, grow the region until its boundaries
     */
    vector<V3DLONG> regionGrowOnPos(V3DLONG _pos_seed,
                                    V3DLONG _threshold_voxelValue,
                                    double _threshold_valueChangeRatio,
                                    V3DLONG _uThreshold_regionSize,
                                    unsigned char *_mask_input,
                                    double _radius) {
      // final segmentation
      vector<V3DLONG> poss_result;
      // voxels being considered
      vector<V3DLONG> poss_growing;
      poss_growing.push_back(_pos_seed);
      poss_result.push_back(_pos_seed);
      V3DLONG count_voxel = 1;

      // voxel has been processed for seed;
      _mask_input[_pos_seed] = 0;
      // store the intensity of the seed voxel
      V3DLONG min_voxelValue = this->Image1D_page[_pos_seed];
      while (true) {
        // growing complete - no more candidate voxels
        if (poss_growing.empty()) {
          return poss_result;
        }
        // take next candidate voxel and get its 3D coordinates
        V3DLONG pos_current = poss_growing.back();
        poss_growing.pop_back();
        vector<V3DLONG> xyz_current = this->index2Coordinate(pos_current);

        // for all possible neighbours inside the image bounds
        for (int j = 0; j < const_count_neighbors; j++) {
          if (((xyz_current[0] + point_neighborRelative[j].x) < 0) ||
              ((xyz_current[0] + point_neighborRelative[j].x) >= this->dim_X) ||
              ((xyz_current[1] + point_neighborRelative[j].y) < 0) ||
              ((xyz_current[1] + point_neighborRelative[j].y) >= this->dim_Y) ||
              ((xyz_current[2] + point_neighborRelative[j].z) < 0) ||
              ((xyz_current[2] + point_neighborRelative[j].z) >= this->dim_Z)) {
            // do nothing because they are invalid
          } else {
            V3DLONG pos_neighbor = pos_current + poss_neighborRelative[j];

            // prevent it from going out of bounds;
            if (this->checkValidity(pos_neighbor)) {
              // only prcoess voxels that haven't yet been processed
              if (_mask_input[pos_neighbor] > 0) {
                // if marker constraint is enabled, skip neighbors that fall
                // outside the sphere
                if (this->applyMarkerConstraint) {
                  double distSq =
                      this->getEuclideanDistance2(_pos_seed, pos_neighbor);
                  if (distSq > _radius * _radius) {
                    continue;
                  }
                }

                V3DLONG value_neighbor = this->Image1D_page[pos_neighbor];
                if ((value_neighbor > _threshold_voxelValue) &&
                    ((min_voxelValue - value_neighbor) <
                     (min_voxelValue * _threshold_valueChangeRatio))) {
                  // value has now been processed
                  _mask_input[pos_neighbor] = 0;
                  poss_growing.push_back(pos_neighbor);
                  poss_result.push_back(pos_neighbor);
                  count_voxel++;

                  // store minimum voxel value
                  if (value_neighbor < min_voxelValue) {
                    min_voxelValue = value_neighbor;
                  }
                  if (count_voxel > (_uThreshold_regionSize + 2))  // too large;
                  {
                    return poss_result;
                  }
                }
              }
            }
          }
        }
      }
    }
#pragma endregion

#pragma region "utility functions"
    static void Image3D2Image1D(int ***Image3D_input,
                                unsigned char *Image1D_output, const int size_X,
                                const int size_Y, const int size_Z) {
      int tmp_value = 0;
      int tmp_idx = 0;
      for (int z = 0; z < size_Z; z++) {
        for (int x = 0; x < size_X; x++) {
          for (int y = 0; y < size_Y; y++) {
            tmp_value = Image3D_input[z][x][y];
            tmp_idx = class_segmentationMain::coordinate2Index(x, y, z, size_X,
                                                               size_X * size_Y);
            Image1D_output[tmp_idx] = tmp_value;
          }
        }
      }
      return;
    }

    void Image1D2Image3D(const unsigned char *Image1D_input,
                         double ***Image3D_output, const int dim_X,
                         const int dim_Y, const int dim_Z) {
      vector<V3DLONG> vct_coordinate;
      int count_page = dim_X * dim_Y * dim_Z;
      for (int i = 0; i < count_page; i++) {
        vct_coordinate = index2Coordinate(i);
        Image3D_output[vct_coordinate[2]][vct_coordinate[0]]
                      [vct_coordinate[1]] = Image1D_input[i];
      }
      return;
    }

    bool checkValidity(V3DLONG idx_input) {
      if ((idx_input >= 0) && (idx_input < this->size_page)) {
        return true;
      } else {
        return false;
      }
    }

    /**
     * @brief - merges all centers of mass into one
     */
    vector<vector<V3DLONG> > mergePossVector(
        vector<vector<V3DLONG> > vctList_input1,
        vector<vector<V3DLONG> > vctList_input2)  // vctList_input2 will be
                                                  // appended to vctList_input1;
    {
      vector<vector<V3DLONG> > vctList_result = vctList_input1;
      vector<V3DLONG> vct_tmp;
      V3DLONG count_region2 = vctList_input2.size();
      for (int i = 0; i < count_region2; i++) {
        vct_tmp = vctList_input2[i];
        vctList_result.push_back(vct_tmp);
      }
      return vctList_result;
    }

    vector<V3DLONG> mergePoss(vector<V3DLONG> poss_input1,
                              vector<V3DLONG> poss_input2) {
      vector<V3DLONG> poss_result = poss_input1;
      V3DLONG count_pos = poss_input2.size();
      for (int i = 0; i < count_pos; i++) {
        poss_result.push_back(poss_input2[i]);
      }
      return poss_result;
    }

    LandmarkList poss2LandMarkList(vector<V3DLONG> vct_index) {
      LandmarkList LandmarkList_result;
      LocationSimple Landmark_tmp;
      for (int i = 0; i < vct_index.size(); i++) {
        Landmark_tmp = index2LandMark(vct_index[i]);
        LandmarkList_result.push_back(Landmark_tmp);
      }
      return LandmarkList_result;
    }

    LocationSimple index2LandMark(V3DLONG idx_Input) {
      vector<V3DLONG> vct_coordinate = index2Coordinate(idx_Input);
      V3DLONG x = vct_coordinate[0] + 1;
      V3DLONG y = vct_coordinate[1] + 1;
      V3DLONG z = vct_coordinate[2] + 1;
      LocationSimple LocationSimple_result(x, y, z);
      return LocationSimple_result;
    }

    vector<V3DLONG> landMarkList2IndexList(LandmarkList LandmarkList_input) {
      vector<V3DLONG> vct_result;
      for (V3DLONG idx_input = 0; idx_input < LandmarkList_input.count();
           idx_input++) {
        vct_result.push_back(landMark2Index(LandmarkList_input.at(idx_input)));
      }
      return vct_result;
    }

    V3DLONG landMark2Index(LocationSimple Landmark_input) {
      float x = 0;
      float y = 0;
      float z = 0;
      Landmark_input.getCoord(x, y, z);
      return (coordinate2Index(x - 1, y - 1, z - 1));
    }

    /**
     * Convert a 1D index to a 3D coordinate
     */
    vector<V3DLONG> index2Coordinate(V3DLONG idx) {
      vector<V3DLONG> vct_result(3, -1);
      vct_result[2] = floor((double)idx / (double)offset_Z);
      vct_result[1] =
          floor((double)(idx - vct_result[2] * offset_Z) / (double)offset_Y);
      vct_result[0] = idx - vct_result[2] * offset_Z - vct_result[1] * offset_Y;
      return vct_result;
    }

    static vector<V3DLONG> index2Coordinate(V3DLONG idx, V3DLONG offset_Y,
                                            V3DLONG offset_Z) {
      vector<V3DLONG> vct_result(3, -1);
      vct_result[2] = floor((double)idx / (double)offset_Z);
      vct_result[1] =
          floor((double)(idx - vct_result[2] * offset_Z) / (double)offset_Y);
      vct_result[0] = idx - vct_result[2] * offset_Z - vct_result[1] * offset_Y;
      return vct_result;
    }

    V3DLONG coordinate2Index(V3DLONG x, V3DLONG y, V3DLONG z) {
      return z * this->offset_Z + y * this->offset_Y + x;
    }

    static V3DLONG coordinate2Index(V3DLONG x, V3DLONG y, V3DLONG z,
                                    V3DLONG offset_Y, V3DLONG offset_Z) {
      return z * offset_Z + y * offset_Y + x;
    }

    /**
     * @brief initialize constants necesssary to perform region growing
     */
    void initializeConstants() {
      this->poss_neighborRelative.clear();
      this->point_neighborRelative.clear();
      this->colors_simpleTable.clear();

      // identify relative neighbours of a voxel
      double3D point_neighbor;
      for (V3DLONG z = -1; z <= 1; z++) {
        for (V3DLONG y = -1; y <= 1; y++) {
          for (V3DLONG x = -1; x <= 1; x++) {
            if (x == 0 && y == 0 && z == 0) {
              // that's itself;
            } else {
              this->poss_neighborRelative.push_back(z * this->offset_Z +
                                                    y * this->offset_Y + x);
              point_neighbor.x = x;
              point_neighbor.y = y;
              point_neighbor.z = z;
              this->point_neighborRelative.push_back(point_neighbor);
            }
          }
        }
      }

      // different predfined colors for visualization
      vector<V3DLONG> color_tmp(3, 0);
      color_tmp[0] = 255;
      color_tmp[1] = 0;
      color_tmp[2] = 0;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 0;
      color_tmp[1] = 255;
      color_tmp[2] = 0;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 0;
      color_tmp[1] = 0;
      color_tmp[2] = 255;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 255;
      color_tmp[1] = 255;
      color_tmp[2] = 0;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 0;
      color_tmp[1] = 255;
      color_tmp[2] = 255;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 255;
      color_tmp[1] = 0;
      color_tmp[2] = 255;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 255;
      color_tmp[1] = 128;
      color_tmp[2] = 0;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 128;
      color_tmp[1] = 255;
      color_tmp[2] = 0;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 0;
      color_tmp[1] = 128;
      color_tmp[2] = 255;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 255;
      color_tmp[1] = 255;
      color_tmp[2] = 128;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 128;
      color_tmp[1] = 255;
      color_tmp[2] = 255;
      this->colors_simpleTable.push_back(color_tmp);
      color_tmp[0] = 255;
      color_tmp[1] = 128;
      color_tmp[2] = 255;
      this->colors_simpleTable.push_back(color_tmp);
    }

    /**
     * @brief - assigns colors to the segmented cell regions
     */
    void possVct2Image1DC(vector<vector<V3DLONG> > possVct_input,
                          unsigned char *Image1D_input) {
      vector<V3DLONG> color_input(3, 0);
      for (int i = 0; i < possVct_input.size(); i++) {
        int idx_color = i % 12;
        color_input[0] = colors_simpleTable[idx_color][0];
        color_input[1] = colors_simpleTable[idx_color][1];
        color_input[2] = colors_simpleTable[idx_color][2];
        poss2Image1DC(Image1D_input, possVct_input[i], color_input);
      }
    }

    void poss2Image1D(vector<V3DLONG> poss_input, unsigned char *Image1D_input,
                      V3DLONG value_input) {
      V3DLONG size_input = poss_input.size();
      for (int i = 0; i < size_input; i++) {
        Image1D_input[poss_input[i]] = value_input;
      }
    }

    void poss2Image1D(vector<V3DLONG> poss_input, V3DLONG *Image1D_input,
                      V3DLONG value_input) {
      V3DLONG size_input = poss_input.size();
      for (int i = 0; i < size_input; i++) {
        Image1D_input[poss_input[i]] = value_input;
      }
    }

    void possVct2Image1D(vector<vector<V3DLONG> > possVct_input,
                         unsigned char *Image1D_input, V3DLONG value_input) {
      V3DLONG count_region = possVct_input.size();
      for (V3DLONG i = 0; i < count_region; i++) {
        poss2Image1D(possVct_input[i], Image1D_input, value_input);
      }
    }

    void poss2Image1DC(unsigned char *Image1D_input, vector<V3DLONG> poss_input,
                       vector<V3DLONG> color_input) {
      for (int i = 0; i < poss_input.size(); i++) {
        if (this->checkValidity(poss_input[i])) {
          Image1D_input[poss_input[i]] = color_input[0];
          Image1D_input[poss_input[i] + this->size_page] = color_input[1];
          Image1D_input[poss_input[i] + this->size_page + this->size_page] =
              color_input[2];
        }
      }
    }

    static void neuronTree2LandmarkList(const NeuronTree &NeuronTree_input,
                                        LandmarkList &LandmarkList_output) {
      LocationSimple LocationSimple_temp(0, 0, 0);
      for (V3DLONG i = 0; i < NeuronTree_input.listNeuron.size(); i++) {
        LocationSimple_temp.x = NeuronTree_input.listNeuron.at(i).x;
        LocationSimple_temp.y = NeuronTree_input.listNeuron.at(i).y;
        LocationSimple_temp.z = NeuronTree_input.listNeuron.at(i).z;
        LandmarkList_output.append(LocationSimple_temp);
      }
    }

    static V3DLONG vctContains(vector<V3DLONG> vct_input, V3DLONG idx_input) {
      for (int i = 0; i < vct_input.size(); i++) {
        if (vct_input[i] == idx_input) {
          return i;
        }
      }
      return -1;
    }
#pragma endregion

#pragma region "sorting and comparison"
    double getMax(vector<double> values_input) {
      double max_result = -INF;
      for (std::vector<double>::iterator it = values_input.begin();
           it != values_input.end(); ++it) {
        if (max_result < *it) {
          max_result = *it;
        }
      }
      return max_result;
    }

    double getMax(vector<V3DLONG> poss_input) {
      double max_result = -INF;
      V3DLONG count_input = poss_input.size();
      for (V3DLONG i = 0; i < count_input; i++) {
        double value_i = this->Image1D_page[poss_input[i]];
        if (max_result < value_i) {
          max_result = value_i;
        }
      }
      return max_result;
    }

    double getMin(vector<double> values_input) {
      double min_result = -INF;
      for (std::vector<double>::iterator it = values_input.begin();
           it != values_input.end(); ++it) {
        if (min_result < *it) {
          min_result = *it;
        }
      }
      return min_result;
    }

    double getMin(vector<V3DLONG> poss_input) {
      double min_result = INF;
      V3DLONG count_input = poss_input.size();
      for (V3DLONG i = 0; i < count_input; i++) {
        double value_i = this->Image1D_page[poss_input[i]];
        if (min_result > value_i) {
          min_result = value_i;
        }
      }
      return min_result;
    }

    vector<V3DLONG> sort(vector<V3DLONG> values_input) {
      V3DLONG count_input = values_input.size();
      vector<V3DLONG> mapping_result;
      for (V3DLONG i = 0; i < count_input; i++) {
        V3DLONG value_i = values_input[i];
        V3DLONG count_greater = 0;
        for (V3DLONG j = 0; j < count_input; j++) {
          V3DLONG value_j = values_input[j];
          if (value_i < value_j) {
            count_greater++;
          }
          if ((value_i == value_j) && (j > i)) {
            count_greater++;
          }
        }
        mapping_result.push_back(count_greater);
      }
      return mapping_result;
    }

    void swap(V3DLONG &x, V3DLONG &y) {
      V3DLONG tmp = x;
      x = y;
      y = tmp;
    }
#pragma endregion

#pragma region "geometry property"
    static vector<V3DLONG> getOffset(const V3DLONG dim_X, const V3DLONG dim_Y,
                                     const V3DLONG dim_Z) {
      vector<V3DLONG> size_result(2, 0);
      size_result[0] = dim_X * dim_Y;
      size_result[1] = dim_X;
      return size_result;
    }

    double getEuclideanDistance2(V3DLONG pos_input1, V3DLONG pos_input2) {
      if ((pos_input1 < 0) || (pos_input2 < 0)) {
        return 0;
      }
      double result = 0;
      vector<V3DLONG> vct_xyz1 = this->index2Coordinate(pos_input1);
      vector<V3DLONG> vct_xyz2 = this->index2Coordinate(pos_input2);
      result += (vct_xyz1[0] - vct_xyz2[0]) * (vct_xyz1[0] - vct_xyz2[0]);
      result += (vct_xyz1[1] - vct_xyz2[1]) * (vct_xyz1[1] - vct_xyz2[1]);
      result += (vct_xyz1[2] - vct_xyz2[2]) * (vct_xyz1[2] - vct_xyz2[2]);
      return result;
    }

    double getEuclideanDistance(V3DLONG pos_input1, V3DLONG pos_input2) {
      if ((pos_input1 < 0) || (pos_input2 < 0)) {
        return 0;
      }
      double result = 0;
      vector<V3DLONG> vct_xyz1 = this->index2Coordinate(pos_input1);
      vector<V3DLONG> vct_xyz2 = this->index2Coordinate(pos_input2);
      result += (vct_xyz1[0] - vct_xyz2[0]) * (vct_xyz1[0] - vct_xyz2[0]);
      result += (vct_xyz1[1] - vct_xyz2[1]) * (vct_xyz1[1] - vct_xyz2[1]);
      result += (vct_xyz1[2] - vct_xyz2[2]) * (vct_xyz1[2] - vct_xyz2[2]);
      result = sqrt(result);
      return result;
    }

    static double getEuclideanDistance(double3D point_input1,
                                       double3D point_input2) {
      double result = 0;
      result +=
          (point_input1.x - point_input2.x) * (point_input1.x - point_input2.x);
      result +=
          (point_input1.y - point_input2.y) * (point_input1.y - point_input2.y);
      result +=
          (point_input1.z - point_input2.z) * (point_input1.z - point_input2.z);
      result = sqrt(result);
      return result;
    }

    void centralizeRegion(const vector<V3DLONG> poss_input,
                          const V3DLONG size_X, const V3DLONG size_Y,
                          const V3DLONG size_Z, const V3DLONG min_X,
                          const V3DLONG min_Y, const V3DLONG min_Z,
                          unsigned char *Image1D_output) {
      V3DLONG pos_voxel = 0;
      V3DLONG pos_centralized = 0;
      V3DLONG x = 0;
      V3DLONG y = 0;
      V3DLONG z = 0;
      vector<V3DLONG> xyz_voxel(0, 0);
      V3DLONG count_voxel = poss_input.size();
      V3DLONG size_region = size_X * size_Y * size_Z;
      for (int i = 0; i < size_region; i++) {
        Image1D_output[i] = 0;
      }
      for (V3DLONG idx_voxel = 0; idx_voxel < count_voxel; idx_voxel++) {
        pos_voxel = poss_input[idx_voxel];
        xyz_voxel = this->index2Coordinate(pos_voxel);
        x = xyz_voxel[0] - min_X;
        y = xyz_voxel[1] - min_Y;
        z = xyz_voxel[2] - min_Z;
        pos_centralized = class_segmentationMain::coordinate2Index(
            x, y, z, size_X, size_X * size_Y);
        Image1D_output[pos_centralized] = (int)this->Image1D_page[pos_voxel];
      }
      return;
    }

    void centralizeRegion(const vector<V3DLONG> poss_input,
                          const V3DLONG size_X, const V3DLONG size_Y,
                          const V3DLONG size_Z, const V3DLONG min_X,
                          const V3DLONG min_Y, const V3DLONG min_Z,
                          double ***Image3D_output) {
      V3DLONG pos_voxel = 0;
      vector<V3DLONG> xyz_voxel(0, 0);
      V3DLONG x = 0;
      V3DLONG y = 0;
      V3DLONG z = 0;
      V3DLONG count_voxel = poss_input.size();
      for (z = 0; z < size_Z; z++) {
        for (y = 0; y < size_Y; y++) {
          for (x = 0; x < size_X; x++) {
            Image3D_output[z][x][y] = 0;
          }
        }
      }
      for (V3DLONG idx_voxel = 0; idx_voxel < count_voxel; idx_voxel++) {
        pos_voxel = poss_input[idx_voxel];
        xyz_voxel = this->index2Coordinate(pos_voxel);
        x = xyz_voxel[0] - min_X;
        y = xyz_voxel[1] - min_Y;
        z = xyz_voxel[2] - min_Z;
        Image3D_output[z][x][y] = (double)this->Image1D_page[pos_voxel];
      }
      return;
    }

    /**
     * @brief Function to calculate the center of mass of a region
     * returns the center as the coordinate in a 1D array
     */
    V3DLONG getCenterByMass(vector<V3DLONG> vct_input) {
      vector<V3DLONG> xyz_voxel;
      V3DLONG x;
      V3DLONG y;
      V3DLONG z;
      double sum_X = 0;
      double sum_Y = 0;
      double sum_Z = 0;
      double sum_mass = 0;
      double value_voxel = 0;
      V3DLONG count_voxel = vct_input.size();
      if (count_voxel < 1) {
        return -1;
      }
      for (int i = 0; i < count_voxel; i++) {
        xyz_voxel = this->index2Coordinate(vct_input[i]);
        value_voxel = this->Image1D_page[vct_input[i]];
        x = xyz_voxel[0];
        y = xyz_voxel[1];
        z = xyz_voxel[2];
        sum_X += (double)x * value_voxel;
        sum_Y += (double)y * value_voxel;
        sum_Z += (double)z * value_voxel;
        sum_mass += value_voxel;
      }
      V3DLONG pos_result = this->coordinate2Index(
          sum_X / sum_mass, sum_Y / sum_mass, sum_Z / sum_mass);
      return pos_result;
    }

    vector<V3DLONG> getBoundBox(
        vector<V3DLONG>
            idxs_input)  // output: min_X, max_X, min_Y, max_Y, min_Z, max_Z;
    {
      V3DLONG x = 0;
      V3DLONG y = 0;
      V3DLONG z = 0;
      V3DLONG max_X = -INF;
      V3DLONG max_Y = -INF;
      V3DLONG max_Z = -INF;
      V3DLONG min_X = INF;
      V3DLONG min_Y = INF;
      V3DLONG min_Z = INF;
      V3DLONG count_voxel = idxs_input.size();
      vector<V3DLONG> xyz_voxel(3, 0);
      vector<V3DLONG> values_result(6, 0);
      V3DLONG idx_tmp;
      for (V3DLONG idx_voxel = 0; idx_voxel < count_voxel; idx_voxel++) {
        idx_tmp = idxs_input[idx_voxel];
        xyz_voxel = this->index2Coordinate(idx_tmp);
        x = xyz_voxel[0];
        y = xyz_voxel[1];
        z = xyz_voxel[2];
        if (x > max_X) {
          max_X = x;
        }
        if (y > max_Y) {
          max_Y = y;
        }
        if (z > max_Z) {
          max_Z = z;
        }
        if (x < min_X) {
          min_X = x;
        }
        if (y < min_Y) {
          min_Y = y;
        }
        if (z < min_Z) {
          min_Z = z;
        }
      }
      values_result[0] = min_X;
      values_result[1] = max_X;
      values_result[2] = min_Y;
      values_result[3] = max_Y;
      values_result[4] = min_Z;
      values_result[5] = max_Z;
      return values_result;
    }

    V3DLONG getMinDimension(
        vector<V3DLONG>
            vct_input)  // input: min_X, max_X, min_Y, max_Y, min_Z, max_Z;
    {
      V3DLONG size_X = vct_input[1] - vct_input[0];
      V3DLONG size_Y = vct_input[3] - vct_input[2];
      V3DLONG size_Z = vct_input[5] - vct_input[4];
      if (size_X < size_Y) {
        size_X = size_Y;
      }
      if (size_X < size_Z) {
        size_X = size_Z;
      }
      return size_X;
    }

    static vector<vector<V3DLONG> > int1D2possVct(
        int *Image1D_label, int count_label, unsigned char *Image1D_image,
        const V3DLONG size_X, const V3DLONG size_Y, const V3DLONG size_Z,
        const V3DLONG min_X, const V3DLONG min_Y, const V3DLONG min_Z,
        const V3DLONG offset_Yglobal, const V3DLONG offset_Zglobal,
        const double min_centerDistance, vector<V3DLONG> &poss_center) {
      int label_voxel = 0;
      V3DLONG pos_voxel = 0;
      double value_voxel = 0;
      vector<V3DLONG> vct_empty(0, 0);
      vector<vector<V3DLONG> > possVct_result;
      vector<vector<V3DLONG> > possVct_resultWithEmpty;
      for (int i = 0; i < count_label; i++) {
        possVct_resultWithEmpty.push_back(vct_empty);
      }
      double3D xyz_zero;
      vector<double3D> mean_center;
      for (int i = 0; i < count_label; i++) {
        mean_center.push_back(xyz_zero);
      }
      vector<int> idxs_remap(count_label, 0);
      V3DLONG label_remap = 0;
      V3DLONG x = 0;
      V3DLONG y = 0;
      V3DLONG z = 0;
      V3DLONG offset_Y = size_X;
      V3DLONG offset_Z = size_X * size_Y;
      vector<double> sums_mass(count_label, 0);
      for (x = 0; x < size_X; x++) {
        for (y = 0; y < size_Y; y++) {
          for (z = 0; z < size_Z; z++) {
            pos_voxel = class_segmentationMain::coordinate2Index(
                x, y, z, offset_Y, offset_Z);
            label_voxel = Image1D_label[pos_voxel];
            if (label_voxel > 0) {
              label_voxel = label_voxel - 1;
              value_voxel = Image1D_image[pos_voxel];
              mean_center[label_voxel].x += value_voxel * x;
              mean_center[label_voxel].y += value_voxel * y;
              mean_center[label_voxel].z += value_voxel * z;
              sums_mass[label_voxel] += value_voxel;
            }
          }
        }
      }
      for (int i = 0; i < count_label; i++) {
        idxs_remap[i] = i;
        if (sums_mass[i] > 0) {
          mean_center[i].x /= sums_mass[i];
          mean_center[i].y /= sums_mass[i];
          mean_center[i].z /= sums_mass[i];
        }
      }
      for (int i = 0; i < count_label; i++) {
        if (sums_mass[i] > 0) {
          for (int j = (i + 1); j < count_label; j++) {
            if (sums_mass[j] > 0) {
              if (class_segmentationMain::getEuclideanDistance(
                      mean_center[i], mean_center[j]) < min_centerDistance) {
                if (idxs_remap[j] == j) {
                  idxs_remap[j] = i;
                }
              }
            }
          }
        }
      }
      mean_center.clear();
      for (int i = 0; i < count_label; i++) {
        mean_center.push_back(xyz_zero);
      }
      fill(sums_mass.begin(), sums_mass.end(), 0);
      for (x = 0; x < size_X; x++) {
        for (y = 0; y < size_Y; y++) {
          for (z = 0; z < size_Z; z++) {
            pos_voxel = class_segmentationMain::coordinate2Index(
                x, y, z, offset_Y, offset_Z);
            label_voxel = Image1D_label[pos_voxel];
            if (label_voxel > 0) {
              label_voxel = label_voxel - 1;
              value_voxel = Image1D_image[pos_voxel];
              label_remap = idxs_remap[label_voxel];
              mean_center[label_remap].x += value_voxel * x;
              mean_center[label_remap].y += value_voxel * y;
              mean_center[label_remap].z += value_voxel * z;
              sums_mass[label_remap] += value_voxel;
              pos_voxel = class_segmentationMain::coordinate2Index(
                  x + min_X, y + min_Y, z + min_Z, offset_Yglobal,
                  offset_Zglobal);
              possVct_resultWithEmpty[label_remap].push_back(pos_voxel);
            }
          }
        }
      }
      for (int i = 0; i < count_label; i++) {
        if (sums_mass[i] > 0) {
          mean_center[i].x /= sums_mass[i];
          mean_center[i].y /= sums_mass[i];
          mean_center[i].z /= sums_mass[i];
        }
      }
      for (int i = 0; i < count_label; i++) {
        if ((sums_mass[i] > 0) && ((possVct_resultWithEmpty[i].size() >
                                    default_threshold_regionSize))) {
          possVct_result.push_back(possVct_resultWithEmpty[i]);
          poss_center.push_back(class_segmentationMain::coordinate2Index(
              mean_center[i].x + min_X, mean_center[i].y + min_Y,
              mean_center[i].z + min_Z, offset_Yglobal, offset_Zglobal));
        }
      }
      return possVct_result;
    }

    static vector<vector<V3DLONG> > int3D2possVct(
        int ***int3D_label, int count_label, double ***Image3D_image,
        const V3DLONG size_X, const V3DLONG size_Y, const V3DLONG size_Z,
        const V3DLONG min_X, const V3DLONG min_Y, const V3DLONG min_Z,
        const V3DLONG offset_Yglobal, const V3DLONG offset_Zglobal,
        vector<V3DLONG> &poss_center) {
      int label_voxel = 0;
      V3DLONG pos_voxel = 0;
      double value_voxel = 0;
      vector<V3DLONG> vct_empty(0, 0);
      vector<vector<V3DLONG> > possVct_result;
      vector<vector<V3DLONG> > possVct_resultWithEmpty;
      for (int i = 0; i < count_label; i++) {
        possVct_resultWithEmpty.push_back(vct_empty);
      }
      double3D xyz_zero;
      vector<double3D> mean_center;
      for (int i = 0; i < count_label; i++) {
        mean_center.push_back(xyz_zero);
      }
      V3DLONG x = 0;
      V3DLONG y = 0;
      V3DLONG z = 0;
      V3DLONG offset_Y = size_X;
      V3DLONG offset_Z = size_X * size_Y;
      vector<double> sums_mass(count_label, 0);
      for (x = 0; x < size_X; x++) {
        for (y = 0; y < size_Y; y++) {
          for (z = 0; z < size_Z; z++) {
            label_voxel = int3D_label[z][x][y];
            if (label_voxel > 0) {
              label_voxel = label_voxel - 1;
              value_voxel = Image3D_image[z][x][y];
              pos_voxel = class_segmentationMain::coordinate2Index(
                  x + min_X, y + min_Y, z + min_Z, offset_Yglobal,
                  offset_Zglobal);
              possVct_resultWithEmpty[label_voxel].push_back(pos_voxel);
              mean_center[label_voxel].x += value_voxel * x;
              mean_center[label_voxel].y += value_voxel * y;
              mean_center[label_voxel].z += value_voxel * z;
              sums_mass[label_voxel] += value_voxel;
            }
          }
        }
      }
      for (int i = 0; i < count_label; i++) {
        if (sums_mass[i] > 0) {
          mean_center[i].x /= sums_mass[i];
          mean_center[i].y /= sums_mass[i];
          mean_center[i].z /= sums_mass[i];
        }
      }
      for (int i = 0; i < count_label; i++) {
        if ((sums_mass[i] > 0) && ((possVct_resultWithEmpty[i].size() >
                                    default_threshold_regionSize))) {
          possVct_result.push_back(possVct_resultWithEmpty[i]);
          poss_center.push_back(class_segmentationMain::coordinate2Index(
              mean_center[i].x + min_X, mean_center[i].y + min_Y,
              mean_center[i].z + min_Z, offset_Yglobal, offset_Zglobal));
        }
      }
      return possVct_result;
    }
#pragma endregion

#pragma region "memoryManagement"
    static double3D ***memory_allocate_double3D3(const int i_size,
                                                 const int j_size,
                                                 const int k_size) {
      double3D ***ptr_result;
      int i, k;
      ptr_result = (double3D ***)calloc(k_size, sizeof(double3D **));
      for (k = 0; k < k_size; k++) {
        ptr_result[k] = (double3D **)calloc(i_size, sizeof(double3D *));
      }
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          ptr_result[k][i] = (double3D *)calloc(j_size, sizeof(double3D));
        }
      }
      return (ptr_result);
    }

    static void memory_free_double3D3(double3D ***ptr_input, const int k_size,
                                      const int i_size) {
      int k, i;
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          free(ptr_input[k][i]);
        }
      }
      for (k = 0; k < k_size; k++) {
        free(ptr_input[k]);
      }
      free(ptr_input);
    }

    static double ***memory_allocate_double3D(int i_size, int j_size,
                                              int k_size) {
      double ***ptr_result;
      int i, k;
      ptr_result = (double ***)calloc(k_size, sizeof(double **));
      for (k = 0; k < k_size; k++) {
        ptr_result[k] = (double **)calloc(i_size, sizeof(double *));
      }
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          ptr_result[k][i] = (double *)calloc(j_size, sizeof(double));
        }
      }
      return (ptr_result);
    }

    static void memory_free_double3D(double ***ptr_input, const int k_size,
                                     const int i_size) {
      int k, i;
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          free(ptr_input[k][i]);
        }
      }
      for (k = 0; k < k_size; k++) {
        free(ptr_input[k]);
      }
      free(ptr_input);
    }

    static unsigned char ***memory_allocate_uchar3D(const int i_size,
                                                    const int j_size,
                                                    const int k_size) {
      unsigned char ***ptr_result;
      int i, k;
      ptr_result = (unsigned char ***)calloc(k_size, sizeof(unsigned char **));
      for (k = 0; k < k_size; k++) {
        ptr_result[k] =
            (unsigned char **)calloc(i_size, sizeof(unsigned char *));
      }
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          ptr_result[k][i] =
              (unsigned char *)calloc(j_size, sizeof(unsigned char));
        }
      }
      return (ptr_result);
    }

    static void memory_free_uchar3D(unsigned char ***ptr_input,
                                    const int k_size, const int i_size) {
      int k, i;
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          free(ptr_input[k][i]);
        }
      }
      for (k = 0; k < k_size; k++) {
        free(ptr_input[k]);
      }
      free(ptr_input);
    }

    static unsigned char **memory_allocate_uchar2D(const V3DLONG i_size,
                                                   const V3DLONG j_size) {
      unsigned char **ptr_result;
      V3DLONG i;
      ptr_result = (unsigned char **)calloc(i_size, sizeof(unsigned char *));
      for (i = 0; i < i_size; i++) {
        ptr_result[i] = (unsigned char *)calloc(j_size, sizeof(unsigned char));
      }
      return (ptr_result);
    }

    static void memory_free_uchar2D(unsigned char **ptr_input,
                                    const V3DLONG i_size) {
      V3DLONG i;
      for (i = 0; i < i_size; i++) free(ptr_input[i]);
      free(ptr_input);
    }

    static long3D ***memory_allocate_int3D3(const int i_size, const int j_size,
                                            const int k_size) {
      long3D ***ptr_result;
      int i, k;
      ptr_result = (long3D ***)calloc(k_size, sizeof(long3D **));
      for (k = 0; k < k_size; k++) {
        ptr_result[k] = (long3D **)calloc(i_size, sizeof(long3D *));
      }
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          ptr_result[k][i] = (long3D *)calloc(j_size, sizeof(long3D));
        }
      }
      return (ptr_result);
    }

    static void memory_free_int3D3(long3D ***ptr_input, const int k_size,
                                   const int i_size) {
      int k, i;
      for (k = 0; k < k_size; k++) {
        for (i = 0; i < i_size; i++) {
          free(ptr_input[k][i]);
        }
      }
      for (k = 0; k < k_size; k++) {
        free(ptr_input[k]);
      }
      free(ptr_input);
    }

    static int **memory_allocate_int2D(const int i_size, const int j_size) {
      int **ptr_result;
      int i;
      ptr_result = (int **)calloc(i_size, sizeof(int *));
      for (i = 0; i < i_size; i++) {
        ptr_result[i] = (int *)calloc(j_size, sizeof(int));
      }
      return (ptr_result);
    }

    static void memory_free_int2D(int **ptr_input, const int i_size) {
      int i;

      for (i = 0; i < i_size; i++) free(ptr_input[i]);
      free(ptr_input);
    }

    static double **memory_allocate_double2D(const int i_size,
                                             const int j_size) {
      double **ptr_result;
      int i;
      ptr_result = (double **)calloc(i_size, sizeof(double *));
      for (i = 0; i < i_size; i++) {
        ptr_result[i] = (double *)calloc(j_size, sizeof(double));
      }
      return (ptr_result);
    }

    static void memory_free_double2D(double **ptr_input, const int i_size) {
      int i;

      for (i = 0; i < i_size; i++) free(ptr_input[i]);
      free(ptr_input);
    }

    static V3DLONG *memory_allocate_int1D(const V3DLONG i_size) {
      V3DLONG *ptr_result;
      ptr_result = (V3DLONG *)calloc(i_size, sizeof(V3DLONG));
      return (ptr_result);
    }

    static void memory_free_int1D(int *ptr_input, const int i_size) {
      free(ptr_input);
    }

    static void memory_free_double2D(bool **ptr_input, const int i_size) {
      int i;

      for (i = 0; i < i_size; i++) free(ptr_input[i]);
      free(ptr_input);
    }

    static unsigned char *memory_allocate_uchar1D(const int i_size) {
      unsigned char *ptr_result;
      ptr_result = (unsigned char *)calloc(i_size, sizeof(unsigned char));
      return (ptr_result);
    }

    static void memory_free_uchar1D(unsigned char *ptr_input) {
      free(ptr_input);
    }

    static int ***memory_allocate_int3D(int i_size, int j_size, int k_size) {
      int ***ptr_result;
      int i, k;
      ptr_result = (int ***)calloc(k_size, sizeof(int **));
      for (k = 0; k < k_size; k++)
        ptr_result[k] = (int **)calloc(i_size, sizeof(int *));
      for (k = 0; k < k_size; k++)
        for (i = 0; i < i_size; i++)
          ptr_result[k][i] = (int *)calloc(j_size, sizeof(int));
      return (ptr_result);
    }

    static void memory_free_int3D(int ***ptr_input, int k_size, int i_size) {
      int k, i;
      for (k = 0; k < k_size; k++)
        for (i = 0; i < i_size; i++) free(ptr_input[k][i]);
      for (k = 0; k < k_size; k++) free(ptr_input[k]);
      free(ptr_input);
    }
#pragma endregion

#pragma region "otsu thresholding"
    /**
     * @brief Global Otsu threshold: use the histogram of the entire image
     */
    V3DLONG globalOtsuThreshold() {
      int hist[256] = {0};
      for (V3DLONG i = 0; i < size_page; i++) {
        hist[Image1D_page[i]]++;
      }
      int total = size_page;
      double sum = 0;
      for (int t = 0; t < 256; t++) {
        sum += t * hist[t];
      }
      double sumB = 0;
      int wB = 0;
      double varMax = 0;
      int threshold = 0;
      for (int t = 0; t < 256; t++) {
        wB += hist[t];
        if (wB == 0) continue;
        int wF = total - wB;
        if (wF == 0) break;
        sumB += t * hist[t];
        double mB = sumB / wB;
        double mF = (sum - sumB) / wF;
        double varBetween = (double)wB * wF * (mB - mF) * (mB - mF);
        if (varBetween > varMax) {
          varMax = varBetween;
          threshold = t;
        }
      }
      return threshold;
    }

    /**
     * @brief Local Otsu threshold: compute a threshold using only voxels
     * within a cubic region of the given radius around a given landmark.
     */
    V3DLONG localOtsuThreshold(V3DLONG landmarkIndex, V3DLONG radius) {
      vector<V3DLONG> coord = index2Coordinate(landmarkIndex);
      V3DLONG cx = coord[0], cy = coord[1], cz = coord[2];
      int hist[256] = {0};
      V3DLONG x_start = max((V3DLONG)0, cx - radius);
      V3DLONG x_end = min(dim_X - 1, cx + radius);
      V3DLONG y_start = max((V3DLONG)0, cy - radius);
      V3DLONG y_end = min(dim_Y - 1, cy + radius);
      V3DLONG z_start = max((V3DLONG)0, cz - radius);
      V3DLONG z_end = min(dim_Z - 1, cz + radius);
      int count_pixels = 0;
      for (V3DLONG z = z_start; z <= z_end; z++) {
        for (V3DLONG y = y_start; y <= y_end; y++) {
          for (V3DLONG x = x_start; x <= x_end; x++) {
            V3DLONG idx = coordinate2Index(x, y, z);
            hist[Image1D_page[idx]]++;
            count_pixels++;
          }
        }
      }
      if (count_pixels == 0) return default_threshold_global;  // fallback value
      double sum = 0;
      for (int t = 0; t < 256; t++) {
        sum += t * hist[t];
      }
      double sumB = 0;
      int wB = 0;
      double varMax = 0;
      int threshold = 0;
      for (int t = 0; t < 256; t++) {
        wB += hist[t];
        if (wB == 0) continue;
        int wF = count_pixels - wB;
        if (wF == 0) break;
        sumB += t * hist[t];
        double mB = sumB / wB;
        double mF = (sum - sumB) / wF;
        double varBetween = (double)wB * wF * (mB - mF) * (mB - mF);
        if (varBetween > varMax) {
          varMax = varBetween;
          threshold = t;
        }
      }
      return threshold;
    }
#pragma endregion

#pragma region "smoothing and filtering"
    /**
     * @brief Apply a median filter to the image with a given radius.
     */
    void filter_Median(V3DLONG radius) {
      if (radius < 1) return;

      // Allocate a temporary output buffer.
      unsigned char *Image1D_output = memory_allocate_uchar1D(this->size_page);

      // Process each slice.
      for (V3DLONG iz = 0; iz < this->dim_Z; iz++) {
        // Print progress.
        cout << "\r median filter, " << (double)(iz + 1) * 100.0 / this->dim_Z
             << "% completed;" << flush;

        V3DLONG offsetZ = iz * this->offset_Z;
        for (V3DLONG iy = 0; iy < this->dim_Y; iy++) {
          V3DLONG offsetY = iy * this->offset_Y;
          for (V3DLONG ix = 0; ix < this->dim_X; ix++) {
            // Compute window bounds (clamped to image boundaries)
            V3DLONG xb = (ix >= radius) ? ix - radius : 0;
            V3DLONG xe =
                (ix + radius < this->dim_X) ? ix + radius : this->dim_X - 1;
            V3DLONG yb = (iy >= radius) ? iy - radius : 0;
            V3DLONG ye =
                (iy + radius < this->dim_Y) ? iy + radius : this->dim_Y - 1;
            V3DLONG zb = (iz >= radius) ? iz - radius : 0;
            V3DLONG ze =
                (iz + radius < this->dim_Z) ? iz + radius : this->dim_Z - 1;

            // Build a histogram of the intensities in the current window.
            int hist[256] = {0};
            int count = 0;
            for (V3DLONG k = zb; k <= ze; k++) {
              V3DLONG offsetK = k * this->offset_Z;
              for (V3DLONG j = yb; j <= ye; j++) {
                V3DLONG offsetJ = j * this->offset_Y;
                for (V3DLONG i = xb; i <= xe; i++) {
                  unsigned char val = this->Image1D_page[offsetK + offsetJ + i];
                  hist[val]++;
                  count++;
                }
              }
            }

            // Find the median value by accumulating histogram counts.
            int mid = count / 2;
            int sum = 0, median = 0;
            for (int v = 0; v < 256; v++) {
              sum += hist[v];
              if (sum > mid) {
                median = v;
                break;
              }
            }

            // Set the output voxel.
            Image1D_output[offsetZ + offsetY + ix] = (unsigned char)median;
          }
        }
      }

      // Copy the filtered result back to the main image array.
      memcpy(this->Image1D_page, Image1D_output,
             this->size_page * sizeof(unsigned char));
      memory_free_uchar1D(Image1D_output);
    }

    static void smooth_GVFkernal(double ***Image3D_input,
                                 int count_smoothIteration, int dim_X,
                                 int dim_Y, int dim_Z) {
      int i, x, y, z;
      double ***Image3D_update;
      Image3D_update = memory_allocate_double3D(dim_X, dim_Y, dim_Z);
      for (i = 0; i < count_smoothIteration; i++) {
        for (z = 0; z < dim_Z; z++)
          for (x = 0; x < dim_X; x++)
            for (y = 0; y < dim_Y; y++) Image3D_update[z][x][y] = 0;

        for (z = 1; z < dim_Z - 1; z++)
          for (x = 1; x < dim_X - 1; x++)
            for (y = 1; y < dim_Y - 1; y++)
              Image3D_update[z][x][y] =
                  0.4 * Image3D_input[z][x][y] +
                  0.1 *
                      (Image3D_input[z - 1][x][y] + Image3D_input[z + 1][x][y] +
                       Image3D_input[z][x - 1][y] + Image3D_input[z][x + 1][y] +
                       Image3D_input[z][x][y - 1] + Image3D_input[z][x][y + 1]);

        for (z = 0; z < dim_Z; z++)
          for (x = 0; x < dim_X; x++)
            for (y = 0; y < dim_Y; y++)
              Image3D_input[z][x][y] = Image3D_update[z][x][y];
      }
      memory_free_double3D(Image3D_update, dim_Z, dim_X);
      return;
    }

    void smooth_GVFkernal(V3DLONG count_smoothRadius) {
      unsigned char *Image1D_update = memory_allocate_uchar1D(this->size_page);
      V3DLONG pos_voxel;
      V3DLONG pos_neighbor1;
      V3DLONG pos_neighbor2;
      V3DLONG pos_neighbor3;
      V3DLONG pos_neighbor4;
      V3DLONG pos_neighbor5;
      V3DLONG pos_neighbor6;
      for (int i = 0; i < count_smoothRadius; i++) {
        memset(Image1D_update, 0, this->size_page);
        for (int z = 1; z < this->dim_Z - 1; z++) {
          for (int x = 1; x < this->dim_X - 1; x++) {
            for (int y = 1; y < this->dim_Y - 1; y++) {
              pos_voxel = this->coordinate2Index(x, y, z);
              pos_neighbor1 = this->coordinate2Index(x, y, z - 1);
              pos_neighbor2 = this->coordinate2Index(x, y, z + 1);
              pos_neighbor3 = this->coordinate2Index(x, y - 1, z);
              pos_neighbor4 = this->coordinate2Index(x, y + 1, z);
              pos_neighbor5 = this->coordinate2Index(x - 1, y, z);
              pos_neighbor6 = this->coordinate2Index(x + 1, y, z);
              Image1D_update[pos_voxel] =
                  0.4 * this->Image1D_page[pos_voxel] +
                  0.1 * (this->Image1D_page[pos_neighbor1] +
                         this->Image1D_page[pos_neighbor2] +
                         this->Image1D_page[pos_neighbor3] +
                         this->Image1D_page[pos_neighbor4] +
                         this->Image1D_page[pos_neighbor5] +
                         this->Image1D_page[pos_neighbor6]);
            }
          }
        }
        for (int i = 0; i < this->size_page; i++) {
          this->Image1D_page[i] = Image1D_update[i];
        }
      }
      return;
    }
#pragma endregion

#pragma region "shapeStat"
    vector<vector<double> > getShapeStat(V3DLONG x, V3DLONG y, V3DLONG z,
                                         V3DLONG value_radius) {
      vector<vector<double> > valuesVct_result;
      double value_PC1 = 0;
      double value_PC2 = 0;
      double value_PC3 = 0;
      vector<double> values_PC1;
      vector<double> values_PC2;
      vector<double> values_PC3;
      double size_step = (double)(value_radius - 2) / 3.0;
      V3DLONG rr = 0;
      for (int i = 1; i <= 4; i++) {
        rr = 2 + size_step * (i - 1);
        // bool is_valid = false;
        /*if(getPCA(this->Image3D_page, this->dim_X, this->dim_Y, this->dim_Z,
        x , y, z, rr, rr, rr,value_PC1, value_PC2, value_PC3, this->idx_shape,
        false))
        {
                is_valid = true;
        }*/
        getPCA(this->Image3D_page, this->dim_X, this->dim_Y, this->dim_Z, x, y,
               z, rr, rr, rr, value_PC1, value_PC2, value_PC3, this->idx_shape,
               false, false);
        values_PC1.push_back(value_PC1 / rr);
        values_PC2.push_back(value_PC2 / rr);
        values_PC3.push_back(value_PC3 / rr);
      }
      valuesVct_result.push_back(values_PC1);
      valuesVct_result.push_back(values_PC2);
      valuesVct_result.push_back(values_PC3);
      /*double value_Lscore = exp( -(
         (value_PC1-value_PC2)*(value_PC1-value_PC2) +
         (value_PC2-value_PC3)*(value_PC2-value_PC3) +
         (value_PC1-value_PC3)*(value_PC1-value_PC3) ) / (value_PC1*value_PC1
         + value_PC2*value_PC2 + value_PC3*value_PC3) );*/
      /*double value_linear =
      (value_PC1-value_PC2)/(value_PC1+value_PC2+value_PC3); double
      value_planar
      = 2.0*(value_PC2-value_PC3)/(value_PC1+value_PC2+value_PC3); double
      value_sphere = 3.0*value_PC3/(value_PC1+value_PC2+value_PC3);*/
      // vct_result.push_back(value_Lscore);
      // vct_result.push_back(value_linear);
      // vct_result.push_back(value_planar);
      // vct_result.push_back(value_sphere);
      return valuesVct_result;
    }

    /**
     * @brief get the shape analysis depending on whether it is a sphere or a
     * cube
     */
    template <class T>
    bool getPCA(
        T ***img3d, V3DLONG sx, V3DLONG sy, V3DLONG sz, V3DLONG x0, V3DLONG y0,
        V3DLONG z0, V3DLONG rx, V3DLONG ry, V3DLONG rz, double &pc1,
        double &pc2, double &pc3, int wintype = 0,
        bool b_disp_CoM_etc = true,      // b_disp_CoM_etc is the display option
                                         // for center of mass )
        bool b_normalize_score = false)  // if the score if normalized with
                                         // respect to the window size
    {
      if (wintype == 0)
        return getPCA_cube(img3d, sx, sy, sz, x0, y0, z0, rx, ry, rz, pc1, pc2,
                           pc3, b_disp_CoM_etc, b_normalize_score);
      else  // wintype==1
        return getPCA_sphere(img3d, sx, sy, sz, x0, y0, z0, rx, ry, rz, pc1,
                             pc2, pc3, b_disp_CoM_etc, b_normalize_score);
    }

    /**
     * @brief get the PCA for a sphereical shape
     */
    template <class T>
    bool getPCA_sphere(
        T ***img3d, V3DLONG sx, V3DLONG sy, V3DLONG sz, V3DLONG x0, V3DLONG y0,
        V3DLONG z0, V3DLONG rx, V3DLONG ry, V3DLONG rz, double &pc1,
        double &pc2, double &pc3,
        bool b_disp_CoM_etc = true,  // b_disp_CoM_etc is the display option
                                     // for center of mass )
        bool b_normalize_score = false) {
      if (!img3d || sx <= 0 || sy <= 0 || sz <= 0 || x0 < 0 || x0 >= sx ||
          y0 < 0 || y0 >= sy || z0 < 0 || z0 >= sz || rx < 0 || ry < 0 ||
          rz < 0)
        return false;

      // get max radius
      V3DLONG maxrr = (rx > ry) ? rx : ry;
      maxrr = (maxrr > rz) ? maxrr : rz;

      // get the boundary

      V3DLONG xb = x0 - rx;
      if (xb < 0)
        xb = 0;
      else if (xb >= sx)
        xb = sx - 1;
      V3DLONG xe = x0 + rx;
      if (xe < 0)
        xe = 0;
      else if (xe >= sx)
        xe = sx - 1;
      V3DLONG yb = y0 - ry;
      if (yb < 0)
        yb = 0;
      else if (yb >= sy)
        yb = sy - 1;
      V3DLONG ye = y0 + ry;
      if (ye < 0)
        ye = 0;
      else if (ye >= sy)
        ye = sy - 1;
      V3DLONG zb = z0 - rz;
      if (zb < 0)
        zb = 0;
      else if (zb >= sz)
        zb = sz - 1;
      V3DLONG ze = z0 + rz;
      if (ze < 0)
        ze = 0;
      else if (ze >= sz)
        ze = sz - 1;

      V3DLONG i, j, k;
      double w;

      // first get the center of mass
      double x2, y2, z2;
      double rx2 = double(rx + 1) * (rx + 1), ry2 = (double)(ry + 1) * (ry + 1),
             rz2 =
                 (double)(rz + 1) *
                 (rz + 1);  //+1 because later need to do use it for radius cmp
      double tmpd;
      double xm = 0, ym = 0, zm = 0, s = 0, mv = 0, n = 0;
      for (k = zb; k <= ze; k++) {
        z2 = k - z0;
        z2 *= z2;
        for (j = yb; j <= ye; j++) {
          y2 = j - y0;
          y2 *= y2;
          tmpd = y2 / ry2 + z2 / rz2;
          if (tmpd > 1.0) continue;

          for (i = xb; i <= xe; i++) {
            x2 = i - x0;
            x2 *= x2;
            if (x2 / rx2 + tmpd > 1.0) continue;

            w = double(img3d[k][j][i]);
            xm += w * i;
            ym += w * j;
            zm += w * k;
            s += w;
            n = n + 1;
          }
        }
      }
      if (s > 0) {
        xm /= s;
        ym /= s;
        zm /= s;
        mv = s / n;
        // if (b_disp_CoM_etc)
        //{
        // printf("center of mass is (xm, ym, zm) = %5.3f, %5.3f,
        // %5.3f\n",xm,ym,zm);
        // }

      } else {
        // printf("Sum of window pixels equals or is smaller than 0. The
        // window is not valid or some other problems in the data. Do
        // nothing.\n");
        return false;
      }

      // get the covariance. Note that the center of mass must be in the
      // ellpsoid

      double cc11 = 0, cc12 = 0, cc13 = 0, cc22 = 0, cc23 = 0, cc33 = 0;
      double dfx, dfy, dfz;
      for (k = zb; k <= ze; k++) {
        z2 = k - z0;
        z2 *= z2;

        dfz = double(k) - zm;
        if (b_normalize_score) dfz /= maxrr;

        for (j = yb; j <= ye; j++) {
          y2 = j - y0;
          y2 *= y2;
          tmpd = y2 / ry2 + z2 / rz2;
          if (tmpd > 1.0) continue;

          dfy = double(j) - ym;
          if (b_normalize_score) dfy /= maxrr;

          for (i = xb; i <= xe; i++) {
            x2 = i - x0;
            x2 *= x2;
            if (x2 / rx2 + tmpd > 1.0) continue;

            dfx = double(i) - xm;
            if (b_normalize_score) dfx /= maxrr;

            //                w = img3d[k][j][i]; //140128
            w = img3d[k][j][i] - mv;
            if (w < 0) w = 0;  // 140128 try the new formula

            cc11 += w * dfx * dfx;
            cc12 += w * dfx * dfy;
            cc13 += w * dfx * dfz;
            cc22 += w * dfy * dfy;
            cc23 += w * dfy * dfz;
            cc33 += w * dfz * dfz;
          }
        }
      }

      cc11 /= s;
      cc12 /= s;
      cc13 /= s;
      cc22 /= s;
      cc23 /= s;
      cc33 /= s;
      // if (b_disp_CoM_etc)
      // printf("convariance value (c11,c12,c13,c22,c23,c33) = %5.3f, %5.3f,
      // %5.3f, %5.3f, %5.3f, %5.3f\n",cc11, cc12, cc13, cc22, cc23, cc33);

      // now get the eigen vectors and eigen values

      try {
        // then find the eigen vector
        SymmetricMatrix Cov_Matrix(3);
        Cov_Matrix.Row(1) << cc11;
        Cov_Matrix.Row(2) << cc12 << cc22;
        Cov_Matrix.Row(3) << cc13 << cc23 << cc33;

        DiagonalMatrix DD;
        Matrix VV;
        EigenValues(Cov_Matrix, DD, VV);

        // output the result
        pc1 = DD(3);
        pc2 = DD(2);
        pc3 = DD(1);
      } catch (...) {
        pc1 = VAL_INVALID;
        pc2 = VAL_INVALID;
        pc3 = VAL_INVALID;
      }

      return true;
    }

    /**
     * @brief get the PCA for a cube shape
     */
    template <class T>
    bool getPCA_cube(
        T ***img3d, V3DLONG sx, V3DLONG sy, V3DLONG sz, V3DLONG x0, V3DLONG y0,
        V3DLONG z0, V3DLONG rx, V3DLONG ry, V3DLONG rz, double &pc1,
        double &pc2, double &pc3,
        bool b_disp_CoM_etc =
            true,  // b_disp_CoM_etc is the display option for center of mass
        bool b_normalize_score = false) {
      if (!img3d || sx <= 0 || sy <= 0 || sz <= 0 || x0 < 0 || x0 >= sx ||
          y0 < 0 || y0 >= sy || z0 < 0 || z0 >= sz || rx < 0 || ry < 0 ||
          rz < 0)
        return false;

      // get max radius
      V3DLONG maxrr = (rx > ry) ? rx : ry;
      maxrr = (maxrr > rz) ? maxrr : rz;

      // get the boundary

      V3DLONG xb = x0 - rx;
      if (xb < 0)
        xb = 0;
      else if (xb >= sx)
        xb = sx - 1;
      V3DLONG xe = x0 + rx;
      if (xe < 0)
        xe = 0;
      else if (xe >= sx)
        xe = sx - 1;
      V3DLONG yb = y0 - ry;
      if (yb < 0)
        yb = 0;
      else if (yb >= sy)
        yb = sy - 1;
      V3DLONG ye = y0 + ry;
      if (ye < 0)
        ye = 0;
      else if (ye >= sy)
        ye = sy - 1;
      V3DLONG zb = z0 - rz;
      if (zb < 0)
        zb = 0;
      else if (zb >= sz)
        zb = sz - 1;
      V3DLONG ze = z0 + rz;
      if (ze < 0)
        ze = 0;
      else if (ze >= sz)
        ze = sz - 1;

      V3DLONG i, j, k;
      double w;

      // first get the center of mass
      double xm = 0, ym = 0, zm = 0, s = 0, mv = 0;
      for (k = zb; k <= ze; k++) {
        for (j = yb; j <= ye; j++) {
          for (i = xb; i <= xe; i++) {
            w = double(img3d[k][j][i]);
            xm += w * i;
            ym += w * j;
            zm += w * k;
            s += w;
          }
        }
      }

      if (s > 0) {
        xm /= s;
        ym /= s;
        zm /= s;
        mv = s / (double(ze - zb + 1) * (ye - yb + 1) * (xe - xb + 1));
        if (b_disp_CoM_etc)
          printf("center of mass is (xm, ym, zm) = %5.3f, %5.3f, %5.3f\n", xm,
                 ym, zm);
      } else {
        printf(
            "Sum of window pixels equals or is smaller than 0. The window is "
            "not valid or some other problems in the data. Do nothing.\n");
        return false;
      }

      // get the covariance

      double cc11 = 0, cc12 = 0, cc13 = 0, cc22 = 0, cc23 = 0, cc33 = 0;
      double dfx, dfy, dfz;
      for (k = zb; k <= ze; k++) {
        dfz = double(k) - zm;
        if (b_normalize_score) dfz /= maxrr;
        for (j = yb; j <= ye; j++) {
          dfy = double(j) - ym;
          if (b_normalize_score) dfy /= maxrr;
          for (i = xb; i <= xe; i++) {
            dfx = double(i) - xm;
            if (b_normalize_score) dfx /= maxrr;

            //                w = img3d[k][j][i]; //140128
            w = img3d[k][j][i] - mv;
            if (w < 0) w = 0;  // 140128 try the new formula

            cc11 += w * dfx * dfx;
            cc12 += w * dfx * dfy;
            cc13 += w * dfx * dfz;
            cc22 += w * dfy * dfy;
            cc23 += w * dfy * dfz;
            cc33 += w * dfz * dfz;
          }
        }
      }

      cc11 /= s;
      cc12 /= s;
      cc13 /= s;
      cc22 /= s;
      cc23 /= s;
      cc33 /= s;
      if (b_disp_CoM_etc)
        printf(
            "convariance value (c11,c12,c13,c22,c23,c33) = %5.3f, %5.3f, "
            "%5.3f, %5.3f, %5.3f, %5.3f\n",
            cc11, cc12, cc13, cc22, cc23, cc33);

      // now get the eigen vectors and eigen values

      try {
        // then find the eigen vector
        SymmetricMatrix Cov_Matrix(3);
        Cov_Matrix.Row(1) << cc11;
        Cov_Matrix.Row(2) << cc12 << cc22;
        Cov_Matrix.Row(3) << cc13 << cc23 << cc33;

        DiagonalMatrix DD;
        Matrix VV;
        EigenValues(Cov_Matrix, DD, VV);

        // output the result
        pc1 = DD(3);
        pc2 = DD(2);
        pc3 = DD(1);
      } catch (...) {
        pc1 = VAL_INVALID;
        pc2 = VAL_INVALID;
        pc3 = VAL_INVALID;
      }

      return true;
    }

#pragma endregion
  };
#pragma endregion

  class_segmentationMain class_segmentationMain1;

#pragma region "interface"
  /**
   * @brief main method that is called when the plugin is clicked
   */
  bool interface_run(V3DPluginCallback2 &_V3DPluginCallback2_currentCallback,
                     QWidget *_QWidget_parent) {
    // generic checks to make sure that the plugin can run
    v3dhandle v3dhandle_currentWindow =
        _V3DPluginCallback2_currentCallback.currentImageWindow();
    if (!v3dhandle_currentWindow) {
      v3d_msg(
          "You have not loaded any image or the image is corrupted, program "
          "canceled!");
      return false;
    }
    Image4DSimple *Image4DSimple_current =
        _V3DPluginCallback2_currentCallback.getImage(v3dhandle_currentWindow);
    if (!Image4DSimple_current) {
      v3d_msg(
          "You have not loaded any image or the image is corrupted, program "
          "canceled!");
      return false;
    }
    V3DLONG count_totalBytes = Image4DSimple_current->getTotalBytes();
    if (count_totalBytes < 1) {
      v3d_msg(
          "You have not loaded any image or the image is corrupted, program "
          "canceled!");
      return false;
    }
    unsigned char *Image1D_current = Image4DSimple_current->getRawData();
    QString name_currentWindow =
        _V3DPluginCallback2_currentCallback.getImageName(
            v3dhandle_currentWindow);

    // get name of the image
    QString fileName = Image4DSimple_current->getFileName();

    bool isTeraFly = false;
    if (fileName.startsWith("ID")) {
      isTeraFly = true;
    }

    // modify name if necessary for TeraFly
    fileName = modifyFileNameForTeraFly(fileName);

    // get image and landmarks
    V3DLONG dim_X = Image4DSimple_current->getXDim();
    V3DLONG dim_Y = Image4DSimple_current->getYDim();
    V3DLONG dim_Z = Image4DSimple_current->getZDim();
    V3DLONG dim_C = Image4DSimple_current->getCDim();
    V3DLONG size_image = dim_X * dim_Y * dim_Z * dim_C;
    LandmarkList LandmarkList_userDefined =
        _V3DPluginCallback2_currentCallback.getLandmark(
            v3dhandle_currentWindow);
    V3DLONG count_userDefinedLandmarkList = LandmarkList_userDefined.count();

    // can also get SWC files if they are present - in our case, they are not
    // needed
    QList<NeuronTree> *SWCList_current =
        _V3DPluginCallback2_currentCallback.getHandleNeuronTrees_3DGlobalViewer(
            v3dhandle_currentWindow);
    V3DLONG count_SWCList = 0;
    if (SWCList_current) {
      count_SWCList = SWCList_current->count();
    }

    // check to make sure that landmarks are defined
    LandmarkList LandmarkList_current;
    V3DLONG count_currentLandmarkList = -1;

    // cases where there are swc files (doesn't apply to us)
    if ((count_SWCList < 1) && (count_userDefinedLandmarkList < 1)) {
      v3d_msg(
          "You have not defined any landmarks or swc structure to run the "
          "segmenation, program canceled!");
      return false;
    } else if ((count_SWCList > 0) && (count_userDefinedLandmarkList > 0)) {
      LandmarkList_current = LandmarkList_userDefined;
      class_segmentationMain::neuronTree2LandmarkList(SWCList_current->first(),
                                                      LandmarkList_current);
      count_currentLandmarkList = LandmarkList_current.count();
    } else if ((count_SWCList > 0) && (count_userDefinedLandmarkList < 1)) {
      class_segmentationMain::neuronTree2LandmarkList(SWCList_current->first(),
                                                      LandmarkList_current);
      count_currentLandmarkList = LandmarkList_current.count();
    }
    // this is triggered by our code - we have landmarks
    if (count_userDefinedLandmarkList > 0) {
      LandmarkList_current = LandmarkList_userDefined;
      count_currentLandmarkList = LandmarkList_current.count();
    }

    // give the user the dialog
    dialogRun dialogRun1(_V3DPluginCallback2_currentCallback, _QWidget_parent,
                         dim_C);

    bool is_success = false;

    /*if (this->class_segmentationMain1.is_initialized) //temporary solution
    for the "parameter window not popped up" problem;
    {
            is_success =
    this->class_segmentationMain1.control_run(this->class_segmentationMain1.Image1D_page,
    this->class_segmentationMain1.dim_X, this->class_segmentationMain1.dim_Y,
    this->class_segmentationMain1.dim_Z,
    this->class_segmentationMain1.idx_channel, LandmarkList_current,
    this->class_segmentationMain1.idx_shape,
                    this->class_segmentationMain1.threshold_deltaShapeStat,
    this->class_segmentationMain1.multiplier_thresholdRegionSize,
                    this->class_segmentationMain1.multiplier_uThresholdRegionSize,
    this->class_segmentationMain1.name_currentWindow,
                    this->class_segmentationMain1.max_movment1,
    this->class_segmentationMain1.max_movment1);
    }*/
    // else
    {
      if (dialogRun1.exec() != QDialog::Accepted) {
        return false;
      }
      // Set the median filtering flag from the dialog
      this->class_segmentationMain1.applyMedianFiltering =
          dialogRun1.applyMedianFiltering;
      // Set the median filtering radius from the dialog
      this->class_segmentationMain1.medianFilteringRadius =
          dialogRun1.medianFilteringRadius;
      // Set the marker flag from the dialog
      this->class_segmentationMain1.applyMarkerConstraint =
          dialogRun1.applyMarkerConstraint;
      // Set the manual thresholding flag from the dialog
      this->class_segmentationMain1.manualThresholding =
          dialogRun1.manualThresholding;
      int idx_shape;  // get shape paramters;
      if (dialogRun1.shape_type_selection == sphere) {
        idx_shape = 1;
      } else if (dialogRun1.shape_type_selection == cube) {
        idx_shape = 0;
      }
      // call control run method to do segmentation
      // In interface_run (or wherever control_run is called), pass the
      // selected mode: (Assume dialogRun1.segmentationMode is set from the
      // combo box.)
      // call control_run method to do segmentation
      // For example:

      is_success = this->class_segmentationMain1.control_run(
          Image1D_current, dim_X, dim_Y, dim_Z,
          dialogRun1.channel_idx_selection, LandmarkList_current, idx_shape,
          dialogRun1.shape_para_delta,
          dialogRun1.shape_multiplier_thresholdRegionSize,
          dialogRun1.shape_multiplier_uThresholdRegionSize, name_currentWindow,
          dialogRun1.exemplar_maxMovement1, dialogRun1.exemplar_maxMovement2,
          fileName, dialogRun1.segmentationMode);

      // Then update the original window with the modified landmarks:
      _V3DPluginCallback2_currentCallback.setLandmark(v3dhandle_currentWindow,
                                                      LandmarkList_current);
      _V3DPluginCallback2_currentCallback.updateImageWindow(
          v3dhandle_currentWindow);
    }
    // if the segmentation is successful, display the results
    QString name_result = "Result";
    if (is_success) {
      // visualizationImage1D(this->class_segmentationMain1.Image1D_exemplar,
      // this->class_segmentationMain1.dim_X,
      // this->class_segmentationMain1.dim_Y,
      // this->class_segmentationMain1.dim_Z, 3,
      // _V3DPluginCallback2_currentCallback, "Exemplar");

      // visualization of result
      //   visualizationImage1D(
      //       this->class_segmentationMain1.Image1D_segmentationResult,
      //       this->class_segmentationMain1.dim_X,
      //       this->class_segmentationMain1.dim_Y,
      //       this->class_segmentationMain1.dim_Z, 3,
      //       _V3DPluginCallback2_currentCallback, name_result);
      //   // visualizationImage1D(this->class_segmentationMain1.Image1D_mask,
      //   // this->class_segmentationMain1.dim_X,
      //   // this->class_segmentationMain1.dim_Y,
      //   // this->class_segmentationMain1.dim_Z, 1,
      //   // _V3DPluginCallback2_currentCallback, "Mask");
      //   v3dhandleList v3dhandleList_current =
      //       _V3DPluginCallback2_currentCallback.getImageWindowList();
      //   V3DLONG count_v3dhandle = v3dhandleList_current.size();

      //   // QString name_exemplar = "Exemplar";
      //   for (V3DLONG i = 0; i < count_v3dhandle; i++) {
      //     if (_V3DPluginCallback2_currentCallback
      //             .getImageName(v3dhandleList_current[i])
      //             .contains(this->class_segmentationMain1.name_currentWindow))
      //             {
      //       _V3DPluginCallback2_currentCallback.setLandmark(
      //           v3dhandleList_current[i],
      //           this->class_segmentationMain1.LandmarkList_segmentationResult);
      //       _V3DPluginCallback2_currentCallback.updateImageWindow(
      //           v3dhandleList_current[i]);
      //       break;
      //     }
      //     // if
      //     //
      //     (_V3DPluginCallback2_currentCallback.getImageName(v3dhandleList_current[i]).contains(name_result))
      //     //{
      //     //_V3DPluginCallback2_currentCallback.setLandmark(v3dhandleList_current[i],
      //     // this->class_segmentationMain1.LandmarkList_exemplar);
      //     //_V3DPluginCallback2_currentCallback.updateImageWindow(v3dhandleList_current[i]);
      //     //}
      //     /*if
      //     (_V3DPluginCallback2_currentCallback.getImageName(v3dhandleList_current[i]).contains(name_exemplar))
      //     {
      //             _V3DPluginCallback2_currentCallback.setLandmark(v3dhandleList_current[i],
      //     this->class_segmentationMain1.LandmarkList_exemplar);
      //             _V3DPluginCallback2_currentCallback.updateImageWindow(v3dhandleList_current[i]);
      //     }*/
      //   }
      // temporary solution for Haru's request;
      ofstream ofstream_log;

      time_t rawtime;
      struct tm *timeinfo;
      time(&rawtime);
      timeinfo = localtime(&rawtime);
      char buffer[80];
      strftime(buffer, 80, "%d_%m_%Y_%I_%M_%S", timeinfo);

      stringstream tt;
      tt << "C:\\segmentationResult_" << buffer << ".csv";

      ofstream_log.open(tt.str().c_str());

      V3DLONG count_segments =
          this->class_segmentationMain1.possVct_segmentationResult.size();
      for (V3DLONG i = 0; i < count_segments; i++) {
        V3DLONG count_poss =
            this->class_segmentationMain1.possVct_segmentationResult[i].size();
        for (V3DLONG j = 0; j < (count_poss - 1); j++) {
          ofstream_log
              << this->class_segmentationMain1.possVct_segmentationResult[i][j]
              << ",";
        }
        ofstream_log << this->class_segmentationMain1
                            .possVct_segmentationResult[i][count_poss - 1]
                     << endl;
      }
      ofstream_log.close();

      // new code for saving a binary TIFF
      V3DLONG size_page = this->class_segmentationMain1.dim_X *
                          this->class_segmentationMain1.dim_Y *
                          this->class_segmentationMain1.dim_Z;

      // compute gradient of image
      unsigned char *gradientImage = new unsigned char[size_page];
      sobel3D(this->class_segmentationMain1.binarySegImage, gradientImage,
              this->class_segmentationMain1.dim_X,
              this->class_segmentationMain1.dim_Y,
              this->class_segmentationMain1.dim_Z);

      // overlay
      overlay2(_V3DPluginCallback2_currentCallback, _QWidget_parent,
               this->class_segmentationMain1.binarySegImage, gradientImage);

      // save binary image
      // QString savePath = QFileDialog::getSaveFileName(
      //     _QWidget_parent, "Save binary segmented image", "",
      //     "TIFF Files (*.tiff *.tif)");
      // if (!savePath.isEmpty()) {
      //   V3DLONG outSZ[4] = {this->class_segmentationMain1.dim_X,
      //                       this->class_segmentationMain1.dim_Y,
      //                       this->class_segmentationMain1.dim_Z, 1};
      //   simple_saveimage_wrapper(_V3DPluginCallback2_currentCallback,
      //                            savePath.toStdString().c_str(),
      //                            binarySegImage, outSZ, 1);
      //   v3d_msg("Binary segmented image saved.");
      // }

      // Automatically save binary segmented image to current directory.
      QString savePath = fileName + "_binary_segmentation.tif";
      V3DLONG outSZ[4] = {this->class_segmentationMain1.dim_X,
                          this->class_segmentationMain1.dim_Y,
                          this->class_segmentationMain1.dim_Z, 1};
      simple_saveimage_wrapper(
          _V3DPluginCallback2_currentCallback, savePath.toStdString().c_str(),
          this->class_segmentationMain1.binarySegImage, outSZ, 1);

      // save original image if we are usingt TeraFly
      if (isTeraFly) {
        savePath = fileName + "original_image.tif";
        simple_saveimage_wrapper(_V3DPluginCallback2_currentCallback,
                                 savePath.toStdString().c_str(),
                                 Image1D_current, outSZ, 1);
      }
      if (this->class_segmentationMain1.errorOccurred) {
        v3d_msg(QString("Some cells were not properly segmented, which may "
                        "give a poor probability model of cell shape. Check "
                        "debuggin log for details. Plugin files saved to %1.")
                    .arg(fileName));
      } else {
        v3d_msg(QString("Plugin files saved to %1.").arg(fileName));
      }
      delete[] this->class_segmentationMain1.binarySegImage;

      return true;
    } else {
      v3dhandleList v3dhandleList_current =
          _V3DPluginCallback2_currentCallback.getImageWindowList();
      V3DLONG count_v3dhandle = v3dhandleList_current.size();
      bool is_foundResultWindow = false;
      for (V3DLONG i = 0; i < count_v3dhandle; i++) {
        if (_V3DPluginCallback2_currentCallback
                .getImageName(v3dhandleList_current[i])
                .contains(name_result)) {
          LandmarkList LandmarkList_empty;
          _V3DPluginCallback2_currentCallback.setLandmark(
              v3dhandleList_current[i], LandmarkList_empty);
          _V3DPluginCallback2_currentCallback.updateImageWindow(
              v3dhandleList_current[i]);
          is_foundResultWindow = true;
          break;
        }
      }
      if (!is_foundResultWindow) {
        for (V3DLONG i = 0; i < count_v3dhandle; i++) {
          if (_V3DPluginCallback2_currentCallback
                  .getImageName(v3dhandleList_current[i])
                  .contains(this->class_segmentationMain1.name_currentWindow)) {
            LandmarkList LandmarkList_empty;
            _V3DPluginCallback2_currentCallback.setLandmark(
                v3dhandleList_current[i], LandmarkList_empty);
            _V3DPluginCallback2_currentCallback.updateImageWindow(
                v3dhandleList_current[i]);
            break;
          }
        }
      }
      v3d_msg("Warning: no exemplar defined, please re-select the exemplars!");
      return false;
    }
  }

  /**
   * @brief visualization of the image
   */
  void visualizationImage1D(
      unsigned char *Image1D_input, V3DLONG dim_X, V3DLONG dim_Y, V3DLONG dim_Z,
      int dim_C, V3DPluginCallback2 &_V3DPluginCallback2_currentCallback,
      QString string_windowName) {
    V3DLONG size_page = dim_X * dim_Y * dim_Z * dim_C;
    unsigned char *Image1D_tmp =
        class_segmentationMain::memory_allocate_uchar1D(size_page);
    for (V3DLONG i = 0; i < size_page; i++) {
      Image1D_tmp[i] = Image1D_input[i];
    }
    Image4DSimple Image4DSimple_temp;
    Image4DSimple_temp.setData(Image1D_tmp, dim_X, dim_Y, dim_Z, dim_C,
                               V3D_UINT8);

    v3dhandleList v3dhandleList_current =
        _V3DPluginCallback2_currentCallback.getImageWindowList();
    V3DLONG count_v3dhandle = v3dhandleList_current.size();
    bool is_found = false;
    for (V3DLONG i = 0; i < count_v3dhandle; i++) {
      if (_V3DPluginCallback2_currentCallback
              .getImageName(v3dhandleList_current[i])
              .contains(string_windowName)) {
        _V3DPluginCallback2_currentCallback.setImage(v3dhandleList_current[i],
                                                     &Image4DSimple_temp);
        _V3DPluginCallback2_currentCallback.updateImageWindow(
            v3dhandleList_current[i]);
        is_found = true;
        break;
      }
    }
    if (!is_found) {
      v3dhandle v3dhandle_main =
          _V3DPluginCallback2_currentCallback.newImageWindow();
      _V3DPluginCallback2_currentCallback.setImage(v3dhandle_main,
                                                   &Image4DSimple_temp);
      _V3DPluginCallback2_currentCallback.setImageName(v3dhandle_main,
                                                       string_windowName);
      _V3DPluginCallback2_currentCallback.updateImageWindow(v3dhandle_main);
    }
  }
#pragma endregion

#pragma region "overlay"
  /**
   * @brief overlay the binary segmentation onto the original image
   */
  void overlay(V3DPluginCallback2 &callback, QWidget *parent,
               unsigned char *binarySegImage) {
    v3dhandle curwin = callback.currentImageWindow();

    Image4DSimple *p4DImage = callback.getImage(curwin);
    unsigned char *data = p4DImage->getRawData();
    unsigned char *newData = new unsigned char[p4DImage->getTotalBytes() * 2];

    memcpy(newData, data, p4DImage->getTotalBytes());
    memcpy(newData + p4DImage->getTotalBytes(), binarySegImage,
           p4DImage->getTotalBytes());

    Image4DSimple *newImage = new Image4DSimple;
    newImage->setData(newData, p4DImage->getXDim(), p4DImage->getYDim(),
                      p4DImage->getZDim(), p4DImage->getCDim() * 2, V3D_UINT8);

    callback.setImage(curwin, newImage);
  }

  /**
   * @brief overlay two images onto the original image. Used to overlay
   * the segmentation and the gradient image.
   */
  void overlay2(V3DPluginCallback2 &callback, QWidget *parent,
                unsigned char *binarySegImage, unsigned char *gradientImage) {
    v3dhandle curwin = callback.currentImageWindow();
    v3dhandle newwin = callback.newImageWindow();

    Image4DSimple *p4DImage = callback.getImage(curwin);
    unsigned char *data = p4DImage->getRawData();
    unsigned char *newData = new unsigned char[p4DImage->getTotalBytes() * 3];

    memcpy(newData, data, p4DImage->getTotalBytes());
    memcpy(newData + p4DImage->getTotalBytes(), binarySegImage,
           p4DImage->getTotalBytes());
    memcpy(newData + p4DImage->getTotalBytes() * 2, gradientImage,
           p4DImage->getTotalBytes());

    Image4DSimple *newImage = new Image4DSimple;
    newImage->setData(newData, p4DImage->getXDim(), p4DImage->getYDim(),
                      p4DImage->getZDim(), p4DImage->getCDim() * 3, V3D_UINT8);

    callback.setImage(newwin, newImage);
  }

  /**
   * @brief Compute 3D Sobel filter
   */
  void sobel3D(unsigned char *data, unsigned char *out, V3DLONG dim_X,
               V3DLONG dim_Y, V3DLONG dim_Z) {
    std::vector<cv::Mat> gradX(dim_Z), gradY(dim_Z);

    // Compute Sobel for x and y for each slice
    for (int k = 0; k < dim_Z; k++) {
      cv::Mat slice(dim_Y, dim_X, CV_8U, data + k * dim_X * dim_Y);
      cv::Mat gx, gy;
      cv::Sobel(slice, gx, CV_16S, 1, 0, 3);
      cv::Sobel(slice, gy, CV_16S, 0, 1, 3);
      gradX[k] = gx.clone();
      gradY[k] = gy.clone();
    }

    // Compute gradient magnitude for each voxel
    for (int k = 0; k < dim_Z; k++) {
      for (int j = 0; j < dim_Y; j++) {
        for (int i = 0; i < dim_X; i++) {
          // Get Sobel derivatives in x and y
          short sx = gradX[k].at<short>(j, i);
          short sy = gradY[k].at<short>(j, i);

          // Compute derivative in z using central difference
          int idx = k * dim_X * dim_Y + j * dim_X + i;
          int center = data[idx];
          int prev =
              (k == 0) ? center : data[(k - 1) * dim_X * dim_Y + j * dim_X + i];
          int next = (k == dim_Z - 1)
                         ? center
                         : data[(k + 1) * dim_X * dim_Y + j * dim_X + i];
          short sz = static_cast<short>((next - prev) / 2);

          // Gradient magnitude (using Euclidean norm)
          int mag = static_cast<int>(std::sqrt(sx * sx + sy * sy + sz * sz));
          if (mag > 255) mag = 255;
          out[idx] = static_cast<unsigned char>(mag);
        }
      }
    }
  }
};
#pragma endregion
#endif
