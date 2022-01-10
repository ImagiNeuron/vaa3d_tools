#ifndef __CLASSIFICATION_H__
#define __CLASSIFICATION_H__

#include "../../../MK/DeepNeuron/Deep_Neuron_plugin.h"
#include "prediction_caffe_plugin.h"

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <caffe/caffe.hpp>

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace caffe; // NOLINT(build/namespaces)
using std::string;

/* Pair (label, confidence) representing a prediction. */
typedef std::pair<string, float> Prediction;
class Classifier {
public:
    Classifier(const string& model_file,
               const string& trained_file,
               const string& mean_file);
    std::vector<std::vector<float> > Predict(const std::vector<cv::Mat>& imgs);
    std::vector<float> Predict_3D(const std::vector<cv::Mat>& imgs);

    std::vector< float > PredictBatch(const vector< cv::Mat > imgs);

    std::vector<std::vector<float> > extractFeature_siamese(const std::vector<cv::Mat>& imgs);


private:
    void SetMean(const string& mean_file);
    void WrapInputLayer(std::vector<cv::Mat>* input_channels, int n);
    void Preprocess(const cv::Mat& img,
                    std::vector<cv::Mat>* input_channels);
    void WrapBatchInputLayer(std::vector<std::vector<cv::Mat> > *input_batch);
    void PreprocessBatch(const vector<cv::Mat> imgs,
                                          std::vector< std::vector<cv::Mat> >* input_batch);

private:
    boost::shared_ptr<Net<float> > net_; // Add boost qualifier to avoid ambiguity in C++ 11, MK, Dec, 2017
    cv::Size input_geometry_;
    int num_channels_;
    cv::Mat mean_;
    int batch_size_ ;
};

QStringList importSeriesFileList_addnumbersort(const QString & curFilePath);
NeuronTree DL_eliminate_swc(NeuronTree nt,QList <ImageMarker> marklist);
NeuronTree remove_swc(NeuronTree nt,double length);
std::vector<std::vector<float> > batch_detection(unsigned char * & data1d,Classifier classifier, int N, int M, int P, int Sxy);
std::vector<std::vector<float> > batch_detection_ref(unsigned char * & data1d,unsigned char * & data1d_ref,Classifier classifier, int N, int M, int P, int Sxy);
QList <ImageMarker> batch_deletion(unsigned char * & data1d,Classifier classifier, QList <ImageMarker> input_marker, int N, int M, int P);
//void connect_swc(NeuronTree nt,QList<NeuronSWC>& newNeuron, double disThr,double angThr);
void gaussian_filter(unsigned char* data1d,V3DLONG *in_sz, unsigned int Wx,unsigned int Wy,unsigned int Wz,unsigned int c,double sigma,float* &outimg);
#endif

