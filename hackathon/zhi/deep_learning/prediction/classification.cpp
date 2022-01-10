#include "classification.h"

#include <caffe/caffe.hpp>

#include "../../../MK/DeepNeuron/Deep_Neuron_plugin.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#if  defined(Q_OS_MAC)
    #include <opencv2/imgcodecs/imgcodecs.hpp>
#endif

#include <algorithm>
#include <iosfwd>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <cstdlib>
#include <string>
#include "../../../../released_plugins/v3d_plugins/istitch/y_imglib.h"

using namespace std;

//#define PI 3.14159265359
//#define NTDIS(a,b) (sqrt(((a).x-(b).x)*((a).x-(b).x)+((a).y-(b).y)*((a).y-(b).y)+((a).z-(b).z)*((a).z-(b).z)))
//#define NTDOT(a,b) ((a).x*(b).x+(a).y*(b).y+(a).z*(b).z)
//#define angle(a,b,c) (acos((((b).x-(a).x)*((c).x-(a).x)+((b).y-(a).y)*((c).y-(a).y)+((b).z-(a).z)*((c).z-(a).z))/(NTDIS(a,b)*NTDIS(a,c)))*180.0/PI)

#ifndef MAX_DOUBLE
#define MAX_DOUBLE 1.79768e+308        //actual: 1.79769e+308
#endif

template <class T> T pow2(T a)
{
    return a*a;

}


Classifier::Classifier(const string& model_file,
                       const string& trained_file,
                       const string& mean_file) {

    /* Load the network. */
    batch_size_=1000;
    new Net<float>(model_file, TEST);
    net_.reset(new Net<float>(model_file, TEST));
    net_->CopyTrainedLayersFrom(trained_file);
    CHECK_EQ(net_->num_inputs(), 1) << "Network should have exactly one input.";
    CHECK_EQ(net_->num_outputs(), 1) << "Network should have exactly one output.";
    Blob<float>* input_layer = net_->input_blobs()[0];
    num_channels_ = input_layer->channels();
    CHECK(num_channels_ == 3 || num_channels_ == 1)
            << "Input layer should have 1 or 3 channels.";
    input_geometry_ = cv::Size(input_layer->width(), input_layer->height());
    /* Load the binaryproto mean file. */
    if(!mean_file.empty()) SetMean(mean_file);
    Caffe::set_mode(Caffe::GPU);

#if  defined(Q_OS_MAC)
    Caffe::set_mode(Caffe::CPU);
#endif

#if defined(Q_OS_WIN) // MK, 2018, Jan. Automatic GPU/CPU selection on Windows
    char* CUDA = getenv("CUDA_PATH");
    if (CUDA == NULL) Caffe::set_mode(Caffe::CPU);
    else Caffe::set_mode(Caffe::GPU);
#endif
    /* Load labels. */
}


/* Load the mean file in binaryproto format. */
void Classifier::SetMean(const string& mean_file) {
    BlobProto blob_proto;
    ReadProtoFromBinaryFileOrDie(mean_file.c_str(), &blob_proto);
    /* Convert from BlobProto to Blob<float> */
    Blob<float> mean_blob;
    mean_blob.FromProto(blob_proto);
    CHECK_EQ(mean_blob.channels(), num_channels_)
            << "Number of channels of mean file doesn't match input layer.";
    /* The format of the mean file is planar 32-bit float BGR or grayscale. */
    std::vector<cv::Mat> channels;
    float* data = mean_blob.mutable_cpu_data();
    for (int i = 0; i < num_channels_; ++i) {
        /* Extract an individual channel. */
        cv::Mat channel(mean_blob.height(), mean_blob.width(), CV_32FC1, data);
        channels.push_back(channel);
        data += mean_blob.height() * mean_blob.width();
    }
    /* Merge the separate channels into a single image. */
    cv::Mat mean;
    cv::merge(channels, mean);
    mean_ = mean;
//    cv::imwrite("/local4/Data/IVSCC_test/comparison/Caffe_testing_3nd/train/test_cpp/mean.tif",mean_);
//    /* Compute the global mean pixel value and create a mean image
//* filled with this value. */
//    cv::Scalar channel_mean = cv::mean(mean);
//    mean_ = cv::Mat(input_geometry_, mean.type(), channel_mean);
//   // cv::imwrite("/local4/Data/IVSCC_test/comparison/Caffe_testing_3nd/train/test_cpp/mean.tif",mean_);
}
std::vector<std::vector<float> > Classifier::Predict(const std::vector<cv::Mat>& imgs) {
    Blob<float>* input_layer = net_->input_blobs()[0];
    input_layer->Reshape(imgs.size(), num_channels_,
                         input_geometry_.height, input_geometry_.width);
    /* Forward dimension change to all layers. */
    net_->Reshape();
    for (int i = 0; i < imgs.size(); ++i) {
        std::vector<cv::Mat> input_channels;
        WrapInputLayer(&input_channels, i);
        Preprocess(imgs[i], &input_channels);
    }
    net_->ForwardPrefilled();
    std::vector<std::vector<float> > outputs;
    Blob<float>* output_layer = net_->output_blobs()[0];
    for (int i = 0; i < output_layer->num(); ++i) {
        const float* begin = output_layer->cpu_data() + i * output_layer->channels();
        const float* end = begin + output_layer->channels();
        /* Copy the output layer to a std::vector */
        outputs.push_back(std::vector<float>(begin, end));
    }
    return outputs;
}

std::vector<float> Classifier::Predict_3D(const std::vector<cv::Mat>& imgs) {
    Blob<float>* input_layer = net_->input_blobs()[0];
    input_layer->Reshape(1, num_channels_,imgs.size(),
                         input_geometry_.height, input_geometry_.width);
    /* Forward dimension change to all layers. */
    net_->Reshape();
    for (int i = 0; i < imgs.size(); ++i) {
        std::vector<cv::Mat> input_channels;
        WrapInputLayer(&input_channels, i);
        Preprocess(imgs[i], &input_channels);
    }
    net_->ForwardPrefilled();
    Blob<float>* output_layer = net_->output_blobs()[0];
    const float* begin = output_layer->cpu_data()+imgs.size()*input_geometry_.height* input_geometry_.width;
    const float* end = begin + imgs.size()*input_geometry_.height* input_geometry_.width;
    return std::vector<float>(begin, end);
}

std::vector<std::vector<float> > Classifier::extractFeature_siamese(const std::vector<cv::Mat>& imgs) {
    Blob<float>* input_layer = net_->input_blobs()[0];
    input_layer->Reshape(imgs.size(), num_channels_,
                         input_geometry_.height, input_geometry_.width);
    /* Forward dimension change to all layers. */
    net_->Reshape();
    for (int i = 0; i < imgs.size(); ++i) {
        std::vector<cv::Mat> input_channels;
        WrapInputLayer(&input_channels, i);
        Preprocess(imgs[i], &input_channels);
    }
    net_->ForwardPrefilled();

    std::vector<std::vector<float> > outputs;
    Blob<float>* output_layer = net_->output_blobs()[0];
    for (int i = 0; i < output_layer->num(); ++i) {
        const float* begin = output_layer->cpu_data() + i * output_layer->channels();
        const float* end = begin + output_layer->channels();
        /* Copy the output layer to a std::vector */
        outputs.push_back(std::vector<float>(begin, end));
    }
    return outputs;

}
/* Wrap the input layer of the network in separate cv::Mat objects
* (one per channel). This way we save one memcpy operation and we
* don't need to rely on cudaMemcpy2D. The last preprocessing
* operation will write the separate channels directly to the input
* layer. */
void Classifier::WrapInputLayer(std::vector<cv::Mat>* input_channels, int n) {
    Blob<float>* input_layer = net_->input_blobs()[0];
    int width = input_layer->width();
    int height = input_layer->height();
    int channels = input_layer->channels();
    float* input_data = input_layer->mutable_cpu_data() + n * width * height * channels;
    for (int i = 0; i < channels; ++i) {
        cv::Mat channel(height, width, CV_32FC1, input_data);
        input_channels->push_back(channel);
        input_data += width * height;
    }
}
void Classifier::Preprocess(const cv::Mat& img,
                            std::vector<cv::Mat>* input_channels) {
    /* Convert the input image to the input image format of the network. */
    cv::Mat sample;
    if (img.channels() == 3 && num_channels_ == 1)
        cv::cvtColor(img, sample, cv::COLOR_BGR2GRAY);
    else if (img.channels() == 4 && num_channels_ == 1)
        cv::cvtColor(img, sample, cv::COLOR_BGRA2GRAY);
    else if (img.channels() == 4 && num_channels_ == 3)
        cv::cvtColor(img, sample, cv::COLOR_BGRA2BGR);
    else if (img.channels() == 1 && num_channels_ == 3)
        cv::cvtColor(img, sample, cv::COLOR_GRAY2BGR);
    else
        sample = img;
    cv::Mat sample_resized;
    if (sample.size() != input_geometry_)
        cv::resize(sample, sample_resized, input_geometry_);
    else
        sample_resized = sample;
    cv::Mat sample_float;
    if (num_channels_ == 3)
        sample_resized.convertTo(sample_float, CV_32FC3);
    else
        sample_resized.convertTo(sample_float, CV_32FC1);
    cv::Mat sample_normalized;
    if(!mean_.empty())
        cv::subtract(sample_float, mean_, sample_normalized);
    else
        sample_normalized = sample_float;// * 0.00390625;
    /* This operation will write the separate BGR planes directly to the
* input layer of the network because it is wrapped by the cv::Mat
* objects in input_channels. */
    cv::split(sample_normalized, *input_channels);
    /*
CHECK(reinterpret_cast<float*>(input_channels->at(0).data)
== net_->input_blobs()[0]->cpu_data())
<< "Input channels are not wrapping the input layer of the network.";
*/
}


std::vector< float > Classifier::PredictBatch(const vector< cv::Mat > imgs) {
    Blob<float>* input_layer = net_->input_blobs()[0];
    input_layer->Reshape(batch_size_, num_channels_,
                         input_geometry_.height,
                         input_geometry_.width);
    /* Forward dimension change to all layers. */
    net_->Reshape();
    std::vector< std::vector<cv::Mat> > input_batch;
    WrapBatchInputLayer(&input_batch);
    PreprocessBatch(imgs, &input_batch);
    net_->ForwardPrefilled();
    /* Copy the output layer to a std::vector */
    Blob<float>* output_layer = net_->output_blobs()[0];
    const float* begin = output_layer->cpu_data();
    const float* end = begin + output_layer->channels()*imgs.size();
    return std::vector<float>(begin, end);
}
void Classifier::WrapBatchInputLayer(std::vector<std::vector<cv::Mat> > *input_batch){
    Blob<float>* input_layer = net_->input_blobs()[0];
    int width = input_layer->width();
    int height = input_layer->height();
    int num = input_layer->num();
    float* input_data = input_layer->mutable_cpu_data();
    for ( int j = 0; j < num; j++){
        vector<cv::Mat> input_channels;
        for (int i = 0; i < input_layer->channels(); ++i){
            cv::Mat channel(height, width, CV_32FC1, input_data);
            input_channels.push_back(channel);
            input_data += width * height;
        }
        input_batch -> push_back(vector<cv::Mat>(input_channels));
    }
    cv::imshow("bla", input_batch->at(1).at(0));
    cv::waitKey(1);
}
void Classifier::PreprocessBatch(const vector<cv::Mat> imgs,
                                      std::vector< std::vector<cv::Mat> >* input_batch){
    for (int i = 0 ; i < imgs.size(); i++){
        cv::Mat img = imgs[i];
        std::vector<cv::Mat> *input_channels = &(input_batch->at(i));
        /* Convert the input image to the input image format of the network. */
        cv::Mat sample;
        if (img.channels() == 3 && num_channels_ == 1)
            cv::cvtColor(img, sample, CV_BGR2GRAY);
        else if (img.channels() == 4 && num_channels_ == 1)
            cv::cvtColor(img, sample, CV_BGRA2GRAY);
        else if (img.channels() == 4 && num_channels_ == 3)
            cv::cvtColor(img, sample, CV_BGRA2BGR);
        else if (img.channels() == 1 && num_channels_ == 3)
            cv::cvtColor(img, sample, CV_GRAY2BGR);
        else
            sample = img;
        cv::Mat sample_resized;
        if (sample.size() != input_geometry_)
            cv::resize(sample, sample_resized, input_geometry_);
        else
            sample_resized = sample;
        cv::Mat sample_float;
        if (num_channels_ == 3)
            sample_resized.convertTo(sample_float, CV_32FC3);
        else
            sample_resized.convertTo(sample_float, CV_32FC1);
        cv::Mat sample_normalized;
        cv::subtract(sample_float, mean_, sample_normalized);
        /* This operation will write the separate BGR planes directly to the
* input layer of the network because it is wrapped by the cv::Mat
* objects in input_channels. */
        cv::split(sample_normalized, *input_channels);
        // CHECK(reinterpret_cast<float*>(input_channels->at(0).data)
        // == net_->input_blobs()[0]->cpu_data())
        // << "Input channels are not wrapping the input layer of the network.";
    }
}

QStringList importSeriesFileList_addnumbersort(const QString & curFilePath)
{
    QStringList myList;
    myList.clear();

    // get the image files namelist in the directory
    QStringList imgSuffix;
    imgSuffix<<"*.tif"<<"*.raw"<<"*.v3draw"<<"*.lsm"
            <<"*.TIF"<<"*.RAW"<<"*.V3DRAW"<<"*.LSM"
           <<"*.jpg"<<"*.JPG";

    QDir dir(curFilePath);
    if (!dir.exists())
    {
        qWarning("Cannot find the directory");
        return myList;
    }

    foreach (QString file, dir.entryList(imgSuffix, QDir::Files, QDir::Name))
    {
        myList += QFileInfo(dir, file).absoluteFilePath();
    }

    // print filenames
    foreach (QString qs, myList)  qDebug() << qs;

    return myList;
}

NeuronTree DL_eliminate_swc(NeuronTree nt,QList <ImageMarker> marklist)
{
    NeuronTree nt_prunned;
    QList <NeuronSWC> listNeuron;
    QHash <int, int>  hashNeuron;
    listNeuron.clear();
    hashNeuron.clear();
    NeuronSWC S;

    for (V3DLONG i=0;i<nt.listNeuron.size();i++)
    {
        bool flag = false;
        for(V3DLONG j=0; j<marklist.size();j++)
        {
            double dis = sqrt(pow2(nt.listNeuron.at(i).x - marklist.at(j).x) + pow2(nt.listNeuron.at(i).y - marklist.at(j).y) + pow2(nt.listNeuron.at(i).z - marklist.at(j).z));
            if(dis < 1.0)
            {
                flag = true;
                break;
            }
        }
        if(!flag)
        {
            NeuronSWC curr = nt.listNeuron.at(i);
            S.n 	= curr.n;
            S.type 	= curr.type;
            S.x 	= curr.x;
            S.y 	= curr.y;
            S.z 	= curr.z;
            S.r 	= curr.r;
            S.pn 	= curr.pn;
            listNeuron.append(S);
            hashNeuron.insert(S.n, listNeuron.size()-1);
        }
    }

    nt_prunned.n = -1;
    nt_prunned.on = true;
    nt_prunned.listNeuron = listNeuron;
    nt_prunned.hashNeuron = hashNeuron;

    return nt_prunned;
}

NeuronTree remove_swc(NeuronTree nt,double length)
{
    QList<NeuronSWC> list = nt.listNeuron;
    V3DLONG *flag = new V3DLONG[list.size()];
    V3DLONG counter = 0;
    V3DLONG root_ID = 0;
    for (V3DLONG i=1;i<list.size();i++)
    {
        if(list.at(i).parent > 0)
        {
            counter++;
        }else
        {
            for(V3DLONG j=root_ID; j<=root_ID+counter;j++)
            {
               if(counter <= length)
                   flag[j] = -1;
               else
                   flag[j] = 1;
            }
            counter = 0;
            root_ID = i;
        }
    }

    NeuronTree nt_kept;
    QList <NeuronSWC> listNeuron;
    QHash <int, int>  hashNeuron;
    listNeuron.clear();
    hashNeuron.clear();

    //set node

    NeuronSWC S;
    for (int i=0;i<list.size();i++)
    {
        if(flag[i] == 1)
        {
             NeuronSWC curr = list.at(i);
             S.n 	= curr.n;
             S.type 	= curr.type;
             S.x 	= curr.x;
             S.y 	= curr.y;
             S.z 	= curr.z;
             S.r 	= curr.r;
             S.pn 	= curr.pn;
             listNeuron.append(S);
             hashNeuron.insert(S.n, listNeuron.size()-1);
        }

   }
    nt_kept.n = -1;
    nt_kept.on = true;
    nt_kept.listNeuron = listNeuron;
    nt_kept.hashNeuron = hashNeuron;
    return nt_kept;
}

std::vector<std::vector<float> > batch_detection(unsigned char * & data1d,Classifier classifier, int N, int M, int P, int Sxy)
{
    std::vector<cv::Mat> imgs;
    int Wx = 30, Wy = 30, Wz = 1;
    int Sz = (int)Sxy;

    V3DLONG num_patches = 0;
    std::vector<std::vector<float> > outputs_overall;
    std::vector<std::vector<float> > outputs;
    for(V3DLONG iz = 0; iz < P; iz = iz+Sz)
    {
        for(V3DLONG iy = Sxy; iy < M; iy = iy+Sxy)
        {
            for(V3DLONG ix = Sxy; ix < N; ix = ix+Sxy)
            {
                    V3DLONG xb = ix-Wx; if(xb<0) xb = 0;if(xb>=N-1) xb = N-1;
                    V3DLONG xe = ix+Wx; if(xe>=N-1) xe = N-1;
                    V3DLONG yb = iy-Wy; if(yb<0) yb = 0;if(yb>=M-1) yb = M-1;
                    V3DLONG ye = iy+Wy; if(ye>=M-1) ye = M-1;
                    V3DLONG zb = iz-Wz; if(zb<0) zb = 0;if(zb>=P-1) zb = P-1;
                    V3DLONG ze = iz+Wz; if(ze>=P-1) ze = P-1;

                    V3DLONG im_cropped_sz[4];
                    im_cropped_sz[0] = xe - xb + 1;
                    im_cropped_sz[1] = ye - yb + 1;
                    im_cropped_sz[2] = 1;
                    im_cropped_sz[3] = 1;

                    unsigned char *im_cropped = 0;

                    V3DLONG pagesz = im_cropped_sz[0]* im_cropped_sz[1]* im_cropped_sz[2]*im_cropped_sz[3];
                    try {im_cropped = new unsigned char [pagesz];}
                    catch(...)  {v3d_msg("cannot allocate memory for im_cropped."); return outputs_overall;}
                    memset(im_cropped, 0, sizeof(unsigned char)*pagesz);

                    for(V3DLONG iiz = zb; iiz <= ze; iiz++)
                    {
                        V3DLONG offsetk = iiz*M*N;
                        V3DLONG j = 0;
                        for(V3DLONG iiy = yb; iiy <= ye; iiy++)
                        {
                            V3DLONG offsetj = iiy*N;
                            for(V3DLONG iix = xb; iix <= xe; iix++)
                            {
                                if(data1d[offsetk + offsetj + iix] >= im_cropped[j])
                                    im_cropped[j] = data1d[offsetk + offsetj + iix];
                                j++;
                            }
                        }
                    }
                    cv::Mat img(im_cropped_sz[1], im_cropped_sz[0], CV_8UC1, im_cropped);
                    imgs.push_back(img);

                    if(num_patches >= 5000)
                    {
                        outputs = classifier.Predict(imgs);
                        for(V3DLONG d = 0; d<outputs.size();d++)
                            outputs_overall.push_back(outputs[d]);
                        outputs.clear();
                        imgs.clear();
                        num_patches = 0;
                    }else
                        num_patches++;
            }
        }
    }

    if(imgs.size()>0)
    {
        outputs = classifier.Predict(imgs);
        for(V3DLONG d = 0; d<outputs.size();d++)
            outputs_overall.push_back(outputs[d]);
    }

    imgs.clear();
    outputs.clear();

    return outputs_overall;
}

std::vector<std::vector<float> > batch_detection_ref(unsigned char * & data1d,unsigned char * & data1d_ref,Classifier classifier, int N, int M, int P, int Sxy)
{
    std::vector<cv::Mat> imgs;
    int Wx = 30, Wy = 30, Wz = 1;
    int Sz = (int)Sxy;

    V3DLONG num_patches = 0;
    std::vector<std::vector<float> > outputs_overall;
    std::vector<std::vector<float> > outputs;
    for(V3DLONG iz = 0; iz < P; iz = iz+Sz)
    {
        for(V3DLONG iy = Sxy; iy < M; iy = iy+Sxy)
        {
            for(V3DLONG ix = Sxy; ix < N; ix = ix+Sxy)
            {
                if(data1d_ref[iy*N +ix]>0)
                {
                    V3DLONG xb = ix-Wx; if(xb<0) xb = 0;if(xb>=N-1) xb = N-1;
                    V3DLONG xe = ix+Wx; if(xe>=N-1) xe = N-1;
                    V3DLONG yb = iy-Wy; if(yb<0) yb = 0;if(yb>=M-1) yb = M-1;
                    V3DLONG ye = iy+Wy; if(ye>=M-1) ye = M-1;
                    V3DLONG zb = iz-Wz; if(zb<0) zb = 0;if(zb>=P-1) zb = P-1;
                    V3DLONG ze = iz+Wz; if(ze>=P-1) ze = P-1;

                    V3DLONG im_cropped_sz[4];
                    im_cropped_sz[0] = xe - xb + 1;
                    im_cropped_sz[1] = ye - yb + 1;
                    im_cropped_sz[2] = 1;
                    im_cropped_sz[3] = 1;

                    unsigned char *im_cropped = 0;

                    V3DLONG pagesz = im_cropped_sz[0]* im_cropped_sz[1]* im_cropped_sz[2]*im_cropped_sz[3];
                    try {im_cropped = new unsigned char [pagesz];}
                    catch(...)  {v3d_msg("cannot allocate memory for im_cropped."); return outputs_overall;}
                    memset(im_cropped, 0, sizeof(unsigned char)*pagesz);

                    for(V3DLONG iiz = zb; iiz <= ze; iiz++)
                    {
                        V3DLONG offsetk = iiz*M*N;
                        V3DLONG j = 0;
                        for(V3DLONG iiy = yb; iiy <= ye; iiy++)
                        {
                            V3DLONG offsetj = iiy*N;
                            for(V3DLONG iix = xb; iix <= xe; iix++)
                            {
                                if(data1d[offsetk + offsetj + iix] >= im_cropped[j])
                                    im_cropped[j] = data1d[offsetk + offsetj + iix];
                                j++;
                            }
                        }
                    }
                    cv::Mat img(im_cropped_sz[1], im_cropped_sz[0], CV_8UC1, im_cropped);
                    imgs.push_back(img);

                    if(num_patches >= 5000)
                    {
                        outputs = classifier.Predict(imgs);
                        for(V3DLONG d = 0; d<outputs.size();d++)
                            outputs_overall.push_back(outputs[d]);
                        outputs.clear();
                        imgs.clear();
                        num_patches = 0;
                    }else
                        num_patches++;
                }
            }
        }
    }

    if(imgs.size()>0)
    {
        outputs = classifier.Predict(imgs);
        for(V3DLONG d = 0; d<outputs.size();d++)
            outputs_overall.push_back(outputs[d]);
    }

    imgs.clear();
    outputs.clear();

    return outputs_overall;
}

QList <ImageMarker> batch_deletion(unsigned char * & data1d,Classifier classifier, QList <ImageMarker> input_markerlist, int N, int M, int P)
{
    QList <ImageMarker> output_marklist;
    unsigned int Wx=30, Wy=30, Wz=15;
    std::vector<cv::Mat> imgs;
    V3DLONG num_patches = 0;
    std::vector<std::vector<float> > outputs_overall;
    std::vector<std::vector<float> > outputs;
    for (V3DLONG i=0;i<input_markerlist.size();i++)
    {
        V3DLONG tmpx = input_markerlist.at(i).x;
        V3DLONG tmpy = input_markerlist.at(i).y;
        V3DLONG tmpz = input_markerlist.at(i).z;

        V3DLONG xb = tmpx-1-Wx; if(xb<0) xb = 0;if(xb>=N-1) xb = N-1;
        V3DLONG xe = tmpx-1+Wx; if(xe>=N-1) xe = N-1;
        V3DLONG yb = tmpy-1-Wy; if(yb<0) yb = 0;if(yb>=M-1) yb = M-1;
        V3DLONG ye = tmpy-1+Wy; if(ye>=M-1) ye = M-1;
        V3DLONG zb = tmpz-1-Wz; if(zb<0) zb = 0;if(zb>=P-1) zb = P-1;
        V3DLONG ze = tmpz-1+Wz; if(ze>=P-1) ze = P-1;

        V3DLONG im_cropped_sz[4];
        im_cropped_sz[0] = xe - xb + 1;
        im_cropped_sz[1] = ye - yb + 1;
        im_cropped_sz[2] = 1;
        im_cropped_sz[3] = 1;

        unsigned char *im_cropped = 0;

        V3DLONG pagesz = im_cropped_sz[0]* im_cropped_sz[1]* im_cropped_sz[2]*im_cropped_sz[3];
        try {im_cropped = new unsigned char [pagesz];}
        catch(...)  {v3d_msg("cannot allocate memory for im_cropped."); return output_marklist;}
        memset(im_cropped, 0, sizeof(unsigned char)*pagesz);

        for(V3DLONG iz = zb; iz <= ze; iz++)
        {
            V3DLONG offsetk = iz*M*N;
            V3DLONG j = 0;
            for(V3DLONG iy = yb; iy <= ye; iy++)
            {
                V3DLONG offsetj = iy*N;
                for(V3DLONG ix = xb; ix <= xe; ix++)
                {
                    if(data1d[offsetk + offsetj + ix] >= im_cropped[j])
                        im_cropped[j] = data1d[offsetk + offsetj + ix];
                    j++;
                }
            }
        }

        cv::Mat img(im_cropped_sz[1], im_cropped_sz[0], CV_8UC1, im_cropped);
        imgs.push_back(img);

        if(num_patches >=10000)
        {
            outputs = classifier.Predict(imgs);
            for(V3DLONG d = 0; d<outputs.size();d++)
                outputs_overall.push_back(outputs[d]);
            outputs.clear();
            imgs.clear();
            num_patches = 0;
        }else
            num_patches++;
    }

    if(imgs.size()>0)
    {
        outputs = classifier.Predict(imgs);
        for(V3DLONG d = 0; d<outputs.size();d++)
            outputs_overall.push_back(outputs[d]);
    }

    for (V3DLONG j=0;j<input_markerlist.size();j++)
    {
        std::vector<float> output = outputs_overall[j];
        if(output.at(0) < output.at(1))
        {
            input_markerlist[j].radius = output.at(1);
            output_marklist.push_back(input_markerlist.at(j));
        }
    }
    imgs.clear();
    outputs_overall.clear();
    outputs.clear();
    return output_marklist;
}

void gaussian_filter(unsigned char* data1d,V3DLONG *in_sz, unsigned int Wx,unsigned int Wy,unsigned int Wz,unsigned int c,double sigma,float* &outimg)
{
    if (!data1d || !in_sz || in_sz[0]<=0 || in_sz[1]<=0 || in_sz[2]<=0 || in_sz[3]<=0 || outimg)
    {
        v3d_msg("Invalid parameters to gaussian_filter().", 0);
        return;
    }

    if (outimg)
    {
        v3d_msg("Warning: you have supplied an non-empty output image pointer. This program will force to free it now. But you may want to double check.");
        delete []outimg;
        outimg = 0;
    }

     // for filter kernel
     double sigma_s2 = 0.5/(sigma*sigma); // 1/(2*sigma*sigma)
     double pi_sigma = 1.0/(sqrt(2*3.1415926)*sigma); // 1.0/(sqrt(2*pi)*sigma)

     float min_val = INF, max_val = 0;

     V3DLONG N = in_sz[0];
     V3DLONG M = in_sz[1];
     V3DLONG P = in_sz[2];
     V3DLONG sc = in_sz[3];
     V3DLONG pagesz = N*M*P;

     //filtering
     V3DLONG offset_init = (c-1)*pagesz;

     //declare temporary pointer
     float *pImage = new float [pagesz];
     if (!pImage)
     {
          printf("Fail to allocate memory.\n");
          return;
     }
     else
     {
          for(V3DLONG i=0; i<pagesz; i++)
               pImage[i] = data1d[i + offset_init];  //first channel data (red in V3D, green in ImageJ)
     }
       //Filtering
     //
     //   Filtering along x
     if(N<2)
     {
          //do nothing
     }
     else
     {
          //create Gaussian kernel
          float  *WeightsX = 0;
          WeightsX = new float [Wx];
          if (!WeightsX)
               return;

          float Half = (float)(Wx-1)/2.0;

          // Gaussian filter equation:
          // http://en.wikipedia.org/wiki/Gaussian_blur
       //   for (unsigned int Weight = 0; Weight < Half; ++Weight)
       //   {
       //        const float  x = Half* float (Weight) / float (Half);
      //         WeightsX[(int)Half - Weight] = WeightsX[(int)Half + Weight] = pi_sigma * exp(-x * x *sigma_s2); // Corresponding symmetric WeightsX
      //    }

          for (unsigned int Weight = 0; Weight <= Half; ++Weight)
          {
              const float  x = float(Weight)-Half;
              WeightsX[Weight] = WeightsX[Wx-Weight-1] = pi_sigma * exp(-(x * x *sigma_s2)); // Corresponding symmetric WeightsX
          }


          double k = 0.;
          for (unsigned int Weight = 0; Weight < Wx; ++Weight)
               k += WeightsX[Weight];

          for (unsigned int Weight = 0; Weight < Wx; ++Weight)
               WeightsX[Weight] /= k;

         printf("\n x dierction");

         for (unsigned int Weight = 0; Weight < Wx; ++Weight)
             printf("/n%f",WeightsX[Weight]);

          //   Allocate 1-D extension array
          float  *extension_bufferX = 0;
          extension_bufferX = new float [N + (Wx<<1)];

          unsigned int offset = Wx>>1;

          //	along x
          const float  *extStop = extension_bufferX + N + offset;

          for(V3DLONG iz = 0; iz < P; iz++)
          {
               for(V3DLONG iy = 0; iy < M; iy++)
               {
                    float  *extIter = extension_bufferX + Wx;
                    for(V3DLONG ix = 0; ix < N; ix++)
                    {
                         *(extIter++) = pImage[iz*M*N + iy*N + ix];
                    }

                    //   Extend image
                    const float  *const stop_line = extension_bufferX - 1;
                    float  *extLeft = extension_bufferX + Wx - 1;
                    const float  *arrLeft = extLeft + 2;
                    float  *extRight = extLeft + N + 1;
                    const float  *arrRight = extRight - 2;

                    while (extLeft > stop_line)
                    {
                         *(extLeft--) = *(arrLeft++);
                         *(extRight++) = *(arrRight--);

                    }

                    //	Filtering
                    extIter = extension_bufferX + offset;

                    float  *resIter = &(pImage[iz*M*N + iy*N]);

                    while (extIter < extStop)
                    {
                         double sum = 0.;
                         const float  *weightIter = WeightsX;
                         const float  *const End = WeightsX + Wx;
                         const float * arrIter = extIter;
                         while (weightIter < End)
                              sum += *(weightIter++) * float (*(arrIter++));
                         extIter++;
                         *(resIter++) = sum;

                         //for rescale
                         if(max_val<*arrIter) max_val = *arrIter;
                         if(min_val>*arrIter) min_val = *arrIter;


                    }

               }
          }
          //de-alloc
           if (WeightsX) {delete []WeightsX; WeightsX=0;}
           if (extension_bufferX) {delete []extension_bufferX; extension_bufferX=0;}

     }

     //   Filtering along y
     if(M<2)
     {
          //do nothing
     }
     else
     {
          //create Gaussian kernel
          float  *WeightsY = 0;
          WeightsY = new float [Wy];
          if (!WeightsY)
               return;

          float Half = (float)(Wy-1)/2.0;

          // Gaussian filter equation:
          // http://en.wikipedia.org/wiki/Gaussian_blur
         /* for (unsigned int Weight = 0; Weight < Half; ++Weight)
          {
               const float  y = Half* float (Weight) / float (Half);
               WeightsY[(int)Half - Weight] = WeightsY[(int)Half + Weight] = pi_sigma * exp(-y * y *sigma_s2); // Corresponding symmetric WeightsY
          }*/

          for (unsigned int Weight = 0; Weight <= Half; ++Weight)
          {
              const float  y = float(Weight)-Half;
              WeightsY[Weight] = WeightsY[Wy-Weight-1] = pi_sigma * exp(-(y * y *sigma_s2)); // Corresponding symmetric WeightsY
          }


          double k = 0.;
          for (unsigned int Weight = 0; Weight < Wy; ++Weight)
               k += WeightsY[Weight];

          for (unsigned int Weight = 0; Weight < Wy; ++Weight)
               WeightsY[Weight] /= k;

          //	along y
          float  *extension_bufferY = 0;
          extension_bufferY = new float [M + (Wy<<1)];

          unsigned int offset = Wy>>1;
          const float *extStop = extension_bufferY + M + offset;

          for(V3DLONG iz = 0; iz < P; iz++)
          {
               for(V3DLONG ix = 0; ix < N; ix++)
               {
                    float  *extIter = extension_bufferY + Wy;
                    for(V3DLONG iy = 0; iy < M; iy++)
                    {
                         *(extIter++) = pImage[iz*M*N + iy*N + ix];
                    }

                    //   Extend image
                    const float  *const stop_line = extension_bufferY - 1;
                    float  *extLeft = extension_bufferY + Wy - 1;
                    const float  *arrLeft = extLeft + 2;
                    float  *extRight = extLeft + M + 1;
                    const float  *arrRight = extRight - 2;

                    while (extLeft > stop_line)
                    {
                         *(extLeft--) = *(arrLeft++);
                         *(extRight++) = *(arrRight--);
                    }

                    //	Filtering
                    extIter = extension_bufferY + offset;

                    float  *resIter = &(pImage[iz*M*N + ix]);

                    while (extIter < extStop)
                    {
                         double sum = 0.;
                         const float  *weightIter = WeightsY;
                         const float  *const End = WeightsY + Wy;
                         const float * arrIter = extIter;
                         while (weightIter < End)
                              sum += *(weightIter++) * float (*(arrIter++));
                         extIter++;
                         *resIter = sum;
                         resIter += N;

                         //for rescale
                         if(max_val<*arrIter) max_val = *arrIter;
                         if(min_val>*arrIter) min_val = *arrIter;


                    }

               }
          }

          //de-alloc
          if (WeightsY) {delete []WeightsY; WeightsY=0;}
          if (extension_bufferY) {delete []extension_bufferY; extension_bufferY=0;}


     }

     //  Filtering  along z
     if(P<2)
     {
          //do nothing
     }
     else
     {
          //create Gaussian kernel
          float  *WeightsZ = 0;
          WeightsZ = new float [Wz];
          if (!WeightsZ)
               return;

          float Half = (float)(Wz-1)/2.0;

         /* for (unsigned int Weight = 1; Weight < Half; ++Weight)
          {
               const float  z = Half * float (Weight) / Half;
               WeightsZ[(int)Half - Weight] = WeightsZ[(int)Half + Weight] = pi_sigma * exp(-z * z * sigma_s2) ; // Corresponding symmetric WeightsZ
          }*/

          for (unsigned int Weight = 0; Weight <= Half; ++Weight)
          {
              const float  z = float(Weight)-Half;
              WeightsZ[Weight] = WeightsZ[Wz-Weight-1] = pi_sigma * exp(-(z * z *sigma_s2)); // Corresponding symmetric WeightsZ
          }


          double k = 0.;
          for (unsigned int Weight = 0; Weight < Wz; ++Weight)
               k += WeightsZ[Weight];

          for (unsigned int Weight = 0; Weight < Wz; ++Weight)
               WeightsZ[Weight] /= k;

          //	along z
          float  *extension_bufferZ = 0;
          extension_bufferZ = new float [P + (Wz<<1)];

          unsigned int offset = Wz>>1;
          const float *extStop = extension_bufferZ + P + offset;

          for(V3DLONG iy = 0; iy < M; iy++)
          {
               for(V3DLONG ix = 0; ix < N; ix++)
               {

                    float  *extIter = extension_bufferZ + Wz;
                    for(V3DLONG iz = 0; iz < P; iz++)
                    {
                         *(extIter++) = pImage[iz*M*N + iy*N + ix];
                    }

                    //   Extend image
                    const float  *const stop_line = extension_bufferZ - 1;
                    float  *extLeft = extension_bufferZ + Wz - 1;
                    const float  *arrLeft = extLeft + 2;
                    float  *extRight = extLeft + P + 1;
                    const float  *arrRight = extRight - 2;

                    while (extLeft > stop_line)
                    {
                         *(extLeft--) = *(arrLeft++);
                         *(extRight++) = *(arrRight--);
                    }

                    //	Filtering
                    extIter = extension_bufferZ + offset;

                    float  *resIter = &(pImage[iy*N + ix]);

                    while (extIter < extStop)
                    {
                         double sum = 0.;
                         const float  *weightIter = WeightsZ;
                         const float  *const End = WeightsZ + Wz;
                         const float * arrIter = extIter;
                         while (weightIter < End)
                              sum += *(weightIter++) * float (*(arrIter++));
                         extIter++;
                         *resIter = sum;
                         resIter += M*N;

                         //for rescale
                         if(max_val<*arrIter) max_val = *arrIter;
                         if(min_val>*arrIter) min_val = *arrIter;

                    }

               }
          }

          //de-alloc
          if (WeightsZ) {delete []WeightsZ; WeightsZ=0;}
          if (extension_bufferZ) {delete []extension_bufferZ; extension_bufferZ=0;}


     }

    outimg = pImage;


    return;
}



//void connect_swc(NeuronTree nt,QList<NeuronSWC>& newNeuron, double disThr,double angThr)
//{
//    //rescale neurons
//    QList<XYZ> scaledXYZ;
//    for(V3DLONG i=0; i<nt.listNeuron.size(); i++){
//        XYZ S;
//        S.x = nt.listNeuron.at(i).x*1;
//        S.y = nt.listNeuron.at(i).y*1;
//        S.z = nt.listNeuron.at(i).z*4;
//        scaledXYZ.append(S);
//    }

//    qDebug()<<"search for components and tips";
//    //initialize tree components and get all tips
//    QList<V3DLONG> cand;
//    QList<XYZ> canddir;
//    QVector<int> childNum(nt.listNeuron.size(), 0);
//    QVector<int> connNum(nt.listNeuron.size(), 0);
//    QList<V3DLONG> components;
//    QList<V3DLONG> pList;
//    V3DLONG curid=0;
//    for(V3DLONG i=0; i<nt.listNeuron.size(); i++){
//        if(nt.listNeuron.at(i).pn<0){
//            connNum[i]--; //root that only have 1 child will also be a dead end
//            components.append(curid); curid++;
//            pList.append(-1);
//        }else{
//            V3DLONG pid = nt.hashNeuron.value(nt.listNeuron.at(i).pn);
//            childNum[pid]++;
//            connNum[pid]++;
//            components.append(-1);
//            pList.append(pid);
//        }
//    }

//    qDebug()<<"components searching";
//    //connected component
//    for(V3DLONG cid=0; cid<curid; cid++){
//        QStack<V3DLONG> pstack;
//        V3DLONG chid;
//        //recursively search for child and mark them as the same component
//        pstack.push(components.indexOf(cid));
//        while(!pstack.isEmpty()){
//            V3DLONG pid=pstack.pop();
//            chid = -1;
//            chid = pList.indexOf(pid,chid+1);
//            while(chid>=0){
//                pstack.push(chid);
//                components[chid]=cid;
//                chid=pList.indexOf(pid,chid+1);
//            }
//        }
//    }

//    qDebug()<<"tips searching";
//    vector< pair<int,int> > tip_pair;
//    //get tips
//    for(V3DLONG i=0; i<childNum.size(); i++){
//        if(connNum.at(i)<1){
//            cand.append(i);
//            //get direction
//            V3DLONG id=i;
//            V3DLONG sid;
//            if(childNum[id]==1){ //single child root
//                sid = pList.indexOf(id);
//            }else{ //tips
//                sid = pList[id];
//            }
//            tip_pair.push_back(std::make_pair(id,sid));
//        }
//    }

//    qDebug()<<connNum.size()<<":"<<childNum.size()<<":"<<cand.size();

//    qDebug()<<"match tips";

//    //match tips
//    multimap<double, QVector<V3DLONG> > connMap;
//    for(V3DLONG tid=0; tid<cand.size(); tid++){
//        V3DLONG tidx=cand.at(tid);
//        V3DLONG mvid=-1, mtid=-1;
//        for(V3DLONG cid=0; cid<curid; cid++){
//            if(cid==components.at(cand[tid])) continue;
//            double mvdis=disThr, mtdis=disThr;
//            V3DLONG id=components.indexOf(cid);
//            while(id>=0){

//                double dis=NTDIS(scaledXYZ.at(tidx),scaledXYZ.at(id));
//                if(dis<mvdis){
//                    mvdis=dis;
//                    mvid=id;
//                }
//                if(dis<mtdis){
//                    if(connNum.at(id)<1){//tips
//                        V3DLONG tmpid=cand.indexOf(id);
//                        double local_ang1 = angle(nt.listNeuron.at(tip_pair[tid].first),nt.listNeuron.at(tip_pair[tid].second),nt.listNeuron.at(tip_pair[tmpid].first));
//                        double local_ang2 = angle(nt.listNeuron.at(tip_pair[tmpid].first),nt.listNeuron.at(tip_pair[tmpid].second),nt.listNeuron.at(tip_pair[tid].first));
//                        if(local_ang1 >= angThr && local_ang2 >= angThr){
//                            mtdis=dis;
//                            mtid=id;
//                        }
//                    }
//                }
//                id=components.indexOf(cid, id+1);
//            }

//            if(mtid>=0){
//                QVector<V3DLONG> tmp;
//                tmp.append(tidx); tmp.append(mtid);
//                connMap.insert(pair<double, QVector<V3DLONG> >(mtdis,tmp));
//            }
//        }
//    }

//    qDebug()<<"connecting tips";
//    //find the best solution for connecting tips
//    QMap<V3DLONG, QVector<V3DLONG> > connectPairs;
//    for(multimap<double, QVector<V3DLONG> >::iterator iter=connMap.begin(); iter!=connMap.end(); iter++){
//        if(components.at(iter->second.at(0))==components.at(iter->second.at(1))) //already connected
//            continue;
//        if(connectPairs.contains(iter->second.at(0))){
//            connectPairs[iter->second.at(0)].append(iter->second.at(1));
//        }else{
//            QVector<V3DLONG> tmp; tmp.append(iter->second.at(1));
//            connectPairs.insert(iter->second.at(0),tmp);
//        }
//        if(connectPairs.contains(iter->second.at(1))){
//            connectPairs[iter->second.at(1)].append(iter->second.at(0));
//        }else{
//            QVector<V3DLONG> tmp; tmp.append(iter->second.at(0));
//            connectPairs.insert(iter->second.at(1),tmp);
//        }
//        V3DLONG cid_0=components.at(iter->second.at(0));
//        V3DLONG cid_1=components.at(iter->second.at(1));
//        V3DLONG tmpid=components.indexOf(cid_1);
//        while(tmpid>=0){
//            components[tmpid]=cid_0;
//            tmpid=components.indexOf(cid_1,tmpid+1);
//        }
//    }

//    qDebug()<<"reconstruct neuron tree";
//    //reconstruct tree
//    QVector<V3DLONG> newid(nt.listNeuron.size(), -1);
//    QVector<V3DLONG> newpn(nt.listNeuron.size(), -1); //id starts from 1, -1: not touched, 0: touched but overlap with parent
//    curid=1;
//    int rootID = -1;
//    int rootidx=nt.hashNeuron.value(rootID);
//    if(nt.listNeuron[rootidx].n != rootID)
//        rootidx=-1;
//    QVector<V3DLONG> prinode;
//    if(rootidx!=-1){
//        prinode.push_back(rootidx);
//    }
//    for(V3DLONG i=0; i<nt.listNeuron.size(); i++){
//        if(nt.listNeuron[i].parent==-1){
//            prinode.push_back(i);
//        }
//    }
//    V3DLONG i=0;
//    V3DLONG priIdx=0;
//    while(1){
//        if(priIdx<prinode.size()){
//            i=prinode[priIdx];
//            priIdx++;
//        }else if(priIdx==prinode.size()){
//            i=0;
//            priIdx++;
//        }else{
//            i++;
//            if(i>=nt.listNeuron.size())
//                break;
//        }
//        if(newid[i]>0) continue;
//        QQueue<V3DLONG> pqueue; pqueue.clear();
//        pqueue.enqueue(i);
//        newid[i]=curid++;
//        while(!pqueue.isEmpty()){
//            //add current node to the listNeuron
//            V3DLONG oid=pqueue.dequeue();

//            if(newid[oid]>0){
//                NeuronSWC tmpNeuron;
//                tmpNeuron.n = newid[oid];
//                tmpNeuron.x = nt.listNeuron.at(oid).x;
//                tmpNeuron.y = nt.listNeuron.at(oid).y;
//                tmpNeuron.z = nt.listNeuron.at(oid).z;
//                tmpNeuron.type = nt.listNeuron.at(oid).type;
//                tmpNeuron.r = nt.listNeuron.at(oid).r;
//                tmpNeuron.fea_val = nt.listNeuron.at(oid).fea_val;
//                tmpNeuron.pn = newpn.at(oid);
//                newNeuron.append(tmpNeuron);
//            }

//            //add current node's children/parent/new-neighbor to the stack
//            //parent
//            if(nt.listNeuron.at(oid).pn>=0){
//                V3DLONG opid = nt.hashNeuron.value(nt.listNeuron.at(oid).pn);
//                if(newid.at(opid)<0){
//                    pqueue.enqueue(opid);
//                    newpn[opid]=newid[oid];
//                    newid[opid]=curid++;
//                }
//            }
//            //child
//            V3DLONG tmpid=pList.indexOf(oid);
//            while(tmpid>=0){
//                if(newid.at(tmpid)<0){
//                    pqueue.enqueue(tmpid);
//                    newpn[tmpid]=newid[oid];
//                    newid[tmpid]=curid++;
//                }
//                tmpid=pList.indexOf(oid,tmpid+1);
//            }
//            //new-neighbor
//            if(connectPairs.contains(oid)){
//                for(V3DLONG j=0; j<connectPairs[oid].size(); j++){
//                    V3DLONG onid=connectPairs[oid].at(j);
//                    if(newid.at(onid)<0){
//                        pqueue.enqueue(onid);
//                        newpn[onid]=newid[oid];
//                        newid[onid]=curid++;
//                    }
//                }
//            }
//        }
//    }
//}
