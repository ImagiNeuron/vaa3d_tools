#ifndef __REFINE_SWC_FUNCS_H__
#define __REFINE_SWC_FUNCS_H__

#include <QtGui>
#include <v3d_interface.h>
#include <vector>

using namespace std;

struct Point;
struct Point
{
    double x,y,z,r;
    V3DLONG type;
    Point* p;
    V3DLONG childNum;
    V3DLONG level,seg_id;
    QList<float> fea_val;
};

typedef vector<Point*> Segment;
typedef vector<Point*> Tree;

NeuronTree refineSWCTerafly(V3DPluginCallback2 &callback,QString fname_img, NeuronTree nt);
NeuronTree breakSWC(NeuronTree nt, double dist);
NeuronTree smoothSWCTerafly(V3DPluginCallback2 &callback,QString fname_img, NeuronTree nt, QList<ImageMarker>& break_points);
double distTwoSegs(NeuronTree nt, Segment * seg);
void break_path(Segment * seg, double step, QList<ImageMarker>& break_points);
NeuronTree breakSWC(NeuronTree input, double step, QList<ImageMarker>& break_points);

NeuronTree initialSWCTerafly(V3DPluginCallback2 &callback,QString fname_img, NeuronTree nt, double ds_rate);
void resample_path(Segment * seg, double step);
NeuronTree resample(NeuronTree input, double step);


// added by DZC 14Mar2019
NeuronTree breakSWC_with_points(NeuronTree input, QList<ImageMarker>& break_points);
void break_path_with_points(Segment * seg, QList<ImageMarker>& break_points);
QList <ImageMarker> break_points_reselect(QList<ImageMarker> & break_points, NeuronTree nt);
NeuronTree neuronlist_2_neurontree(QList<NeuronSWC> neuronlist);



#endif

