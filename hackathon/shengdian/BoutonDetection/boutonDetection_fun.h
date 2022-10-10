#ifndef BOUTONDETECTION_FUN_H
#define BOUTONDETECTION_FUN_H

#include <QtGui>
#include <v3d_interface.h>
#include "basic_surf_objs.h"
#include "color_xyz.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include "neuron_format_converter.h"
#include "v_neuronswc.h"
#include <functional>
#include <QHash>
#include "volimg_proc.h"
#include <fstream>
#include "FL_coordDefinition.h"
#include "FL_defType.h"
#include "neuronbranchtree.h"
#include "../../SGuo/imPreProcess/helpFunc.h"
#include "../../SGuo/imPreProcess/bilateral_filter.h"

using namespace std;

#define BoutonType 1
#define BoutonSWCNodeType 5
#define MinRefineAngle 60.0
#define PI 3.1415926
#define NeuronTreeRefine 0
#define NodeRefine 1
#define RefineAllinOne 2
#define RadiusEstimation 3

struct Bouton_Color_Basic
{
    int color_ele_r;
    int color_ele_g;
    int color_ele_b;
    Bouton_Color_Basic() {color_ele_r=color_ele_g=0;color_ele_b=255;}
};
struct Bouton_Color_List
{
    QList <Bouton_Color_Basic> listcolor;
    bool on;
    Bouton_Color_List()
    {
        listcolor.clear();
        on=true;
        init();
    }
    void init()
    {
        int colorR[10]={120,200,0,200,0,220,188,180,250,120};
        int colorG[10]={200,20,20,0,200,200,200,94,200,100};
        int colorB[10]={200,0,200,200,200,0,20,37,120,120};
        for(int i=0;i<10;i++)
        {
            Bouton_Color_Basic tmp;
            tmp.color_ele_r=colorR[i];
            tmp.color_ele_g=colorG[i];
            tmp.color_ele_b=colorB[i];
            listcolor.append(tmp);
        }
    }
};
struct AxonalBouton
{
    float x,y,z;
    float ccfx,ccfy,ccfz;
    int nodetype,broder,btype/*0:not botuon;1: bouton;2: tip bouton;3: terminaux bouton*/;
    Bouton_Color_Basic color;
    float path_dist_to_soma,euler_dist_to_soma;
    float nodeR,ccfR,br_r_mean,br_r_std;
    int intensity, br_intensity_mean,br_intensity_std;
    float density,vol_density;
    long region_id;
    AxonalBouton() {
        x=y=z=ccfx=ccfy=ccfz=0.0;
        path_dist_to_soma=euler_dist_to_soma=0.0;
        nodetype=2;btype=0;
        broder=1;
        intensity=br_intensity_mean=br_intensity_std=0;
        density=vol_density=nodeR=br_r_mean=br_r_std=ccfR=0.0;
        color.color_ele_r=color.color_ele_g=0;color.color_ele_b=255;
        region_id=0;
    }
    void init_bouton(NeuronSWC s){
        this->x=s.x;this->y=s.y;this->z=s.z;
        this->nodetype=s.type;
        this->broder=s.level;
        this->nodeR=s.radius;
        int sfsize=s.fea_val.size();
//        cout<<"fea_val size="<<sfsize<<endl;
        switch (sfsize) {
        case 0:
            this->ccfx=this->ccfy=this->ccfz=0.0;
            break;
        case 3:
            this->ccfx=s.fea_val.at(0);
            this->ccfy=s.fea_val.at(1);
            this->ccfz=s.fea_val.at(2);
            break;
        case 4:
            this->ccfx=s.fea_val.at(0);
            this->ccfy=s.fea_val.at(1);
            this->ccfz=s.fea_val.at(2);
            this->btype=s.fea_val.at(3);
            break;
        case 10:
            this->ccfx=s.fea_val.at(0);
            this->ccfy=s.fea_val.at(1);
            this->ccfz=s.fea_val.at(2);
            this->btype=s.fea_val.at(3);
            this->ccfR=s.fea_val.at(4);
            this->br_r_mean=s.fea_val.at(5);
            this->br_r_std=s.fea_val.at(6);
            this->intensity=s.fea_val.at(7);
            this->br_intensity_mean=s.fea_val.at(8);
            this->br_intensity_std=s.fea_val.at(9);
            break;
        case 14:
            this->ccfx=s.fea_val.at(0);
            this->ccfy=s.fea_val.at(1);
            this->ccfz=s.fea_val.at(2);
            this->btype=s.fea_val.at(3);
            this->ccfR=s.fea_val.at(4);
            this->br_r_mean=s.fea_val.at(5);
            this->br_r_std=s.fea_val.at(6);
            this->intensity=s.fea_val.at(7);
            this->br_intensity_mean=s.fea_val.at(8);
            this->br_intensity_std=s.fea_val.at(9);
            this->density=s.fea_val.at(10);
            this->vol_density=s.fea_val.at(11);
            this->path_dist_to_soma=s.fea_val.at(12);
            this->euler_dist_to_soma=s.fea_val.at(13);
            break;
        case 15:
            this->ccfx=s.fea_val.at(0);
            this->ccfy=s.fea_val.at(1);
            this->ccfz=s.fea_val.at(2);
            this->btype=s.fea_val.at(3);
            this->ccfR=s.fea_val.at(4);
            this->br_r_mean=s.fea_val.at(5);
            this->br_r_std=s.fea_val.at(6);
            this->intensity=s.fea_val.at(7);
            this->br_intensity_mean=s.fea_val.at(8);
            this->br_intensity_std=s.fea_val.at(9);
            this->density=s.fea_val.at(10);
            this->vol_density=s.fea_val.at(11);
            this->path_dist_to_soma=s.fea_val.at(12);
            this->euler_dist_to_soma=s.fea_val.at(13);
            this->region_id=s.fea_val.at(14);
            break;
        default:
            break;
        }
        if(this->ccfR==0)
            this->ccfR=this->nodeR;
    }
    void reset(){
        x=y=z=ccfx=ccfy=ccfz=0.0;
        path_dist_to_soma=euler_dist_to_soma=0.0;
        nodetype=2;btype=0;
        broder=1;
        intensity=br_intensity_mean=br_intensity_std=0;
        density=vol_density=nodeR=br_r_mean=br_r_std=ccfR=0.0;
        color.color_ele_r=color.color_ele_g=0;color.color_ele_b=255;
        region_id=0;
    }
    NeuronSWC out_to_NeuronSWC(){
        NeuronSWC s; s.fea_val.clear();
        s.x=this->x; s.y=this->y; s.z=this->z;
        s.type=this->nodetype;
        s.r=this->nodeR;
        s.level=this->broder;
        s.fea_val.append(this->ccfx);
        s.fea_val.append(this->ccfy);
        s.fea_val.append(this->ccfz);
        s.fea_val.append(this->btype);
        s.fea_val.append(this->ccfR);
        s.fea_val.append(this->br_r_mean);
        s.fea_val.append(this->br_r_std);
        s.fea_val.append(this->intensity);
        s.fea_val.append(this->br_intensity_mean);
        s.fea_val.append(this->br_intensity_std);
        s.fea_val.append(this->density);
        s.fea_val.append(this->vol_density);
        s.fea_val.append(this->path_dist_to_soma);
        s.fea_val.append(this->euler_dist_to_soma);
        return s;
    }
    void out_to_NeuronSWC(NeuronSWC& s){
        s.x=this->x; s.y=this->y; s.z=this->z;
        s.type=this->nodetype;
        s.r=this->nodeR;
        s.level=this->broder;
        s.fea_val.clear();
        s.fea_val.append(this->ccfx);
        s.fea_val.append(this->ccfy);
        s.fea_val.append(this->ccfz);
        s.fea_val.append(this->btype);//index=3
        s.fea_val.append(this->ccfR);
        s.fea_val.append(this->br_r_mean);
        s.fea_val.append(this->br_r_std);
        s.fea_val.append(this->intensity);
        s.fea_val.append(this->br_intensity_mean);
        s.fea_val.append(this->br_intensity_std);
        s.fea_val.append(this->density);//index=10
        s.fea_val.append(this->vol_density);//index=11
        s.fea_val.append(this->path_dist_to_soma);
        s.fea_val.append(this->euler_dist_to_soma);
    }
    CellAPO out_to_APO(){
        CellAPO s;
        switch (this->btype) {
        case 1:
            s.color.r=0; s.color.g=255; s.color.b=0;
            break;
        case 2:
            s.color.r=0; s.color.g=255; s.color.b=0;
            break;
        case 3:
            s.color.r=0; s.color.g=0; s.color.b=255;
            break;
        default:
            return s;
        }
        s.x=this->x; s.y=this->y; s.z=this->z;
        s.volsize=this->nodeR;
        s.comment="bouton sites";
        s.intensity=this->intensity;
        return s;
    }
    ImageMarker out_to_ImageMarker(bool ccf=true){
        ImageMarker s;
        switch (this->btype) {
        case 1:
            s.color.r=0; s.color.g=255; s.color.b=0;
            break;
        case 2:
            s.color.r=0; s.color.g=255; s.color.b=0;
            break;
        case 3:
            s.color.r=0; s.color.g=0; s.color.b=255;
            break;
        default:
            return s;
        }
        if(ccf){
            s.x=this->ccfx; s.y=this->ccfy; s.z=this->ccfz;
            s.radius=this->ccfR;
        }
        else{
            s.x=this->x; s.y=this->y; s.z=this->z;
            s.radius=this->nodeR;
        }
        s.x+=(float)1.0;
        s.y+=(float)1.0;
        s.z+=(float)1.0;
        s.comment="bouton sites";
        if(this->btype==BoutonType||this->btype==BoutonType+1)
            s.comment="En-passant-bouton";
        else if(this->btype==BoutonType+2)
            s.comment="Terminaux-bouton";
        s.shape=1;
        return s;
    }
};
/*preprocess*/
void neuron_checking_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);
void preprocess_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);
QList<NeuronTree> nt_2_trees(NeuronTree& nt);
NeuronTree preprocess_simple(NeuronTree nt);//for swc in an image block

void postprocess_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);

/*refinement: 1. mean-shift; 2. node_refine;3.line_refine*/
void refinement_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output,bool in_terafly=true);
void refinement_terafly_fun(V3DPluginCallback2 &callback,string imgPath, NeuronTree& nt,int method_code=0,int refine_radius=8,long half_block_size=128,int nodeRefine_radius=2);
void refinement_Image_fun(V3DPluginCallback2 &callback,string imgPath, NeuronTree& nt,int method_code=0,int refine_radius=8,int nodeRefine_radius=2);
NeuronSWC calc_mean_shift_center(unsigned char * & inimg1d,NeuronSWC snode,V3DLONG sz_image[], double bkg_thre=40, int windowradius=3);

//NeuronSWC calc_mean_shift_center_v4(unsigned char * & inimg1d,NeuronSWC snode,V3DLONG sz_image[], double &bkg_thre,int &windowradius_pid,int windowradius,int windowradius_limit);
//NeuronSWC calc_mean_shift_center_v5(unsigned char * & inimg1d,NeuronSWC snode,V3DLONG sz_image[], double &bkg_thre,int &windowradius_pid,int windowradius,int windowradius_limit);
NeuronSWC calc_mean_shift_center_v4(unsigned char * & inimg1d,NeuronSWC snode,V3DLONG sz_image[], double &bkg_thre,int &windowradius_pid,int windowradius=3,int windowradius_limit=15);
NeuronSWC calc_mean_shift_center_v5(unsigned char * & inimg1d,NeuronSWC snode,V3DLONG sz_image[], double &bkg_thre,int &windowradius_pid,int windowradius=3,int windowradius_limit=15);

NeuronSWC nodeRefine(unsigned char * & inimg1d,NeuronSWC s,V3DLONG * sz,int neighbor_size=2);
double getAngleofNodeVector(NeuronSWC n0,NeuronSWC n1,NeuronSWC n2);
NeuronSWC lineRefine(unsigned char * & inimg1d,V3DLONG * sz,NeuronSWC snode,NeuronSWC spnode, int sqhere_radius=5,int searching_line_radius=2);

/*bouton detection*/
void boutonDetection_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output,bool in_terafly);
void swc_profile_terafly_fun(V3DPluginCallback2 &callback,string imgPath, NeuronTree& nt,long block_size=128,int bkg_bias=15,double upfactor=2.0);
void swc_profile_image_fun(V3DPluginCallback2 &callback,string inimg_file, NeuronTree& nt,int bkg_bias=15,double upfactor=2.0);
void boutonFilter_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);

/*filter*/
QList <AxonalBouton> boutonFilter_fun(NeuronTree nt,double radius_delta=1.5,double intensity_delta=1.5,double AXON_BACKBONE_RADIUS=2);
void map_bouton_2_neuronTree(NeuronTree& nt_bouton,QList <AxonalBouton>  bouton_sites);
void swc_profile_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output,bool in_terafly);
std::vector<double> get_sorted_fea_of_seg(V_NeuronSWC inseg,bool radiusfea=false);
std::vector<double> mean_and_std_seg_fea(std::vector<double> input);
std::vector<int> peaks_in_seg(std::vector<double> input,int isRadius_fea=0,float delta=0.5);

/*feature computation*/
void bouton_feature_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);
void boutonType_label(NeuronTree& nt,bool ccf_domain=true);
void rendering_different_bouton(NeuronTree& nt, int type_bias=BoutonSWCNodeType);
void boutonVolDesity_computing(NeuronTree& nt,float vol_r,float xy_res=0.3,bool ccf_domain=true);
void boutonDesity_computing(NeuronTree& nt,float xy_res=0.3,bool ccf_domain=true);
void dist_to_soma(NeuronTree& nt,float xy_res=0.3,bool ccf_domain=true,bool bouton_only=true);
/*pruning*/
void boutonswc_pruning_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);
NeuronTree boutonSWC_internode_pruning(NeuronTree nt,float pruning_dist=1.0,bool ccf_domain=false);
void nearBouton_pruning(NeuronTree& nt,float pruning_dist=5.0,bool ccf_domain=false);
NeuronTree tipNode_pruning(NeuronTree nt, float pruning_dist=1.0,bool ccf_domain=false);
void sparseBouton_pruning(NeuronTree& nt,float pruning_dist,int pruning_num,bool ccf_domain=false);
bool smooth_radius(NeuronTree& nt,int win=5,bool notbouton=true);
/*file io*/
void bouton_file_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);
void boutonswc_to_ccf(NeuronTree& nt,float scale=1.0);
QList<CellAPO> bouton_to_apo(NeuronTree nt);
QList<ImageMarker> bouton_to_imageMarker(NeuronTree nt,bool ccf=true);
void featureTable(const QString &filename,NeuronTree nt,float *res);
void getBoutonMIP(V3DPluginCallback2 &callback, unsigned char *& inimg1d,V3DLONG in_sz[],QString outpath);
void getBoutonBlock(V3DPluginCallback2 &callback, string imgPath,NeuronTree nt,QString outpath,
                    int getbtype=BoutonType,int crop_half_size=16,bool mip_flag=false,int mask_size=0);
void getBoutonBlock_inImg(V3DPluginCallback2 &callback,string inimg_file,QList <CellAPO> apolist,string outpath,int block_size=16);
void maskImg(unsigned char *&inimg1d,unsigned char * & im_transfer,long in_sz[], NeuronTree nt, int maskRadius=8);

void ccf_profile_dofunc(V3DPluginCallback2 & callback, const V3DPluginArgList & input,V3DPluginArgList & output);// for processing registered swc
void scale_registered_swc(NeuronTree& nt,float xshift_pixels=20.0,float scale_xyz=25.0);
void merge_raw_swc_onto_reg(NeuronTree nt_raw,NeuronTree& nt_registered, float rscale=1.0);
void merge_registered_swc_onto_raw(NeuronTree& nt_raw,NeuronTree nt_registered);

/*radius estimation*/
double radiusEstimation(unsigned char *&inimg1d, long in_zz[], NeuronSWC s, double dfactor, double bkg_thresh);
double radiusEstimation(V3DPluginCallback2 &callback,unsigned char * & inimg1d,V3DLONG in_zz[4], NeuronSWC s,double dfactor,double bkg_thresh);
double getNodeRadius(unsigned char * & inimg1d,V3DLONG in_sz[4], NeuronSWC s,double bkg_thresh=40);
double getNodeRadius_XY(unsigned char *&inimg1d, long in_sz[], NeuronSWC s, double bkg_thresh,int z_half_win_size);
double getNodeRadius_XY(unsigned char *&inimg1d, long in_sz[], NeuronSWC s, double bkg_thresh);
void getNTRadius_XY(unsigned char *&inimg1d, long in_sz[], NeuronTree& nt, double bkg_thresh=40);

/*Image processing*/
bool enhanceImage(unsigned char * & data1d,unsigned char * & dst,V3DLONG *mysz,bool biilateral_filter=false);
void adaptiveThresholding(unsigned char * & data1d,unsigned char * & dst,V3DLONG *mysz);
bool upsampleImage(unsigned char * & inimg1d,unsigned char * & outimg1d,V3DLONG *szin, V3DLONG *szout, double *dfactor);
template <class T> bool upsample3dvol(T *outdata, T *indata, V3DLONG *szout, V3DLONG *szin, double *dfactor);
template <class T> bool interpolate_coord_linear(T * interpolatedVal, Coord3D *c, V3DLONG numCoord,
                       T *** templateVol3d, V3DLONG tsz0, V3DLONG tsz1, V3DLONG tsz2,
                       V3DLONG tlow0, V3DLONG tup0, V3DLONG tlow1, V3DLONG tup1, V3DLONG tlow2, V3DLONG tup2);
template <class T> bool interpolate_coord_linear(T * interpolatedVal, Coord3D *c, V3DLONG numCoord,
                                                 T * templateVol3d, V3DLONG tsz0, V3DLONG tsz1, V3DLONG tsz2,
                                                 V3DLONG tlow0, V3DLONG tup0, V3DLONG tlow1, V3DLONG tup1, V3DLONG tlow2, V3DLONG tup2);

/*other func:*/
bool teraImage_swc_crop(V3DPluginCallback2 &callback, string inimg, string inswc, string inapo,QString save_path, int cropx, int cropy, int cropz);
QList <CellAPO> nt_2_multi_centers(NeuronTree nt,float xs=256.0,float ys=256.0,float zs=256.0);
void printHelp();
#endif // BOUTONDETECTION_FUN_H
