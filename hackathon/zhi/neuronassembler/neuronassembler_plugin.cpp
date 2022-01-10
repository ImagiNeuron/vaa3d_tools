/* neuronassembler_plugin.cpp
 * This is a test plugin, you can use it as a demo.
 * 2014-11-09 : by Zhi Zhou
 */
 
#include "v3d_message.h"
#include <vector>
#include "neuronassembler_plugin.h"
#include "../../../released_plugins/v3d_plugins/istitch/y_imglib.h"
#include "../../../released_plugins/v3d_plugins/neurontracing_vn2/app2/my_surf_objs.h"

#include "stackutil.h"
#include "../APP2_large_scale/readRawfile_func.h"
#include <boost/lexical_cast.hpp>


using namespace std;
Q_EXPORT_PLUGIN2(neuronassembler, neuronassembler);

struct root_node
{
    V3DLONG root_x;
    V3DLONG root_y;
    V3DLONG root_z;
    V3DLONG tc_index;
    V3DLONG ref_index;

    V3DLONG start[3];
    V3DLONG end[3];
    V3DLONG parent;
    QString tilename;

    int direction; //1 left, 2 right, 3 up, and 4 down
    struct root_node* next;

};

struct NA_PARA
{
    int  bkg_thresh;
    int  channel;
    int block_size;
    int tracing_method;
    int root_1st[3];

    //NeuTube
    int merge;

    //TReMap
    int is_gsdt;
    int is_break_accept;
    double length_thresh;
    int  cnn_type;
    double SR_ratio;
    int  b_256cube;
    int b_RadiusFrom2D;
    int mip_plane;

    //MST
    int win_size;

    //APP2
    int region_number;
    QString tcfilename,inimg_file,rawfilename,markerfilename;
};

struct Point;
struct Point
{
    double x,y,z,r;
    V3DLONG type;
    Point* p;
    V3DLONG childNum;
};

typedef vector<Point*> Segment;
typedef vector<Point*> Tree;

// group images blending function
template <class SDATATYPE>
int region_groupfusing(SDATATYPE *pVImg, Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim, unsigned char *relative1d,
                       V3DLONG vx, V3DLONG vy, V3DLONG vz, V3DLONG vc, V3DLONG rx, V3DLONG ry, V3DLONG rz, V3DLONG rc,
                       V3DLONG tile2vi_zs, V3DLONG tile2vi_ys, V3DLONG tile2vi_xs,
                       V3DLONG z_start, V3DLONG z_end, V3DLONG y_start, V3DLONG y_end, V3DLONG x_start, V3DLONG x_end, V3DLONG *start)
{

    //
    SDATATYPE *prelative = (SDATATYPE *)relative1d;

    if(x_end<x_start || y_end<y_start || z_end<z_start)
        return false;

    // update virtual image pVImg
    V3DLONG offset_volume_v = vx*vy*vz;
    V3DLONG offset_volume_r = rx*ry*rz;

    V3DLONG offset_pagesz_v = vx*vy;
    V3DLONG offset_pagesz_r = rx*ry;

    for(V3DLONG c=0; c<rc; c++)
    {
        V3DLONG o_c = c*offset_volume_v;
        V3DLONG o_r_c = c*offset_volume_r;
        for(V3DLONG k=z_start; k<z_end; k++)
        {
            V3DLONG o_k = o_c + (k-start[2])*offset_pagesz_v;
            V3DLONG o_r_k = o_r_c + (k-tile2vi_zs)*offset_pagesz_r;

            for(V3DLONG j=y_start; j<y_end; j++)
            {
                V3DLONG o_j = o_k + (j-start[1])*vx;
                V3DLONG o_r_j = o_r_k + (j-tile2vi_ys)*rx;

                for(V3DLONG i=x_start; i<x_end; i++)
                {
                    V3DLONG idx = o_j + i-start[0];
                    V3DLONG idx_r = o_r_j + (i-tile2vi_xs);

                    if(pVImg[idx])
                    {
                        pVImg[idx] = 0.5*(pVImg[idx] + prelative[idx_r]); // Avg. Intensity
                    }
                    else
                    {
                        pVImg[idx] = prelative[idx_r];
                    }
                }
            }
        }
    }

    return true;
}


bool assembler_tc(V3DPluginCallback2 &callback, QWidget *parent,NA_PARA &p,bool bmenu);
bool assembler_raw(V3DPluginCallback2 &callback, QWidget *parent,NA_PARA &p,bool bmenu);

void save_region(V3DPluginCallback2 &callback, V3DLONG *start, V3DLONG *end, QString &tcfile,QString &region_name,
                 Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim, NA_PARA &p);
NeuronTree eliminate(NeuronTree input, double length);
 
QStringList neuronassembler::menulist() const
{
	return QStringList() 
		<<tr("trace_tc")
		<<tr("trace_raw")
		<<tr("about");
}

QStringList neuronassembler::funclist() const
{
	return QStringList()
		<<tr("trace_tc")
		<<tr("trace_raw")
		<<tr("help");
}

void neuronassembler::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("trace_tc"))
	{
        NA_PARA P;
        bool bmenu = true;

        NeuronAssemblerDialog dialog(callback, parent);
        if (!dialog.image)
            return;
        if(dialog.listLandmarks.count() ==0)
            return;

        if (dialog.exec()!=QDialog::Accepted)
            return;
        P.inimg_file = dialog.image->getFileName();
        P.tcfilename = dialog.tcfilename;

        if(P.tcfilename.isEmpty())
        {
            v3d_msg("Please select the tc file.");
            return;
        }

        if(dialog.tracing_method ==0)
        {
            v3d_msg("Please select a tracing method.");
            return;
        }


        P.bkg_thresh = dialog.bkg_thresh;
        P.channel = dialog.channel;
        P.block_size = dialog.block_size;
        P.tracing_method = dialog.tracing_method;
        LocationSimple tmpLocation(0,0,0);
        tmpLocation = dialog.listLandmarks.at(0);
        tmpLocation.getCoord(P.root_1st[0],P.root_1st[1],P.root_1st[2]);

        switch(P.tracing_method)
        {
            case 2: P.merge = dialog.merge; break;
            case 4: P.mip_plane = dialog.mip_plane;P.b_256cube = dialog.b_256cube;P.is_gsdt = dialog.is_gsdt;P.is_break_accept = dialog.is_break_accept;P.length_thresh = dialog.length_thresh;break;
            case 5: P.b_256cube = dialog.b_256cube;
           // case 6: P.win_size = dialog.win_size;
        }

        assembler_tc(callback,parent,P,bmenu);
	}
	else if (menu_name == tr("trace_raw"))
	{
        NA_PARA P;
        bool bmenu = true;

        NeuronAssemblerDialog_raw dialog(callback, parent);

        if (dialog.exec()!=QDialog::Accepted)
            return;

        if(dialog.rawfilename.isEmpty())
        {
            v3d_msg("Please select the image file.");
            return;
        }

        if(dialog.markerfilename.isEmpty())
        {
            v3d_msg("Please select the marker file.");
            return;
        }

        if(dialog.tracing_method ==0)
        {
            v3d_msg("Please select a tracing method.");
            return;
        }



        vector<MyMarker> file_inmarkers;
        file_inmarkers = readMarker_file(string(qPrintable(dialog.markerfilename)));

        P.bkg_thresh = dialog.bkg_thresh;
        P.channel = dialog.channel;
        P.block_size = dialog.block_size;
        P.tracing_method = dialog.tracing_method;
        P.inimg_file = dialog.rawfilename;

        P.root_1st[0] = file_inmarkers[0].x;
        P.root_1st[1] = file_inmarkers[0].y;
        P.root_1st[2] = file_inmarkers[0].z;

        switch(P.tracing_method)
        {
            case 2: P.merge = dialog.merge; break;
            case 4: P.mip_plane = dialog.mip_plane;P.b_256cube = dialog.b_256cube;P.is_gsdt = dialog.is_gsdt;P.is_break_accept = dialog.is_break_accept;P.length_thresh = dialog.length_thresh;break;
            case 5: P.b_256cube = dialog.b_256cube;
          //  case 6: P.win_size = dialog.win_size;

        }
        assembler_raw(callback,parent,P,bmenu);
	}
	else
	{
		v3d_msg(tr("This is a test plugin, you can use it as a demo.. "
			"Developed by Zhi Zhou, 2014-11-09"));
	}
}

bool neuronassembler::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
	if (func_name == tr("trace_tc"))
	{
        NA_PARA P;
        bool bmenu = false;

        vector<char*> * pinfiles = (input.size() >= 1) ? (vector<char*> *) input[0].p : 0;
        vector<char*> * pparas = (input.size() >= 2) ? (vector<char*> *) input[1].p : 0;
        vector<char*> infiles = (pinfiles != 0) ? * pinfiles : vector<char*>();
        vector<char*> paras = (pparas != 0) ? * pparas : vector<char*>();

        int k=0;
        QString inmarker_file = paras.empty() ? "" : paras[k]; if(inmarker_file == "NULL") inmarker_file = ""; k++;
        if(inmarker_file.isEmpty())
        {
            cerr<<"Need a marker file"<<endl;
            return false;
        }
        P.tcfilename  = paras.empty() ? "" : paras[k]; k++;


        if(infiles.empty())
        {
            vector<MyMarker> file_inmarkers;
            file_inmarkers = readMarker_file(string(qPrintable(inmarker_file)));
            int soma_x = file_inmarkers[0].x;
            int soma_y = file_inmarkers[0].y;

            Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim;

            if( !vim.y_load( P.tcfilename.toStdString()) )
            {
                printf("Wrong stitching configuration file to be load!\n");
                return false;
            }

            QString curFilePath = QFileInfo( P.tcfilename).path();
            curFilePath.append("/");
            for(V3DLONG ii=0; ii<vim.number_tiles; ii++)
            {
                if(soma_x > vim.lut[ii].start_pos[0] && soma_y > vim.lut[ii].start_pos[1] && soma_x < vim.lut[ii].end_pos[0] && soma_y < vim.lut[ii].end_pos[1])
                {
                    P.inimg_file = curFilePath + QString(vim.lut[ii].fn_img.c_str());
                    v3d_msg(P.inimg_file,0);
                    inmarker_file =  P.inimg_file + ".marker";
                    break;
                }

            }
        }
        else
        {

            P.inimg_file = infiles[0];
        }

        vector<MyMarker> file_inmarkers;
        file_inmarkers = readMarker_file(string(qPrintable(inmarker_file)));
        P.root_1st[0] = file_inmarkers[0].x;
        P.root_1st[1] = file_inmarkers[0].y;
        P.root_1st[2] = file_inmarkers[0].z;

        P.tracing_method = (paras.size() >= k+1) ? atoi(paras[k]) : 1;  k++;
        P.channel = (paras.size() >= k+1) ? atoi(paras[k]) : 1;  k++;
        P.bkg_thresh = (paras.size() >= k+1) ? atof(paras[k]) : 10; k++;
        P.block_size = (paras.size() >= k+1) ? atoi(paras[k]) : 1024; k++;

        switch(P.tracing_method)
        {
            case 2: P.merge = (paras.size() >= k+1) ? atoi(paras[k]) : 1; k++; break;
            case 4: P.mip_plane = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.b_256cube = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.is_gsdt = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.is_break_accept = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.length_thresh = (paras.size() >= k+1) ? atoi(paras[k]) : 20;  k++;
                    break;
            case 5: P.b_256cube = (paras.size() >= k+1) ? atoi(paras[k]) : 0; k++;
                    P.region_number = (paras.size() >= k+1) ? atoi(paras[k]) : 6000;  k++;
                    break;
            case 6: P.b_256cube = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.is_gsdt = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.is_break_accept = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.length_thresh = (paras.size() >= k+1) ? atoi(paras[k]) : 20;  k++;
                    P.region_number = (paras.size() >= k+1) ? atoi(paras[k]) : 6000;  k++;
                    break;
        }

        assembler_tc(callback,parent,P,bmenu);
	}
	else if (func_name == tr("trace_raw"))
	{
        NA_PARA P;
        bool bmenu = false;

        vector<char*> * pinfiles = (input.size() >= 1) ? (vector<char*> *) input[0].p : 0;
        vector<char*> * pparas = (input.size() >= 2) ? (vector<char*> *) input[1].p : 0;
        vector<char*> infiles = (pinfiles != 0) ? * pinfiles : vector<char*>();
        vector<char*> paras = (pparas != 0) ? * pparas : vector<char*>();

        if(infiles.empty())
        {
            cerr<<"Need input image"<<endl;
            return false;
        }

        P.inimg_file = infiles[0];
        int k=0;
        QString inmarker_file = paras.empty() ? "" : paras[k]; if(inmarker_file == "NULL") inmarker_file = ""; k++;
        if(inmarker_file.isEmpty())
        {
            cerr<<"Need a marker file"<<endl;
            return false;
        }
        else
        {
            vector<MyMarker> file_inmarkers;
            file_inmarkers = readMarker_file(string(qPrintable(inmarker_file)));
            P.root_1st[0] = file_inmarkers[0].x;
            P.root_1st[1] = file_inmarkers[0].y;
            P.root_1st[2] = file_inmarkers[0].z;
        }

        P.tracing_method = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
        P.channel = (paras.size() >= k+1) ? atoi(paras[k]) : 1;  k++;
        P.bkg_thresh = (paras.size() >= k+1) ? atoi(paras[k]) : 10; k++;
        P.block_size = (paras.size() >= k+1) ? atof(paras[k]) : 1024; k++;


        switch(P.tracing_method)
        {
            case 2: P.merge = (paras.size() >= k+1) ? atoi(paras[k]) : 1; k++; break;
            case 4: P.mip_plane = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;;
                    P.b_256cube = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;;
                    P.is_gsdt = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;;
                    P.is_break_accept = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;;
                    P.length_thresh = (paras.size() >= k+1) ? atoi(paras[k]) : 20;  k++;;
                    break;
            case 5: P.b_256cube = (paras.size() >= k+1) ? atoi(paras[k]) : 0; k++;
                    P.region_number = (paras.size() >= k+1) ? atoi(paras[k]) : 6000;  k++;
                    break;
            case 6: P.b_256cube = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.is_gsdt = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.is_break_accept = (paras.size() >= k+1) ? atoi(paras[k]) : 0;  k++;
                    P.length_thresh = (paras.size() >= k+1) ? atoi(paras[k]) : 20;  k++;
                    P.region_number = (paras.size() >= k+1) ? atoi(paras[k]) : 6000;  k++;
                    break;

        }

        assembler_raw(callback,parent,P,bmenu);
	}
	else if (func_name == tr("help"))
	{
        printf("\n**** Usage of Neuron Assembler ****\n");
        printf("vaa3d -x plugin_name -f trace_tc -i <inimg_file> -p <inmarker_file> <tc_file> <tracing_method> <channel> <bkg_thresh> <block size>\n");
        printf("inimg_file       Should be 8 bit image\n");
        printf("inmarker_file    Please specify the path of the marker file\n");
        printf("tc_file          Please specify the path of the tc file\n");

        printf("tracing_method   Tracing method.0 for MOST tracing, 1 for NeuTube tracing, 2 for Farsight snake tracing, (default 0).\n");
        printf("channel          Data channel for tracing. Start from 1 (default 1).\n");
        printf("bkg_thresh       Default 10 (is specified as -1 then auto-thresolding)\n");
        printf("block size       Default 1024\n");

        printf("outswc_file      Will be named automatically based on the input image file name, so you don't have to specify it.\n\n");


        printf("vaa3d -x plugin_name -f trace_raw -i <inimg_file> -p <inmarker_file> <tracing_method> <channel> <bkg_thresh> <block size>\n");
        printf("inimg_file       Should be 8 bit image\n");
        printf("inmarker_file    Please specify the path of the marker file\n");

        printf("tracing_method   Tracing method.0 for MOST tracing, 1 for NeuTube tracing, 2 for Farsight snake tracing, (default 0).\n");
        printf("channel          Data channel for tracing. Start from 1 (default 1).\n");
        printf("bkg_thresh       Default 10 (is specified as -1 then auto-thresolding)\n");
        printf("block size       Default 1024\n");

        printf("outswc_file      Will be named automatically based on the input image file name, so you don't have to specify it.\n");

	}
	else return false;

	return true;
}


bool assembler_tc(V3DPluginCallback2 &callback, QWidget *parent,NA_PARA &P,bool bmenu)
{
    QElapsedTimer timer1;
    timer1.start();

    QString tcfile = P.tcfilename;
    QString fileOpenName = P.inimg_file;

    Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim;

    if( !vim.y_load(tcfile.toStdString()) )
    {
        printf("Wrong stitching configuration file to be load!\n");
        return false;
    }


    int soma_index = -1;
    for(V3DLONG ii=0; ii<vim.number_tiles; ii++)
    {
        QString curPath = QFileInfo(tcfile).path();
        QString curtilename = curPath.append("/").append(QString(vim.lut[ii].fn_img.c_str()));
        if(curtilename == fileOpenName)
        {
            soma_index = ii;
            break;
        }
    }

    if(soma_index == -1 )
    {
        printf("Can not find the selected image in the tc file!\n");
        return false;
    }

    V3DLONG start[3], end[3];

    start[0] = vim.lut[soma_index].start_pos[0] - vim.min_vim[0] + P.root_1st[0] - int(P.block_size/2);
    start[1] = vim.lut[soma_index].start_pos[1] - vim.min_vim[1] + P.root_1st[1] - int(P.block_size/2);
    start[2] = 0;

    end[0] = start[0] + P.block_size - 1;
    end[1] = start[1] + P.block_size - 1;
    end[2] = vim.sz[2] - 1;


    struct root_node *head = new root_node[1];
    struct root_node *walker;

    QString startingpos="", tmps;
    tmps.setNum(start[0]).prepend("x"); startingpos += tmps;
    tmps.setNum(start[1]).prepend("_y"); startingpos += tmps;
    QString region_name = startingpos + ".raw";

    head->root_x = int(P.block_size/2);
    head->root_y = int(P.block_size/2);
    head->root_z = vim.lut[soma_index].start_pos[2] - vim.min_vim[2] + P.root_1st[2] ;

    head->tc_index = soma_index;
    head->tilename = QFileInfo(tcfile).path().append("/tmp/").append(QString(region_name));
    head->ref_index = -1;

    head->start[0] = start[0];
    head->start[1] = start[1];
    head->start[2] = start[2];
    head->end[0] = end[0];
    head->end[1] = end[1];
    head->end[2] = end[2];

    head->parent = -1;
    head->direction = 0;
    head->next = NULL;
    walker = head;

    QString rootposstr="", tmps2;
    tmps2.setNum(int(P.root_1st[0]+0.5)).prepend("_x"); rootposstr += tmps2;
    tmps2.setNum(int(P.root_1st[1]+0.5)).prepend("_y"); rootposstr += tmps2;
    tmps2.setNum(int(P.root_1st[2]+0.5)).prepend("_z"); rootposstr += tmps2;

    QString finalswcfilename;

    switch(P.tracing_method){
    case 1: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_MOST.swc"; break;
    case 2: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_NeuTube.swc"; break;
    case 3: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_Snake.swc"; break;
    case 4: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_TReMap.swc"; break;
    case 5: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_APP1.swc"; break;
    case 6: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_APP2.swc"; break;

    }

    QString tmpfolder = QFileInfo(tcfile).path()+("/tmp");
    system(qPrintable(QString("mkdir %1").arg(tmpfolder.toStdString().c_str())));
    if(tmpfolder.isEmpty())
    {

        printf("Can not create a tmp folder!\n");
        return false;
    }


    int swc_type = 2;
    while(walker != NULL)
    {

        ifstream ifs_image(walker->tilename.toStdString().c_str());
        if(!ifs_image)
         save_region(callback,walker->start,walker->end,tcfile,walker->tilename,vim, P);
        else
        {
            walker = walker->next;
            continue;
        }


        QString swcfilename;
        switch(P.tracing_method){
        case 1: swcfilename =  walker->tilename + QString("_MOST.swc"); break;
        case 2: swcfilename =  walker->tilename + QString("_neutube.swc"); break;
        case 3: swcfilename =  walker->tilename + QString("_snake.swc"); break;
        case 4: swcfilename =  walker->tilename + QString("_XY_3D_TreMap.swc"); break;
        case 5: swcfilename =  walker->tilename + QString("_region_APP1.swc"); break;
        case 6: swcfilename =  walker->tilename + QString("_region_APP2.swc"); break;
        }


        V3DPluginArgItem arg;
        V3DPluginArgList input;
        V3DPluginArgList output;

        QString full_plugin_name;
        QString func_name;

        arg.type = "random";std::vector<char*> arg_input;
        std:: string fileName_Qstring(walker->tilename.toStdString());char* fileName_string =  new char[fileName_Qstring.length() + 1]; strcpy(fileName_string, fileName_Qstring.c_str());
        arg_input.push_back(fileName_string);
        arg.p = (void *) & arg_input; input<< arg;

        char* char_swcout =  new char[swcfilename.length() + 1];strcpy(char_swcout, swcfilename.toStdString().c_str());
        arg.type = "random";std::vector<char*> arg_output;arg_output.push_back(char_swcout); arg.p = (void *) & arg_output; output<< arg;


        char channel = '0' + P.channel;
        string T_background = boost::lexical_cast<string>(P.bkg_thresh);
        char* Th =  new char[T_background.length() + 1];
        strcpy(Th, T_background.c_str());

        //NeuTube
        string merge_string = boost::lexical_cast<string>(P.merge);char* merge_char =  new char[T_background.length() + 1]; strcpy(merge_char, merge_string.c_str());

        //TReMap
        string mip_plane_string = boost::lexical_cast<string>(P.mip_plane);char* mip_plane_char =  new char[mip_plane_string.length() + 1]; strcpy(mip_plane_char, mip_plane_string.c_str());
        string b_256cube_string = boost::lexical_cast<string>(P.b_256cube);char* b_256cube_char =  new char[b_256cube_string.length() + 1]; strcpy(b_256cube_char, b_256cube_string.c_str());
        string is_gsdt_string = boost::lexical_cast<string>(P.is_gsdt);char* is_gsdt_char =  new char[is_gsdt_string.length() + 1]; strcpy(is_gsdt_char, is_gsdt_string.c_str());
        string is_break_accept_string = boost::lexical_cast<string>(P.is_break_accept);char* is_break_accept_char =  new char[is_break_accept_string.length() + 1]; strcpy(is_break_accept_char, is_break_accept_string.c_str());
        string length_thresh_string = boost::lexical_cast<string>(P.length_thresh);char* length_thresh_char =  new char[length_thresh_string.length() + 1]; strcpy(length_thresh_char, length_thresh_string.c_str());

        //MST
        string win_size_string = boost::lexical_cast<string>(P.win_size);char* win_size_char =  new char[win_size_string.length() + 1]; strcpy(win_size_char, win_size_string.c_str());

        //APP1
       // if(P.tracing_method == 5) channel = '0' + P.channel - 1;

        //APP2
        string region_number_string = boost::lexical_cast<string>(P.region_number);char* region_number_char =  new char[region_number_string.length() + 1]; strcpy(region_number_char, region_number_string.c_str());


        arg.type = "random";
        std::vector<char*> arg_para;

        switch(P.tracing_method)
        {
            case 1: arg_para.push_back(&channel);arg_para.push_back(Th);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "mostVesselTracer"; func_name = "MOST_trace";break;
            case 2: arg_para.push_back(&channel); arg_para.push_back(merge_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "neuTube"; func_name = "neutube_trace";break;
            case 3: arg_para.push_back(&channel);arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "snake"; func_name = "snake_trace";break;
            case 4: arg_para.push_back(mip_plane_char);arg_para.push_back(&channel);arg_para.push_back(Th);arg_para.push_back(b_256cube_char);arg_para.push_back(is_gsdt_char);arg_para.push_back(is_break_accept_char);arg_para.push_back(length_thresh_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "TReMap"; func_name = "trace_mip";break;
            case 5: arg_para.push_back(&channel);arg_para.push_back(Th);arg_para.push_back(b_256cube_char);arg_para.push_back(region_number_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "Region_Neuron2"; func_name = "trace_app1";break;
            case 6: arg_para.push_back(&channel);arg_para.push_back(Th);arg_para.push_back(b_256cube_char);arg_para.push_back(is_gsdt_char);arg_para.push_back(is_break_accept_char);arg_para.push_back(length_thresh_char);arg_para.push_back(region_number_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "Region_Neuron2"; func_name = "trace_app2";break;
//            case 6: arg_para.push_back(&channel);arg_para.push_back(win_size_char);
//                    arg.p = (void *) & arg_para; input << arg;
//                    full_plugin_name = "MST_tracing"; func_name = "trace_mst";break;
        }
       // arg.type = "random";std::vector<char*> arg_output;arg_output.push_back(fileName_string); arg.p = (void *) & arg_output; output<< arg;

        if(!callback.callPluginFunc(full_plugin_name,func_name,input,output))
        {

             printf("Can not find the tracing plugin!\n");
             return false;
        }




//        #if  defined(Q_OS_LINUX)
//            QString cmd_method;
//            switch(P.tracing_method){
//            case 0: cmd_method = QString("%1/vaa3d -x mostVesselTracer -f MOST_trace -i %2 -p %3 %4").arg(getAppPath().toStdString().c_str()).arg(walker->tilename.toStdString().c_str()).arg(P.channel).arg(P.bkg_thresh); break;
//            case 1: cmd_method = QString("%1/vaa3d -x neuTube -f neutube_trace -i %2 -p %3 1").arg(getAppPath().toStdString().c_str()).arg(walker->tilename.toStdString().c_str()).arg(P.channel); break;
//            case 2: cmd_method = QString("%1/vaa3d -x snake -f snake_trace -i %2 -p %3").arg(getAppPath().toStdString().c_str()).arg(walker->tilename.toStdString().c_str()).arg(P.channel); break;
//             }
//            system(qPrintable(cmd_method));


//        #elif defined(Q_OS_MAC)
//            QString cmd_method;
//            switch(P.tracing_method){
//            case 0: cmd_method = QString("%1/vaa3d64.app/Contents/MacOS/vaa3d64 -x mostVesselTracer -f MOST_trace -i %2 -p %3 %4").arg(getAppPath().toStdString().c_str()).arg(walker->tilename.toStdString().c_str()).arg(P.channel).arg(P.bkg_thresh); break;
//            case 1: cmd_method = QString("%1/vaa3d64.app/Contents/MacOS/vaa3d64 -x neuTube -f neutube_trace -i %2 -p %3 1").arg(getAppPath().toStdString().c_str()).arg(walker->tilename.toStdString().c_str()).arg(P.channel); break;
//            case 2: cmd_method = QString("%1/vaa3d64.app/Contents/MacOS/vaa3d64 -x snake -f snake_trace -i %2 -p %3").arg(getAppPath().toStdString().c_str()).arg(walker->tilename.toStdString().c_str()).arg(P.channel); break;
//             }
//            system(qPrintable(cmd_method));
//        #else
//                 v3d_msg("The OS is not Linux or Mac. Do nothing.",bmenu);
//                 return false;
//        #endif

        ifstream ifs_swc(swcfilename.toStdString().c_str());
        if(!ifs_swc)
        {
            walker = walker->next;
            continue;
        }

        NeuronTree nt = readSWC_file(swcfilename);
        if(nt.listNeuron.empty())
        {
            remove(swcfilename.toStdString().c_str());
            walker = walker->next;

            continue;
        }


        NeuronTree nt_pruned = eliminate(nt,3.0);

        struct root_node *walker_inside;
        struct root_node *newNode;
        walker_inside = head;
        while(walker_inside->next != NULL)
        {
             walker_inside = walker_inside->next;
        }


        bool left = 0, right = 0,up = 0,down = 0;
        for(V3DLONG i = 0; i < nt_pruned.listNeuron.size(); i++)
        {
            NeuronSWC curr = nt_pruned.listNeuron.at(i);
            newNode =  new root_node[1];
            if(curr.x < 0.05* P.block_size && left == 0 && walker->direction !=2)
            {

                if(walker->end[0] == vim.sz[0] - 1 || walker->start[0] == 0)       continue;
                newNode->start[0] = walker->start[0] - P.block_size;
                newNode->start[1] = walker->start[1];
                newNode->start[2] = walker->start[2];
                newNode->end[0] = newNode->start[0] + P.block_size - 1;;
                newNode->end[1] = walker->end[1];
                newNode->end[2] = walker->end[2];

                if( newNode->start[0] < 0)  newNode->start[0] = 0;
                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 2;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;

                left = 1;
            }
            else if(curr.x > 0.95 * P.block_size && right == 0 && walker->direction !=1)
            {
                if(walker->end[0] == vim.sz[0] - 1 || walker->start[0] == 0)       continue;

                newNode->start[0] = walker->start[0] + P.block_size;
                newNode->start[1] = walker->start[1];
                newNode->start[2] = walker->start[2];

                newNode->end[0] = newNode->start[0] + P.block_size - 1;
                newNode->end[1] = walker->end[1];
                newNode->end[2] = walker->end[2];

                if( newNode->end[0] > vim.sz[0] - 1)  newNode->end[0] = vim.sz[0] - 1;

                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 1;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;
                right = 1;


            }
            else if(curr.y < 0.05* P.block_size && up == 0 && walker->direction !=4)
            {
                if(walker->end[1] == vim.sz[1] - 1 || walker->start[1] == 0)       continue;

                newNode->start[0] = walker->start[0];
                newNode->start[1] = walker->start[1] - P.block_size;
                newNode->start[2] = walker->start[2];
                newNode->end[0] = walker->end[0];
                newNode->end[1] = newNode->start[1]+ P.block_size - 1;
                newNode->end[2] = walker->end[2];

                if(newNode->start[1] < 0) newNode->start[1] = 0;

                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 4;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;
                up = 1;

            }
            else if(curr.y > 0.95 * P.block_size && down == 0 && walker->direction !=3)
            {
                if(walker->end[1] == vim.sz[1] - 1 || walker->start[1] == 0)       continue;

                newNode->start[0] = walker->start[0];
                newNode->start[1] = walker->start[1]  + P.block_size;
                newNode->start[2] = walker->start[2];
                newNode->end[0] = walker->end[0];
                newNode->end[1] = newNode->start[1] + P.block_size - 1;
                newNode->end[2] = walker->end[2];

                if( newNode->end[1] > vim.sz[1] - 1)  newNode->end[1] = vim.sz[1] - 1;

                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 3;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;
                down = 1;
            }

        }

        vector<MyMarker*> temp_out_swc = readSWC_file(swcfilename.toStdString());

        ifstream ifs(finalswcfilename.toStdString().c_str());
        if(walker->ref_index == -1)
        {
            if(ifs)
               remove(finalswcfilename.toStdString().c_str());
            for(int j = 0; j < temp_out_swc.size(); j++)
            {
                temp_out_swc[j]->x = temp_out_swc[j]->x + walker->start[0];
                temp_out_swc[j]->y = temp_out_swc[j]->y + walker->start[1];
                temp_out_swc[j]->z = temp_out_swc[j]->z + walker->start[2];
                temp_out_swc[j]->type = swc_type;
            }

           saveSWC_file(finalswcfilename.toStdString(), temp_out_swc);
        }
        else
        {
            vector<MyMarker*> final_out_swc = readSWC_file(finalswcfilename.toStdString());
            vector<MyMarker*> final_out_swc_updated =  final_out_swc;
           // temp_out_swc[0]->parent = final_out_swc[walker->parent];

            for(int j = 0; j < temp_out_swc.size(); j++)
            {
                temp_out_swc[j]->x = temp_out_swc[j]->x + walker->start[0];
                temp_out_swc[j]->y = temp_out_swc[j]->y + walker->start[1];
                temp_out_swc[j]->z = temp_out_swc[j]->z + walker->start[2];
                temp_out_swc[j]->type = swc_type;
                final_out_swc_updated.push_back(temp_out_swc[j]);
            }
           saveSWC_file(finalswcfilename.toStdString(), final_out_swc_updated);

        }

        swc_type++;
        walker = walker->next;

    }

    v3d_msg(QString("Now you can drag and drop the generated swc fle [%1] into Vaa3D.").arg(finalswcfilename),bmenu);

    system(qPrintable(QString("rm -rf %1").arg(tmpfolder.toStdString().c_str())));
    return true;

}

bool assembler_raw(V3DPluginCallback2 &callback, QWidget *parent,NA_PARA &P,bool bmenu)
{
    QString fileOpenName = P.inimg_file;

    unsigned char * datald = 0;
    V3DLONG *in_zz = 0;
    V3DLONG *in_sz = 0;

    int datatype;
    if (!loadRawRegion(const_cast<char *>(P.inimg_file.toStdString().c_str()), datald, in_zz, in_sz,datatype,0,0,0,1,1,1))
    {
        return false;
    }

    if(datald) {delete []datald; datald = 0;}

    V3DLONG start[3], end[3];

    start[0] = P.root_1st[0] - int(P.block_size/2);
    start[1] = P.root_1st[1] - int(P.block_size/2);
    start[2] = 0;

    end[0] = start[0] + P.block_size - 1;
    end[1] = start[1] + P.block_size - 1;
    end[2] = in_zz[2] - 1;


    struct root_node *head = new root_node[1];
    struct root_node *walker;

    QString startingpos="", tmps;
    tmps.setNum(start[0]).prepend("x"); startingpos += tmps;
    tmps.setNum(start[1]).prepend("_y"); startingpos += tmps;
    QString region_name = startingpos + ".raw";

    head->root_x = int(P.block_size/2);
    head->root_y = int(P.block_size/2);
    head->root_z = P.root_1st[2] ;

    head->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
    head->ref_index = -1;

    head->start[0] = start[0];
    head->start[1] = start[1];
    head->start[2] = start[2];
    head->end[0] = end[0];
    head->end[1] = end[1];
    head->end[2] = end[2];

    head->parent = -1;
    head->direction = 0;
    head->next = NULL;
    walker = head;

    QString rootposstr="", tmps2;
    tmps2.setNum(int(P.root_1st[0]+0.5)).prepend("_x"); rootposstr += tmps2;
    tmps2.setNum(int(P.root_1st[1]+0.5)).prepend("_y"); rootposstr += tmps2;
    tmps2.setNum(int(P.root_1st[2]+0.5)).prepend("_z"); rootposstr += tmps2;

    QString finalswcfilename;

    switch(P.tracing_method){
        case 1: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_MOST.swc"; break;
        case 2: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_NeuTube.swc"; break;
        case 3: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_Snake.swc"; break;
        case 4: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_TReMap.swc"; break;
        case 5: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_APP1.swc"; break;
        case 6: finalswcfilename = fileOpenName + rootposstr + "_NeuronAssembler_APP2.swc"; break;

    }

    QString tmpfolder = QFileInfo(fileOpenName).path()+("/tmp");
    system(qPrintable(QString("mkdir %1").arg(tmpfolder.toStdString().c_str())));
    if(tmpfolder.isEmpty())
    {

        printf("Can not create a tmp folder!\n");
        return false;
    }


    int swc_type = 2;
    while(walker != NULL)
    {

        ifstream ifs_image(walker->tilename.toStdString().c_str());
        if(!ifs_image)
        {
            unsigned char * datald = 0;
            V3DLONG *in_zz = 0;
            V3DLONG *in_sz = 0;

            int datatype;
            if (!loadRawRegion(const_cast<char *>(P.inimg_file.toStdString().c_str()), datald, in_zz, in_sz,datatype,walker->start[0],walker->start[1],walker->start[2],walker->end[0]+1,walker->end[1]+1,walker->end[2]+1))
            {
                printf("can not load the region");
                if(datald) {delete []datald; datald = 0;}
                return false;
            }

            if(P.tracing_method == 2 || P.tracing_method ==3 || P.tracing_method ==6)
            {
                for(V3DLONG i = 0; i < in_sz[0]*in_sz[1]*in_sz[2]; i++)
                {
                    if(datald[i] < P.bkg_thresh)
                        datald[i] = 0;

                }
            }
            simple_saveimage_wrapper(callback, walker->tilename.toStdString().c_str(),  (unsigned char *)datald, in_sz, V3D_UINT8);
            if(datald) {delete []datald; datald = 0;}
        }
        else
        {
            walker = walker->next;
            continue;
        }

        QString swcfilename;
        switch(P.tracing_method){
        case 1: swcfilename =  walker->tilename + QString("_MOST.swc"); break;
        case 2: swcfilename =  walker->tilename + QString("_neutube.swc"); break;
        case 3: swcfilename =  walker->tilename + QString("_snake.swc"); break;
        case 4: swcfilename =  walker->tilename + QString("_XY_3D_TreMap.swc"); break;
        case 5: swcfilename =  walker->tilename + QString("_region_APP1.swc"); break;
        case 6: swcfilename =  walker->tilename + QString("_region_APP2.swc"); break;


        }

        V3DPluginArgItem arg;
        V3DPluginArgList input;
        V3DPluginArgList output;

        QString full_plugin_name;
        QString func_name;

        arg.type = "random";std::vector<char*> arg_input;
        std:: string fileName_Qstring(walker->tilename.toStdString());char* fileName_string =  new char[fileName_Qstring.length() + 1]; strcpy(fileName_string, fileName_Qstring.c_str());
        arg_input.push_back(fileName_string);
        arg.p = (void *) & arg_input; input<< arg;

        char* char_swcout =  new char[swcfilename.length() + 1];strcpy(char_swcout, swcfilename.toStdString().c_str());
        arg.type = "random";std::vector<char*> arg_output;arg_output.push_back(char_swcout); arg.p = (void *) & arg_output; output<< arg;


        char channel = '0' + P.channel;
        string T_background = boost::lexical_cast<string>(P.bkg_thresh);
        char* Th =  new char[T_background.length() + 1];
        strcpy(Th, T_background.c_str());

        //NeuTube
        string merge_string = boost::lexical_cast<string>(P.merge);char* merge_char =  new char[merge_string.length() + 1]; strcpy(merge_char, merge_string.c_str());

        //TReMap
        string mip_plane_string = boost::lexical_cast<string>(P.mip_plane);char* mip_plane_char =  new char[mip_plane_string.length() + 1]; strcpy(mip_plane_char, mip_plane_string.c_str());
        string b_256cube_string = boost::lexical_cast<string>(P.b_256cube);char* b_256cube_char =  new char[b_256cube_string.length() + 1]; strcpy(b_256cube_char, b_256cube_string.c_str());
        string is_gsdt_string = boost::lexical_cast<string>(P.is_gsdt);char* is_gsdt_char =  new char[is_gsdt_string.length() + 1]; strcpy(is_gsdt_char, is_gsdt_string.c_str());
        string is_break_accept_string = boost::lexical_cast<string>(P.is_break_accept);char* is_break_accept_char =  new char[is_break_accept_string.length() + 1]; strcpy(is_break_accept_char, is_break_accept_string.c_str());
        string length_thresh_string = boost::lexical_cast<string>(P.length_thresh);char* length_thresh_char =  new char[length_thresh_string.length() + 1]; strcpy(length_thresh_char, length_thresh_string.c_str());

        //MST
        string win_size_string = boost::lexical_cast<string>(P.win_size);char* win_size_char =  new char[win_size_string.length() + 1]; strcpy(win_size_char, win_size_string.c_str());

        //APP1
       // if(P.tracing_method == 5) channel = '0' + P.channel - 1;

        //APP2
        string region_number_string = boost::lexical_cast<string>(P.region_number);char* region_number_char =  new char[region_number_string.length() + 1]; strcpy(region_number_char, region_number_string.c_str());


        arg.type = "random";
        std::vector<char*> arg_para;

        switch(P.tracing_method)
        {
            case 1: arg_para.push_back(&channel);arg_para.push_back(Th);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "mostVesselTracer"; func_name = "MOST_trace";break;
            case 2: arg_para.push_back(&channel); arg_para.push_back(merge_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "neuTube"; func_name = "neutube_trace";break;
            case 3: arg_para.push_back(&channel);arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "snake"; func_name = "snake_trace";break;
            case 4: arg_para.push_back(mip_plane_char);arg_para.push_back(&channel);arg_para.push_back(Th);arg_para.push_back(b_256cube_char);arg_para.push_back(is_gsdt_char);arg_para.push_back(is_break_accept_char);arg_para.push_back(length_thresh_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "TReMap"; func_name = "trace_mip";break;
            case 5: arg_para.push_back(&channel);arg_para.push_back(Th);arg_para.push_back(b_256cube_char);arg_para.push_back(region_number_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "Region_Neuron2"; func_name = "trace_app1";break;
            case 6: arg_para.push_back(&channel);arg_para.push_back(Th);arg_para.push_back(b_256cube_char);arg_para.push_back(is_gsdt_char);arg_para.push_back(is_break_accept_char);arg_para.push_back(length_thresh_char);arg_para.push_back(region_number_char);
                    arg.p = (void *) & arg_para; input << arg;
                    full_plugin_name = "Region_Neuron2"; func_name = "trace_app2";break;
//            case 6: arg_para.push_back(&channel);arg_para.push_back(win_size_char);
//                    arg.p = (void *) & arg_para; input << arg;
//                    full_plugin_name = "MST_tracing"; func_name = "trace_mst";break;
        }

        if(!callback.callPluginFunc(full_plugin_name,func_name,input,output))
        {

             printf("Can not find the tracing plugin!\n");
             return false;
        }


        ifstream ifs_swc(swcfilename.toStdString().c_str());
        if(!ifs_swc)
        {
            walker = walker->next;
            continue;
        }

        NeuronTree nt = readSWC_file(swcfilename);
        if(nt.listNeuron.empty())
        {
            remove(swcfilename.toStdString().c_str());
            walker = walker->next;

            continue;
        }


        NeuronTree nt_pruned = eliminate(nt,3.0);

        struct root_node *walker_inside;
        struct root_node *newNode;
        walker_inside = head;
        while(walker_inside->next != NULL)
        {
             walker_inside = walker_inside->next;
        }


        bool left = 0, right = 0,up = 0,down = 0;
        for(V3DLONG i = 0; i < nt_pruned.listNeuron.size(); i++)
        {
            NeuronSWC curr = nt_pruned.listNeuron.at(i);
            newNode =  new root_node[1];
            if(curr.x < 0.05* P.block_size && left == 0 && walker->direction !=2)
            {

                if(walker->end[0] == in_zz[0] - 1|| walker->start[0] == 0)       continue;

                newNode->start[0] = walker->start[0] - P.block_size;
                newNode->start[1] = walker->start[1];
                newNode->start[2] = walker->start[2];
                newNode->end[0] = newNode->start[0] + P.block_size - 1;;
                newNode->end[1] = walker->end[1];
                newNode->end[2] = walker->end[2];

                if( newNode->start[0] < 0)  newNode->start[0] = 0;
                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 2;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;

                left = 1;
            }
            else if(curr.x > 0.95 * P.block_size && right == 0 && walker->direction !=1)
            {
                if(walker->end[0] == in_zz[0] - 1|| walker->start[0] == 0)      continue;

                newNode->start[0] = walker->start[0] + P.block_size;
                newNode->start[1] = walker->start[1];
                newNode->start[2] = walker->start[2];

                newNode->end[0] = newNode->start[0] + P.block_size - 1;
                newNode->end[1] = walker->end[1];
                newNode->end[2] = walker->end[2];

                if( newNode->end[0] > in_zz[0] - 1)  newNode->end[0] = in_zz[0] - 1;

                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 1;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;
                right = 1;


            }
            else if(curr.y < 0.05* P.block_size && up == 0 && walker->direction !=4)
            {
                if(walker->end[1] == in_zz[1] - 1 || walker->start[1] == 0)       continue;

                newNode->start[0] = walker->start[0];
                newNode->start[1] = walker->start[1] - P.block_size;
                newNode->start[2] = walker->start[2];
                newNode->end[0] = walker->end[0];
                newNode->end[1] = newNode->start[1]+ P.block_size - 1;
                newNode->end[2] = walker->end[2];

                if(newNode->start[1] < 0) newNode->start[1] = 0;

                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 4;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;
                up = 1;

            }
            else if(curr.y > 0.95 * P.block_size && down == 0 && walker->direction !=3)
            {
                if(walker->end[1] == in_zz[1] - 1 || walker->start[1] == 0)       continue;

                newNode->start[0] = walker->start[0];
                newNode->start[1] = walker->start[1]  + P.block_size;
                newNode->start[2] = walker->start[2];
                newNode->end[0] = walker->end[0];
                newNode->end[1] = newNode->start[1] + P.block_size - 1;
                newNode->end[2] = walker->end[2];

                if( newNode->end[1] > in_zz[1] - 1)  newNode->end[1] = in_zz[1] - 1;

                QString startingpos="", tmps;
                tmps.setNum(newNode->start[0]).prepend("x"); startingpos += tmps;
                tmps.setNum(newNode->start[1]).prepend("_y"); startingpos += tmps;
                QString region_name = startingpos + ".raw";

                newNode->tilename = QFileInfo(fileOpenName).path().append("/tmp/").append(QString(region_name));
                newNode->ref_index = walker->tc_index;

                newNode->direction = 3;
                newNode->next = NULL;
                walker_inside->next = newNode;
                walker_inside = walker_inside->next;
                down = 1;
            }

        }

        vector<MyMarker*> temp_out_swc = readSWC_file(swcfilename.toStdString());

        ifstream ifs(finalswcfilename.toStdString().c_str());
        if(walker->ref_index == -1)
        {
            if(ifs)
               remove(finalswcfilename.toStdString().c_str());
            for(int j = 0; j < temp_out_swc.size(); j++)
            {
                temp_out_swc[j]->x = temp_out_swc[j]->x + walker->start[0];
                temp_out_swc[j]->y = temp_out_swc[j]->y + walker->start[1];
                temp_out_swc[j]->z = temp_out_swc[j]->z + walker->start[2];
                temp_out_swc[j]->type = swc_type;
            }

           saveSWC_file(finalswcfilename.toStdString(), temp_out_swc);
        }
        else
        {
            vector<MyMarker*> final_out_swc = readSWC_file(finalswcfilename.toStdString());
            vector<MyMarker*> final_out_swc_updated =  final_out_swc;
           // temp_out_swc[0]->parent = final_out_swc[walker->parent];

            for(int j = 0; j < temp_out_swc.size(); j++)
            {
                temp_out_swc[j]->x = temp_out_swc[j]->x + walker->start[0];
                temp_out_swc[j]->y = temp_out_swc[j]->y + walker->start[1];
                temp_out_swc[j]->z = temp_out_swc[j]->z + walker->start[2];
                temp_out_swc[j]->type = swc_type;
                final_out_swc_updated.push_back(temp_out_swc[j]);
            }
           saveSWC_file(finalswcfilename.toStdString(), final_out_swc_updated);

        }

        swc_type++;
        walker = walker->next;

    }


    v3d_msg(QString("Now you can drag and drop the generated swc fle [%1] into Vaa3D.").arg(finalswcfilename),bmenu);

    system(qPrintable(QString("rm -rf %1").arg(tmpfolder.toStdString().c_str())));

    return true;

}
void save_region(V3DPluginCallback2 &callback, V3DLONG *start, V3DLONG *end, QString &tcfile, QString &region_name,
                 Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim, NA_PARA &P)
{

    //virtual image
    V3DLONG vx, vy, vz, vc;

    vx = end[0] - start[0] + 1;
    vy = end[1] - start[1] + 1;
    vz = end[2] - start[2] + 1;
    vc = vim.sz[3];

    V3DLONG pagesz_vim = vx*vy*vz*vc;

    // flu bird algorithm
    bitset<3> lut_ss, lut_se, lut_es, lut_ee;

    //
    V3DLONG x_s = start[0] + vim.min_vim[0];
    V3DLONG y_s = start[1] + vim.min_vim[1];
    V3DLONG z_s = start[2] + vim.min_vim[2];

    V3DLONG x_e = end[0] + vim.min_vim[0];
    V3DLONG y_e = end[1] + vim.min_vim[1];
    V3DLONG z_e = end[2] + vim.min_vim[2];
    printf("%d, %d, ,%d, %d, %d, %d\n\n\n\n\n",x_s, y_s, z_s, x_e,y_e,z_e);

    bool flag_init = true;
    unsigned char *pVImg_UINT8 = NULL;

    QString curFilePath = QFileInfo(tcfile).path();
    curFilePath.append("/");

    for(V3DLONG ii=0; ii<vim.number_tiles; ii++)
    {
        // init
        lut_ss.reset();
        lut_se.reset();
        lut_es.reset();
        lut_ee.reset();

        //
        if(x_s < vim.lut[ii].start_pos[0]) lut_ss[1] = 1; // r  0 l
        if(y_s < vim.lut[ii].start_pos[1]) lut_ss[0] = 1; // d  0 u
        if(z_s < vim.lut[ii].start_pos[2]) lut_ss[2] = 1; // b  0 f

        if(x_e < vim.lut[ii].start_pos[0]) lut_se[1] = 1; // r  0 l
        if(y_e < vim.lut[ii].start_pos[1]) lut_se[0] = 1; // d  0 u
        if(z_e < vim.lut[ii].start_pos[2]) lut_se[2] = 1; // b  0 f

        if(x_s < vim.lut[ii].end_pos[0]) lut_es[1] = 1; // r  0 l
        if(y_s < vim.lut[ii].end_pos[1]) lut_es[0] = 1; // d  0 u
        if(z_s < vim.lut[ii].end_pos[2]) lut_es[2] = 1; // b  0 f

        if(x_e < vim.lut[ii].end_pos[0]) lut_ee[1] = 1; // r  0 l
        if(y_e < vim.lut[ii].end_pos[1]) lut_ee[0] = 1; // d  0 u
        if(z_e < vim.lut[ii].end_pos[2]) lut_ee[2] = 1; // b  0 f

        // copy data
        if( (!lut_ss.any() && lut_ee.any()) || (lut_es.any() && !lut_ee.any()) || (lut_ss.any() && !lut_se.any()) )
        {
            //
            cout << "satisfied image: "<< vim.lut[ii].fn_img << endl;

            // loading relative image files
            V3DLONG *sz_relative = 0;
            int datatype_relative = 0;
            unsigned char* relative1d = 0;

            QString curPath = curFilePath;

            string fn = curPath.append( QString(vim.lut[ii].fn_img.c_str()) ).toStdString();

            qDebug()<<"testing..."<<curFilePath<< fn.c_str();


            if(flag_init)
            {
                try
                {
                    pVImg_UINT8 = new unsigned char [pagesz_vim];
                }
                catch (...)
                {
                    printf("Fail to allocate memory.\n");
                    return;
                }

                // init
                memset(pVImg_UINT8, 0, pagesz_vim*sizeof(unsigned char));

                flag_init = false;

            }

            if (loadImage(const_cast<char *>(fn.c_str()), relative1d, sz_relative, datatype_relative)!=true)
            {
                fprintf (stderr, "Error happens in reading the subject file [%s]. Exit. \n",vim.tilesList.at(ii).fn_image.c_str());
                continue;
            }
            V3DLONG rx=sz_relative[0], ry=sz_relative[1], rz=sz_relative[2], rc=sz_relative[3];


            //
            V3DLONG tile2vi_xs = vim.lut[ii].start_pos[0]-vim.min_vim[0];
            V3DLONG tile2vi_xe = vim.lut[ii].end_pos[0]-vim.min_vim[0];
            V3DLONG tile2vi_ys = vim.lut[ii].start_pos[1]-vim.min_vim[1];
            V3DLONG tile2vi_ye = vim.lut[ii].end_pos[1]-vim.min_vim[1];
            V3DLONG tile2vi_zs = vim.lut[ii].start_pos[2]-vim.min_vim[2];
            V3DLONG tile2vi_ze = vim.lut[ii].end_pos[2]-vim.min_vim[2];

            V3DLONG x_start = (start[0] > tile2vi_xs) ? start[0] : tile2vi_xs;
            V3DLONG x_end = (end[0] < tile2vi_xe) ? end[0] : tile2vi_xe;
            V3DLONG y_start = (start[1] > tile2vi_ys) ? start[1] : tile2vi_ys;
            V3DLONG y_end = (end[1] < tile2vi_ye) ? end[1] : tile2vi_ye;
            V3DLONG z_start = (start[2] > tile2vi_zs) ? start[2] : tile2vi_zs;
            V3DLONG z_end = (end[2] < tile2vi_ze) ? end[2] : tile2vi_ze;

            x_end++;
            y_end++;
            z_end++;

            //
            cout << x_start << " " << x_end << " " << y_start << " " << y_end << " " << z_start << " " << z_end << endl;
            region_groupfusing<unsigned char>(pVImg_UINT8, vim, relative1d,
                                              vx, vy, vz, vc, rx, ry, rz, rc,
                                              tile2vi_zs, tile2vi_ys, tile2vi_xs,
                                              z_start, z_end, y_start, y_end, x_start, x_end, start);

            //de-alloc
            if(relative1d) {delete []relative1d; relative1d=0;}
            if(sz_relative) {delete []sz_relative; sz_relative=0;}
        }

    }

    V3DLONG in_sz[4];
    in_sz[0] = vx;
    in_sz[1] = vy;
    in_sz[2] = vz;
    in_sz[3] = vc;

    if(P.tracing_method == 2 || P.tracing_method ==3 || P.tracing_method ==6)
    {
        for(V3DLONG i = 0; i < vx*vy*vz*vc; i++)
        {
            if(pVImg_UINT8[i] < P.bkg_thresh)
                pVImg_UINT8[i] = 0;
        }
    }

    simple_saveimage_wrapper(callback, region_name.toStdString().c_str(),  (unsigned char *)pVImg_UINT8, in_sz, V3D_UINT8);
    if(pVImg_UINT8) {delete []pVImg_UINT8; pVImg_UINT8=0;}

    return;
}

NeuronTree eliminate(NeuronTree input, double length)
{
    NeuronTree result;
    V3DLONG siz = input.listNeuron.size();
    Tree tree;
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = input.listNeuron[i];
        Point* pt = new Point;
        pt->x = s.x;
        pt->y = s.y;
        pt->z = s.z;
        pt->r = s.r;
        pt ->type = s.type;
        pt->p = NULL;
        pt->childNum = 0;
        tree.push_back(pt);
    }
    for (V3DLONG i=0;i<siz;i++)
    {
        if (input.listNeuron[i].pn<0) continue;
        V3DLONG pid = input.hashNeuron.value(input.listNeuron[i].pn);
        tree[i]->p = tree[pid];
        tree[pid]->childNum++;
    }
//	printf("tree constructed.\n");
    vector<Segment*> seg_list;
    for (V3DLONG i=0;i<siz;i++)
    {
        if (tree[i]->childNum!=1)//tip or branch point
        {
            Segment* seg = new Segment;
            Point* cur = tree[i];
            do
            {
                seg->push_back(cur);
                cur = cur->p;
            }
            while(cur && cur->childNum==1);
            seg_list.push_back(seg);
        }
    }


    tree.clear();
    map<Point*, V3DLONG> index_map;
    for (V3DLONG i=0;i<seg_list.size();i++)
    {
        if(seg_list[i]->size()>length)// || space_distace > length)
        {
            for (V3DLONG j=0;j<seg_list[i]->size();j++)
            {
                tree.push_back(seg_list[i]->at(j));
                index_map.insert(pair<Point*, V3DLONG>(seg_list[i]->at(j), tree.size()-1));
            }
        }
    }
    for (V3DLONG i=0;i<tree.size();i++)
    {
        NeuronSWC S;
        Point* p = tree[i];
        S.n = i+1;
        if (p->p==NULL || index_map[p->p] ==0) S.pn = -1;
        else
            S.pn = index_map[p->p]+1;
        if (p->p==p) printf("There is loop in the tree!\n");
        S.x = p->x;
        S.y = p->y;
        S.z = p->z;
        S.r = p->r;
        S.type = p->type;
        result.listNeuron.push_back(S);
    }
    for (V3DLONG i=0;i<tree.size();i++)
    {
        if (tree[i]) {delete tree[i]; tree[i]=NULL;}
    }
    for (V3DLONG j=0;j<seg_list.size();j++)
        if (seg_list[j]) {delete seg_list[j]; seg_list[j] = NULL;}
    for (V3DLONG i=0;i<result.listNeuron.size();i++)
    {
        result.hashNeuron.insert(result.listNeuron[i].n, i);
    }
    return result;
}

