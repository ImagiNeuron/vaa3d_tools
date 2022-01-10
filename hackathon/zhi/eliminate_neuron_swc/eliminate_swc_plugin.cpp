/* eliminate_swc_plugin.cpp
 * This plugin is used to eliminate small segments in swc file
 * 2013-12-16 : by Zhi Zhou
 */
 
#include "v3d_message.h"
#include <vector>
#include "eliminate_swc_plugin.h"
#include "basic_surf_objs.h"
#include <iostream>
#include <map>
#include "../../../released_plugins/v3d_plugins/istitch/y_imglib.h"
#include "my_surf_objs.h"
#include "../../../released_plugins/v3d_plugins/bigneuron_zz_neurontracing_TReMAP/smooth_curve.h"
#include "../../../released_plugins/v3d_plugins/neuron_radius/hierarchy_prune.h"



using namespace std;

#define b_VERBOSE_PRINT 0
#define getParent(n,nt) ((nt).listNeuron.at(n).pn<0)?(1000000000):((nt).hashNeuron.value((nt).listNeuron.at(n).pn))
#define DEFINE_NBYTE2G \
  V3DLONG nBytes2G = (V3DLONG(1024)*V3DLONG(1024)*V3DLONG(1024)-1)*V3DLONG(2);

typedef int BIT32_UNIT;


// Open a series of inputs
QStringList importSeriesFileList_addnumbersort(const QString & curFilePath)
{
    QStringList myList;
    myList.clear();

    // get the image files namelist in the directory
    QStringList imgSuffix;
    imgSuffix<<"*.tif"<<"*.raw"<<"*.v3draw"<<"*.lsm"
            <<"*.TIF"<<"*.RAW"<<"*.V3DRAW"<<"*.LSM";

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

QStringList importSWCFileList_addnumbersort(const QString & curFilePath)
{
    QStringList myList;
    myList.clear();

    // get the image files namelist in the directory
    QStringList imgSuffix;
    imgSuffix<<"*.swc"<<"*.eswc"<<"*.SWC"<<"*.ESWC";

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

bool saveSWC_file_TreMap(string swc_file, vector<MyMarker*> & outmarkers)
{
    if(swc_file.find_last_of(".dot") == swc_file.size() - 1) return saveDot_file(swc_file, outmarkers);

    cout<<"marker num = "<<outmarkers.size()<<", save swc file to "<<swc_file<<endl;
    map<MyMarker*, int> ind;
    ofstream ofs(swc_file.c_str());

    if(ofs.fail())
    {
        cout<<"open swc file error"<<endl;
        return false;
    }
    ofs<<"#name "<<"TreMap_Tracing"<<endl;
    ofs<<"#comment "<<endl;

    ofs<<"##n,type,x,y,z,radius,parent"<<endl;
    for(int i = 0; i < outmarkers.size(); i++) ind[outmarkers[i]] = i+1;

    for(int i = 0; i < outmarkers.size(); i++)
    {
        MyMarker * marker = outmarkers[i];
        int parent_id;
        if(marker->parent == 0) parent_id = -1;
        else parent_id = ind[marker->parent];
        if(parent_id == 0)  parent_id = -1;
        ofs<<i+1<<" "<<marker->type<<" "<<marker->x<<" "<<marker->y<<" "<<marker->z<<" "<<marker->radius<<" "<<parent_id<<endl;
    }
    ofs.close();
    return true;
}


Q_EXPORT_PLUGIN2(eliminate_swc, eliminate_swc);
 
QStringList eliminate_swc::menulist() const
{
	return QStringList() 
        <<tr("eliminate_swc")
        <<tr("smooth_swc")
        <<tr("combine_swc_group")
        <<tr("combine_swc_pair")
        <<tr("align_swc")
        <<tr("prun_swc")
        <<tr("z_section to tiles")
        <<tr("generate qsub list")
        <<tr("3D images to tiles")
		<<tr("about");
}

QStringList eliminate_swc::funclist() const
{
	return QStringList()
        <<tr("smooth_swc")
        <<tr("combine_swc")
		<<tr("help");
}

NeuronTree eliminate(NeuronTree input, double length);
NeuronTree eliminate_overlap(NeuronTree target, NeuronTree subject, double length);
void combineSWC_group(V3DPluginCallback2 &callback, QWidget *parent);
void combineSWC_group_dupcheck(V3DPluginCallback2 &callback, QWidget *parent);
void combineSWC_pair(V3DPluginCallback2 &callback, QWidget *parent);
void prunSWC(V3DPluginCallback2 &callback, QWidget *parent);
void zsectionsTotiles(V3DPluginCallback2 &callback, QWidget *parent);
void qsublist(V3DPluginCallback2 &callback, QWidget *parent);
void threeDimageTotiles(V3DPluginCallback2 &callback, QWidget *parent);
NeuronTree smooth_swc(NeuronTree input, double length);



char checkMachineEndian();
void swap2bytes(void *targetp);
void swap4bytes(void *targetp);


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

void extend_path(Segment * seg);
bool export_list2file(QList<NeuronSWC> & lN, QString fileSaveName, QString fileOpenName)
{
    QFile file(fileSaveName);
    if (!file.open(QIODevice::WriteOnly|QIODevice::Text))
        return false;
    QTextStream myfile(&file);
    myfile<<"# generated by Vaa3D Plugin eliminate_swc"<<endl;
    myfile<<"# source file(s): "<<fileOpenName<<endl;
    myfile<<"# id,type,x,y,z,r,pid"<<endl;
    for (V3DLONG i=0;i<lN.size();i++)
        myfile << lN.at(i).n <<" " << lN.at(i).type << " "<< lN.at(i).x <<" "<<lN.at(i).y << " "<< lN.at(i).z << " "<< lN.at(i).r << " " <<lN.at(i).pn << "\n";

    file.close();
    cout<<"swc file "<<fileSaveName.toStdString()<<" has been generated, size: "<<lN.size()<<endl;
    return true;
};

int loadRaw2Stack_2byte(char * filename, unsigned char * & img, V3DLONG * & sz, int & datatype);

void eliminate_swc::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
    if (menu_name == tr("eliminate_swc"))
	{
        QString fileOpenName;
        fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Open File"),
                "",
                QObject::tr("Supported file (*.swc *.eswc)"
                    ";;Neuron structure	(*.swc)"
                    ";;Extended neuron structure (*.eswc)"
                    ));
        if(fileOpenName.isEmpty())
            return;
        double length = 0;
        NeuronTree nt;
        if (fileOpenName.toUpper().endsWith(".SWC") || fileOpenName.toUpper().endsWith(".ESWC"))
        {
            bool ok;
            nt = readSWC_file(fileOpenName);
            length = QInputDialog::getDouble(parent, "Please specify the minimum segment length","segment length:",1,0,2147483647,0.1,&ok);
            if (!ok)
                return;
        }

        NeuronTree result = eliminate(nt, length);

        QString fileDefaultName = fileOpenName+QString("thresholded.swc");
        //write new SWC to file
        QString fileSaveName = QFileDialog::getSaveFileName(0, QObject::tr("Save File"),
                fileDefaultName,
                QObject::tr("Supported file (*.swc)"
                    ";;Neuron structure	(*.swc)"
                    ));
        if (!export_list2file(result.listNeuron,fileSaveName,fileOpenName))
        {
            v3d_msg("fail to write the output swc file.");
            return;
        }


	}
    else if (menu_name == tr("smooth_swc"))
    {
        QString fileOpenName;
        fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Open File"),
                "",
                QObject::tr("Supported file (*.swc *.eswc)"
                    ";;Neuron structure	(*.swc)"
                    ";;Extended neuron structure (*.eswc)"
                    ));
        if(fileOpenName.isEmpty())
            return;
        double length = 0;
        vector<MyMarker*> inswc;
        if (fileOpenName.toUpper().endsWith(".SWC") || fileOpenName.toUpper().endsWith(".ESWC"))
        {
            bool ok;
            inswc = readSWC_file(fileOpenName.toStdString());

            length = QInputDialog::getDouble(parent, "Please specify the smooth step size","step size:",1,0,2147483647,0.1,&ok);
            if (!ok)
                return;
        }

        unsigned char* inimg1d = 0;
        vector<HierarchySegment*> topo_segs;
        swc2topo_segs(inswc, topo_segs, 1, inimg1d, 0, 0, 0);

        cout<<"Smooth the final curve"<<endl;
        for(int i = 0; i < topo_segs.size(); i++)
        {
            HierarchySegment * seg = topo_segs[i];
            MyMarker * leaf_marker = seg->leaf_marker;
            MyMarker * root_marker = seg->root_marker;
            vector<MyMarker*> seg_markers;
            MyMarker * p = leaf_marker;
            while(p != root_marker)
            {
                seg_markers.push_back(p);
                p = p->parent;
            }
            seg_markers.push_back(root_marker);
            smooth_curve_Z(seg_markers, length);
        }
        inswc.clear();

        QString outswc_file = fileOpenName + "_smoothed.swc";
        topo_segs2swc(topo_segs, inswc, 0);
        saveSWC_file(outswc_file.toStdString(), inswc);
        for(int i = 0; i < inswc.size(); i++) delete inswc[i];

    }
    else if(menu_name == tr("combine_swc_group"))
    {
         combineSWC_group_dupcheck(callback, parent);
    }
    else if(menu_name == tr("combine_swc_pair"))
    {
        combineSWC_pair(callback, parent);

    }
    else if(menu_name == tr("prun_swc"))
    {
        prunSWC(callback, parent);

    }
    else if(menu_name == tr("z_section to tiles"))
    {
        zsectionsTotiles(callback,parent);
    }
    else if(menu_name == tr("generate qsub list"))
    {
        qsublist(callback,parent);
    }
    else if(menu_name == tr("3D images to tiles"))
    {
        threeDimageTotiles(callback,parent);
    }
    else if(menu_name == tr("align_swc"))
    {
        QString m_InputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory including all swc files "),
                                              QDir::currentPath(),
                                              QFileDialog::ShowDirsOnly);

        QStringList swcList = importSWCFileList_addnumbersort(m_InputfolderName);

        double xoriginal = 0;
        double yoriginal = 0;
        double a11, a12, a21,a22,xshift,yshift;

        for(int i = 0; i < swcList.size(); i++)
        {
            QString fileOpenName = swcList.at(i);
            NeuronTree nt = readSWC_file(fileOpenName);

            switch (i)
            {
            case 0: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 2238 - xoriginal,yshift = 1226 - yoriginal; break;
            case 1: a11 = 0.9969223686848976; a12 = -0.07839509433435808; a21 = 0.07839509433435808; a22 = 0.9969223686848976; xshift = -511.4940779028416 - xoriginal,yshift = 241.60709307748505- yoriginal; break;
            case 2: a11 = 0.9992834654465653; a12 = 0.03784911736232148; a21 = -0.03784911736232148; a22 = 0.9992834654465653; xshift = 3625.9847459171847 - xoriginal,yshift = 1801.3516938690236- yoriginal; break;

            case 3: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 3712 - xoriginal,yshift = 1184 - yoriginal; break;
            case 4: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 3917; yshift = 692; break;

            case 5: a11 = 0.9997030084369902; a12 = 0.024369959417943652; a21 = -0.024369959417943652; a22 = 0.9997030084369902; xshift = 8011.313822690522; yshift = 998.8719300918069; break;
            case 6: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 13612 - xoriginal,yshift = 2084 - yoriginal; break;
            case 7: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 8916 - xoriginal,yshift = 1025 - yoriginal; break;
            case 8: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 12708 - xoriginal,yshift = 162 - yoriginal; break;
            case 9: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 13638 - xoriginal,yshift = 162 - yoriginal; break;

            case 10: a11 = 0.999966760546087; a12 = -0.00815339211405012; a21 = 0.00815339211405012; a22 = 0.999966760546087; xshift = 3473.58427039758; yshift = 3087.6355858839374; break;
            case 11: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 8386 - xoriginal,yshift = 2229 - yoriginal; break;

            case 12: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 11752; yshift = 2048; break;
            case 13: a11 = 1; a12 = 0; a21 = 0; a22 = 1; xshift = 1086 - xoriginal,yshift = 912 - yoriginal; break;

            }

            NeuronTree nt_aligned;

            QList <NeuronSWC> listNeuron;
            QHash <int, int>  hashNeuron;
            listNeuron.clear();
            hashNeuron.clear();

            QList<NeuronSWC> list = nt.listNeuron;
            NeuronSWC S;
            for (int i=0;i<list.size();i++)
            {
                NeuronSWC curr = list.at(i);
                float x_old = curr.x;
                float y_old = curr.y;
                S.x = x_old*a11 + y_old*a21 + xshift;
                S.y = x_old*a12 + y_old*a22 + yshift;
                S.z = curr.z;
                S.n 	= curr.n;
                S.type 	= curr.type;
                S.r 	= curr.r;
                S.pn 	= curr.pn;
                listNeuron.append(S);
                hashNeuron.insert(S.n, listNeuron.size()-1);
            }

            nt_aligned.n = -1;
            nt_aligned.on = true;
            nt_aligned.listNeuron = listNeuron;
            nt_aligned.hashNeuron = hashNeuron;

            float min = 2000000;
            float max = -2000000;
            float offsecZ;
            switch (i)
            {
                case 0: offsecZ = -1144; break;
                case 1: offsecZ = -858; break;
                case 2: offsecZ = -286; break;
                case 3: offsecZ = 0; break;
                case 4: offsecZ = 286; break;
                case 5: offsecZ = 572; break;
                case 6: offsecZ = 572; break;
                case 7: offsecZ = 572; break;
                case 8: offsecZ = 572; break;
                case 9: offsecZ = 572; break;
                case 10: offsecZ = 858; break;
                case 11: offsecZ = 858; break;
                case 12: offsecZ = 858; break;
                case 13: offsecZ = 1144; break;

            }


            NeuronTree nt_norm;
            listNeuron.clear();
            hashNeuron.clear();

            list = nt_aligned.listNeuron;

            for (int i=0;i<list.size();i++)
            {
                NeuronSWC curr = list.at(i);
                if(curr.z > max ) max = curr.z;
                if(curr.z < min ) min = curr.z;
            }

           // NeuronSWC S;
            for (int i=0;i<list.size();i++)
            {
                NeuronSWC curr = list.at(i);
                S.x = curr.x;
                S.y = curr.y;
                S.z = offsecZ + (285 * (curr.z - min) /(max - min));
                S.n 	= curr.n;
                S.type 	= curr.type;
                S.r 	= curr.r;
                S.pn 	= curr.pn;
                listNeuron.append(S);
                hashNeuron.insert(S.n, listNeuron.size()-1);
            }


            nt_norm.n = -1;
            nt_norm.on = true;
            nt_norm.listNeuron = listNeuron;
            nt_norm.hashNeuron = hashNeuron;

            //write new SWC to file
            QString fileSaveName =  fileOpenName+QString("_aligned.swc");

            if (!export_list2file(nt_norm.listNeuron,fileSaveName,fileOpenName))
            {
                v3d_msg("fail to write the output swc file.");
                return;
            }
        }
    }
    else
	{
		v3d_msg(tr("This plugin is used to eliminate small segments in swc file. "
			"Developed by Zhi Zhou, 2013-12-16"));
	}
}

bool eliminate_swc::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
    if (func_name == tr("smooth_swc"))
	{
        cout<<"Welcome to resample_swc"<<endl;
        vector<char*>* inlist = (vector<char*>*)(input.at(0).p);
        vector<char*>* outlist = NULL;
        vector<char*>* paralist = NULL;

        if(input.size() != 2)
        {
            printf("Please specify both input file and step length parameter.\n");
            return false;
        }
        paralist = (vector<char*>*)(input.at(1).p);
        if (paralist->size()!=1)
        {
            printf("Please specify only two parameters.\n");
            return false;
        }
        double length = atof(paralist->at(0));

        QString fileOpenName = QString(inlist->at(0));
        QString fileSaveName;
        if (output.size()==0)
        {
            printf("No outputfile specified.\n");
            fileSaveName = fileOpenName + QString("_Z_T%1.swc").arg(length);
        }
        else if (output.size()==1)
        {
            outlist = (vector<char*>*)(output.at(0).p);
            fileSaveName = QString(outlist->at(0));
        }
        else
        {
            printf("You have specified more than 1 output file.\n");
            return false;
        }

        vector<MyMarker*> inswc;
        if (fileOpenName.toUpper().endsWith(".SWC") || fileOpenName.toUpper().endsWith(".ESWC"))
            inswc = readSWC_file(fileOpenName.toStdString());

        unsigned char* inimg1d = 0;
        vector<HierarchySegment*> topo_segs;
        swc2topo_segs(inswc, topo_segs, 1, inimg1d, 0, 0, 0);

        cout<<"Smooth the final curve"<<endl;
        for(int i = 0; i < topo_segs.size(); i++)
        {
            HierarchySegment * seg = topo_segs[i];
            MyMarker * leaf_marker = seg->leaf_marker;
            MyMarker * root_marker = seg->root_marker;
            vector<MyMarker*> seg_markers;
            MyMarker * p = leaf_marker;
            while(p != root_marker)
            {
                seg_markers.push_back(p);
                p = p->parent;
            }
            seg_markers.push_back(root_marker);
            smooth_curve_Z(seg_markers, length);
           // smooth_curve_XY(seg_markers, length);
        }
        inswc.clear();

        topo_segs2swc(topo_segs, inswc, 0);
        saveSWC_file(fileSaveName.toStdString(), inswc);
        for(int i = 0; i < inswc.size(); i++) delete inswc[i];

        return true;
    }else if (func_name == tr("combine_swc"))
    {
        vector<char*>* inlist = (vector<char*>*)(input.at(0).p);
        vector<char*>* outlist = NULL;
        QString targetname = QString(inlist->at(0));
        outlist = (vector<char*>*)(output.at(0).p);
        QString subjectname = QString(outlist->at(0));

        NeuronTree nt_target = readSWC_file(targetname);
        NeuronTree nt_subject = readSWC_file(subjectname);

        NeuronTree result = nt_subject;
        V3DLONG target_size = nt_target.listNeuron.size();
        for(V3DLONG i = 0; i < result.listNeuron.size();i++)
        {
           NeuronSWC S = result.listNeuron[i];
           int flag_prun = 0;
           for(int jj = 0; jj < target_size;jj++)
           {
               NeuronSWC S2 = nt_target.listNeuron[jj];
               int dis_prun = sqrt(pow(S.x - S2.x,2) + pow(S.y - S2.y,2) + pow(S.z-S2.z,2));
               if( dis_prun < 10)
               {
                   flag_prun = 1;
                   break;
               }

           }
           if(flag_prun == 0)
           {
               // S.parent = -1;//S.parent + target_size;
               // S.n = S.n + target_size;
               S.parent = -1;
               S.n = S.n + target_size;
               nt_target.listNeuron.push_back(S);
           }

        }

        QString fileDefaultName = targetname+QString("_combined.swc");
        export_list2file(nt_target.listNeuron,fileDefaultName,fileDefaultName);
    }
	else if (func_name == tr("help"))
	{
		v3d_msg("To be implemented.");
	}
	else return false;

	return true;
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
//	printf("segment list constructed.\n");
  /*  for (V3DLONG i=0;i<seg_list.size();i++)
    {
        if(seg_list[i]->size() > length)
            extend_path(seg_list[i]);
    }*/

//	printf("resample done.\n");
    tree.clear();
    map<Point*, V3DLONG> index_map;
    for (V3DLONG i=0;i<seg_list.size();i++)
    {
        //double space_distace = sqrt(pow(seg_list[i]->at(0)-seg_list[i]->at((seg_list[i]->size())),2.0));
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

void extend_path(Segment * seg)
{

    Point* last_pt = seg->at(0);
    Point* last_pt2 = seg->at(1);
    Segment seg_r = *seg;
    if(last_pt->z > 45.0)
    {
        for(int i = 0; i <200;i++)
        {
           Point* pt = new Point;;
           pt->x = 2*last_pt->x - last_pt2->x;
           pt->y = 2*last_pt->y - last_pt2->y;
           pt->z = 2*last_pt->z - last_pt2->z;
           pt ->type = 3;
           pt->p = last_pt;
           pt->r = last_pt->r;
           seg_r.push_back(pt);
           last_pt2 = last_pt;
           last_pt = pt;
        }
        printf("end point is %.2f\n",last_pt->z);
    }
    *seg = seg_r;
}

void combineSWC_group(V3DPluginCallback2 &callback, QWidget *parent)
{

    QString m_InputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory including all swc files "),
                                          QDir::currentPath(),
                                          QFileDialog::ShowDirsOnly);

    QStringList swcList = importSWCFileList_addnumbersort(m_InputfolderName);

    vector<MyMarker*> outswc;
    for(int i = 0; i < swcList.size(); i++)
    {

        QString curPathSWC = swcList.at(i);
        vector<MyMarker*> inputswc = readSWC_file(curPathSWC.toStdString());

        for(V3DLONG d = 0; d < inputswc.size(); d++)
        {
            outswc.push_back(inputswc[d]);
        }

    }
    QString swc_combined = m_InputfolderName + "/combined.swc";
    saveSWC_file(swc_combined.toStdString().c_str(), outswc);
}

//add by: Hanbo 2014/12/5
void combineSWC_group_dupcheck(V3DPluginCallback2 &callback, QWidget *parent)
{

    QString m_InputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory including all swc files "),
                                          QDir::currentPath(),
                                          QFileDialog::ShowDirsOnly);

    QStringList swcList = importSWCFileList_addnumbersort(m_InputfolderName);

    bool ok;
    double disthr=QInputDialog::getDouble(parent,"group combine swc",QString::number(swcList.size()) + " files found. Please specificy a distance threshold for checking duplicate branches.",50,0,2147483647,1,&ok);
    if(!ok)
        return;
    int count = 0;
    vector<MyMarker*> outswc;
    multimap<int, vector<MyMarker*>, std::greater<int> > allinswc;
    multimap<int, QString, std::greater<int> > allinname;
    for(int i = 0; i < swcList.size(); i++)
    {
        QString curPathSWC = swcList.at(i);
        vector<MyMarker*> inputswc = readSWC_file(curPathSWC.toStdString());

        allinswc.insert(pair<int, vector<MyMarker*> >(inputswc.size(), inputswc));
        allinname.insert(pair<int, QString >(inputswc.size(), curPathSWC));
    }

    multimap<int, vector<MyMarker*>, std::greater<int> >::iterator iter_swc;
    multimap<int, QString, std::greater<int> >::iterator iter_name=allinname.begin();
    for(iter_swc = allinswc.begin(); iter_swc!=allinswc.end(); iter_swc++)
    {
        double dis = 1e16;
        if(outswc.size()>0){
            dis=getHDisBetweenTwoMarkers(iter_swc->second, outswc);
            printf("file:%s; distance:%f\n",iter_name->second.toStdString().c_str(),dis);
        }

        if(dis>disthr){
            count++;
            for(V3DLONG d = 0; d < iter_swc->second.size(); d++){
                outswc.push_back(iter_swc->second[d]);
            }
        }

        iter_name++;
    }


    QString swc_combined = m_InputfolderName + "/combined.swc";
    saveSWC_file(swc_combined.toStdString().c_str(), outswc);
    v3d_msg(QString::number(count) + " out of " + QString::number(swcList.size()) + " swc files were combined to combined.swc");
}

void combineSWC_pair(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString targetname;
    targetname = QFileDialog::getOpenFileName(0, QObject::tr("Open Target File"),
            "",
            QObject::tr("Supported file (*.swc *.eswc)"
                ";;Neuron structure	(*.swc)"
                ";;Extended neuron structure (*.eswc)"
                ));
    if(targetname.isEmpty())
        return;
    NeuronTree nt_target;
    if (targetname.toUpper().endsWith(".SWC") || targetname.toUpper().endsWith(".ESWC"))
    {
        nt_target = readSWC_file(targetname);

    }

    QString subjectname;
    subjectname = QFileDialog::getOpenFileName(0, QObject::tr("Open Subject File"),
            "",
            QObject::tr("Supported file (*.swc *.eswc)"
                ";;Neuron structure	(*.swc)"
                ";;Extended neuron structure (*.eswc)"
                ));
    if(subjectname.isEmpty())
        return;
    NeuronTree nt_subject;
    if (subjectname.toUpper().endsWith(".SWC") || subjectname.toUpper().endsWith(".ESWC"))
    {
        nt_subject = readSWC_file(subjectname);
    }

  //  NeuronTree result = eliminate_overlap(nt_target,nt_subject,5.0);
    NeuronTree result = nt_subject;
    V3DLONG target_size = nt_target.listNeuron.size();
    for(V3DLONG i = 0; i < result.listNeuron.size();i++)
    {
      //  V3DLONG target_size = nt_target.listNeuron.size();

        NeuronSWC S = result.listNeuron[i];
//        if(S.parent >0) S.parent = S.parent + target_size;
//        S.n = S.n + target_size;
//        nt_target.listNeuron.push_back(S);

        int flag_prun = 0;
        for(int jj = 0; jj < target_size;jj++)
        {
             NeuronSWC S2 = nt_target.listNeuron[jj];
            int dis_prun = sqrt(pow(S.x - S2.x,2) + pow(S.y - S2.y,2) + pow(S.z-S2.z,2));
            if( dis_prun < 10)
            {
                flag_prun = 1;
                break;
            }

        }
        if(flag_prun == 0)
        {
           // S.parent = -1;//S.parent + target_size;
           // S.n = S.n + target_size;
            S.parent = -1;
            S.n = S.n + target_size;
            nt_target.listNeuron.push_back(S);
        }

    }

    QString fileDefaultName = targetname+QString("_combined.swc");
    //write new SWC to file
    QString fileSaveName = QFileDialog::getSaveFileName(0, QObject::tr("Save File"),
            fileDefaultName,
            QObject::tr("Supported file (*.swc)"
                ";;Neuron structure	(*.swc)"
                ));
    if (!export_list2file(nt_target.listNeuron,fileSaveName,fileDefaultName))
    {
        v3d_msg("fail to write the output swc file.");
        return;
    }

}

NeuronTree eliminate_overlap(NeuronTree target_swc, NeuronTree subject_swc, double length)
{
    NeuronTree result;
    V3DLONG siz = subject_swc.listNeuron.size();
    Tree tree;
    for (V3DLONG i=0;i<siz;i++)
    {
        NeuronSWC s = subject_swc.listNeuron[i];
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
        if (subject_swc.listNeuron[i].pn<0) continue;
        V3DLONG pid = subject_swc.hashNeuron.value(subject_swc.listNeuron[i].pn);
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

    V3DLONG sz[4];
    sz[0] = 4080/3;
    sz[1] = 13297/3;
    sz[2] = 210/3;
    V3DLONG pagesz = sz[0]*sz[1]*sz[2];

    char *target_map = new char[pagesz];
    if(!target_map) return result;
    for(V3DLONG i =0; i < pagesz; i++) target_map[i] = 0;
    for(V3DLONG i = 0; i <target_swc.listNeuron.size();i++)
    {

        NeuronSWC S = target_swc.listNeuron[i];
        V3DLONG dz = (int)(S.z/length);
        V3DLONG dy = (int)(S.y/length);
        V3DLONG dx = (int)(S.x/length);
        target_map[dz*sz[0]*sz[1] + dy*sz[0] + dx] = 1;
    }


    for (V3DLONG i=0;i<seg_list.size();i++)
    {

        bool flag = 1;
        for (V3DLONG j=0;j<seg_list[i]->size();j++)
        {
            Point* node = seg_list[i]->at(j);
            V3DLONG dz = (int)(node->z/length);
            V3DLONG dy = (int)(node->y/length);
            V3DLONG dx = (int)(node->x/length);

            if(target_map[dz*sz[0]*sz[1] + dy*sz[0] + dx] ==1)
            {
               flag = 0;
               break;
            }
        }
        if(flag)
        {
            for (V3DLONG j=0;j<seg_list[i]->size();j++)
            {
                tree.push_back(seg_list[i]->at(j));
                index_map.insert(pair<Point*, V3DLONG>(seg_list[i]->at(j), tree.size()-1));
            }
        }
    }

    if(target_map) {delete target_map; target_map = 0;}
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

void prunSWC(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString fileOpenName;
    fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Open File"),
            "",
            QObject::tr("Supported file (*.swc *.eswc)"
                ";;Neuron structure	(*.swc)"
                ";;Extended neuron structure (*.eswc)"
                ));
    if(fileOpenName.isEmpty())
        return;
    double length = 0;
    NeuronTree nt;
    if (fileOpenName.toUpper().endsWith(".SWC") || fileOpenName.toUpper().endsWith(".ESWC"))
    {
        bool ok;
        nt = readSWC_file(fileOpenName);
        length = QInputDialog::getDouble(parent, "Please specify the minimum segment length","segment length:",1,0,2147483647,0.1,&ok);
        if (!ok)
            return;
    }

    QVector<QVector<V3DLONG> > childs;


    V3DLONG neuronNum = nt.listNeuron.size();
    childs = QVector< QVector<V3DLONG> >(neuronNum, QVector<V3DLONG>() );

    for (V3DLONG i=0;i<neuronNum;i++)
    {
        V3DLONG par = nt.listNeuron[i].pn;
        if (par<0) continue;
        childs[nt.hashNeuron.value(par)].push_back(i);
    }


    V3DLONG *AAA = new V3DLONG[neuronNum];
    for (V3DLONG i=0;i<neuronNum;i++)
    {
        AAA[i] = 1;
    }

    QList<NeuronSWC> list = nt.listNeuron;
   for (int i=0;i<list.size();i++)
    {


        if (childs[i].size()==0)
        {
            int index_tip = 0;
            int parent_tip = getParent(i,nt);
            while(childs[parent_tip].size()<2)
            {

                parent_tip = getParent(parent_tip,nt);
                index_tip++;
            }
            if(index_tip < length)
            {
                AAA[i] = -1;

                int parent_tip = getParent(i,nt);
                while(childs[parent_tip].size()<2)
                {
                    AAA[parent_tip] = -1;

                    parent_tip = getParent(parent_tip,nt);
                }
            }

        }

    }

   //NeutronTree structure
   NeuronTree marker_MST;
   QList <NeuronSWC> listNeuron;
   QHash <int, int>  hashNeuron;
   listNeuron.clear();
   hashNeuron.clear();

   //set node

   NeuronSWC S;
   listNeuron.append(S);
   hashNeuron.insert(S.n, listNeuron.size()-1);


   for (int i=0;i<list.size();i++)
   {
       if(AAA[i] == 1)
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
   marker_MST.n = -1;
   marker_MST.on = true;
   marker_MST.listNeuron = listNeuron;
   marker_MST.hashNeuron = hashNeuron;

   QString outfilename = "tips.swc";

   writeSWC_file(outfilename,marker_MST);
}

void zsectionsTotiles(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString m_InputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory including all images "),
                                          QDir::currentPath(),
                                          QFileDialog::ShowDirsOnly);

    QString m_OutputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory to save all tiles "),
                                          QDir::currentPath(),
                                          QFileDialog::ShowDirsOnly);

    QStringList imgList = importSeriesFileList_addnumbersort(m_InputfolderName);

    Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim;

    V3DLONG count=0;
    foreach (QString img_str, imgList)
    {
        V3DLONG offset[3];
        offset[0]=0; offset[1]=0; offset[2]=0;

        indexed_t<V3DLONG, REAL> idx_t(offset);

        idx_t.n = count;
        idx_t.ref_n = 0; // init with default values
        idx_t.fn_image = img_str.toStdString();
        idx_t.score = 0;

        vim.tilesList.push_back(idx_t);
        count++;
    }

    int NTILES  = vim.tilesList.size();


    unsigned char * data1d = 0;
    V3DLONG in_sz[4];
    int datatype;

    if (!simple_loadimage_wrapper(callback,const_cast<char *>(vim.tilesList.at(0).fn_image.c_str()), data1d, in_sz, datatype))
    {
        fprintf (stderr, "Error happens in reading the subject file [%0]. Exit. \n",vim.tilesList.at(0).fn_image.c_str());
        return;
    }

    int Ws = 1024;
    V3DLONG N = in_sz[0];
    V3DLONG M = in_sz[1];
    V3DLONG P = NTILES;
    V3DLONG pagesz = N*M*1;

    if(data1d) {delete []data1d; data1d=0;}

    for(V3DLONG iz = 0; iz < NTILES; iz = iz + 100)
    {
        V3DLONG zb = iz;
        V3DLONG ze = iz + 100 - 1; if(ze>=NTILES-1) ze = NTILES-1;

        unsigned char *sub_image=0;
        V3DLONG sub_image_sz = N*M*(ze-zb+1);
        sub_image = new unsigned char [sub_image_sz];
        for(V3DLONG i = 0; i < sub_image_sz; i++)
            sub_image[i] = 0;

        V3DLONG j = 0;
        for(int ii = zb; ii < ze + 1; ii++)
        {
            unsigned char * data1d = 0;
            V3DLONG in_sz[4];
            int datatype;

            if (!simple_loadimage_wrapper(callback,const_cast<char *>(vim.tilesList.at(ii).fn_image.c_str()), data1d, in_sz, datatype))
            {
                fprintf (stderr, "Error happens in reading the subject file [%0]. Exit. \n",vim.tilesList.at(ii).fn_image.c_str());
                return;
            }

            for(V3DLONG i = 0; i < pagesz; i++)
            {
                sub_image[j] = 255 - data1d[i];
                j++;
            }
            if(data1d) {delete []data1d; data1d=0;}

        }

        for(V3DLONG iy = 0; iy < M; iy = iy+Ws-Ws/10)
        {
            V3DLONG yb = iy;
            V3DLONG ye = iy+Ws-1; if(ye>=M-1) ye = M-1;

            for(V3DLONG ix = 0; ix < N; ix = ix+Ws-Ws/10)
            {
                V3DLONG xb = ix;
                V3DLONG xe = ix+Ws-1; if(xe>=N-1) xe = N-1;

                unsigned char *blockarea=0;
                V3DLONG blockpagesz = (xe-xb+1)*(ye-yb+1)*(ze-zb+1);
                blockarea = new unsigned char [blockpagesz];
                for(V3DLONG i = 0; i < blockpagesz; i++)
                    blockarea[i] = 0;

                V3DLONG i = 0;
                for(V3DLONG iz = 0; iz < ze-zb+1; iz++)
                {
                    V3DLONG offsetk = iz*M*N;
                    for(V3DLONG iy = yb; iy < ye+1; iy++)
                    {
                        V3DLONG offsetj = iy*N;
                        for(V3DLONG ix = xb; ix < xe+1; ix++)
                        {

                            blockarea[i] = sub_image[offsetk + offsetj + ix];
                            i++;
                        }
                    }
                }

                V3DLONG block_sz[4];
                block_sz[0] = xe-xb+1; block_sz[1] = ye-yb+1; block_sz[2] = (ze-zb+1); block_sz[3] = 1;



                QString outputTile(m_OutputfolderName);
                outputTile.append(QString("/x_%1_%2_y_%3_%4_z_%5_%6.raw").arg(xb).arg(xe).arg(yb).arg(ye).arg(zb).arg(ze));
                simple_saveimage_wrapper(callback, outputTile.toStdString().c_str(), (unsigned char *)blockarea, block_sz, 1);
                if(blockarea) {delete []blockarea; blockarea=0;}

            }
        }

        if(sub_image) {delete []sub_image; sub_image=0;}


    }

    V3DLONG tilenum = (floor(N/(0.9*Ws))+1.0)*(floor(M/(0.9*Ws))+1.0);
    QString tc_name(m_OutputfolderName);
    tc_name.append("/stitched_image.tc");

    ofstream myfile;
    myfile.open (tc_name.toStdString().c_str(),ios::out | ios::app );
    myfile << "# thumbnail file \n";
    myfile << "NULL \n\n";
    myfile << "# tiles \n";
    myfile << tilenum << " \n\n";
    myfile << "# dimensions (XYZC) \n";
    myfile << N << " " << M << " " << P << " " << 1 << " ";
    myfile << "\n\n";
    myfile << "# origin (XYZ) \n";
    myfile << "0.000000 0.000000 0.000000 \n\n";
    myfile << "# resolution (XYZ) \n";
    myfile << "1.000000 1.000000 1.000000 \n\n";
    myfile << "# image coordinates look up table \n";
    myfile.close();

    for(V3DLONG iy = 0; iy < M; iy = iy+Ws-Ws/10)
    {
        V3DLONG yb = iy;
        V3DLONG ye = iy+Ws-1; if(ye>=M-1) ye = M-1;

        for(V3DLONG ix = 0; ix < N; ix = ix+Ws-Ws/10)
        {
            V3DLONG xb = ix;
            V3DLONG xe = ix+Ws-1; if(xe>=N-1) xe = N-1;

            unsigned char *tilearea=0;
            V3DLONG tilepagesz = (xe-xb+1)*(ye-yb+1)*P;
            tilearea = new unsigned char [tilepagesz];
            for(V3DLONG i = 0; i < tilepagesz; i++)
                tilearea[i] = 0;

            V3DLONG tilearea_sz[4];
            tilearea_sz[0] = xe-xb+1; tilearea_sz[1] = ye-yb+1; tilearea_sz[2] = P; tilearea_sz[3] = 1;
            V3DLONG i = 0;

            for(V3DLONG iz = 0; iz < NTILES; iz = iz + 100)
            {
                V3DLONG zb = iz;
                V3DLONG ze = iz + 100 - 1; if(ze>=NTILES-1) ze = NTILES-1;

                QString inputTile(m_OutputfolderName);
                inputTile.append(QString("/x_%1_%2_y_%3_%4_z_%5_%6.raw").arg(xb).arg(xe).arg(yb).arg(ye).arg(zb).arg(ze));

                unsigned char * sub_data1d = 0;
                V3DLONG in_sz_sub[4];
                int datatype;

                if (!simple_loadimage_wrapper(callback,inputTile.toStdString().c_str(), sub_data1d, in_sz_sub, datatype))
                {
                    return;
                }
                V3DLONG inputTilepagesz = in_sz_sub[0]*in_sz_sub[1]*in_sz_sub[2];
                for(V3DLONG j = 0; j < inputTilepagesz; j++)
                {
                    tilearea[i] = sub_data1d[j];
                    i++;
                }
                if(sub_data1d) {delete []sub_data1d; sub_data1d = 0;}
                remove(inputTile.toStdString().c_str());
            }


            QString outputTile(m_OutputfolderName);
            outputTile.append(QString("/x_%1_%2_y_%3_%4.raw").arg(xb).arg(xe).arg(yb).arg(ye));
            simple_saveimage_wrapper(callback, outputTile.toStdString().c_str(), (unsigned char *)tilearea, tilearea_sz, 1);

            myfile.open (tc_name.toStdString().c_str(),ios::out | ios::app );
            QString outputilefull;
            outputilefull.append(QString("x_%1_%2_y_%3_%4.raw").arg(xb).arg(xe).arg(yb).arg(ye));
            outputilefull.append(QString("   ( %1, %2, 0) ( %3, %4, %5)").arg(xb).arg(yb).arg(xe).arg(ye).arg(P-1));
            myfile << outputilefull.toStdString();
            myfile << "\n";
            myfile.close();

            if(tilearea) {delete []tilearea; tilearea =0;}

        }
    }

    myfile.open (tc_name.toStdString().c_str(),ios::out | ios::app );
    myfile << "\n# MST LUT\n";
    myfile.close();

}

void qsublist(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString m_InputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory including all images "),
                                                                   QDir::currentPath(),
                                                                   QFileDialog::ShowDirsOnly);

    QString m_OutputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory to save all tiles "),
                                                                   QDir::currentPath(),
                                                                   QFileDialog::ShowDirsOnly);

    QStringList imgList = importSeriesFileList_addnumbersort(m_InputfolderName);

    Y_VIM<REAL, V3DLONG, indexed_t<V3DLONG, REAL>, LUT<V3DLONG> > vim;

    V3DLONG count=0;
    foreach (QString img_str, imgList)
    {
        V3DLONG offset[3];
        offset[0]=0; offset[1]=0; offset[2]=0;

        indexed_t<V3DLONG, REAL> idx_t(offset);

        idx_t.n = count;
        idx_t.ref_n = 0; // init with default values
        idx_t.fn_image = img_str.toStdString();
        idx_t.score = 0;

        vim.tilesList.push_back(idx_t);
        count++;
    }

    int NTILES  = vim.tilesList.size();

    for(int ii = 0; ii < NTILES; ii++)
    {
        ofstream myfile;
        QString filename(const_cast<char *>(vim.tilesList.at(ii).fn_image.c_str()));
        filename.append(".qsub");
        myfile.open (filename.toStdString().c_str(),ios::out | ios::app );
        myfile << "## There are several queues available.  Please check with \!ITsupport to verify which queue you should use \n";
        myfile << "#PBS -q mindscope \n";
        myfile << "# declare that your job will use no more than 4gb of memory _at_peak_ \n";
        myfile << "#PBS -l vmem=20g\n";
        myfile << "# Allow up to 36hr of walltime.  Default is 12 hours\n";
        myfile << "#PBS -l walltime=12:00:00\n";
        myfile << "# request just one core on the host\n";
        myfile << "#PBS -l ncpus=1\n";
        myfile << "# Give your job a descriptive name. This is visible in qstat and other job reports.  Also serves as the default basename for log files\n";
        myfile << "#PBS -N a" << ii <<"\n";
        myfile << "# should torque automatically re-run the job on error?\n";
        myfile << "#PBS -r n\n";
        myfile << "# merge STDOUT into STDERR\n";
        myfile << "#PBS -j oe\n";
        myfile << "# location for stderr/stdout log files _after_ job completion\n";
        myfile << "#PBS -o /data/mat/zhi/a" << ii << ".out\n";
        myfile << "#PBS -o /data/mat/zhi/a" << ii << "_error.out\n\n\n";
        myfile << "# send email on job error\n";
        myfile << "#PBS -m a\n";
        myfile << "export DISPLAY=:$RANDOM\n";
        myfile << "Xvfb $DISPLAY -auth /dev/null &\n";
        myfile << "export LD_PRELOAD=/usr/lib64/libstdc++.so.6\n";
        myfile << "cd /data/mat/zhi/VAA3D/vaa3d_redhat_fedora_ubuntu_64bit_v2.868/\n";
        myfile << "./start_vaa3d.sh -x multiscaleEnhancement -f adaptive_auto -i " << vim.tilesList.at(ii).fn_image.c_str() << " -o " << vim.tilesList.at(ii).fn_image.c_str() <<"_enhanced.raw -p 2 1 1 0 0 \n";
        myfile << "kill %1\n";


        myfile.close();


      }
}



void threeDimageTotiles(V3DPluginCallback2 &callback, QWidget *parent)
{

    QString fileOpenName;
    fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Open File"),
            "",
            QObject::tr("Supported file (*.raw *.RAW *.V3DRAW *.v3draw)"
                ));
    if(fileOpenName.isEmpty())
        return;

    unsigned char * datald = 0;
    V3DLONG *in_zz = 0;

    int datatype;

    if (!loadRaw2Stack_2byte(const_cast<char *>(fileOpenName.toStdString().c_str()), datald, in_zz, datatype))
    {
        return;
    }
    printf("%d %d %d\n\n",in_zz[0],in_zz[1],in_zz[2]);

    Image4DSimple * new4DImage = new Image4DSimple();
    new4DImage->setData((unsigned char *)datald, 50, 50, in_zz[2], 1, V3D_UINT8);
    v3dhandle newwin = callback.newImageWindow();
    callback.setImage(newwin, new4DImage);
    callback.setImageName(newwin, "Cropped image");
    callback.updateImageWindow(newwin);



    /*v3dhandle curwin = callback.currentImageWindow();
    if (!curwin)
    {
        v3d_msg("You don't have any image open in the main window.");
        return;
    }

    Image4DSimple* p4DImage = callback.getImage(curwin);

    if (!p4DImage)
    {
        v3d_msg("The image pointer is invalid. Ensure your data is valid and try again!");
        return;
    }


    QString m_OutputfolderName = QFileDialog::getExistingDirectory(parent, QObject::tr("Choose the directory to save all tiles "),
                                                                   QDir::currentPath(),
                                                                   QFileDialog::ShowDirsOnly);

    unsigned char* data1d = p4DImage->getRawData();
    V3DLONG pagesz = p4DImage->getTotalUnitNumberPerChannel();

    V3DLONG N = p4DImage->getXDim();
    V3DLONG M = p4DImage->getYDim();
    V3DLONG P = p4DImage->getZDim();
    V3DLONG sc = p4DImage->getCDim();
    ImagePixelType pixeltype = p4DImage->getDatatype();

    int Ws = 1024;


    V3DLONG in_sz[4];
    in_sz[0] = N; in_sz[1] = M; in_sz[2] = P; in_sz[3] = 1;

    V3DLONG tilenum = (floor(N/(0.9*Ws))+1.0)*(floor(M/(0.9*Ws))+1.0);

    QString tc_name = m_OutputfolderName;
    tc_name.append("/stitched_image.tc");

    ofstream myfile;
    myfile.open (tc_name.toStdString().c_str(),ios::out | ios::app );
    myfile << "# thumbnail file \n";
    myfile << "NULL \n\n";
    myfile << "# tiles \n";
    myfile << tilenum << " \n\n";
    myfile << "# dimensions (XYZC) \n";
    myfile << N << " " << M << " " << P << " " << 1 << " ";
    myfile << "\n\n";
    myfile << "# origin (XYZ) \n";
    myfile << "0.000000 0.000000 0.000000 \n\n";
    myfile << "# resolution (XYZ) \n";
    myfile << "1.000000 1.000000 1.000000 \n\n";
    myfile << "# image coordinates look up table \n";
    myfile.close();

    for(V3DLONG iy = 0; iy < M; iy = iy+Ws-Ws/10)
    {
        V3DLONG yb = iy;
        V3DLONG ye = iy+Ws-1; if(ye>=M-1) ye = M-1;

        for(V3DLONG ix = 0; ix < N; ix = ix+Ws-Ws/10)
        {
            V3DLONG xb = ix;
            V3DLONG xe = ix+Ws-1; if(xe>=N-1) xe = N-1;

            unsigned char *blockarea=0;
            V3DLONG blockpagesz = (xe-xb+1)*(ye-yb+1)*P;
            blockarea = new unsigned char [blockpagesz];
            int i = 0;
            for(V3DLONG iz = 0; iz < P; iz++)
            {
                V3DLONG offsetk = iz*M*N;
                for(V3DLONG iy = yb; iy < ye+1; iy++)
                {
                    V3DLONG offsetj = iy*N;
                    for(V3DLONG ix = xb; ix < xe+1; ix++)
                    {

                        blockarea[i] = data1d[offsetk + offsetj + ix];
                        i++;
                    }
                }
            }

            V3DLONG block_sz[4];
            block_sz[0] = xe-xb+1; block_sz[1] = ye-yb+1; block_sz[2] = P; block_sz[3] = 1;

            QString outputTile = m_OutputfolderName;
            outputTile.append(QString("/x_%1_%2_y_%3_%4.raw").arg(xb).arg(xe).arg(yb).arg(ye));
            simple_saveimage_wrapper(callback, outputTile.toStdString().c_str(), (unsigned char *)blockarea, block_sz, 1);

            myfile.open (tc_name.toStdString().c_str(),ios::out | ios::app );
            QString outputilefull;
            outputilefull.append(QString("x_%1_%2_y_%3_%4.raw").arg(xb).arg(xe).arg(yb).arg(ye));
            outputilefull.append(QString("   ( %1, %2, 0) ( %3, %4, %5)").arg(xb).arg(yb).arg(xe).arg(ye).arg(P-1));
            myfile << outputilefull.toStdString();
            myfile << "\n";
            myfile.close();
            if(blockarea) {delete []blockarea; blockarea =0;}
        }

    }
    myfile.open (tc_name.toStdString().c_str(),ios::out | ios::app );
    myfile << "\n# MST LUT\n";
    myfile.close();

    v3d_msg("done!");*/
    return;
}

int loadRaw2Stack_2byte(char * filename, unsigned char * & img, V3DLONG * & sz, int & datatype)
{
    FILE * fid = fopen(filename, "rb");
    if (!fid)
    {
        return 0;
    }

    fseek (fid, 0, SEEK_END);
    V3DLONG fileSize = ftell(fid);
    rewind(fid);
/*
#endif
*/
    /* Read header */

    char formatkey[] = "raw_image_stack_by_hpeng";
    V3DLONG lenkey = strlen(formatkey);
    char * keyread = new char [lenkey+1];
    if (!keyread)
        return 0;

    V3DLONG nread = fread(keyread, 1, lenkey, fid);
    if (nread!=lenkey) {
        if (keyread) {delete []keyread; keyread=0;}
        return 0;
    }

    keyread[lenkey] = '\0';

    V3DLONG i;
    if (strcmp(formatkey, keyread)) /* is non-zero then the two strings are different */
    {
        if (keyread) {delete []keyread; keyread=0;}
        return 0;
    }

    char endianCodeData;
    int dummy = fread(&endianCodeData, 1, 1, fid);
    printf("The data endian code is [%c]\n", endianCodeData);
    if (endianCodeData!='B' && endianCodeData!='L')
    {
        if (keyread) {delete []keyread; keyread=0;}
        return 0;    }

    char endianCodeMachine;
    endianCodeMachine = checkMachineEndian();
    printf("The machine endian code is [%c]\n", endianCodeMachine);
    if (endianCodeMachine!='B' && endianCodeMachine!='L')
    {
        if (keyread) {delete []keyread; keyread=0;}
        return 0;
       }

    int b_swap = (endianCodeMachine==endianCodeData)?0:1;

    short int dcode = 0;
    dummy = fread(&dcode, 2, 1, fid); /* because I have already checked the file size to be bigger than the header, no need to check the number of actual bytes read. */
    if (b_swap)
        swap2bytes((void *)&dcode);

    switch (dcode)
    {
        case 1:
            datatype = 1; /* temporarily I use the same number, which indicates the number of bytes for each data point (pixel). This can be extended in the future. */
            break;

        case 2:
            datatype = 2;
            break;

        case 4:
            datatype = 4;
            break;

        default:
            if (keyread) {delete []keyread; keyread=0;}
            return  0;
    }

    V3DLONG unitSize = datatype; // temporarily I use the same number, which indicates the number of bytes for each data point (pixel). This can be extended in the future.

    BIT32_UNIT mysz[4];
    mysz[0]=mysz[1]=mysz[2]=mysz[3]=0;
    int tmpn=(int)fread(mysz, 4, 4, fid); // because I have already checked the file size to be bigger than the header, no need to check the number of actual bytes read.
    if (tmpn!=4) {
        if (keyread) {delete []keyread; keyread=0;}
        return 0;
    }

    if (b_swap)
    {
        for (i=0;i<4;i++)
        {
            //swap2bytes((void *)(mysz+i));
            if (b_VERBOSE_PRINT)
                printf("mysz raw read unit[%ld]: [%d] ", i, mysz[i]);
            swap4bytes((void *)(mysz+i));
            if (b_VERBOSE_PRINT)
                printf("swap unit: [%d][%0x] \n", mysz[i], mysz[i]);
        }
    }
    if (sz) {delete []sz; sz=0;}
    sz = new V3DLONG [4]; // reallocate the memory if the input parameter is non-null. Note that this requests the input is also an NULL point, the same to img.
    if (!sz)
    {
        if (keyread) {delete []keyread; keyread=0;}
        return 0;
    }

    V3DLONG totalUnit = 1;

    for (i=0;i<4;i++)
    {
        sz[i] = (V3DLONG)mysz[i];
        totalUnit *= sz[i];
    }

    V3DLONG endx = 50, endy = 50, endz = sz[2];
    V3DLONG startx = 0, starty = 0, startz = 0;

    V3DLONG tmpw = endx - startx;
    V3DLONG tmph = endy - starty;
    V3DLONG tmpz = endz - startz;

    V3DLONG head = 4*4+2+1+lenkey; // header_len ?
    V3DLONG pgsz1=sz[2]*sz[1]*sz[0], pgsz2=sz[1]*sz[0], pgsz3=sz[0];
    V3DLONG cn = tmpw*tmph*tmpz;
    V3DLONG kn = tmpw*tmph;
    V3DLONG total = tmpw*tmph*tmpz*sz[3];

    if (img) {delete []img; img=0;}
    V3DLONG totalBytes = V3DLONG(unitSize)*V3DLONG(total);
    try
    {
        img = new unsigned char [totalBytes];
    }
    catch (...)
    {
        if (keyread) {delete []keyread; keyread=0;}
        if (sz) {delete []sz; sz=0;}
        return 0;
    }

    //V3DLONG count=0; // unused
    V3DLONG c,j,k;
    for (c = 0; c < sz[3]; c++)
    {
        for (k = startz; k < endz; k++)
        {
            for (j = starty; j< endy; j++)
            {
                rewind(fid);
                fseek(fid, (long)(head+(c*pgsz1 + k*pgsz2 + j*pgsz3 + startx)*unitSize), SEEK_SET);
                int dummy = ftell(fid);
                dummy = fread(img+(c*cn+(k-startz)*kn + (j-starty)*tmpw)*unitSize,unitSize,tmpw,fid);
            }
        }
    }
    // swap the data bytes if necessary

    if (b_swap==1)
    {
        if (unitSize==2)
        {
            for (i=0;i<total; i++)
            {
                swap2bytes((void *)(img+i*unitSize));
            }
        }
        else if (unitSize==4)
        {
            for (i=0;i<total; i++)
            {
                swap4bytes((void *)(img+i*unitSize));
            }
        }
    }

    // clean and return

    if (keyread) {delete [] keyread; keyread = 0;}
    fclose(fid); //bug fix on 060412


}


char checkMachineEndian()
{
    char e='N'; //for unknown endianness

    V3DLONG int a=0x44332211;
    unsigned char * p = (unsigned char *)&a;
    if ((*p==0x11) && (*(p+1)==0x22) && (*(p+2)==0x33) && (*(p+3)==0x44))
        e = 'L';
    else if ((*p==0x44) && (*(p+1)==0x33) && (*(p+2)==0x22) && (*(p+3)==0x11))
        e = 'B';
    else if ((*p==0x22) && (*(p+1)==0x11) && (*(p+2)==0x44) && (*(p+3)==0x33))
        e = 'M';
    else
        e = 'N';

    //printf("[%c] \n", e);
    return e;
}

void swap2bytes(void *targetp)
{
    unsigned char * tp = (unsigned char *)targetp;
    unsigned char a = *tp;
    *tp = *(tp+1);
    *(tp+1) = a;
}

void swap4bytes(void *targetp)
{
    unsigned char * tp = (unsigned char *)targetp;
    unsigned char a = *tp;
    *tp = *(tp+3);
    *(tp+3) = a;
    a = *(tp+1);
    *(tp+1) = *(tp+2);
    *(tp+2) = a;
}






























