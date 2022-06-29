/* DendritePruning_plugin.cpp
 * This is a test plugin, you can use it as a demo.
 * 2021-4-6 : by XZ
 */
 
#include "v3d_message.h"
#include <vector>
#include "DendritePruning_plugin.h"

#include "branchtree.h"

#include "swccompare.h"

#include "neuronsplit.h"

using namespace std;
Q_EXPORT_PLUGIN2(DendritePruning, DendritePruningPlugin);
 
QStringList DendritePruningPlugin::menulist() const
{
	return QStringList() 
        <<tr("Pruning")
        <<tr("NeuronSplit")
		<<tr("about");
}

QStringList DendritePruningPlugin::funclist() const
{
	return QStringList()
        <<tr("func1")
		<<tr("func2")
		<<tr("help");
}

void DendritePruningPlugin::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
    if (menu_name == tr("Pruning"))
	{
        pruningMenu(callback);
	}
    else if (menu_name == tr("NeuronSplit"))
	{
        splitNeuronTreeManu(callback);
	}
	else
	{
		v3d_msg(tr("This is a test plugin, you can use it as a demo.. "
			"Developed by XZ, 2021-4-6"));
	}
}

bool DendritePruningPlugin::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
	vector<char*> infiles, inparas, outfiles;
	if(input.size() >= 1) infiles = *((vector<char*> *)input.at(0).p);
	if(input.size() >= 2) inparas = *((vector<char*> *)input.at(1).p);
	if(output.size() >= 1) outfiles = *((vector<char*> *)output.at(0).p);

    if (func_name == tr("pruning"))
	{

        QString swcPath = infiles.size()>=1 ? infiles[0] : "";
        QString markerIndexPath = infiles.size()>=2 ? infiles[1] : "";
        QString imgPath = infiles.size()>=3 ? infiles[2] : "";



        double lengthTh = inparas.size()>=1 ? atof(inparas[0]) : 11;
        double linearityTh = inparas.size()>=2 ? atof(inparas[1]) : 1.4;
        double angleTh = inparas.size()>=3 ? atof(inparas[2]) : 160;
        int lamda = inparas.size()>=4 ? atoi(inparas[3]) : 2;
        float lengthMax = inparas.size()>=5 ? atof(inparas[4]) : 850;
        double times = inparas.size()>=6 ? atof(inparas[5]) : 5;

        bool noisyPruning = inparas.size()>=7 ? atoi(inparas[6]) : true;
        bool multiSomaPruning = inparas.size()>=8 ? atoi(inparas[7]) : true;
        bool structurePruning = inparas.size()>=9 ? atoi(inparas[8]) : true;
        bool inflectionPruning = inparas.size()>=10 ? atoi(inparas[9]) : false;

        if(imgPath == ""){
            noisyPruning = false;
            inflectionPruning = false;
        }
        if(markerIndexPath == ""){
            multiSomaPruning = false;
        }
        if(swcPath == ""){
            qDebug()<<"Please specify SWC file!";
            return false;
        }

        unsigned char* pdata = 0;
        int dataType = 1;
        V3DLONG sz[4] = {0,0,0,0};
        if(noisyPruning || inflectionPruning){
            simple_loadimage_wrapper(callback,imgPath.toStdString().c_str(),pdata,sz,dataType);
        }


        NeuronTree t = readSWC_file(swcPath);
        for(int i=0; i<t.listNeuron.size(); ++i){
            t.listNeuron[i].type = 3;
        }
        BranchTree bt = BranchTree();
        bt.initialize(t);

        if(noisyPruning){
            bt.pruningByLength(pdata,sz,lengthTh,linearityTh);
            bt.savePrunedNeuronTree(swcPath + "_outf1.swc");
            bt.update();

            bt.pruningSoma(times);
            bt.savePrunedNeuronTree(swcPath + "_outf2.swc");
            bt.update();
        }

        if(multiSomaPruning){
            if(markerIndexPath.endsWith(".marker")){
                QList<ImageMarker> markers = readMarker_file(markerIndexPath);
                vector<V3DLONG> nIndexList = convertMarkers2SomaNList(t,markers);
                bt.pruningAdjacentSoma3(nIndexList,lamda);
            }else{
                bt.pruningAdjacentSoma3(markerIndexPath,lamda);
            }

            bt.update();
            bt.pruningAdjacentSoma4(lengthMax,lamda);

            bt.savePrunedNeuronTree(swcPath + "_outf3.swc");
            bt.update();
        }

        if(structurePruning){
            bt.pruningCross(angleTh,lengthTh);
            bt.update();
            bt.pruningByLength(lengthTh);
            bt.savePrunedNeuronTree(swcPath + "_outf4.swc");
            bt.update();
        }

        if(inflectionPruning){
            bt.pruningSuspectedBranch(pdata,sz,75,lengthTh);
            bt.update();
            bt.pruningInflectionPoints(pdata,sz,lengthTh,0);
            bt.savePrunedNeuronTree(swcPath + "_outf5.swc");
        }

        if(pdata){
            delete[] pdata;
            pdata = 0;
        }


	}
    else if (func_name == tr("calManualFeature"))
	{
        QString swcPath = infiles[0];
        QString imgPath = infiles[1];
        NeuronTree t = readSWC_file(swcPath);

        unsigned char* pdata = 0;
        int dataType = 1;
        V3DLONG sz[4] = {0,0,0,0};
        simple_loadimage_wrapper(callback,imgPath.toStdString().c_str(),pdata,sz,dataType);

        NeuronTree dendrite;
        for(int i=0; i<t.listNeuron.size(); ++i){
            NeuronSWC s = t.listNeuron[i];
            if(s.type == 1 || s.type == 3 ){
                dendrite.listNeuron.push_back(s);
            }
        }
        dendrite.hashNeuron.clear();
        for(int i=0; i<dendrite.listNeuron.size(); ++i){
            dendrite.hashNeuron.insert(dendrite.listNeuron[i].n,i);
        }

        BranchTree bt = BranchTree();
        ofstream csvFile;
        QString rLevel0Path = "F:\\manual\\rLevel0.csv";
        QString localAnglePath = "F:\\manual\\localAngle.csv";
        csvFile.open(rLevel0Path.toStdString().c_str(),ios::app);
        bt.initialize(dendrite);
//        bt.show(pdata,sz);
        bt.calRlevel0Branches(pdata,sz,csvFile);
        csvFile.close();
        csvFile.open(localAnglePath.toStdString().c_str(),ios::app);
        bt.calBifurcationLocalAngle(csvFile);
        csvFile.close();

        if(pdata){
            delete[] pdata;
            pdata = 0;
        }


	}
    else if (func_name == tr("getSomaMarker")) {
        QString imgPath = infiles.size()>=1 ? infiles[0] : "";
        unsigned char* pdata = 0;
        int dataType = 1;
        V3DLONG sz[4] = {0,0,0,0};
        simple_loadimage_wrapper(callback,imgPath.toStdString().c_str(),pdata,sz,dataType);
        ImageMarker soma = getCenterSoma(pdata,sz);
        QList <ImageMarker> markers;
        markers.push_back(soma);
        QString markerPath = imgPath + ".marker";
        writeMarker_file(markerPath,markers);
        if(pdata){
            delete[] pdata;
            pdata = 0;
        }
    }
    else if (func_name == tr("compareSwc")) {
        QString swcAutoPath = infiles.size()>=1 ? infiles[0] : "";
//        QString swcRemoveNPath = infiles.size()>=2 ? infiles[1] : "";
        QString swcManualPath = infiles.size()>=2 ? infiles[1] : "";
        QString swcPrunedPath = infiles.size()>=3 ? infiles[2] : "";
        QString swcOptimalPath = outfiles.size()>=1 ? outfiles[0] : "";

        float dTh = inparas.size()>=1 ? atof(inparas[0]) : 5;

        QString csvPath = outfiles.size()>=2 ? outfiles[1] : "F:\\manual\\compareResult.csv";

        NeuronTree swcAuto = readSWC_file(swcAutoPath);
        NeuronTree swcManual = readSWC_file(swcManualPath);
        NeuronTree swcPruned = readSWC_file(swcPrunedPath);


        compareResult cr;
        if(!compareSwc(swcAuto,swcManual,swcPruned,cr,dTh,swcOptimalPath)){
            return false;
        }
        ofstream csvFile;

        BranchTree bt_manual = BranchTree();
        bt_manual.initialize(swcManual);
        csvFile.open(csvPath.toStdString().c_str(),ios::app);
        csvFile<<swcAutoPath.toStdString()<<','<<cr.branchAutoSN<<','
              <<cr.branchAutoSL<<','<<cr.branchOptimalSN<<','<<cr.branchOptimalSL<<','
              <<bt_manual.branches.size()<<','
              <<cr.rBranchNumber<<','<<cr.fBranchNumber<<','
             <<cr.negativeBranchNumber<<','<<cr.aBranchNumber<<','
            <<cr.rLength<<','<<cr.fLength<<','<<cr.negativeLength<<','<<cr.alength<<','
           <<cr.otherSoma<<endl;
        csvFile.close();


    }else if (func_name == tr("compareSwc2")) {
        QString swcAutoPath = infiles.size()>=1 ? infiles[0] : "";
        QString swcPrunedPath = infiles.size()>=2 ? infiles[1] : "";
        QString swcOptimalPath = infiles.size()>=3 ? infiles[2] : "";

        QString csvPath = outfiles.size()>=1 ? outfiles[0] : "F:\\manual\\compareResult.csv";

        NeuronTree swcAuto = readSWC_file(swcAutoPath);
        NeuronTree swcPruned = readSWC_file(swcPrunedPath);
        NeuronTree swcOptimal = readSWC_file(swcOptimalPath);

        qDebug()<<"----read swc end----";

        compareResult cr;
        if(!compareSwc3(swcAuto,swcPruned,swcOptimal,cr)){
            return false;
        }
        ofstream csvFile;

        csvFile.open(csvPath.toStdString().c_str(),ios::app);
        csvFile<<swcAutoPath.toStdString()<<','<<cr.branchAutoSN<<','
              <<cr.branchAutoSL<<','<<cr.branchOptimalSN<<','<<cr.branchOptimalSL<<','
              <<cr.rBranchNumber<<','<<cr.fBranchNumber<<','
             <<cr.negativeBranchNumber<<','<<cr.aBranchNumber<<','
            <<cr.rLength<<','<<cr.fLength<<','<<cr.negativeLength<<','<<cr.alength<<endl;
        csvFile.close();


    }
    else if (func_name == tr("shiftSwc")) {
        QString swcManualPath = infiles.size()>=1 ? infiles[0] : "";
        QString imgPath = infiles.size()>=2 ? infiles[1] : "";
        QString swcOutPath = outfiles.size()>=1 ? outfiles[0] : "";

        unsigned char* pdata = 0;
        int dataType = 1;
        V3DLONG sz[4] = {0,0,0,0};
        simple_loadimage_wrapper(callback,imgPath.toStdString().c_str(),pdata,sz,dataType);

        NeuronTree swcManual = readSWC_file(swcManualPath);
        shiftSwc(swcManual,pdata,sz);
        writeSWC_file(swcOutPath,swcManual);
        if(pdata){
            delete[] pdata;
            pdata = 0;
        }
    }
    else if (func_name == tr("getBranchNumber")) {
        vector<NeuronTree> trees;
        QString swcAutoPath = "";
        for(int i=0; i<infiles.size(); ++i){
            QString swcPath = infiles[i];
            if(i == 0)
                swcAutoPath = swcPath;
            NeuronTree nt = readSWC_file(swcPath);
            trees.push_back(nt);
        }
        ofstream csvFile;
        QString csvPath = "F:\\manual\\branchNumber.csv";
        csvFile.open(csvPath.toStdString().c_str(),ios::app);
        csvFile<<swcAutoPath.toStdString();
        getBranchNumber(trees,csvFile);
        csvFile.close();

    }
    else if (func_name == tr("getd")) {
        QString swcAutoPath = infiles.size()>=1 ? infiles[0] : "";
        QString swcManualPath = infiles.size()>=2 ? infiles[1] : "";

        NeuronTree swcAuto = readSWC_file(swcAutoPath);
        NeuronTree swcManual = readSWC_file(swcManualPath);

        ofstream csvFile;
        QString csvPath = "F:\\manual\\d.csv";
        csvFile.open(csvPath.toStdString().c_str(),ios::app);
        compareSwc2(swcAuto,swcManual,csvFile);
        csvFile.close();

    }
    else if (func_name == tr("convertMarkers2SomaNList")) {
        QString swcPath = infiles.size()>=1 ? infiles[0] : "";
        QString markerPath = infiles.size()>=2 ? infiles[1] : "";
        QString nIndexPath = outfiles.size()>=1 ? outfiles[0]: "";

        NeuronTree nt = readSWC_file(swcPath);
        QList<ImageMarker> markers = readMarker_file(markerPath);
        vector<V3DLONG> nIndexList = convertMarkers2SomaNList(nt,markers);
        ofstream o;
        o.open(nIndexPath.toStdString().c_str());
        for(V3DLONG n : nIndexList){
            o<<n<<endl;
        }
        o.close();
    }
	else if (func_name == tr("help"))
	{
        cout<<"Usage : v3d -x dllname -f pruning -i <swcPath> <imagePath> -p <d> <angleTh> <somaRth> <times> <inflectionPruning>"<<endl;
        cout<<endl;
        cout<<"d       the minimum pruning length   "<<endl;
        cout<<"angleTh the threshold of angle         "<<endl;
        cout<<"somaRth the radius of other soma    "<<endl;
        cout<<"times   the parameter in the soma pruning step    "<<endl;
        cout<<"inflectionPruning  select inflection pruning or not   "<<endl;
        cout<<endl;
	}
	else return false;

	return true;
}

