#define FNUM 22

#include "Nfmain.h"
#include "../neuron_editing/global_feature_compute.h"
#include "customary_structs/vaa3d_neurontoolbox_para.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include "../sort_neuron_swc/openSWCDialog.h"
#include <QMessageBox>

using namespace std;

void nf_main(const V3DPluginArgList & input, V3DPluginArgList & output)
{
	vector<char*>* inlist = (vector<char*>*)(input.at(0).p);
	if (inlist->size()!=1)
	{
		cerr<<"Error in input list. You must input .ano file or a list of .swc file name"<<endl;
		return;
	}

	QStringList nameList;
	QString qs_name(inlist->at(0));
	qs_name = qs_name.simplified();
	int neuronNum;
	vector<NeuronTree> nt_list;

	if (qs_name.toUpper().endsWith(".ANO"))
	{
		cout<<"reading a linker file..."<<endl;;
		P_ObjectFileType linker_object;
		if (!loadAnoFile(QString(qs_name),linker_object))
		{
			cerr<<"Error in reading the linker file."<<endl;
			return;
		}
		nameList = linker_object.swc_file_list;
		neuronNum = nameList.size();
		for (V3DLONG i=0;i<neuronNum;i++)
		{
			NeuronTree tmp = readSWC_file(nameList.at(i));
			nt_list.push_back(tmp);
		}
	}
	else if (qs_name.toUpper().endsWith(".SWC") || qs_name.toUpper().endsWith(".ESWC"))
	{
		cout<<"reading a list of swc file names."<<endl;
		nameList = qs_name.split(" ");
		neuronNum = nameList.size();
		for (V3DLONG i=0;i<neuronNum;i++)
		{
			NeuronTree tmp = readSWC_file(nameList.at(i));
			nt_list.push_back(tmp);
		}
	}

	for (int i=0;i<neuronNum;i++)
	{
		NeuronTree nt = nt_list[i];
		cout<<"\n--------------Neuron #"<<(i+1)<<"----------------\n";
		double * features = new double[FNUM];
		computeFeature(nt, features);
		printFeature(features);
		if (features) {delete []features; features = NULL;}
	}
}
void nf_main_infolder(const V3DPluginArgList &input, V3DPluginArgList &output)
{
//    vector<char*>* inlist = (vector<char*>*)(input.at(0).p);

    vector<char*> infiles, inparas, outfiles;
    if(input.size() >= 1) infiles = *((vector<char*> *)input.at(0).p);
    if(input.size() >= 2) inparas = *((vector<char*> *)input.at(1).p);
    if(output.size() >= 1) outfiles = *((vector<char*> *)output.at(0).p);

    if (outfiles.size()!=1)
    {
        cerr<<"Input Error. out Error. You must input an out file"<<endl;
        return;
    }
    QString outfile=outfiles[0];
    if (infiles.size()!=1)
    {
        cerr<<"Input Error. You must input one and only one folder path"<<endl;
        return;
    }
    string inputpath = infiles[0];
    QFileInfo inputInfo(QString::fromStdString(inputpath));

    if(!inputInfo.isDir())
    {
        cerr<<"Error in input. You must input an existed folder path"<<endl;
        return;
    }
    QStringList nameFilters;
    nameFilters<<"*.swc";
    nameFilters<<"*.eswc";
    nameFilters<<"*.ESWC";
    nameFilters<<"*.SWC";

    QDir indir(QString::fromStdString(inputpath));
    QStringList nameList = indir.entryList(nameFilters,QDir::Files|QDir::Readable, QDir::Name);
    vector<NeuronTree> nt_list;
    if(nameList.size()==0)
    {
        cerr<<"Error in input. This is an empty folder or I can't find swc and eswc files in it"<<endl;
        return;
    }
    for(V3DLONG i=0;i<nameList.size();i++)
    {
        QString thisname=QString::fromStdString(inputpath)+"/"+nameList.at(i);
        NeuronTree tmp = readSWC_file(thisname);
        nt_list.push_back(tmp);
    }

    ofstream csvOutFile(outfile.toStdString().c_str());
    if(!csvOutFile.is_open()){
         cerr<<"out Error: cannot open file to save"<<endl;
         return;
    }
    csvOutFile<<"Name,Nodes,SomaSurface,Stems,Bifurcations,Branches,Tips,OverallWidth,OverallHeight,OverallDepth";
    csvOutFile<<",AverageDiameter,Length,Surface,Volume,MaxEuclideanDistance,MaxPathDistance,MaxBranchOrder";
    csvOutFile<<",AverageContraction,AverageFragmentation,AverageParent-daughterRatio,AverageBifurcationAngleLocal,AverageBifurcationAngleRemote,HausdorffDimension"<<endl;

    for (V3DLONG i=0;i<nameList.size();i++)
    {
        NeuronTree nt = nt_list[i];
        cout<<"\n----------Neuron #"<<(i+1)<<",Total #"<<nameList.size()<<"-----------\n";
        double * features = new double[FNUM];
        computeFeature(nt, features);
        //save features to file

        QString thisname=nameList[i];
        csvOutFile<<thisname.toStdString();
        for (int i=0;i<FNUM;i++)
        {
            csvOutFile<<","<<features[i];
        }
        csvOutFile<<endl;

        if (features) {delete []features; features = NULL;}
    }
    csvOutFile.close();
}

void nf_main(V3DPluginCallback2 &callback, QWidget *parent)
{
    OpenSWCDialog * openDlg = new OpenSWCDialog(0, &callback);
    if (!openDlg->exec())
        return;

    NeuronTree nt = openDlg->nt;
	double * features = new double[FNUM];
	computeFeature(nt,features);
	QMessageBox infoBox;
	infoBox.setText("Global features of the neuron:");
	infoBox.setInformativeText(QString("<pre><font size='4'>"
				"number of nodes                  : %1<br>"
				"soma surface                     : %2<br>"
				"number of stems                  : %3<br>"
				"number of bifurcations           : %4<br>"
				"number of branches               : %5<br>"
				"number of tips                   : %6<br>"
				"overall width                    : %7<br>"
				"overall height                   : %8<br>"
				"overall depth                    : %9<br>"
				"average diameter                 : %10<br>"
				"total length                     : %11<br>"
				"total surface                    : %12<br>"
				"total volume                     : %13<br>"
				"max euclidean distance           : %14<br>"
				"max path distance                : %15<br>"
				"max branch order                 : %16<br>"
				"average contraction              : %17<br>"
				"average fragmentation            : %18<br>"
				"average parent-daughter ratio    : %19<br>"
				"average bifurcation angle local  : %20<br>"
				"average bifurcation angle remote : %21<br>"
				"Hausdorff dimension              : %22</font></pre>")
				.arg(features[0])
				.arg(features[1])
				.arg(features[2])
				.arg(features[3])
				.arg(features[4])
				.arg(features[5])
				.arg(features[6])
				.arg(features[7])
				.arg(features[8])
				.arg(features[9])
				.arg(features[10])
				.arg(features[11])
				.arg(features[12])
				.arg(features[13])
				.arg(features[14])
				.arg(features[15])
				.arg(features[16])
				.arg(features[17])
				.arg(features[18])
				.arg(features[19])
				.arg(features[20])
				.arg(features[21]));
	infoBox.exec();


	if (features) {delete []features; features = NULL;}

	delete openDlg; // MK, Oct, 2017, free up dialog pointer to memory violation.
}

void nf_main_batch(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString inputSwcFolder = QFileDialog::getExistingDirectory(parent,
                                                           QString(QObject::tr("Choose the directory that including all swcs."))
                                                           );
    if(inputSwcFolder.size()==0){
        v3d_msg("Empty input folder.\n Let the developer know if you see this message.");
        return;
    }
    if(!inputSwcFolder.endsWith("/")){inputSwcFolder = inputSwcFolder + "/";}

    QString CsvName;
    CsvName = QFileDialog::getOpenFileName(0,QObject::tr("Open empty csv file(output result)."),"",QObject::tr("Supported file(*.csv)"));
    if(CsvName.isEmpty())
        return;

//    QString outputFolder = QFileDialog::getExistingDirectory(parent,
//                                                              QString(QObject::tr("Choose the output file folder")));

//    if((outputFolder.size()>0) && (!outputFolder.endsWith("/"))){outputFolder = outputFolder + "/";}
//    if(outputFolder.size()==0){outputFolder = inputSwcFolder;}
    QStringList nameFilters;
    nameFilters<<"*.swc";
    nameFilters<<"*.eswc";
    nameFilters<<"*.ESWC";
    nameFilters<<"*.SWC";

    QDir inputDir(inputSwcFolder);
    QStringList nameList = inputDir.entryList(nameFilters,QDir::Files|QDir::Readable, QDir::Name);
    vector<NeuronTree> nt_list;
    if(nameList.size()==0)
    {
        cerr<<"Error in input. This is an empty folder or I can't find swc and eswc files in it"<<endl;
        return;
    }
    for(V3DLONG i=0;i<nameList.size();i++)
    {
        QString thisname=inputSwcFolder+"/"+nameList.at(i);
        NeuronTree tmp = readSWC_file(thisname);
        nt_list.push_back(tmp);
    }



    ofstream csvOutFile(CsvName.toStdString().c_str());
    if(!csvOutFile.is_open()){
         cerr<<"out Error: cannot open file to save"<<endl;
         return;
    }
    csvOutFile<<"Name,Nodes,SomaSurface,Stems,Bifurcations,Branches,Tips,OverallWidth,OverallHeight,OverallDepth";
    csvOutFile<<",AverageDiameter,Length,Surface,Volume,MaxEuclideanDistance,MaxPathDistance,MaxBranchOrder";
    csvOutFile<<",AverageContraction,AverageFragmentation,AverageParent-daughterRatio,AverageBifurcationAngleLocal,AverageBifurcationAngleRemote,HausdorffDimension"<<endl;

    for (V3DLONG i=0;i<nameList.size();i++)
    {
        NeuronTree nt = nt_list[i];
        cout<<"\n----------Neuron #"<<(i+1)<<",Total #"<<nameList.size()<<"-----------\n";
        double * features = new double[FNUM];
        computeFeature(nt, features);
        //save features to file

        QString thisname=nameList[i];
        csvOutFile<<thisname.toStdString();
        for (int i=0;i<FNUM;i++)
        {
            csvOutFile<<","<<features[i];
        }
        csvOutFile<<endl;

        if (features) {delete []features; features = NULL;}
    }
    csvOutFile.close();
    return;


}
// Added by zll for global neuron batch in gui.07.22.2022
void nf_first_main(V3DPluginCallback2 &callback, QWidget *parent)
{
    OpenSWCDialog * openDlg = new OpenSWCDialog(0, &callback);
    if (!openDlg->exec())
        return;

    NeuronTree nt = openDlg->nt;
    if(nt.listNeuron[0].pn >=0)
    {
        v3d_msg("Please sort the swc file first");
        return;
    }
    QList<NeuronSWC> list = nt.listNeuron;
    for (int i=0;i<list.size();i++)
    {
        if(i>0 && nt.listNeuron[i].pn < 0)
        {
            nt.listNeuron.erase(nt.listNeuron.begin()+i,nt.listNeuron.end());
            break;
        }
    }

    double * features = new double[FNUM];
    computeFeature(nt,features);
    QMessageBox infoBox;
    infoBox.setText("Global features of the neuron:");
    infoBox.setInformativeText(QString("<pre><font size='4'>"
                "number of nodes                  : %1<br>"
                "soma surface                     : %2<br>"
                "number of stems                  : %3<br>"
                "number of bifurcations           : %4<br>"
                "number of branches               : %5<br>"
                "number of tips                   : %6<br>"
                "overall width                    : %7<br>"
                "overall height                   : %8<br>"
                "overall depth                    : %9<br>"
                "average diameter                 : %10<br>"
                "total length                     : %11<br>"
                "total surface                    : %12<br>"
                "total volume                     : %13<br>"
                "max euclidean distance           : %14<br>"
                "max path distance                : %15<br>"
                "max branch order                 : %16<br>"
                "average contraction              : %17<br>"
                "average fragmentation            : %18<br>"
                "average parent-daughter ratio    : %19<br>"
                "average bifurcation angle local  : %20<br>"
                "average bifurcation angle remote : %21<br>"
                "Hausdorff dimension              : %22</font></pre>")
                .arg(features[0])
                .arg(features[1])
                .arg(features[2])
                .arg(features[3])
                .arg(features[4])
                .arg(features[5])
                .arg(features[6])
                .arg(features[7])
                .arg(features[8])
                .arg(features[9])
                .arg(features[10])
                .arg(features[11])
                .arg(features[12])
                .arg(features[13])
                .arg(features[14])
                .arg(features[15])
                .arg(features[16])
                .arg(features[17])
                .arg(features[18])
                .arg(features[19])
                .arg(features[20])
                .arg(features[21]));
    infoBox.exec();


    if (features) {delete []features; features = NULL;}
	delete openDlg; // MK, Oct, 2017, free up dialog pointer to memory violation.
}

void printFeature(double * features)
{

	for (int i=0;i<FNUM;i++)
	{
		switch (i)
		{
			case 0:
				cout<<"Number of Nodes:";
				break;
			case 1:
				cout<<"Soma Surface:\t";
				break;
			case 2:
				cout<<"Number of Stems:";
				break;
			case 3:
				cout<<"Number of Bifurcatons:";
				break;
			case 4:
				cout<<"Number of Branches:";
				break;
			case 5:
				cout<<"Number of Tips:\t";
				break;
			case 6:
				cout<<"Overall Width:\t";
				break;
			case 7:
				cout<<"Overall Height:\t";
				break;
			case 8:
				cout<<"Overall Depth:\t";
				break;
			case 9:
				cout<<"Average Diameter:";
				break;
			case 10:
				cout<<"Total Length:\t";
				break;
			case 11:
				cout<<"Total Surface:\t";
				break;
			case 12:
				cout<<"Total Volume:\t";
				break;
			case 13:
				cout<<"Max Euclidean Distance:";
				break;
			case 14:
				cout<<"Max Path Distance:\t\t";
				break;
			case 15:
				cout<<"Max Branch Order:\t\t";
				break;
			case 16:
				cout<<"Average Contraction:\t\t";
				break;
			case 17:
				cout<<"Average Fragmentation:\t\t";
				break;
			case 18:
				cout<<"Average Parent-daughter Ratio:\t";
				break;
			case 19:
				cout<<"Average Bifurcation Angle Local:";
				break;
			case 20:
				cout<<"Average Bifurcation Angle Remote:";
				break;
			case 21:
				cout<<"Hausdorff Dimension:\t\t";
		}
		cout<<"\t"<<features[i]<<endl;
	}
}

void nf_toolbox(const V3DPluginArgList & input)
{
	vaa3d_neurontoolbox_paras * paras = (vaa3d_neurontoolbox_paras *)(input.at(0).p);
	NeuronTree nt = paras->nt;
	QString fileOpenName = nt.file;
	
	double * features = new double[FNUM];
	computeFeature(nt,features);
	QMessageBox infoBox;
	infoBox.setText("Global features of the neuron:");
	infoBox.setInformativeText(QString("<pre><font size='4'>"
				"number of nodes                  : %1<br>"
				"soma surface                     : %2<br>"
				"number of stems                  : %3<br>"
				"number of bifurcations           : %4<br>"
				"number of branches               : %5<br>"
				"number of tips                   : %6<br>"
				"overall width                    : %7<br>"
				"overall height                   : %8<br>"
				"overall depth                    : %9<br>"
				"average diameter                 : %10<br>"
				"total length                     : %11<br>"
				"total surface                    : %12<br>"
				"total volume                     : %13<br>"
				"max euclidean distance           : %14<br>"
				"max path distance                : %15<br>"
				"max branch order                 : %16<br>"
				"average contraction              : %17<br>"
				"average fragmentation            : %18<br>"
				"average parent-daughter ratio    : %19<br>"
				"average bifurcation angle local  : %20<br>"
				"average bifurcation angle remote : %21<br>"
				"Hausdorff dimension              : %22</font></pre>")
				.arg(features[0])
				.arg(features[1])
				.arg(features[2])
				.arg(features[3])
				.arg(features[4])
				.arg(features[5])
				.arg(features[6])
				.arg(features[7])
				.arg(features[8])
				.arg(features[9])
				.arg(features[10])
				.arg(features[11])
				.arg(features[12])
				.arg(features[13])
				.arg(features[14])
				.arg(features[15])
				.arg(features[16])
				.arg(features[17])
				.arg(features[18])
				.arg(features[19])
				.arg(features[20])
				.arg(features[21]));
	infoBox.exec();


	if (features) {delete []features; features = NULL;}

}
