//last changed by Hanchuan Peng, 2013-08-06. fix a return value bug and also add some v3d_msg calls to indicate where the code should be added

#ifndef __CREATE_PLUGIN_H__
#define __CREATE_PLUGIN_H__
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

string toupper(string instr)
{
	int n = instr.size();
	string new_str = instr;
	for(int i = 0; i < n; i++)
	{
		new_str[i] = (char) ::toupper(instr.at(i));
	}
	return new_str;
}

struct PluginTemplate
{
	string PLUGIN_NAME;
	string PLUGIN_CLASS;
	string PLUGIN_DESCRIPTION;
	string PLUGIN_DATE;
	string PLUGIN_AUTHOR;
	string PLUGIN_GUI;                              // plugin gui file name,     test_gui.h
	string PLUGIN_HEADER;                           // plugin header file name , test_plugin.h
	string PLUGIN_CPP;
	string FUNC_HEADER;
	string FUNC_CPP;
	string PRO_FILE;

	string WINDOW_TITLE;
	string VAA3D_PATH;
	bool DOFUNC;
    vector<string> PARA_NAME;
    vector<string> PARA_TYPE;
    vector<string> PARA_VALUE;

	vector<string> MAINFUNCS;
	vector<string> SYSINVOKES;

    string OUTPUTSWC;
    string FUNC_NAME;
    string TRACINGPLUGIN_NAME;
};


void create_plugin_pro(PluginTemplate & pt)
{
	cout<<"create "<<pt.PRO_FILE<<" ... "<<endl;
	ofstream ofs(pt.PRO_FILE.c_str());
	ofs<<""<<endl;
	ofs<<"TEMPLATE\t= lib"<<endl;
	ofs<<"CONFIG\t+= qt plugin warn_off"<<endl;
	ofs<<"#CONFIG\t+= x86_64"<<endl;
	ofs<<"VAA3DPATH = "<<pt.VAA3D_PATH<<endl;
    ofs<<"INCLUDEPATH\t+= $$VAA3DPATH/basic_c_fun"<<endl;
    ofs<<"INCLUDEPATH\t+= $$VAA3DPATH/common_lib/include"<<endl;
	ofs<<endl;
	ofs<<"HEADERS\t+= "<<pt.PLUGIN_HEADER<<endl;
	if(pt.PLUGIN_GUI != "") ofs<<"HEADERS\t+= "<<pt.PLUGIN_GUI<<endl; 
 //   ofs<<"SOURCES\t+= $$VAA3DPATH/../../vaa3d_tools/hackathon/zhi/APP2_large_scale/readrawfile_func.h"<<endl;

    ofs<<""<<endl;

	ofs<<"SOURCES\t+= "<<pt.PLUGIN_CPP<<endl;
    ofs<<"SOURCES\t+= $$VAA3DPATH/basic_c_fun/v3d_message.cpp"<<endl;
    ofs<<"SOURCES\t+= $$VAA3DPATH/basic_c_fun/basic_surf_objs.cpp"<<endl;
    ofs<<"SOURCES\t+= $$VAA3DPATH/../released_plugins_more/v3d_plugins/neurontracing_vn2/app2/my_surf_objs.cpp"<<endl;

    ofs<<"SOURCES\t+= $$VAA3DPATH/basic_c_fun/stackutil.cpp"<<endl;
    ofs<<"SOURCES\t+= $$VAA3DPATH/basic_c_fun/mg_utilities.cpp"<<endl;
    ofs<<"SOURCES\t+= $$VAA3DPATH/basic_c_fun/mg_image_lib.cpp"<<endl;
    ofs<<"SOURCES\t+= $$VAA3DPATH/../../vaa3d_tools/hackathon/zhi/APP2_large_scale/readrawfile_func.cpp"<<endl;
	ofs<<""<<endl;

    ofs<<"LIBS\t+= -lm -L$$VAA3DPATH/common_lib/lib -lv3dtiff"<<endl;
    ofs<<"LIBS\t+= -lpthread"<<endl;
    ofs<<"LIBS\t+= -lv3dfftw3f -lv3dfftw3f_threads"<<endl;

    ofs<<"TARGET\t= $$qtLibraryTarget("<<pt.PLUGIN_NAME<<")"<<endl;
    ofs<<"DESTDIR\t= $$VAA3DPATH/../bin/plugins/neuron_tracing/"<<pt.PLUGIN_NAME<<"/"<<endl;
	ofs.close();
}

void create_plugin_cpp(PluginTemplate & pt)
{
	cout<<"create "<<pt.PLUGIN_CPP<<" ... "<<endl;
	ofstream ofs(pt.PLUGIN_CPP.c_str());
    string line;
    string template_path = pt.VAA3D_PATH + "/../../vaa3d_tools/hackathon/zhi/neuronassembler_plugin_creator/neuronassembler_template.cpp";
    ifstream templatefile (template_path.c_str());
    if (templatefile.is_open())
    {
        for(int i = 0; i< 7; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"#include \""<<pt.PLUGIN_HEADER<<"\""<<endl;
        ofs<<"#include \""<<pt.VAA3D_PATH<< "/../../vaa3d_tools/hackathon/zhi/APP2_large_scale/readRawfile_func.h"<<"\""<<endl;
        ofs<<"#include \""<<pt.VAA3D_PATH<< "/../released_plugins_more/v3d_plugins/istitch/y_imglib.h"<<"\""<<endl;
        ofs<<"#include \""<<pt.VAA3D_PATH<< "/../released_plugins_more/v3d_plugins/neurontracing_vn2/app2/my_surf_objs.h"<<"\""<<endl;
        ofs<<"#include \""<<pt.VAA3D_PATH<< "/../../vaa3d_tools/hackathon/zhi/neuronassembler_plugin_creator/sort_swc.h"<<"\""<<endl;

        ofs<<"Q_EXPORT_PLUGIN2("<<pt.PLUGIN_NAME<<", "<<pt.PLUGIN_CLASS<<");"<<endl;

        for(int i = 7;i< 31; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }


        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                     ofs<<"\tQString "<<pt.PARA_NAME.at(i)<<";"<<endl;
                }
                else
                {
                    ofs<<"\t"<<pt.PARA_TYPE.at(i)<<" "<<pt.PARA_NAME.at(i)<<";"<<endl;
                }
            }

        }

        for(int i = 31;i< 114; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }


        ofs<<"QStringList "<<pt.PLUGIN_CLASS<<"::menulist() const"<<endl;

        for(int i = 114;i< 122; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"QStringList "<<pt.PLUGIN_CLASS<<"::funclist() const"<<endl;

        for(int i = 123;i< 130; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"void "<<pt.PLUGIN_CLASS<<"::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)"<<endl;
        for(int i = 130;i< 155; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                  ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = dialog."<<pt.PARA_NAME.at(i)<<";"<<endl;
            }
        }

        for(int i = 155;i< 194; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }


        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                  ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = dialog."<<pt.PARA_NAME.at(i)<<";"<<endl;
            }
        }

        for(int i = 194;i< 202; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\tv3d_msg(tr(\""<<pt.PLUGIN_DESCRIPTION<<". \"\n\t\t\t\"Developed by "<<pt.PLUGIN_AUTHOR<<", "<<pt.PLUGIN_DATE<<"\"));"<<endl;
        ofs<<"\t}"<<endl;
        ofs<<"}"<<endl;
        ofs<<""<<endl;

        ofs<<"bool "<<pt.PLUGIN_CLASS<<"::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)"<<endl;

        for(int i = 202;i< 274; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                    ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = (paras.size() >= k+1) ? paras[k]:"<< pt.PARA_VALUE.at(i) <<"; k++;"<<endl;
                }
                else if(pt.PARA_TYPE.at(i) == "double")
                {
                    ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = (paras.size() >= k+1) ? atof(paras[k]):"<< pt.PARA_VALUE.at(i) <<"; k++;"<<endl;
                }
                else
                    ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = (paras.size() >= k+1) ? atoi(paras[k]):"<< pt.PARA_VALUE.at(i) <<"; k++;"<<endl;

            }
        }

        for(int i = 274;i< 313; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                    ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = (paras.size() >= k+1) ? paras[k]:"<< pt.PARA_VALUE.at(i) <<"; k++;"<<endl;
                }
                else if(pt.PARA_TYPE.at(i) == "double")
                {
                    ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = (paras.size() >= k+1) ? atof(paras[k]):"<< pt.PARA_VALUE.at(i) <<"; k++;"<<endl;
                }
                else
                    ofs<<"\t\tP."<<pt.PARA_NAME.at(i)<<" = (paras.size() >= k+1) ? atoi(paras[k]):"<< pt.PARA_VALUE.at(i) <<"; k++;"<<endl;

            }
        }

        for(int i = 313;i< 319; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\tprintf(\"vaa3d -x " << pt.PLUGIN_NAME <<" -f trace_tc -i <inimg_file> -p <inmarker_file> <tc_file> <block size>";
        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                ofs<<" <" << pt.PARA_NAME.at(i) <<">";
            }
        }

        ofs<<"\\n\");\n";
        ofs<<"\t\tprintf(\"inimg_file\t\tShould be 8 bit image\\n\");"<<endl;
        ofs<<"\t\tprintf(\"inmarker_file\t\tPlease specify the path of the marker file\\n\");"<<endl;
        ofs<<"\t\tprintf(\"tc_file\t\t\tPlease specify the path of the tc file\\n\");"<<endl;
        ofs<<"\t\tprintf(\"block size\t\tDefault 1024\\n\");"<<endl;
        ofs<<"\n";
        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_NAME.at(i).size() > 7)
                     ofs<<"\t\tprintf(\"" << pt.PARA_NAME.at(i) <<"\t\tRequired by the tracing algorithm. Default value is "<<pt.PARA_VALUE.at(i)<<"\\n\");"<<endl;
                else
                    ofs<<"\t\tprintf(\"" << pt.PARA_NAME.at(i) <<"\t\t\tRequired by the tracing algorithm. Default value is "<<pt.PARA_VALUE.at(i)<<"\\n\");"<<endl;

            }
        }
        ofs<<"\n";
        ofs<<"\t\tprintf(\"outswc_file\t\tWill be named automatically based on the input image file name, so you don't have to specify it.\\n\\n\");\n\n"<<endl;


        ofs<<"\t\tprintf(\"vaa3d -x " << pt.PLUGIN_NAME <<" -f trace_raw -i <inimg_file> -p <inmarker_file> <block size> <tracing_entire_image>";
        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                ofs<<" <" << pt.PARA_NAME.at(i) <<">";
            }
        }

        ofs<<"\\n\");\n";
        ofs<<"\t\tprintf(\"inimg_file\t\tShould be 8 bit v3draw/raw image\\n\");"<<endl;
        ofs<<"\t\tprintf(\"inmarker_file\t\tPlease specify the path of the marker file, Default value is NULL\\n\");"<<endl;
        ofs<<"\t\tprintf(\"block size\t\tDefault 1024\\n\");"<<endl;
        ofs<<"\t\tprintf(\"tracing_entire_image\tYES:1, NO:0. Default value is 0\\n\");"<<endl;
        ofs<<"\n";
        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_NAME.at(i).size() > 7)
                     ofs<<"\t\tprintf(\"" << pt.PARA_NAME.at(i) <<"\t\tRequired by the tracing algorithm. Default value is "<<pt.PARA_VALUE.at(i)<<"\\n\");"<<endl;
                else
                    ofs<<"\t\tprintf(\"" << pt.PARA_NAME.at(i) <<"\t\t\tRequired by the tracing algorithm. Default value is "<<pt.PARA_VALUE.at(i)<<"\\n\");"<<endl;

            }
        }
        ofs<<"\n";
        ofs<<"\t\tprintf(\"outswc_file\t\tWill be named automatically based on the input image file name, so you don't have to specify it.\\n\\n\");\n\n"<<endl;



        for(int i = 319;i< 430; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\tQString finalswcfilename = fileOpenName + rootposstr + \"_"<< pt.PLUGIN_NAME <<".swc\";"<<endl;
        ofs<<"\tQString tmpfolder = QFileInfo(tcfile).path()+(\"/tmp_"<< pt.PLUGIN_NAME <<"\");"<<endl;
        ofs<<"\thead->tilename = QFileInfo(tcfile).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;


        for(int i = 430;i< 455; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\tQString swcfilename =  walker->tilename + QString(\""<< pt.OUTPUTSWC <<"\");"<<endl;

        for(int i = 455;i< 475; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                    ofs<<"\t\tstring S_"<<pt.PARA_NAME.at(i) <<" = "<<"P."<<pt.PARA_NAME.at(i) <<".toStdString();"<<endl;
                }
                else
                {
                    ofs<<"\t\tstring S_"<<pt.PARA_NAME.at(i) <<" = boost::lexical_cast<string>(P." << pt.PARA_NAME.at(i) <<");"<<endl;
                }

                ofs<<"\t\tchar* C_"<<pt.PARA_NAME.at(i) <<" = new char[S_"<<pt.PARA_NAME.at(i) <<".length() + 1];"<<endl;
                ofs<<"\t\tstrcpy(C_"<<pt.PARA_NAME.at(i) <<",S_"<<pt.PARA_NAME.at(i) <<".c_str());"<<endl;
                ofs<<"\t\targ_para.push_back(C_"<<pt.PARA_NAME.at(i) <<");\n"<<endl;
            }
        }

        ofs<<"\t\tfull_plugin_name = \""<<pt.TRACINGPLUGIN_NAME <<"\";  func_name =  \"" <<pt.FUNC_NAME<<"\";"<<endl;

        for(int i = 475;i< 542; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;

        for(int i = 542;i< 569; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;

        for(int i = 569;i< 596; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;

        for(int i = 596;i< 622; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;

        for(int i = 622;i< 746; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }




        ofs<<"\tQString finalswcfilename = fileOpenName + rootposstr + \"_"<< pt.PLUGIN_NAME <<".swc\";"<<endl;
        ofs<<"\tQString tmpfolder = QFileInfo(fileOpenName).path()+(\"/tmp_"<< pt.PLUGIN_NAME <<"\");"<<endl;
        ofs<<"\thead->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;


        for(int i = 746;i< 785; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\tQString swcfilename =  walker->tilename + QString(\""<< pt.OUTPUTSWC <<"\");"<<endl;


        for(int i = 785;i< 805; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                    ofs<<"\t\tstring S_"<<pt.PARA_NAME.at(i) <<" = "<<"P."<<pt.PARA_NAME.at(i) <<".toStdString();"<<endl;
                }
                else
                {
                    ofs<<"\t\tstring S_"<<pt.PARA_NAME.at(i) <<" = boost::lexical_cast<string>(P." << pt.PARA_NAME.at(i) <<");"<<endl;
                }

                ofs<<"\t\tchar* C_"<<pt.PARA_NAME.at(i) <<" = new char[S_"<<pt.PARA_NAME.at(i) <<".length() + 1];"<<endl;
                ofs<<"\t\tstrcpy(C_"<<pt.PARA_NAME.at(i) <<",S_"<<pt.PARA_NAME.at(i) <<".c_str());"<<endl;
                ofs<<"\t\targ_para.push_back(C_"<<pt.PARA_NAME.at(i) <<");\n"<<endl;
            }
        }

        ofs<<"\t\tfull_plugin_name = \""<<pt.TRACINGPLUGIN_NAME <<"\";  func_name =  \"" <<pt.FUNC_NAME<<"\";"<<endl;

        for(int i = 805;i< 872; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;


        for(int i = 872;i< 899; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;


        for(int i = 899;i< 926; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }


        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;


        for(int i = 926;i< 952; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        ofs<<"\t\t\t\tnewNode->tilename = QFileInfo(fileOpenName).path().append(\"/tmp_"<< pt.PLUGIN_NAME <<"/\").append(QString(region_name));"<<endl;


        for(int i = 952;i< 1250; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }


        templatefile.close();

    }

//	ofs<<"/* "<<pt.PLUGIN_CPP<<endl;
//	ofs<<" * "<<pt.PLUGIN_DESCRIPTION<<endl;
//	ofs<<" * "<<pt.PLUGIN_DATE<<" : by "<<pt.PLUGIN_AUTHOR<<endl;
//	ofs<<" */"<<endl;
//	ofs<<" "<<endl;
//	ofs<<"#include \"v3d_message.h\""<<endl;
//	ofs<<"#include <vector>"<<endl;
//	ofs<<"#include \""<<pt.PLUGIN_HEADER<<"\""<<endl;
//	ofs<<"using namespace std;"<<endl;
//	ofs<<"Q_EXPORT_PLUGIN2("<<pt.PLUGIN_NAME<<", "<<pt.PLUGIN_CLASS<<");"<<endl;
//	ofs<<" "<<endl;
//	ofs<<"QStringList "<<pt.PLUGIN_CLASS<<"::menulist() const"<<endl;
//	ofs<<"{"<<endl;
//	ofs<<"\treturn QStringList() "<<endl;
//	for(int i = 0; i < pt.MENUS.size(); i++) ofs<<"\t\t<<tr(\""<<pt.MENUS[i]<<"\")"<<endl;
//	ofs<<"\t\t<<tr(\"about\");"<<endl;
//	ofs<<"}"<<endl;
//	ofs<<""<<endl;
//	ofs<<"QStringList "<<pt.PLUGIN_CLASS<<"::funclist() const"<<endl;
//	ofs<<"{"<<endl;
//	ofs<<"\treturn QStringList()";
//	for(int i = 0; i < pt.FUNCS.size(); i++) ofs<<endl<<"\t\t<<tr(\""<<pt.FUNCS[i]<<"\")";
//	ofs<<";"<<endl;
//	ofs<<"}"<<endl;
//	ofs<<""<<endl;
//	ofs<<"void "<<pt.PLUGIN_CLASS<<"::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)"<<endl;
//	ofs<<"{"<<endl;
//	ofs<<"\tif (menu_name == tr(\""<<pt.MENUS[0]<<"\"))"<<endl;
//	ofs<<"\t{"<<endl;
//    ofs<<"\t\tv3d_msg(\"To be implemented.\");"<<endl;
//	ofs<<"\t}"<<endl;
//	for(int i = 1; i < pt.MENUS.size(); i++)
//	{
//		ofs<<"\telse if (menu_name == tr(\""<<pt.MENUS[i]<<"\"))"<<endl;
//		ofs<<"\t{"<<endl;
//  	    ofs<<"\t\tv3d_msg(\"To be implemented.\");"<<endl;
//		ofs<<"\t}"<<endl;
//	}
//	ofs<<"\telse"<<endl;
//	ofs<<"\t{"<<endl;
//	ofs<<"\t\tv3d_msg(tr(\""<<pt.PLUGIN_DESCRIPTION<<". \"\n\t\t\t\"Developed by "<<pt.PLUGIN_AUTHOR<<", "<<pt.PLUGIN_DATE<<"\"));"<<endl;
//	ofs<<"\t}"<<endl;
//	ofs<<"}"<<endl;
//	ofs<<""<<endl;

//	ofs<<"bool "<<pt.PLUGIN_CLASS<<"::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)"<<endl;
//	ofs<<"{"<<endl;
//	ofs<<"\tvector<char*> infiles, inparas, outfiles;"<<endl;
//	ofs<<"\tif(input.size() >= 1) infiles = *((vector<char*> *)input.at(0).p);"<<endl;
//	ofs<<"\tif(input.size() >= 2) inparas = *((vector<char*> *)input.at(1).p);"<<endl;
//	ofs<<"\tif(output.size() >= 1) outfiles = *((vector<char*> *)output.at(0).p);"<<endl;
//	ofs<<endl;
//	ofs<<"\tif (func_name == tr(\""<<pt.FUNCS[0]<<"\"))"<<endl;
//	ofs<<"\t{"<<endl;
//	ofs<<"\t\tv3d_msg(\"To be implemented.\");"<<endl;
//	ofs<<"\t}"<<endl;
//	for(int i = 1; i < pt.FUNCS.size(); i++)
//	{
//		ofs<<"\telse if (func_name == tr(\""<<pt.FUNCS[i]<<"\"))"<<endl;
//		ofs<<"\t{"<<endl;
//  	    ofs<<"\t\tv3d_msg(\"To be implemented.\");"<<endl;
//		ofs<<"\t}"<<endl;
//	}
//	ofs<<"\telse return false;"<<endl;
//	ofs<<endl<<"\treturn true;"<<endl;
//	ofs<<"}"<<endl;
//	ofs<<""<<endl;
	ofs.close();
}

void create_plugin_header(PluginTemplate & pt)  // PLUGIN_HEADER
{
	cout<<"create "<<pt.PLUGIN_HEADER<<" ... "<<endl;
	ofstream ofs(pt.PLUGIN_HEADER.c_str());
	ofs<<"/* "<<pt.PLUGIN_HEADER<<endl;
	ofs<<" * "<<pt.PLUGIN_DESCRIPTION<<endl;
	ofs<<" * "<<pt.PLUGIN_DATE<<" : by "<<pt.PLUGIN_AUTHOR<<endl;
	ofs<<" */"<<endl;
	ofs<<" "<<endl;
	ofs<<"#ifndef __"<<toupper(pt.PLUGIN_NAME)<<"_PLUGIN_H__"<<endl;
	ofs<<"#define __"<<toupper(pt.PLUGIN_NAME)<<"_PLUGIN_H__"<<endl;
	ofs<<""<<endl;
	ofs<<"#include <QtGui>"<<endl;
	ofs<<"#include <v3d_interface.h>"<<endl;
	ofs<<""<<endl;
	ofs<<"class "<<pt.PLUGIN_CLASS<<" : public QObject, public V3DPluginInterface2_1"<<endl;
	ofs<<"{"<<endl;
	ofs<<"\tQ_OBJECT"<<endl;
	ofs<<"\tQ_INTERFACES(V3DPluginInterface2_1);"<<endl;
	ofs<<""<<endl;
	ofs<<"public:"<<endl;
	ofs<<"\tfloat getPluginVersion() const {return 1.1f;}"<<endl;
	ofs<<""<<endl;
	ofs<<"\tQStringList menulist() const;"<<endl;
	ofs<<"\tvoid domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);"<<endl;
	ofs<<""<<endl;
	if(!pt.DOFUNC)
	{
		ofs<<"\tQStringList funclist() const {return QStringList();}"<<endl;
		ofs<<"\tbool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent)"<<endl;
		ofs<<"\t{return false;}"<<endl;
	}
	else
	{
		ofs<<"\tQStringList funclist() const ;"<<endl; 
		ofs<<"\tbool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);"<<endl;
	}
	ofs<<"};"<<endl;
	ofs<<""<<endl;
	ofs<<"#endif"<<endl;
	ofs<<""<<endl;

    string line;
    string template_path = pt.VAA3D_PATH + "/../../vaa3d_tools/hackathon/zhi/neuronassembler_plugin_creator/neuronassembler_template.h";
    ifstream templatefile (template_path.c_str());
    if (templatefile.is_open())
    {
        for(int i = 0; i< 37; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        int d = 0;
        if(pt.PARA_NAME.size() >0)
        {
            for(d = 0; d < pt.PARA_NAME.size() ; d++)
            {
                if(pt.PARA_TYPE.at(d) == "string")
                {
                    ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_edit = new QLineEdit();"<<endl;
                    ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_edit->setText(\""<<pt.PARA_VALUE.at(d)<<"\");"<<endl;
                    ofs<<"\t\t\tlayout->addWidget(new QLabel(QObject::tr(\""<<pt.PARA_NAME.at(d)<<":\"))," << d+4 <<",0);"<<endl;
                    ofs<<"\t\t\tlayout->addWidget("<<pt.PARA_NAME.at(d)<<"_edit,"<<d+4<<",1,1,5);\n"<<endl;

                }
                else
                {
                    if(pt.PARA_TYPE.at(d) == "int")
                        ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_box = new QSpinBox();"<<endl;
                    else
                        ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_box = new QDoubleSpinBox();"<<endl;
                    ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_box->setValue("<<pt.PARA_VALUE.at(d)<<");"<<endl;
                    ofs<<"\t\t\tlayout->addWidget(new QLabel(QObject::tr(\""<<pt.PARA_NAME.at(d)<<":\"))," << d+4 <<",0);"<<endl;
                    ofs<<"\t\t\tlayout->addWidget("<<pt.PARA_NAME.at(d)<<"_box,"<<d+4<<",1,1,5);\n"<<endl;

                }
            }

        }

        ofs<<"\t\t\tlayout->addLayout(hbox2,"<<d+4<<",0,3,6);"<<endl;
        ofs<<"\t\t\tsetWindowTitle(QString(\"Vaa3D-NeuronAssembler("<<pt.TRACINGPLUGIN_NAME <<")\"));"<<endl;

        for(int i = 37; i< 46; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) != "string")
                {
                    ofs<<"\t\t\tconnect("<<pt.PARA_NAME.at(i)<<"_box, SIGNAL(valueChanged("<< pt.PARA_TYPE.at(i)<<")), this, SLOT(update()));"<<endl;
                }
            }

        }

        for(int i = 46; i< 60; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                        ofs<<"\t\t\t"<< pt.PARA_NAME.at(i)<<" = "<<pt.PARA_NAME.at(i) <<"_edit->text();"<<endl;
                }
                else
                {
                        ofs<<"\t\t\t"<< pt.PARA_NAME.at(i)<<" = "<<pt.PARA_NAME.at(i) <<"_box->value();"<<endl;

                }
            }

        }
        for(int i = 60; i< 117; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                     ofs<<"\t\tQLineEdit *"<<pt.PARA_NAME.at(i)<<"_edit;"<<endl;
                     ofs<<"\t\tQString "<<pt.PARA_NAME.at(i)<<";"<<endl;
                }
                else
                {
                    if(pt.PARA_TYPE.at(i) == "int")
                        ofs<<"\t\tQSpinBox *"<<pt.PARA_NAME.at(i)<<"_box;"<<endl;
                    else
                        ofs<<"\t\tQDoubleSpinBox *"<<pt.PARA_NAME.at(i)<<"_box;"<<endl;
                    ofs<<"\t\t"<<pt.PARA_TYPE.at(i)<<" "<<pt.PARA_NAME.at(i)<<";"<<endl;
                }
            }

        }

        for(int i = 117; i< 174; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(d = 0; d < pt.PARA_NAME.size() ; d++)
            {
                if(pt.PARA_TYPE.at(d) == "string")
                {
                    ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_edit = new QLineEdit();"<<endl;
                    ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_edit->setText(\""<<pt.PARA_VALUE.at(d)<<"\");"<<endl;
                    ofs<<"\t\t\tlayout->addWidget(new QLabel(QObject::tr(\""<<pt.PARA_NAME.at(d)<<":\"))," << d+2 <<",0);"<<endl;
                    ofs<<"\t\t\tlayout->addWidget("<<pt.PARA_NAME.at(d)<<"_edit,"<<d+2<<",1,1,5);\n"<<endl;

                }
                else
                {
                    if(pt.PARA_TYPE.at(d) == "int")
                        ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_box = new QSpinBox();"<<endl;
                    else
                        ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_box = new QDoubleSpinBox();"<<endl;

                    ofs<<"\t\t\t"<<pt.PARA_NAME.at(d)<<"_box->setValue("<<pt.PARA_VALUE.at(d)<<");"<<endl;
                    ofs<<"\t\t\tlayout->addWidget(new QLabel(QObject::tr(\""<<pt.PARA_NAME.at(d)<<":\"))," << d+2 <<",0);"<<endl;
                    ofs<<"\t\t\tlayout->addWidget("<<pt.PARA_NAME.at(d)<<"_box,"<<d+2<<",1,1,5);\n"<<endl;

                }
            }

        }

        ofs<<"\t\t\tlayout->addLayout(hbox2,"<<d+2<<",0,2,6);"<<endl;
        ofs<<"\t\t\tsetWindowTitle(QString(\"Vaa3D-NeuronAssembler("<<pt.TRACINGPLUGIN_NAME <<")\"));"<<endl;

        for(int i = 174; i< 181; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) != "string")
                {
                    ofs<<"\t\t\tconnect("<<pt.PARA_NAME.at(i)<<"_box, SIGNAL(valueChanged("<< pt.PARA_TYPE.at(i)<<")), this, SLOT(update()));"<<endl;
                }
            }

        }

        for(int i = 181; i< 194; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                        ofs<<"\t\t\t"<< pt.PARA_NAME.at(i)<<" = "<<pt.PARA_NAME.at(i) <<"_edit->text();"<<endl;
                }
                else
                {
                        ofs<<"\t\t\t"<< pt.PARA_NAME.at(i)<<" = "<<pt.PARA_NAME.at(i) <<"_box->value();"<<endl;

                }
            }

        }
        for(int i = 194; i< 225; i++)
        {
            getline (templatefile,line);
            ofs<<line<<endl;
        }

        if(pt.PARA_NAME.size() >0)
        {
            for(int i = 0; i < pt.PARA_NAME.size() ; i++)
            {
                if(pt.PARA_TYPE.at(i) == "string")
                {
                     ofs<<"\t\tQLineEdit *"<<pt.PARA_NAME.at(i)<<"_edit;"<<endl;
                     ofs<<"\t\tQString "<<pt.PARA_NAME.at(i)<<";"<<endl;
                }
                else
                {
                    if(pt.PARA_TYPE.at(i) == "int")
                        ofs<<"\t\tQSpinBox *"<<pt.PARA_NAME.at(i)<<"_box;"<<endl;
                    else
                        ofs<<"\t\tQDoubleSpinBox *"<<pt.PARA_NAME.at(i)<<"_box;"<<endl;
                    ofs<<"\t\t"<<pt.PARA_TYPE.at(i)<<" "<<pt.PARA_NAME.at(i)<<";"<<endl;
                }
            }

        }
        ofs<<"\t};"<<endl;
        templatefile.close();

    }
	ofs.close();
}

void create_plugin_all(PluginTemplate & pt)
{
	create_plugin_header(pt);
	create_plugin_cpp(pt);
	//create_func_header(pt);             // May 11, 2012, no longer export func_header
	//create_func_cpp(pt);                // May 11, 2012, no longer export func_cpp
	create_plugin_pro(pt);
}

bool get_next_string(string &val, istream &is)
{
	int c = is.get(); 
	while(c == ' ' || c == '\t') c = is.get();
	if(c == '"') {getline(is, val, '"'); is.ignore();} 
	else if(c == ')' || c == '\n') return false;
	else is >> val;
	return true;
}


#endif
