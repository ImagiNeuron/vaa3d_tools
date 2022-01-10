/* plugin_creator_func.cpp
 * This plugin is used to produce v3d plugin project from a template file
 * 2012-01-27 : by Hang Xiao
 */

#include <vector>
#include <string>
#include <iostream>

#include <QDir>
#include <v3d_interface.h>
#include "v3d_message.h"
#include "plugin_creator_func.h"
#include "plugin_creator_gui.h"
#include "create_plugin.h"


using namespace std;

#define STRING2VECTSTRING(values, value) \
{ \
	if(value.find('"') != string::npos) \
	{\
		int pos1 = value.find('"'); \
		while(pos1 != string::npos)\
		{\
			int pos2 = value.find('"',pos1+1); \
			string str = value.substr(pos1+1, pos2 - pos1 - 1);\
			values.push_back(str);\
			pos1 = value.find('"', pos2+1);\
		}\
	}\
	else\
	{\
		int pos1 = value.find_first_not_of(' ');\
		while(pos1 != string::npos)\
		{\
			int pos2 = value.find(' ', pos1);\
			string str = value.substr(pos1, pos2-pos1);values.push_back(str);\
			pos1 = value.find_first_not_of(' ', pos2);\
		}\
	}\
} 

const QString title = QObject::tr("Plugin Creator Plugin");


QString getAppPath()
{
    QDir testPluginsDir = QDir(qApp->applicationDirPath());

#if defined(Q_OS_WIN)
    if (testPluginsDir.dirName().toLower() == "debug" || testPluginsDir.dirName().toLower() == "release")
        testPluginsDir.cdUp();
#elif defined(Q_OS_MAC)
    // In a Mac app bundle, plugins directory could be either
    //  a - below the actual executable i.e. v3d.app/Contents/MacOS/plugins/
    //  b - parallel to v3d.app i.e. foo/v3d.app and foo/plugins/
    if (testPluginsDir.dirName() == "MacOS") {
        QDir testUpperPluginsDir = testPluginsDir;
        testUpperPluginsDir.cdUp();
        testUpperPluginsDir.cdUp();
        testUpperPluginsDir.cdUp(); // like foo/plugins next to foo/v3d.app
        if (testUpperPluginsDir.cd("plugins")) testPluginsDir = testUpperPluginsDir;
        testPluginsDir.cdUp();
    }
#endif

    return testPluginsDir.absolutePath();
}


int create_plugin(V3DPluginCallback2 &callback, QWidget *parent)
{
	GuidingDialog dialog(parent);

	if (dialog.exec()!=QDialog::Accepted) return -1;
	dialog.update();

	PluginTemplate pt;
    string plugin_name_new = dialog.plugin_name;
    islower(plugin_name_new[0])? plugin_name_new[0] = toupper(plugin_name_new[0]):plugin_name_new[0] = tolower(plugin_name_new[0]);
    pt.PLUGIN_NAME = "NeuronAssembler_" + plugin_name_new;
    pt.PLUGIN_CLASS = "NeuronAssembler_" + plugin_name_new;
    pt.WINDOW_TITLE = "NeuronAssembler_" + plugin_name_new;
	pt.PLUGIN_DESCRIPTION = dialog.plugin_desp;
	pt.PLUGIN_DATE = dialog.plugin_date;
	pt.PLUGIN_AUTHOR = dialog.plugin_author;
    pt.VAA3D_PATH = getAppPath().toStdString() + "/../v3d_main/";
    pt.OUTPUTSWC = dialog.outputswc_name;
    pt.FUNC_NAME = dialog.function_name;
    pt.TRACINGPLUGIN_NAME = dialog.plugin_name;
    STRING2VECTSTRING(pt.PARA_NAME, dialog.paranamelist);
    STRING2VECTSTRING(pt.PARA_TYPE, dialog.paratypelist);
    STRING2VECTSTRING(pt.PARA_VALUE, dialog.paravaluelist);

	pt.DOFUNC = true;

//    cout<<"menus.size() = "<<pt.MENUS.size()<<endl;
//    cout<<"funcs.size() = "<<pt.FUNCS.size()<<endl;

    string save_folder = dialog.save_folder;
	if(save_folder[0] == '~') {
		save_folder.erase(0,1);
		save_folder = QDir::homePath().toStdString() + save_folder;
	}
    else if (save_folder[1] == ':') //for windows. added by PHC, 20130905
    {
        //do nothing
    }
    else if (save_folder[0] != '/' && save_folder[0] != '.')
        save_folder = "./" + save_folder;

	pt.PLUGIN_HEADER = pt.PLUGIN_NAME + "_plugin.h";
	pt.PLUGIN_CPP =  pt.PLUGIN_NAME + "_plugin.cpp";
	pt.FUNC_HEADER = pt.PLUGIN_NAME + "_func.h";
	pt.FUNC_CPP =  pt.PLUGIN_NAME + "_func.cpp";
	pt.PRO_FILE =  pt.PLUGIN_NAME + ".pro";

	QDir dir(save_folder.c_str()); 
	if(!dir.exists()){QMessageBox::warning(0, QObject::tr("Error"), QObject::tr("Un existing foler : %1").arg(save_folder.c_str())); return 0;}
	else 
    {
        v3d_msg(QString("Files:\n \t%1\n \t%2\n \t%3\n have been saved to directory: [%4]. \n\nYou can go to that folder now and run the following command to build the plugin (assuming you have Qt installed and gcc/make on your computer Mac/Linux/Windows): \n\n>qmake\n>make (or change to nmake -f Makefile.Release for instance for Windows)\n").arg(pt.PLUGIN_HEADER.c_str()).arg(pt.PLUGIN_CPP.c_str()).arg(pt.PRO_FILE.c_str()).arg(save_folder.c_str()));
	}
	QString cur_path = QDir::current().dirName();
	cout<<"current path : "<<QDir::current().dirName().toStdString()<<endl;
	QDir::setCurrent(save_folder.c_str());
	cout<<"current path : "<<QDir::current().dirName().toStdString()<<endl;
	create_plugin_all(pt);
	QDir::setCurrent(cur_path);
	cout<<"current path : "<<QDir::current().dirName().toStdString()<<endl;

	return 1;
}

