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
#include "common_dialog.h"
#include <QMessageBox>
using namespace std;


QString getVaa3dPath()
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
    testPluginsDir.cdUp();
    return testPluginsDir.absolutePath();
}

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

int create_plugin(V3DPluginCallback2 &callback, QWidget *parent)
{
	GuidingDialog dialog(parent);
	if (dialog.exec()!=QDialog::Accepted) return -1;
    dialog.update();

	PluginTemplate pt;
	pt.PLUGIN_NAME = dialog.plugin_name;
	pt.PLUGIN_CLASS = dialog.plugin_class;
	pt.WINDOW_TITLE = dialog.win_title;
	pt.PLUGIN_DESCRIPTION = dialog.plugin_desp;
	pt.PLUGIN_DATE = dialog.plugin_date;
	pt.PLUGIN_AUTHOR = dialog.plugin_author;
	pt.VAA3D_PATH = dialog.vaa3d_path;
	STRING2VECTSTRING(pt.MENUS, dialog.menulst);
	STRING2VECTSTRING(pt.FUNCS, dialog.funclst);
	if(dialog.funclst.find("help") == string::npos) pt.FUNCS.push_back("help");
	pt.DOFUNC = true;

	cout<<"menus.size() = "<<pt.MENUS.size()<<endl;
	cout<<"funcs.size() = "<<pt.FUNCS.size()<<endl;

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


int create_demo1(V3DPluginCallback2 &callback, QWidget *parent)
{
	vector<string> items;
	items.push_back("VAA3D Path");
	items.push_back("Save Path");
	CommonDialog dialog(items);
	dialog.setWindowTitle(QObject::tr("Create Demo1 Project"));
	if(dialog.exec() != QDialog::Accepted) return 1;
	string vaa3d_path = dialog.get_para("VAA3D Path");
	string save_path = dialog.get_para("Save Path");
	produce_demo1(save_path, vaa3d_path);
	QMessageBox::information(0, "Finished", QObject::tr("<pre>The demo project have been saved to %1\n\n"
													    "> cd %1\n"
														"> qmake\n"
														"> make\n</pre>")
			                                         .arg(save_path.c_str()));
}

int create_demo2(V3DPluginCallback2 &callback, QWidget *parent)
{
	vector<string> items;
	items.push_back("Plugin Name");
	items.push_back("Menu Name");
	items.push_back("Func Name");
	items.push_back("VAA3D Path");
	items.push_back("Save Path");
	CommonDialog dialog(items);
	dialog.setWindowTitle(QObject::tr("Create Demo2 Project"));
	dialog.setHelp(QObject::tr("Only one menu name and one func name will be accepted. There should be no empty space in Plugin Name, and Func Name. Menu Name with space will be treated as one menu."));
	dialog.setHistory(true);
	if(dialog.exec() != QDialog::Accepted) return 1;

	string save_path = dialog.get_para("Save Path");
	string vaa3d_path = dialog.get_para("VAA3D Path");
	string plugin_name = dialog.get_para("Plugin Name");
	string menu_name = dialog.get_para("Menu Name");
	string func_name = dialog.get_para("Func Name");

	produce_demo2(save_path, vaa3d_path, plugin_name, menu_name, func_name);
	QMessageBox::information(0, "Finished", QObject::tr("<pre>The demo project have been saved to %1\n\n"
													    "> cd %1\n"
														"> qmake\n"
														"> make\n</pre>")
			                                         .arg(save_path.c_str()));
}

int create_plugin_neuronrec(V3DPluginCallback2 &callback, QWidget *parent)
{
    GuidingDialog dialog(parent);
   // dialog.editor_vaa3d_path->setText(getVaa3dPath());
    dialog.editor_menu_list->setText(QString("tracing_menu"));
    dialog.editor_func_list->setText(QString("tracing_func"));

    if (dialog.exec()!=QDialog::Accepted) return -1;
    dialog.update();

    QFile v3d_interface_file(QString("%1/v3d_main/basic_c_fun/v3d_interface.h").arg(dialog.vaa3d_path.c_str()));
    if(!v3d_interface_file.exists())
    {
        v3d_msg("Vaa3D whole-project path is not correct!");
        return -1;
    }
    PluginTemplate pt;
    pt.PLUGIN_NAME = dialog.plugin_name;
    pt.PLUGIN_CLASS = dialog.plugin_class;
    pt.WINDOW_TITLE = dialog.win_title;
    pt.PLUGIN_DESCRIPTION = dialog.plugin_desp;
    pt.PLUGIN_DATE = dialog.plugin_date;
    pt.PLUGIN_AUTHOR = dialog.plugin_author;
    pt.VAA3D_PATH = dialog.vaa3d_path;
    STRING2VECTSTRING(pt.MENUS, dialog.menulst);
    STRING2VECTSTRING(pt.FUNCS, dialog.funclst);
    if(dialog.funclst.find("help") == string::npos) pt.FUNCS.push_back("help");
    pt.DOFUNC = true;

    cout<<"menus.size() = "<<pt.MENUS.size()<<endl;
    cout<<"funcs.size() = "<<pt.FUNCS.size()<<endl;

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
    create_plugin_neuronrec_all(pt);
    QDir::setCurrent(cur_path);
    cout<<"current path : "<<QDir::current().dirName().toStdString()<<endl;

    return 1;
}

int create_plugin_neuronML(V3DPluginCallback2 &callback, QWidget *parent)
{
    GuidingDialog dialog(parent);
   // dialog.editor_vaa3d_path->setText(getVaa3dPath());
    dialog.editor_menu_list->setText(QString("tracing_menu"));
    dialog.editor_func_list->setText(QString("tracing_func"));

    if (dialog.exec()!=QDialog::Accepted) return -1;
    dialog.update();

    QFile v3d_interface_file(QString("%1/v3d_main/basic_c_fun/v3d_interface.h").arg(dialog.vaa3d_path.c_str()));
    if(!v3d_interface_file.exists())
    {
        v3d_msg("Vaa3D whole-project path is not correct!");
        return -1;
    }
    PluginTemplate pt;
    pt.PLUGIN_NAME = dialog.plugin_name;
    pt.PLUGIN_CLASS = dialog.plugin_class;
    pt.WINDOW_TITLE = dialog.win_title;
    pt.PLUGIN_DESCRIPTION = dialog.plugin_desp;
    pt.PLUGIN_DATE = dialog.plugin_date;
    pt.PLUGIN_AUTHOR = dialog.plugin_author;
    pt.VAA3D_PATH = dialog.vaa3d_path;
    STRING2VECTSTRING(pt.MENUS, dialog.menulst);
    STRING2VECTSTRING(pt.FUNCS, dialog.funclst);
    if(dialog.funclst.find("help") == string::npos) pt.FUNCS.push_back("help");
    pt.DOFUNC = true;

    cout<<"menus.size() = "<<pt.MENUS.size()<<endl;
    cout<<"funcs.size() = "<<pt.FUNCS.size()<<endl;

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
    create_plugin_neuronML_all(pt);
    QDir::setCurrent(cur_path);
    cout<<"current path : "<<QDir::current().dirName().toStdString()<<endl;

    return 1;
}
