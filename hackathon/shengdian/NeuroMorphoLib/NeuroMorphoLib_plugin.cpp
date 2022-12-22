/* NeuroMorphoLib_plugin.cpp
 * functions for processing neuron morphology
 * 2022-8-9 : by SD-Jiang
 */
 
#include "v3d_message.h"
#include <vector>
#include "NeuroMorphoLib_plugin.h"
using namespace std;
Q_EXPORT_PLUGIN2(NeuroMorphoLib, NMorphoPlugin);
 
QStringList NMorphoPlugin::menulist() const
{
	return QStringList() 
		<<tr("None")
		<<tr("about");
}

QStringList NMorphoPlugin::funclist() const
{
	return QStringList()
		<<tr("func1")
		<<tr("help");
}

void NMorphoPlugin::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("None"))
	{
		v3d_msg("To be implemented.");
	}
	else
	{
		v3d_msg(tr("functions for processing neuron morphology. "
			"Developed by SD-Jiang, 2022-8-9"));
	}
}

bool NMorphoPlugin::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
	vector<char*> infiles, inparas, outfiles;
	if(input.size() >= 1) infiles = *((vector<char*> *)input.at(0).p);
	if(input.size() >= 2) inparas = *((vector<char*> *)input.at(1).p);
	if(output.size() >= 1) outfiles = *((vector<char*> *)output.at(0).p);

    if (func_name == tr("lm_feas"))
	{
        return lm_statistic_features(callback,input,output);
	}
	else if (func_name == tr("help"))
	{
		v3d_msg("To be implemented.");
	}
	else return false;

	return true;
}

