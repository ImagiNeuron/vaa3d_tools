#include "v3d_message.h"

#include "neuron_dist_plugin.h"
#include "neuron_dist_func.h"
 
Q_EXPORT_PLUGIN2(neuron_dist, NeuronDistPlugin);
 
QStringList NeuronDistPlugin::menulist() const
{
	return QStringList() 
		<<tr("neuron_weighted_dist")
		<<tr("about");
}

QStringList NeuronDistPlugin::funclist() const
{
	return QStringList()
		<<tr("neuron_weighted_distance")
		<<tr("help");
}

void NeuronDistPlugin::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("neuron_weighted_dist"))
	{
		neuron_dist_io(callback,parent);
	}
	else
	{
        v3d_msg(tr("The plugin to calculate wighted distance between a eswc that encodes feature value (consensus neuron) and a reference neuron (gold standard). "
                   "Distance is defined as the weighted distance among all nearest pairs in two neurons. "
            "Developed by Xiaoxiao Liu, 2016-02-04"));
	}
}

bool NeuronDistPlugin::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
	if (func_name == tr("neuron_weighted_distance"))
	{
		return neuron_dist_io(input, output);
	}
	else if (func_name == tr("help"))
	{
		printHelp();
	}
//	else if (func_name == tr("TOOLBOXneuron_dist"))
//	{
//		neuron_dist_toolbox(input, callback);
//		return true;
//	}
}

