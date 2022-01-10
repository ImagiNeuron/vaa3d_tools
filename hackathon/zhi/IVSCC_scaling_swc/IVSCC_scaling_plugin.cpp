/* IVSCC_scaling_plugin.cpp
 * This is a test plugin, you can use it as a demo.
 * 2014-10-30 : by Zhi Zhou
 */
 
#include "v3d_message.h"
#include <vector>
#include "IVSCC_scaling_plugin.h"
#include <iostream>
#include "../IVSCC_sort_swc/openSWCDialog.h"
#include "ui_SWC_scaling.h"

using namespace std;
Q_EXPORT_PLUGIN2(IVSCC_scaling, IVSCC_scaling);
 
QStringList IVSCC_scaling::menulist() const
{
	return QStringList()
		<<("scale swc")
		<<tr("about");
}

QStringList IVSCC_scaling::funclist() const
{
	return QStringList()
        <<tr("scale")
        <<tr("fmost_scale")
		<<tr("help");
}

void IVSCC_scaling::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("scale swc")) // As per annotation team's request, an UI is added. MK, Oct, 2017.
	{
		swcScalingUI* inputForm = new swcScalingUI(0, &callback);
		inputForm->exec();
		inputForm->~swcScalingUI();
	}
    else if (menu_name == tr("about"))
	{
        v3d_msg(tr("This is a plugin for scaling swc files in IVSCC pipeline. \n"
            "Developed by Zhi Zhou, 2014-10-30\n"
			"Revised by MK, 2017-10-10"));
	}
}

bool IVSCC_scaling::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
	vector<char*> infiles, inparas, outfiles;
	if(input.size() >= 1) infiles = *((vector<char*> *)input.at(0).p);
	if(input.size() >= 2) inparas = *((vector<char*> *)input.at(1).p);
	if(output.size() >= 1) outfiles = *((vector<char*> *)output.at(0).p);

    if (func_name == tr("scale"))
	{
          cout<<"Welcome to IVSCC swc scaling processing plugin"<<endl;
          if(infiles.empty())
          {
              cerr<<"Need input swc file"<<endl;
              return false;
          }

          QString  inswc_file =  infiles[0];
          QString  outswc_file =  outfiles[0];

          double x_scale, y_scale, z_scale, r_scale;
          if (input.size()>=2)
          {

              vector<char*> paras = (*(vector<char*> *)(input.at(1).p));
              if(paras.size() >= 4)
              {
                  cout<<paras.size()<<endl;
                  x_scale = atof(paras.at(0));
                  y_scale = atof(paras.at(1));
                  z_scale = atof(paras.at(2));
                  r_scale = atof(paras.at(3));
              }
              else
              {
                  cerr<<"Need four scaling factors"<<endl;
                  return false;
              }
          }
          else
          {
              cerr<<"Need four scaling factors"<<endl;
              return false;
          }

          cout<<"inswc_file = "<<inswc_file.toStdString().c_str()<<endl;
          cout<<"outswc_file = "<<outswc_file.toStdString().c_str()<<endl;

          NeuronTree nt = readSWC_file(inswc_file);

          NeuronTree nt_scaled;
          QList <NeuronSWC> listNeuron;
          QHash <int, int>  hashNeuron;
          listNeuron.clear();
          hashNeuron.clear();

          NeuronSWC S;
          QList<NeuronSWC> list = nt.listNeuron;

          for(V3DLONG i = 0; i < list.size();i++)
          {
              S.n =  nt.listNeuron.at(i).n;
              S.x = x_scale *  nt.listNeuron.at(i).x;
              S.y = y_scale *  nt.listNeuron.at(i).y;
              S.z = z_scale *  nt.listNeuron.at(i).z;
              S.r = r_scale *  nt.listNeuron.at(i).r;
              S.pn =  nt.listNeuron.at(i).pn;
              S.type =  nt.listNeuron.at(i).type;
              listNeuron.append(S);
              hashNeuron.insert(S.n, listNeuron.size()-1);
          }

          nt_scaled.n = -1;
          nt_scaled.on = true;
          nt_scaled.listNeuron = listNeuron;
          nt_scaled.hashNeuron = hashNeuron;

          writeSWC_file(outswc_file,nt_scaled);
    }else if (func_name == tr("fmost_scale"))
    {
          cout<<"Welcome to fMOST swc scaling processing plugin"<<endl;
          if(infiles.empty())
          {
              cerr<<"Need input swc file"<<endl;
              return false;
          }

          QString  inswc_file =  infiles[0];
          QString  outswc_file =  outfiles[0];

          double x_scale, y_scale, z_scale, r_scale;
          if (input.size()>=2)
          {

              vector<char*> paras = (*(vector<char*> *)(input.at(1).p));
              if(paras.size() >= 4)
              {
                  cout<<paras.size()<<endl;
                  x_scale = atof(paras.at(0));
                  y_scale = atof(paras.at(1));
                  z_scale = atof(paras.at(2));
                  r_scale = atof(paras.at(3));
              }
              else
              {
                  cerr<<"Need four scaling factors"<<endl;
                  return false;
              }
          }
          else
          {
              cerr<<"Need four scaling factors"<<endl;
              return false;
          }

          cout<<"inswc_file = "<<inswc_file.toStdString().c_str()<<endl;
          cout<<"outswc_file = "<<outswc_file.toStdString().c_str()<<endl;

          NeuronTree nt = readSWC_file(inswc_file);
          QList<NeuronSWC> list = nt.listNeuron;
          double max_y = -100000000, max_z = -100000000;
          double min_y = 100000000,min_z = 100000000;

          for(V3DLONG i = 0; i < list.size();i++)
          {
              nt.listNeuron[i].x = x_scale *  nt.listNeuron.at(i).x;
              nt.listNeuron[i].y = y_scale *  nt.listNeuron.at(i).y;
              nt.listNeuron[i].z = z_scale *  nt.listNeuron.at(i).z;
              nt.listNeuron[i].r = r_scale *  nt.listNeuron.at(i).r;

              if(nt.listNeuron[i].y >= max_y) max_y = nt.listNeuron[i].y;
              if(nt.listNeuron[i].y <= min_y) min_y = nt.listNeuron[i].y;

              if(nt.listNeuron[i].z >= max_z) max_z = nt.listNeuron[i].z;
              if(nt.listNeuron[i].z <= min_z) min_z = nt.listNeuron[i].z;


          }
          for(V3DLONG i = 0; i < nt.listNeuron.size();i++)
          {
              nt.listNeuron[i].y = max_y + min_y - nt.listNeuron[i].y;
              nt.listNeuron[i].z = max_z + min_z - nt.listNeuron[i].z;

          }

          writeSWC_file(outswc_file,nt);
    }
	else if (func_name == tr("help"))
	{
        cout<<"Usage : v3d -x dllname -f scale -i <inswc_folder> -o <outswc_file> -p <x_factor> <y_factor> <z_factor> <radius_factor>"<<endl;
        cout<<endl;
	}
	else return false;

	return true;
}

