/*
 *  ada_threshold.cpp
 *
 *  Created by Yang, Jinzhu and Hanchuan Peng, on 11/22/10.
 *  Add dofunc() interface by Jianlong Zhou, 2012-04-18.
 */

#include <iostream>
#include "ada_threshold.h"
#include "v3d_message.h"

#include "stackutil.h"

using namespace std;

//Q_EXPORT_PLUGIN2 ( PluginName, ClassName )
//The value of PluginName should correspond to the TARGET specified in the plugin's project file.
Q_EXPORT_PLUGIN2(threshold, ThPlugin);

template <class T>
void BinaryProcess(T *apsInput, T * aspOutput, V3DLONG iImageWidth, V3DLONG iImageHeight, V3DLONG iImageLayer, V3DLONG h, V3DLONG d)
{
	V3DLONG i, j,k,n,count;
	double t, temp;

	V3DLONG mCount = iImageHeight * iImageWidth;
	for (i=0; i<iImageLayer; i++)
	{
		for (j=0; j<iImageHeight; j++)
		{
			for (k=0; k<iImageWidth; k++)
			{
				V3DLONG curpos = i * mCount + j*iImageWidth + k;
				V3DLONG curpos1 = i* mCount + j*iImageWidth;
				V3DLONG curpos2 = j* iImageWidth + k;
				temp = 0;
				count = 0;
				for(n =1 ; n <= d  ;n++)
				{
					if (k>h*n) {temp += apsInput[curpos1 + k-(h*n)]; count++;}
					if (k+(h*n)< iImageWidth) { temp += apsInput[curpos1 + k+(h*n)]; count++;}
                    if (j>h*n) {temp += apsInput[i* mCount + (j-(h*n))*iImageWidth + k]; count++;}//
					if (j+(h*n)<iImageHeight) {temp += apsInput[i* mCount + (j+(h*n))*iImageWidth + k]; count++;}//
					if (i>(h*n)) {temp += apsInput[(i-(h*n))* mCount + curpos2]; count++;}//
					if (i+(h*n)< iImageLayer) {temp += apsInput[(i+(h*n))* mCount + j* iImageWidth + k ]; count++;}
				}
				t =  apsInput[curpos]-temp/(count);
				aspOutput[curpos]= (t > 0)? t : 0;
			}
		}
	}
}

void thimg(V3DPluginCallback2 &callback, QWidget *parent, int method_code);
bool thimg(V3DPluginCallback2 &callback, const V3DPluginArgList & input, V3DPluginArgList & output);

//plugin funcs
const QString title = "adaptive threshold transform";
QStringList ThPlugin::menulist() const
{
     return QStringList()
          << tr("3D (w/o parameters)")
          << tr("3D (set parameters)")
          << tr("Help");
}

void ThPlugin::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("3D (w/o parameters)"))
	{
          thimg(callback, parent,1 );
     }
	else if (menu_name == tr("3D (set parameters)"))
	{
		thimg(callback, parent, 2 );
	}
	else if (menu_name == tr("Help"))
	{
		v3d_msg("Simple adaptive thresholding: for each voxel, compute a threshold which is the average intensity of its neighboring voxels and then subtract the threshold from the current voxel's intensity value. If the result is <0, then set it as 0. The neighborhood is defined along 6 axial directions in 3D, with N samples of each direction (N -- the 'number of sampling points' in the parameter setting dialog), and M voxels between every nearest pair of samples (M -- the 'sampling interval' in the parameter setting dialog).");
		return;
	}

}

QStringList ThPlugin::funclist() const
{
	return QStringList()
		<<tr("adath")
       << tr("simple_thresholding")
		<<tr("help");
}


bool ThPlugin::dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent)
{
     if (func_name == tr("adath"))
	{
        return thimg(callback, input, output);
	}
     else if(func_name == tr("simple_thresholding"))
     {
         vector<char*> infiles, inparas, outfiles;
         if(input.size() >= 1) infiles = *((vector<char*> *)input.at(0).p);
         if(input.size() >= 2) inparas = *((vector<char*> *)input.at(1).p);
         if(output.size() >= 1) outfiles = *((vector<char*> *)output.at(0).p);

         cout<<"Welcome to simple thresholding"<<endl;
          string inimg_file;
          if(infiles.size()>=1) {
              inimg_file = infiles[0];
              int thre_thre=(inparas.size()>=1)?atoi(inparas[0]):20;
              int toMinValue=(inparas.size()>=2)?atoi(inparas[1]):0;
              int toMaxValue=(inparas.size()>=3)?atoi(inparas[2]):0;

              //get out_image_path and out_swc_path
              string out_image_file=(outfiles.size()>=1)?outfiles[0]:(inimg_file+"_out.v3draw");

              unsigned char * data1d = 0;
              V3DLONG in_sz[4];
              int datatype;
              if(!simple_loadimage_wrapper(callback, (char*)inimg_file.c_str(), data1d, in_sz, datatype))
              {
                  cerr<<"load image "<<inimg_file<<" error!"<<endl;
                  return false;
              }

              V3DLONG iImageWidth = in_sz[0];
              V3DLONG iImageHeight = in_sz[1];
              V3DLONG iImageLayer = in_sz[2];
              V3DLONG pagesz = iImageWidth*iImageHeight*iImageLayer;
              in_sz[3] = 1;

              unsigned char * outImg = 0;
              try {outImg = new unsigned char [pagesz];}
              catch(...)  {v3d_msg("cannot allocate memory for data_blended."); return false;}

              V3DLONG i, j,k;
              double t;

              V3DLONG mCount = iImageHeight * iImageWidth;
              for (i=0; i<iImageLayer; i++)
              {
                  for (j=0; j<iImageHeight; j++)
                  {
                      for (k=0; k<iImageWidth; k++)
                      {
                          V3DLONG curpos = i * mCount + j*iImageWidth + k;
                          t =  data1d[curpos];
                          if(toMaxValue>thre_thre)
                              outImg[curpos]= (t > thre_thre)? toMaxValue:t;
                          else
                              outImg[curpos]= (t < thre_thre)? toMinValue:t;
                      }
                  }
              }

              simple_saveimage_wrapper(callback, (char*)out_image_file.c_str(),(unsigned char *)outImg, in_sz, 1);
              if(data1d) {delete []data1d; data1d =0;}
              if(outImg) {delete []outImg; outImg =0;}
              return true;
          }
          return true;
     }
	else if(func_name == tr("help"))
	{
		cout<<"Usage : v3d -x threshold -f adath -i <inimg_file> -o <outimg_file> -p <h> <d>"<<endl;
		cout<<endl;
		cout<<"h       sampling interval, default 5,"<<endl;
		cout<<"d       number of sampling points, default 3,"<<endl;
		cout<<endl;
		cout<<"e.g. v3d -x threshold -f adath -i input.raw -o output.raw -p 5 3"<<endl;
		cout<<endl;
        cout<<"Usage : v3d -x  <dll>  -f  simple_thresholding -i <inimg_file> -o <outimg_file> -p <bkg_thre> <min_value> <max_value>"<<endl;
        cout<<endl;
        cout<<"bkg_thre: background threshold (default 20),"<<endl;
        cout<<"min_value: pixel value below bkg_thre will set to min_value (default 0),"<<endl;
        cout<<"max_value: pixel value higher than bkg_thre will set to max_value (default 0)"<<endl;
        cout<<"Note: only when max_value is higher than bkg_thre, pixel intensity will set to max_value"<<endl;
        cout<<endl;
        cout<<"e.g. v3d -x threshold -f simple_thresholding -i input.raw -o output.raw -p 15 0 0"<<endl;
        cout<<endl;
		return true;
	}
}

bool thimg(V3DPluginCallback2 &callback, const V3DPluginArgList & input, V3DPluginArgList & output)
{
	cout<<"Welcome to Gaussian filter"<<endl;
	if (input.size()<1 || output.size() != 1) return false;

	V3DLONG h = 5, d = 3;
     if (input.size()>=2)
     {
          vector<char*> paras = (*(vector<char*> *)(input.at(1).p));
          if(paras.size() >= 1) h = atoi(paras.at(0));
          if(paras.size() >= 2) d = atoi(paras.at(1));
	}

	char * inimg_file = ((vector<char*> *)(input.at(0).p))->at(0);
	char * outimg_file = ((vector<char*> *)(output.at(0).p))->at(0);
	cout<<"h = "<<h<<endl;
     cout<<"d = "<<d<<endl;
	cout<<"inimg_file = "<<inimg_file<<endl;
	cout<<"outimg_file = "<<outimg_file<<endl;

    Image4DSimple *subject = callback.loadImage(inimg_file);
    if(!subject || !subject->valid())
    {
         v3d_msg("Fail to load the input image.");
         if (subject) {delete subject; subject=0;}
         return false;
    }

     clock_t start_t = clock(); // record time point
     // =====================================================================>>>>>>>>>
     V3DLONG sz0 = subject->getXDim();
     V3DLONG sz1 = subject->getYDim();
     V3DLONG sz2 = subject->getZDim();
     V3DLONG sz3 = subject->getCDim();
    V3DLONG pagesz_sub = sz0*sz1*sz2;

	//----------------------------------------------------------------------------------------------------------------------------------
	V3DLONG channelsz = sz0*sz1*sz2;
	void *pData=NULL;

	V3DLONG sz_data[4]; sz_data[0]=sz0; sz_data[1]=sz1; sz_data[2]=sz2; sz_data[3]=1;
        switch (subject->getDatatype())
		{
            case V3D_UINT8:
				try
				{
					pData = (void *)(new unsigned char [sz3*channelsz]);
				}
					catch (...)
				{
					v3d_msg("Fail to allocate memory in Distance Transform.",0);
					if (pData) {delete []pData; pData=0;}
					return false;
				}

				{
                    unsigned char * pSubtmp_uint8 = subject->getRawData();

					for (V3DLONG ich=0; ich<sz3; ich++)
						BinaryProcess(pSubtmp_uint8+ich*channelsz, (unsigned char *)pData+ich*channelsz, sz0, sz1, sz2, h, d  );
				}
				break;

            case V3D_UINT16:
				try
				{
					pData = (void *)(new short int [sz3*channelsz]);
				}
					catch (...)
				{
					v3d_msg("Fail to allocate memory in Distance Transform.",0);
					if (pData) {delete []pData; pData=0;}
					return false;
				}

				{
                    short int * pSubtmp_uint16 = (short int *)subject->getRawData();

					for (V3DLONG ich=0; ich<sz3; ich++)
						BinaryProcess(pSubtmp_uint16+ich*channelsz, (short int *)pData+ich*channelsz, sz0, sz1, sz2, h, d );
				}

				break;

            case V3D_FLOAT32:
				try
				{
					pData = (void *)(new float [sz3*channelsz]);
				}
					catch (...)
				{
					v3d_msg("Fail to allocate memory in Distance Transform.", 0);
					if (pData) {delete []pData; pData=0;}
					return false;
				}

				{
                    float * pSubtmp_float32 = (float *)subject->getRawData();

					for (V3DLONG ich=0; ich<sz3; ich++)
						BinaryProcess(pSubtmp_float32+ich*channelsz, (float *)pData+ich*channelsz, sz0, sz1, sz2, h, d );
				}

				break;

			default:
				break;
		}

	//----------------------------------------------------------------------------------------------------------------------------------

	clock_t end_t = clock();
	printf("time eclapse %d s for dist computing!\n", (end_t-start_t)/1000000);

     // =====================================================================<<<<<<<<<
    Image4DSimple outimg;
    outimg.setData((unsigned char *)pData, sz0, sz1, sz2, sz3, subject->getDatatype());

    callback.saveImage(&outimg, outimg_file);
    v3d_msg("Finish saving output file.",0);

     if(subject) {delete []subject; subject=0;}

     return true;
}





void thimg(V3DPluginCallback2 &callback, QWidget *parent, int method_code)
{
	v3dhandle curwin = callback.currentImageWindow();
	V3DLONG h;
	V3DLONG d;
	if (!curwin)
	{
		v3d_msg("You don't have any image open in the main window.");
		return;
	}

	if (method_code == 1)
	{
		h = 5;
		d = 3;
	}
	else
	{
		if( method_code == 2)
		{
			AdaTDialog dialog(callback, parent);
			if (dialog.exec()!=QDialog::Accepted)
			return;
			else
			{
				h = dialog.Ddistance->text().toLong()-1;
				d = dialog.Dnumber->text().toLong()-1;
				printf("d% h,d% d \n ",h,d);
			}
		}
	}

	clock_t start_t = clock(); // record time point

	Image4DSimple* subject = callback.getImage(curwin);
	QString m_InputFileName = callback.getImageName(curwin);

	if (!subject)
	{
		QMessageBox::information(0, title, QObject::tr("No image is open."));
		return;
	}
	Image4DProxy<Image4DSimple> pSub(subject);

	V3DLONG sz0 = subject->getXDim();
    V3DLONG sz1 = subject->getYDim();
    V3DLONG sz2 = subject->getZDim();
	V3DLONG sz3 = subject->getCDim();
	V3DLONG pagesz_sub = sz0*sz1*sz2;

	//----------------------------------------------------------------------------------------------------------------------------------
	V3DLONG channelsz = sz0*sz1*sz2;
	void *pData=NULL;

	V3DLONG sz_data[4]; sz_data[0]=sz0; sz_data[1]=sz1; sz_data[2]=sz2; sz_data[3]=1;
		switch (subject->getDatatype())
		{
			case V3D_UINT8:
				try
				{
					pData = (void *)(new unsigned char [sz3*channelsz]);
				}
					catch (...)
				{
					v3d_msg("Fail to allocate memory in Distance Transform.");
					if (pData) {delete []pData; pData=0;}
					return;
				}

				{
					unsigned char * pSubtmp_uint8 = pSub.begin();

					for (V3DLONG ich=0; ich<sz3; ich++)
						BinaryProcess(pSubtmp_uint8+ich*channelsz, (unsigned char *)pData+ich*channelsz, sz0, sz1, sz2, h, d  );
				}
				break;

			case V3D_UINT16:
				try
				{
					pData = (void *)(new short int [sz3*channelsz]);
				}
					catch (...)
				{
					v3d_msg("Fail to allocate memory in Distance Transform.");
					if (pData) {delete []pData; pData=0;}
					return;
				}

				{
					short int * pSubtmp_uint16 = (short int *)pSub.begin();

					for (V3DLONG ich=0; ich<sz3; ich++)
						BinaryProcess(pSubtmp_uint16+ich*channelsz, (short int *)pData+ich*channelsz, sz0, sz1, sz2, h, d );
				}

				break;

			case V3D_FLOAT32:
				try
				{
					pData = (void *)(new float [sz3*channelsz]);
				}
					catch (...)
				{
					v3d_msg("Fail to allocate memory in Distance Transform.");
					if (pData) {delete []pData; pData=0;}
					return;
				}

				{
					float * pSubtmp_float32 = (float *)pSub.begin();

					for (V3DLONG ich=0; ich<sz3; ich++)
						BinaryProcess(pSubtmp_float32+ich*channelsz, (float *)pData+ich*channelsz, sz0, sz1, sz2, h, d );
				}

				break;

			default:
				break;
		}

	//----------------------------------------------------------------------------------------------------------------------------------

	clock_t end_t = clock();
	printf("time eclapse %d s for dist computing!\n", (end_t-start_t)/1000000);

	Image4DSimple p4DImage;
	p4DImage.setData((unsigned char*)pData, sz0, sz1, sz2, sz3, subject->getDatatype());

	v3dhandle newwin;
	if(QMessageBox::Yes == QMessageBox::question (0, "", QString("Do you want to use the existing window?"), QMessageBox::Yes, QMessageBox::No))
		newwin = callback.currentImageWindow();
	else
		newwin = callback.newImageWindow();

	callback.setImage(newwin, &p4DImage);
	callback.setImageName(newwin, QString("thresholded image"));
	callback.updateImageWindow(newwin);
}

void AdaTDialog::update()
{
	//get current data
	Dn = Dnumber->text().toLong()-1;
	Dh = Ddistance->text().toLong()-1;
		//printf("channel %ld val %d x %ld y %ld z %ld ind %ld \n", c, data1d[ind], nx, ny, nz, ind);
}
