/* movieZCswitch.cpp
 * 2009-09-22: create this program by Yang Yu
 * 2009-09-28. Last edit by Hanchuan Peng. only change the texts in the options
 * 2010-08-06. by hanchuan Peng, to adapt to the patch of the image4dsimple data structure
 * 2011-01-27. by Yang Yu, to support multiple datatype V3D_UINT8, V3D_UINT16, V3D_FLOAT32
 */
 
#include <QtGui>

#include <string>
#include <exception>
#include <iostream>
#include <algorithm>
#include <math.h>


#include "movieZCswitch.h"

using namespace std;

QStringList importSeriesFileList_addnumbersort(const QString & curFile)
{
    QStringList myList;
    myList.clear();

    // get the image files namelist in the directory
    QStringList imgSuffix;
    imgSuffix<<  QString("*.%1").arg(QFileInfo(curFile).completeSuffix().toStdString().c_str());

    QString curFilePath = QFileInfo(curFile).absolutePath();
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
    myList.sort();

    QStringList sortedList, tmpList;

    //-----------------------------------------------------------------------
    // 090731 RZC: fixed numerically sorting file names list, for XFormWidget::importGeneralImgSeries
    //-----------------------------------------------------------------------
    QString fileNameStr, fileNameDigits;	//absolute file name is separated to 2 parts: strings and digits
    QRegularExpression r("(\\d+)");	//find digits
    QMap<V3DLONG, QString> mapList;

    mapList.clear();
    for(V3DLONG i=0; i<myList.size(); ++i)
    {
        fileNameStr = myList.at(i);
        QFileInfo fileFullName(myList.at(i));
        QString fileFullNameStr = fileFullName.completeBaseName();


        //extract the fileNameDigits from fileNameStr
        //e.g. "a9_b2009051801.tif.raw" into "a9_b.tif.raw" & "2009051801"

        V3DLONG pos = 0;
        fileNameDigits = "";
        QRegularExpressionMatch match = r.match(fileFullNameStr);
        if(match.hasMatch()){
            fileNameDigits = match.captured(1);

        }

//        while ((pos = r.indexIn(fileFullNameStr, pos)) != -1)
//        {
//                    fileNameDigits = r.cap(1);
//                    pos += r.matchedLength();
//        }

        if (fileNameDigits.isEmpty()) continue;


        V3DLONG num = fileNameDigits.toULong();
        mapList.insert(num, fileNameStr);
    }
    // must be sorted by QMap
    myList = mapList.values();
    //foreach (QString qs, myList)  qDebug() << qs;

    //no need to pop-out a dialog if no meaningful file has been detected. 131017
    if (myList.isEmpty())
    {
        v3d_msg("It seems no file contains a digit-portion in the file name. Naming convention should be sth like xxx_000001.yyy, xxx_000002.yyy, .... Check your data before importing.");
        return myList;
    }

    // print filenames
    foreach (QString qs, myList)  qDebug() << qs;

    return myList;
}

//Q_EXPORT_PLUGIN2 ( PluginName, ClassName )
//The value of PluginName should correspond to the TARGET specified in the plugin's project file.
//Q_EXPORT_PLUGIN2(movieZCswitch, MovieZCswitchPlugin)

int changeMS(V3DPluginCallback2 &callback, QWidget *parent);
int packInChannel(V3DPluginCallback2 &callback, QWidget *parent);
int import5D(V3DPluginCallback2 &callback, QWidget *parent);


bool changeMS(V3DPluginCallback2 &callback, const V3DPluginArgList & input, V3DPluginArgList & output);



const QString title = "5D Stack Converter";
QStringList MovieZCswitchPlugin::menulist() const
{
    return QStringList() << tr("5D Stack Converter")
                         << tr("pack_in_selected_channels")
                         << tr("import_5D_image")
						 << tr("about this plugin");
}

QStringList MovieZCswitchPlugin::funclist() const
{
    return QStringList()
    <<tr("4D_to_5D")
    <<tr("help");
}

void MovieZCswitchPlugin::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("5D Stack Converter"))
    {
    	changeMS(callback, parent);
    }else if (menu_name == tr("pack_in_channels"))
    {
        packInChannel(callback, parent);
    }
    else if (menu_name == tr("import_5D_image"))
    {
        import5D(callback, parent);
    }
	else if (menu_name == tr("about this plugin"))
	{
		QMessageBox::information(parent, "Version info", 
            QString("5D Stack Converter Plugin Demo %1 (2009-Sep-22) developed by Yang Yu. (Hanchuan Peng Lab, Janelia Research Farm Campus, HHMI)"
            ).arg(getPluginVersion()));
	}
	
}

bool MovieZCswitchPlugin::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback, QWidget * parent)
{
    if (func_name==tr("4D_to_5D"))
    {
        return changeMS(callback, input, output);
    }

    else if (func_name==tr("help"))
    {
        cout<<"Usage : v3d -x 5D_Stack_Converter -f 4D_to_5D -i <inimg_file> -o <outfolder_file> -p <tp>"<<endl;
        cout<<endl;
        cout<<"outfolder_file          The output folder to save all timepoint image stacks"<<endl;
        cout<<"tp                      timepoints, default 1"<<endl;
        cout<<endl;
        return true;
    }
    return false;
}


int changeMS(V3DPluginCallback2 &callback, QWidget *parent)
{
	v3dhandle oldwin = callback.currentImageWindow();
	Image4DSimple* image = callback.getImage(oldwin);
	if (! image)
	{
		QMessageBox::information(0, title, QObject::tr("No image is open."));
		return -1;
	}
	
	ImagePixelType imgdatatype = image->getDatatype();
	
	V3DLONG N = image->getTotalBytes();
	
	Image4DSimple p4DImage;
	
	unsigned char* image1d = image->getRawData();
	V3DLONG sx=image->getXDim(), sy=image->getYDim(), sz=image->getZDim(), sc=image->getCDim();
	
	V3DLONG ts_max = (sz>sc)?sz:sc;
	
	//get Z->C or C->Z command
	bool ok1, ok2;
    int timepoints=0;

	QStringList items;
	items << QObject::tr("4D {XYZ,Color} Stack --> 5D {XYZ,Color,Time} Stack") << QObject::tr("5D {XYZ,Color,Time} Stack --> 4D {XYZ,Color} Stack");

    QString item = QInputDialog::getItem(parent, QObject::tr("Change Movie Stack"),
                                        QObject::tr("Which direction do you want to change:"), items, 0, false, &ok1);

	
    if(ok1)
    {
        timepoints = QInputDialog::getInt(parent, QObject::tr("Set time points"),
                                      QObject::tr("Enter the number of time points:"),
                                      1, 1, ts_max, 1, &ok2);
    }
    else
    {
        return -1;
    }
	
	void *pResult = NULL;
	
	//Get the old image info
	if(!(QString::compare(item, "4D {XYZ,Color} Stack --> 5D {XYZ,Color,Time} Stack")))
	{
		p4DImage.setTimePackType(TIME_PACK_C);

//		V3DLONG imagecount = image->getTDim();
		
		V3DLONG pagesz=sx*sy;
		V3DLONG channelsz=sx*sy*sz;
		
		p4DImage.setTDim(timepoints);
		V3DLONG imagecount = timepoints;
		
		sz /= imagecount;
		
		if(imgdatatype == V3D_UINT8)
		{
			// init
			unsigned char *data1d = NULL;
			
			try
			{
				data1d = new unsigned char [N];
				
				memset(data1d, 0, sizeof(unsigned char)*N);
			}
			catch (...) 
			{
				printf("Fail to allocate memory.\n");
				return -1;
			}
			
			// assign
			unsigned char *pImg = (unsigned char *)image1d;
			for(V3DLONG no=0; no<imagecount; no++)
			{
				for(V3DLONG c=0; c<sc; c++)
				{
					for(V3DLONG k =0; k<sz; k++)
					{
						V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
						V3DLONG offsets = k*pagesz + c*channelsz + no*pagesz*sz;
						for(V3DLONG i=0; i<pagesz; i++)
						{
							data1d[i+offsetc] = pImg[i + offsets];
						}
					}
				}
			}
			
			sc = sc*imagecount;
			
			pResult = data1d; //
		}
		else if(imgdatatype == V3D_UINT16)
		{
			// init
			unsigned short *data1d = NULL;
			
			try
			{
				data1d = new unsigned short [N];
				
				memset(data1d, 0, sizeof(unsigned short)*N);
			}
			catch (...) 
			{
				printf("Fail to allocate memory.\n");
				return -1;
			}
			
			// assign
			unsigned short *pImg = (unsigned short *)image1d;
			for(V3DLONG no=0; no<imagecount; no++)
			{
				for(V3DLONG c=0; c<sc; c++)
				{
					for(V3DLONG k =0; k<sz; k++)
					{
						V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
						V3DLONG offsets = k*pagesz + c*channelsz + no*pagesz*sz;
						for(V3DLONG i=0; i<pagesz; i++)
						{
							data1d[i+offsetc] = pImg[i + offsets];
						}
					}
				}
			}
			
			sc = sc*imagecount;
			
			pResult = data1d; //
		}
		else if(imgdatatype == V3D_FLOAT32)
		{
			// init
			float *data1d = NULL;
			
			try
			{
				data1d = new float [N];
				
				memset(data1d, 0, sizeof(float)*N);
			}
			catch (...) 
			{
				printf("Fail to allocate memory.\n");
				return -1;
			}
			
			// assign
			float *pImg = (float *)image1d;
			for(V3DLONG no=0; no<imagecount; no++)
			{
				for(V3DLONG c=0; c<sc; c++)
				{
					for(V3DLONG k =0; k<sz; k++)
					{
						V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
						V3DLONG offsets = k*pagesz + c*channelsz + no*pagesz*sz;
						for(V3DLONG i=0; i<pagesz; i++)
						{
							data1d[i+offsetc] = pImg[i + offsets];
						}
					}
				}
			}
			
			sc = sc*imagecount;
			
			pResult = data1d; //
		}
		else
		{
			printf("Currently this program only support UINT8, UINT16, and FLOAT32 datatype.\n");
			return -1;
		}
	}
	else if(!(QString::compare(item, "5D {XYZ,Color,Time} Stack --> 4D {XYZ,Color} Stack")))
	{
		p4DImage.setTimePackType(TIME_PACK_Z);
		
		V3DLONG pagesz=sx*sy;
//		V3DLONG channelsz=sx*sy*sz;
//		V3DLONG imagecount = image->getTDim();
		
		V3DLONG imagecount = timepoints;
		
		if(imagecount>sc)
		{
			QMessageBox::information(0, title, QObject::tr("# time points should not be greater than # color channel."));
			return -1;
		}

		sc /= imagecount;

		if(imgdatatype == V3D_UINT8)
		{
			// init
			unsigned char *data1d = NULL;
			
			try
			{
				data1d = new unsigned char [N];
				
				memset(data1d, 0, sizeof(unsigned char)*N);
			}
			catch (...) 
			{
				printf("Fail to allocate memory.\n");
				return -1;
			}
			
			// assign
			unsigned char *pImg = (unsigned char *)image1d;
			for(V3DLONG no=0; no<imagecount; no++)
			{
				for(V3DLONG c=0; c<sc; c++)
				{
					for(V3DLONG k =0; k<sz; k++)
					{
						V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
						V3DLONG offsets = k*pagesz + c*pagesz*sz*imagecount + no*pagesz*sz;
						for(V3DLONG i=0; i<pagesz; i++)
						{
							data1d[i+offsets] = pImg[i + offsetc];
						}
					}
				}
			}
			
			sz = sz*imagecount;
			
			pResult = data1d; //
		}
		else if(imgdatatype == V3D_UINT16)
		{
			// init
			unsigned short *data1d = NULL;
			
			try
			{
				data1d = new unsigned short [N];
				
				memset(data1d, 0, sizeof(unsigned short)*N);
			}
			catch (...) 
			{
				printf("Fail to allocate memory.\n");
				return -1;
			}
			
			// assign
			unsigned short *pImg = (unsigned short *)image1d;
			for(V3DLONG no=0; no<imagecount; no++)
			{
				for(V3DLONG c=0; c<sc; c++)
				{
					for(V3DLONG k =0; k<sz; k++)
					{
						V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
						V3DLONG offsets = k*pagesz + c*pagesz*sz*imagecount + no*pagesz*sz;
						for(V3DLONG i=0; i<pagesz; i++)
						{
							data1d[i+offsets] = pImg[i + offsetc];
						}
					}
				}
			}
			
			sz = sz*imagecount;
			
			pResult = data1d; //
		}
		else if(imgdatatype == V3D_FLOAT32)
		{
			// init
			float *data1d = NULL;
			
			try
			{
				data1d = new float [N];
				
				memset(data1d, 0, sizeof(float)*N);
			}
			catch (...) 
			{
				printf("Fail to allocate memory.\n");
				return -1;
			}
			
			// assign
			float *pImg = (float *)image1d;
			for(V3DLONG no=0; no<imagecount; no++)
			{
				for(V3DLONG c=0; c<sc; c++)
				{
					for(V3DLONG k =0; k<sz; k++)
					{
						V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
						V3DLONG offsets = k*pagesz + c*pagesz*sz*imagecount + no*pagesz*sz;
						for(V3DLONG i=0; i<pagesz; i++)
						{
							data1d[i+offsets] = pImg[i + offsetc];
						}
					}
				}
			}
			
			sz = sz*imagecount;
			
			pResult = data1d; //
		}
		else
		{
			printf("Currently this program only support UINT8, UINT16, and FLOAT32 datatype.\n");
			return -1;
		}

	}
	else
	{
		QMessageBox::information(0, title, QObject::tr("This program only supports time series data. Your current image data type is not supported."));
		return -3;
	}

	// show in v3d
	p4DImage.setData((unsigned char*)pResult, sx,sy,sz,sc, imgdatatype);
	
	v3dhandle newwin = callback.newImageWindow();
	callback.setImage(newwin, &p4DImage);
	callback.setImageName(newwin,  callback.getImageName(oldwin)+"_changed");
	callback.updateImageWindow(newwin);

	return 0;
}

int packInChannel(V3DPluginCallback2 &callback, QWidget *parent)
{
    v3dhandle oldwin = callback.currentImageWindow();
    Image4DSimple* image = callback.getImage(oldwin);
    if (! image)
    {
        QMessageBox::information(0, title, QObject::tr("No image is open."));
        return -1;
    }

    V3DLONG sz=image->getZDim(), sc=image->getCDim();
    int timepoints=sc;

    V3DLONG ts_max = (sz>sc)?sz:sc;
    bool ok;
    timepoints = QInputDialog::getInt(parent, QObject::tr("Set time points"),
                                          QObject::tr("Enter the number of time points:"),
                                          1, 1, ts_max, 1, &ok);
    if(!ok)
        return -1;

    image->setTimePackType(TIME_PACK_C);
    image->setTDim(timepoints);
    callback.updateImageWindow(oldwin);
    callback.open3DWindow(oldwin);
    return 0;
}

int import5D(V3DPluginCallback2 &callback, QWidget *parent)
{
    v3dhandle oldwin = callback.currentImageWindow();
    Image4DSimple* image = callback.getImage(oldwin);
    if (! image)
    {
        QMessageBox::information(0, title, QObject::tr("No image is open."));
        return -1;
    }

    ImagePixelType imgdatatype = image->getDatatype();
    V3DLONG sx=image->getXDim(), sy=image->getYDim(), sz=image->getZDim(), sc=image->getCDim();

    QStringList imgList = importSeriesFileList_addnumbersort(callback.getImageName(oldwin));
    V3DLONG st = imgList.size();

    V3DLONG N = sx*sy*sz*sc*st;

    Image4DSimple p4DImage;
    void *pResult = NULL;

    if(imgdatatype == V3D_UINT8)
    {
        // init
        unsigned char *data5d = NULL;

        try
        {
            data5d = new unsigned char [N];

            memset(data5d, 0, sizeof(unsigned char)*N);
        }
        catch (...)
        {
            printf("Fail to allocate memory.\n");
            return -1;
        }

        // assign
        V3DLONG j = 0;
        for(V3DLONG no=0; no<st; no++)
        {
            unsigned char * data1d = 0;
            V3DLONG in_sz[4];
            int datatype;

            if (!simple_loadimage_wrapper(callback,imgList.at(no).toStdString().c_str(), data1d, in_sz, datatype))
            {
                fprintf (stderr, "Error happens in reading the subject file [%0]. Exit. \n",imgList.at(no).toStdString().c_str());
                return -1;
            }

            for(V3DLONG i=0; i<sx*sy*sz*sc; i++)
            {
                data5d[j] = data1d[i];
                j++;
            }
            if(data1d) {delete []data1d; data1d = 0;}

        }
        pResult = data5d; //
    }
    else if(imgdatatype == V3D_UINT16)
    {
        // init
        unsigned short *data5d = NULL;

        try
        {
            data5d = new unsigned short [N];

            memset(data5d, 0, sizeof(unsigned short)*N);
        }
        catch (...)
        {
            printf("Fail to allocate memory.\n");
            return -1;
        }

        // assign
        V3DLONG j = 0;
        for(V3DLONG no=0; no<st; no++)
        {
            unsigned short * data1d = 0;
            V3DLONG in_sz[4];
            int datatype;

            if (!simple_loadimage_wrapper(callback,imgList.at(no).toStdString().c_str(), (unsigned char * &)data1d, in_sz, datatype))
            {
                fprintf (stderr, "Error happens in reading the subject file [%0]. Exit. \n",imgList.at(no).toStdString().c_str());
                return -1;
            }

            for(V3DLONG i=0; i<sx*sy*sz*sc; i++)
            {
                data5d[j] = data1d[i];
                j++;
            }
            if(data1d) {delete []data1d; data1d = 0;}

        }

        pResult = data5d; //
    }
    else if(imgdatatype == V3D_FLOAT32)
    {
        // init
        float *data5d = NULL;

        try
        {
            data5d = new float [N];

            memset(data5d, 0, sizeof(float)*N);
        }
        catch (...)
        {
            printf("Fail to allocate memory.\n");
            return -1;
        }

        // assign
        V3DLONG j = 0;
        for(V3DLONG no=0; no<st; no++)
        {
            float * data1d = 0;
            V3DLONG in_sz[4];
            int datatype;

            if (!simple_loadimage_wrapper(callback,imgList.at(no).toStdString().c_str(), (unsigned char * &)data1d, in_sz, datatype))
            {
                fprintf (stderr, "Error happens in reading the subject file [%0]. Exit. \n",imgList.at(no).toStdString().c_str());
                return -1;
            }

            for(V3DLONG i=0; i<sx*sy*sz*sc; i++)
            {
                data5d[j] = data1d[i];
                j++;
            }
            if(data1d) {delete []data1d; data1d = 0;}

        }

        pResult = data5d; //
    }
    else
    {
        printf("Currently this program only support UINT8, UINT16, and FLOAT32 datatype.\n");
        return -1;
    }

    p4DImage.setTimePackType(TIME_PACK_C);
    p4DImage.setTDim(st);
    p4DImage.setData((unsigned char*)pResult, sx,sy,sz,sc*st, imgdatatype);

    v3dhandle newwin = callback.newImageWindow();
    callback.setImage(newwin, &p4DImage);
    callback.setImageName(newwin,  callback.getImageName(oldwin)+"_imported");
    callback.open3DWindow(newwin);
    callback.updateImageWindow(newwin);

    return 0;
}


bool changeMS(V3DPluginCallback2 &callback, const V3DPluginArgList & input, V3DPluginArgList & output)
{
    cout<<"Welcome to 5D Stack Converter plugin"<<endl;
    if (output.size() != 1) return false;
    unsigned int timepoints = 1;
    if (input.size()>=1)
    {
        vector<char*> paras = (*(vector<char*> *)(input.at(1).p));
        cout<<paras.size()<<endl;
        if(paras.size() >= 1) timepoints = atoi(paras.at(0));
    }

    char * inimg_file = ((vector<char*> *)(input.at(0).p))->at(0);
    char * outfolder_file = ((vector<char*> *)(output.at(0).p))->at(0);

    cout<<"timepoints = "<<timepoints<<endl;
    cout<<"inimg_file = "<<inimg_file<<endl;
    cout<<"outfolder_file = "<<outfolder_file<<endl;

    unsigned char * image1d = 0;
    V3DLONG in_sz[4];

    int datatype;
    if(!simple_loadimage_wrapper(callback, inimg_file, image1d, in_sz, datatype))
    {
        cerr<<"load image "<<inimg_file<<" error!"<<endl;
        return false;
    }

    V3DLONG sx=in_sz[0], sy = in_sz[1], sz=in_sz[2], sc=in_sz[3];
    V3DLONG ts_max = (sz>sc)?sz:sc;

    if(timepoints < 1 || timepoints > ts_max)
    {
        cerr<<"Invalid timepoints, error!"<<endl;
        return false;
    }

    V3DLONG pagesz=sx*sy;
    V3DLONG channelsz=sx*sy*sz;
    V3DLONG imagecount = timepoints;

    sz /= imagecount;
    V3DLONG N_times = sx*sy*sz*sc;

    V3DLONG in_sz_times[4];
    in_sz_times[0] = sx;
    in_sz_times[1] = sy;
    in_sz_times[2] = sz;
    in_sz_times[3] = sc;

    if(datatype == 1)
    {
        // assign
        unsigned char *pImg = (unsigned char *)image1d;
        for(V3DLONG no=0; no<imagecount; no++)
        {
            V3DLONG d = 0;
            // init
            unsigned char *data1d = NULL;

            try
            {
                data1d = new unsigned char [N_times];

                memset(data1d, 0, sizeof(unsigned char)*N_times);
            }
            catch (...)
            {
                printf("Fail to allocate memory.\n");
                return false;
            }

            for(V3DLONG c=0; c<sc; c++)
            {
                for(V3DLONG k =0; k<sz; k++)
                {
                    // V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
                    V3DLONG offsets = k*pagesz + c*channelsz + no*pagesz*sz;
                    for(V3DLONG i=0; i<pagesz; i++)
                    {
                        data1d[d] = pImg[i + offsets];
                        d++;
                    }
                }
            }
            QString outimg_file = QString(outfolder_file)+ QString("%1.v3draw").arg(no);
            simple_saveimage_wrapper(callback, outimg_file.toStdString().c_str(), (unsigned char *)data1d, in_sz_times, datatype);
            if(data1d) {delete []data1d; data1d=0;}
        }
    }
    else if(datatype == 2)
    {
        // assign
        unsigned short *pImg = (unsigned short *)image1d;
        for(V3DLONG no=0; no<imagecount; no++)
        {
            V3DLONG d = 0;
            // init
            unsigned short *data1d = NULL;

            try
            {
                data1d = new unsigned short [N_times];

                memset(data1d, 0, sizeof(unsigned short)*N_times);
            }
            catch (...)
            {
                printf("Fail to allocate memory.\n");
                return false;
            }

            for(V3DLONG c=0; c<sc; c++)
            {
                for(V3DLONG k =0; k<sz; k++)
                {
                    // V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
                    V3DLONG offsets = k*pagesz + c*channelsz + no*pagesz*sz;
                    for(V3DLONG i=0; i<pagesz; i++)
                    {
                        data1d[d] = pImg[i + offsets];
                        d++;
                    }
                }
            }

            QString outimg_file = QString(outfolder_file)+ QString("%1.v3draw").arg(no);
            simple_saveimage_wrapper(callback, outimg_file.toStdString().c_str(), (unsigned char *)data1d, in_sz_times, datatype);
            if(data1d) {delete []data1d; data1d=0;}
        }
    }
    else if(datatype == 4)
    {
        // init
        // assign
        float *pImg = (float *)image1d;
        for(V3DLONG no=0; no<imagecount; no++)
        {
            V3DLONG d = 0;

            float *data1d = NULL;

            try
            {
                data1d = new float [N_times];

                memset(data1d, 0, sizeof(float)*N_times);
            }
            catch (...)
            {
                printf("Fail to allocate memory.\n");
                return false;
            }
            for(V3DLONG c=0; c<sc; c++)
            {
                for(V3DLONG k =0; k<sz; k++)
                {
                    // V3DLONG offsetc = k*pagesz + c*pagesz*sz + no*pagesz*sz*sc;
                    V3DLONG offsets = k*pagesz + c*channelsz + no*pagesz*sz;
                    for(V3DLONG i=0; i<pagesz; i++)
                    {
                        data1d[d] = pImg[i + offsets];
                        d++;
                    }
                }
            }
            QString outimg_file = QString(outfolder_file)+ QString("%1.v3draw").arg(no);
            simple_saveimage_wrapper(callback, outimg_file.toStdString().c_str(), (unsigned char *)data1d, in_sz_times, datatype);
            if(data1d) {delete []data1d; data1d=0;}
        }
    }
    else
    {
        printf("Currently this program only support UINT8, UINT16, and FLOAT32 datatype.\n");
        return false;
    }

    return true;
}
