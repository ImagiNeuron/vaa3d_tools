/* IVSCC_radius_estimation_plugin.h
 * This is a test plugin, you can use it as a demo.
 * 2015-11-19 : by Zhi Zhou
 */
 
#ifndef __IVSCC_RADIUS_ESTIMATION_PLUGIN_H__
#define __IVSCC_RADIUS_ESTIMATION_PLUGIN_H__

#include <QtGui>
#include <v3d_interface.h>

class IVSCC_radius_estimation : public QObject, public V3DPluginInterface2_1
{
	Q_OBJECT
	Q_INTERFACES(V3DPluginInterface2_1);

public:
	~IVSCC_radius_estimation();

	float getPluginVersion() const {return 1.1f;}

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const ;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);
};

class radiusEstimationDialog : public QDialog
    {
        Q_OBJECT

    public:
        radiusEstimationDialog(V3DPluginCallback2 &cb, QWidget *parent)
        {
            image = 0;

            v3dhandle curwin = cb.currentImageWindow();
            if (!curwin)
            {
                v3d_msg("You don't have any images open in the main window.");
                return;
            }

            image = cb.getImage(curwin);

            if (!image)
            {
                v3d_msg("The image pointer is invalid. Ensure your data is valid and try again!");
                return;
            }

            QGridLayout * layout = new QGridLayout();
            inswc_box = new QLineEdit("");
            outswc_box = new QLineEdit("(optional)");
            bkg_thresh_box = new QDoubleSpinBox();
            bkg_thresh_box->setMinimum(0);
            bkg_thresh_box->setMaximum(255);
            bkg_thresh_box->setValue(5);
            winX_box = new QLineEdit("1000");
            winY_box = new QLineEdit("1000");


            stop_thresh_box = new QDoubleSpinBox();
            stop_thresh_box->setMinimum(0);
            stop_thresh_box->setMaximum(255);
            stop_thresh_box->setValue(5);

            is2d_checker = new QCheckBox("Is 2D radius");
            is2d_checker->setChecked(true);
            pPushButton_openFileDlg = new QPushButton(QObject::tr("..."));

            layout->addWidget(new QLabel("In swc path"), 0, 0, 1, 1);
            layout->addWidget(inswc_box, 0, 1, 1, 3);
            layout->addWidget(pPushButton_openFileDlg,0, 4, 1, 2);
            layout->addWidget(new QLabel("Background Threshold (soma area)"), 1, 0, 1, 1);
            layout->addWidget(bkg_thresh_box, 1, 1, 1, 5);
            layout->addWidget(new QLabel("The size of soma area (X****Y)"), 2, 0, 1, 1);
            layout->addWidget(winX_box, 2, 1, 1, 2);
            layout->addWidget(new QLabel("******"), 2, 3, 1, 1);
             layout->addWidget(winY_box, 2, 4, 1, 2);
            layout->addWidget(new QLabel("Background Threshold (other areas)"), 3, 0, 1, 1);
            layout->addWidget(stop_thresh_box, 3, 1, 1, 5);
            layout->addWidget(new QLabel("Out swc path"), 4, 0, 1, 1);
            layout->addWidget(outswc_box, 4, 1, 1, 5);
            layout->addWidget(is2d_checker, 5, 0, 1, 6);
            QPushButton * ok = new QPushButton("Ok");
            QPushButton * cancel = new QPushButton("Cancel");
            ok->setDefault(true);
            layout->addWidget(ok, 7, 0, 1, 3);
            layout->addWidget(cancel, 7, 3, 1, 3);
            this->setLayout(layout);
            connect(ok, SIGNAL(clicked()), this, SLOT(_slot_start()));
            connect(cancel, SIGNAL(clicked()), this, SLOT(reject()));
            connect(pPushButton_openFileDlg, SIGNAL(clicked()), this, SLOT(_slots_openFileDlg()));
        }

        ~radiusEstimationDialog(){}

        public slots:
        void update()
        {
            inswc_file = inswc_box->text().toStdString();
            if(inswc_file == "")
            {
                v3d_msg("Please select the input swc file.");
                return;
            }
            outswc_file =outswc_box->text().toStdString();
            if(outswc_file == "" || outswc_file == "(optional)") outswc_file = inswc_file + ".out.swc";
            is_2d = is2d_checker->isChecked();
            bkg_thresh = atof(bkg_thresh_box->text().toStdString().c_str());
            stop_thresh = atof(stop_thresh_box->text().toStdString().c_str());
            winX_size = atof(winX_box->text().toStdString().c_str());
            winY_size = atof(winY_box->text().toStdString().c_str());
        }
        void _slots_openFileDlg()
        {
            QString fileOpenName;
            fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Choose swc file:"),
                                                        "",
                                                        QObject::tr("Supported file (*.swc *.eswc)"
                                                            ));
            if(!fileOpenName.isEmpty())
            {
                inswc_box->setText(fileOpenName);
            }
            update();
        }
        void _slot_start()
        {
            update();
            if(inswc_file != "") accept();
        }
    public:
        QLineEdit * inswc_box;
        QLineEdit * outswc_box;
        QDoubleSpinBox * bkg_thresh_box;
        QDoubleSpinBox *  stop_thresh_box;
        QCheckBox * is2d_checker;
        QPushButton *pPushButton_openFileDlg;
        Image4DSimple* image;
        QLineEdit * winX_box;
        QLineEdit * winY_box;


        string inswc_file;
        string outswc_file;
        bool is_2d;
        double bkg_thresh;
        double stop_thresh;
        double winX_size;
        double winY_size;

    };

class radiusEstimationFolderDialog : public QDialog
    {
        Q_OBJECT

    public:
        radiusEstimationFolderDialog(V3DPluginCallback2 &cb, QWidget *parent)
        {

            QGridLayout * layout = new QGridLayout();
            infolder_box = new QLineEdit("");
            inswc_box = new QLineEdit("");
            outswc_box = new QLineEdit("(optional)");
            bkg_thresh_box = new QDoubleSpinBox();
            bkg_thresh_box->setMinimum(0);
            bkg_thresh_box->setMaximum(255);
            bkg_thresh_box->setValue(5);

            stop_thresh_box = new QDoubleSpinBox();
            stop_thresh_box->setMinimum(0);
            stop_thresh_box->setMaximum(255);
            stop_thresh_box->setValue(5);

            is2d_checker = new QCheckBox("Is 2D radius");
            is2d_checker->setChecked(true);
            pPushButton_openFileDlg = new QPushButton(QObject::tr("..."));
            pPushButton_openImageDlg = new QPushButton(QObject::tr("..."));

            layout->addWidget(new QLabel("Input IVSCC image folder path"), 0, 0, 1, 1);
            layout->addWidget(infolder_box, 0, 1, 1, 3);
            layout->addWidget(pPushButton_openImageDlg,0, 4, 1, 2);
            layout->addWidget(new QLabel("Input swc path"), 1, 0, 1, 1);
            layout->addWidget(inswc_box, 1, 1, 1, 3);
            layout->addWidget(pPushButton_openFileDlg,1, 4, 1, 2);
            layout->addWidget(new QLabel("Background Threshold (soma area)"), 2, 0, 1, 1);
            layout->addWidget(bkg_thresh_box, 2, 1, 1, 5);
            layout->addWidget(new QLabel("Stopping Threshold (other areas)"), 3, 0, 1, 1);
            layout->addWidget(stop_thresh_box, 3, 1, 1, 5);
            layout->addWidget(new QLabel("Out swc path"), 4, 0, 1, 1);
            layout->addWidget(outswc_box, 4, 1, 1, 5);
            layout->addWidget(is2d_checker, 5, 0, 1, 6);
            QPushButton * ok = new QPushButton("Ok");
            QPushButton * cancel = new QPushButton("Cancel");
            ok->setDefault(true);
            layout->addWidget(ok, 7, 0, 1, 3);
            layout->addWidget(cancel, 7, 3, 1, 3);
            this->setLayout(layout);
            connect(ok, SIGNAL(clicked()), this, SLOT(_slot_start()));
            connect(cancel, SIGNAL(clicked()), this, SLOT(reject()));
            connect(pPushButton_openFileDlg, SIGNAL(clicked()), this, SLOT(_slots_openFileDlg()));
            connect(pPushButton_openImageDlg, SIGNAL(clicked()), this, SLOT(_slots_openImageDlg()));

        }

        ~radiusEstimationFolderDialog(){}

        public slots:
        void update()
        {
            inswc_file = inswc_box->text().toStdString();
            image_folder = infolder_box->text().toStdString();

            outswc_file =outswc_box->text().toStdString();
            if(outswc_file == "" || outswc_file == "(optional)") outswc_file = inswc_file + ".out.swc";
            is_2d = is2d_checker->isChecked();
            bkg_thresh = atof(bkg_thresh_box->text().toStdString().c_str());
            stop_thresh = atof(stop_thresh_box->text().toStdString().c_str());
        }
        void _slots_openFileDlg()
        {
            QString fileOpenName;
            fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Choose swc file:"),
                                                        "",
                                                        QObject::tr("Supported file (*.swc *.eswc)"
                                                            ));
            if(!fileOpenName.isEmpty())
            {
                inswc_box->setText(fileOpenName);
            }
            update();
        }
        void _slots_openImageDlg()
        {
            QString fileOpenName;
            fileOpenName = QFileDialog::getExistingDirectory(0, QObject::tr("Choose the directory including all IVSCC images "),
                                                             QDir::currentPath(),
                                                             QFileDialog::ShowDirsOnly);
            if(!fileOpenName.isEmpty())
            {
                infolder_box->setText(fileOpenName);
            }
            update();
        }
        void _slot_start()
        {
            update();
            if(inswc_file == "" || image_folder == "")
            {
                v3d_msg("Please select the input swc file and input IVSCC image folder.");
                return;
            }
            else
                accept();
        }
    public:
        QLineEdit * inswc_box;
        QLineEdit * outswc_box;
        QDoubleSpinBox * bkg_thresh_box;
        QDoubleSpinBox * stop_thresh_box;
        QCheckBox * is2d_checker;
        QPushButton *pPushButton_openFileDlg;
        QPushButton *pPushButton_openImageDlg;

        QLineEdit * infolder_box;

        string inswc_file;
        string outswc_file;
        string image_folder;
        bool is_2d;
        double bkg_thresh;
        double stop_thresh;
    };

#endif

