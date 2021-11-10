/* ZMovieMaker_plugin.h
 * This plugin can be used to generate a smooth movie by several points
 * 2013-11-21 : by Zhi Zhou
 */
 
#ifndef __ZMOVIEMAKER_PLUGIN_H__
#define __ZMOVIEMAKER_PLUGIN_H__

#include <QtGui>

#include "v3d_interface.h"
#include <QComboBox>
#include <QGridLayout>
#include <QLabel>
#include <QSpinBox>
#include <QListWidget>
#include <QListWidgetItem>

class ZMovieMaker : public QObject, public V3DPluginInterface2_1
{
	Q_OBJECT
	Q_INTERFACES(V3DPluginInterface2_1);
    Q_PLUGIN_METADATA(IID"com.janelia.v3d.V3DPluginInterface/2.1")

public:
    float getPluginVersion() const {return 0.96f;}

	QStringList menulist() const;
	void domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent);

	QStringList funclist() const ;
	bool dofunc(const QString &func_name, const V3DPluginArgList &input, V3DPluginArgList &output, V3DPluginCallback2 &callback, QWidget *parent);
};


class MyComboBox : public QComboBox
{
    Q_OBJECT

public:
    V3DPluginCallback2 * m_v3d;
    MyComboBox(V3DPluginCallback2 * ini_v3d) {m_v3d = ini_v3d;}

    void enterEvent(QEnterEvent * event);

public slots:
    void updateList();
};

class controlPanel: public QDialog
{
    Q_OBJECT

public:
    controlPanel(V3DPluginCallback2 &v3d, QWidget *parent);
    ~controlPanel();

public:
    v3dhandle curwin;
    V3dR_MainWindow *surface_win;
    View3DControl *view;
    static controlPanel *panel;


    QList <V3dR_MainWindow *> list_3dviewer;
    v3dhandleList list_triview;

    V3DPluginCallback2 & m_v3d;

    QGridLayout *gridLayout;
    QListWidget *list_anchors;
    QSpinBox* box_SampleRate;

    MyComboBox* combo_surface;
    QLabel* label_surface;

    bool saveAnchorFile(QString filename);

private slots:
    void _slot_record();
    void _slot_preview();
    void _slot_show_item(QListWidgetItem *item);
    void _slot_show();
    void _slot_delete();
    void _slot_save();
    void _slot_load();
    void _slot_up();
    void _slot_down();
    void _slot_batch();
    void _slot_snapshot();

};

bool _saveAnchorFile(QString filename, QStringList ParaLists, bool b_append);


#endif

