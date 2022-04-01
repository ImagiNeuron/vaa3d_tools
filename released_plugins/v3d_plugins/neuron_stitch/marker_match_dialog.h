#ifndef MARKER_MATCH_DIALOG_H
#define MARKER_MATCH_DIALOG_H

#include <QDialog>
#include <v3d_interface.h>
#include "neuron_stitch_func.h"
#include <QGridLayout>
#include <QComboBox>
#include <QListWidget>
//#include <QDoubleSpinBox>
//#include <QCheckBox>


class marker_match_dialog: public QDialog
{
    Q_OBJECT

public:
    marker_match_dialog(V3DPluginCallback2 * cb, LandmarkList * mList);

private:
    LandmarkList * mList;
    QList<int> m0, m1;
    QList<int> mm0, mm1;

public:
    QGridLayout *gridLayout;
    QListWidget * list_pairs;
    QPushButton * btn_add;
    QPushButton * btn_free;
    QPushButton * btn_done;
    QComboBox * combobox_n0;
    QComboBox * combobox_n1;

private slots:
    void add();
    void free();
};

#endif // MARKER_MATCH_DIALOG_H
