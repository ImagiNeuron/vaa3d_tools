#ifndef PAINT_DIALOG_H
#define PAINT_DIALOG_H

#include <QtGui>
#include <QDialog>
#include "scribblearea.h"
#include <QToolBar>
#include <QVBoxLayout>
#include <QToolButton>
#include <v3d_interface.h>
#include <QSpinBox>
#include <QPlainTextEdit>
class Paint_Dialog : public QDialog
{
    Q_OBJECT
public:
    explicit Paint_Dialog(V3DPluginCallback2 *callback, QWidget *parent = 0);

public:
    ScribbleArea* paintarea;

signals:

public slots:
    bool load();
    bool saveimage();
    void penColor();
    void penWidth();
    void zdisplay(int);
    void fetch();
    bool pushback();
    void clearimage();
    void zoomin();
    void zoomout();
    bool saveFile();
    void dosavemenu();
    void doopenmenu();
    void help();
    bool save_label();

private:
    void create();
    bool maybeSave();
    unsigned char * datacopy(unsigned char *data,long size);
    void savezimage(int z);
    void convert2UINT8(unsigned short *pre1d, unsigned char *pPost, V3DLONG imsz);
    void convert2UINT8(float *pre1d, unsigned char *pPost, V3DLONG imsz);
    void resetdata();
    void closeEvent(QCloseEvent *event);
    QMenu *savemenu,*openmenu;
    QToolBar *tool;
    void createsavemenu();
    void createopenmenu();

    QString fileName;
    V3DPluginCallback2 *callback;
    QPlainTextEdit *edit;
    QSpinBox *spin;
    V3DLONG sz_img[4];
    int intype;
    int datasource; //loaded data:datasource=1; fetched data:datasource=2;
    int newWidth;
    v3dhandle curwin;
    bool zoominflag;
    int previousz;
    unsigned char *paint_1DC,*backupdata,*image1Dc_in;

protected:
    bool eventFilter(QObject *obj, QEvent *event);


};

#endif // PAINT_DIALOG_H
