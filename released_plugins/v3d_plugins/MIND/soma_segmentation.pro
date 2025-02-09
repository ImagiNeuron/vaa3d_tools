
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
QT += widgets
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external

INCLUDEPATH += $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH += $$VAA3DPATH/v3d_main/v3d
INCLUDEPATH += $$VAA3DPATH/v3d_main/jba/newmat11
INCLUDEPATH += $$VAA3DPATH/v3d_main/jba/c++
INCLUDEPATH += $$VAA3DPATH/v3d_main/common_lib/include

HEADERS	+= soma_segmentation_plugin.h
HEADERS += cellSegmentation_plugin.h
SOURCES	+= soma_segmentation_plugin.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES += $$VAA3DPATH/v3d_main/neuron_editing/v_neuronswc.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat1.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat2.cpp 
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat3.cpp 
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat4.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat5.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat6.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat7.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat8.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmatex.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/bandmat.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/submat.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/myexcept.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/cholesky.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/evalue.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/fft.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/hholder.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/jacobi.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newfft.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/sort.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/svd.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/nm_misc.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmatrm.cpp
SOURCES += $$VAA3DPATH/v3d_main/jba/newmat11/newmat9.cpp

TARGET	= $$qtLibraryTarget(soma_segmentation)
DESTDIR	= $$VAA3DPATH/bin/plugins/MIND/soma_segmentation/
