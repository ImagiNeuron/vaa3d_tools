
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
QT += widgets
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/common_lib/include

HEADERS	+= soma_segmentation_plugin.h
HEADERS += ResolutionDialog.h
SOURCES	+= soma_segmentation_plugin.cpp
SOURCES += ResolutionDialog.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp

TARGET	= $$qtLibraryTarget(soma_segmentation)
DESTDIR	= $$VAA3DPATH/bin/plugins/MIND/soma_segmentation/
