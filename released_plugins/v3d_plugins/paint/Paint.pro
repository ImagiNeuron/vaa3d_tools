
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH =  ../../../../v3d_external
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
QT += widgets
QT += printsupport
HEADERS	+= Paint_plugin.h \
    scribblearea.h \
   # mainwindow.h \
   # mydialog.h \
    paint_dialog.h
SOURCES	+= Paint_plugin.cpp \
    scribblearea.cpp \
   # mainwindow.cpp \
   # mydialog.cpp \
    paint_dialog.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp

TARGET	= $$qtLibraryTarget(paint)
DESTDIR	= $$VAA3DPATH/bin/plugins/misc/paint/

RESOURCES += \
    paint.qrc


