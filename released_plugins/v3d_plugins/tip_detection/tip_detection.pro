
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external/v3d_main
INCLUDEPATH	+= $$VAA3DPATH/basic_c_fun
QT += widgets
HEADERS	+= tip_detection_plugin.h
SOURCES	+= tip_detection_plugin.cpp
SOURCES	+= $$VAA3DPATH/basic_c_fun/v3d_message.cpp

TARGET	= $$qtLibraryTarget(tip_detection)
DESTDIR	= $$VAA3DPATH/bin/plugins/tip_detection/
