
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external/v3d_main
INCLUDEPATH	+= $$VAA3DPATH/basic_c_fun
QT += widgets
HEADERS	+= demo_opensurfacefile_plugin.h
SOURCES	+= demo_opensurfacefile_plugin.cpp
SOURCES	+= $$VAA3DPATH/basic_c_fun/v3d_message.cpp

TARGET	= $$qtLibraryTarget(demo_opensurfacefile)
DESTDIR	= $$VAA3DPATH/../bin/plugins/Vaa3D_PluginInterface_Demos/demo_opensurfacefile/
