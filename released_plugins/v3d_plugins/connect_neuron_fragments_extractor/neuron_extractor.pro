
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
QT += widgets
HEADERS	+= neuron_extractor_plugin.h
SOURCES	+= neuron_extractor_plugin.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp

TARGET	= $$qtLibraryTarget(neuron_fragment_extractor)
DESTDIR	= $$VAA3DPATH/bin/plugins/neuron_utilities/neuron_fragment_extractor/
