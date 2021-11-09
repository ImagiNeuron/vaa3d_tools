
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external/v3d_main
INCLUDEPATH	+= $$VAA3DPATH/basic_c_fun
INCLUDEPATH += $$VAA3DPATH/common_lib/include
QT += widgets
HEADERS	+= inter_node_pruning_plugin.h
HEADERS      += ../sort_neuron_swc/openSWCDialog.h

SOURCES	+= inter_node_pruning_plugin.cpp
SOURCES	+= $$VAA3DPATH/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/basic_c_fun/basic_surf_objs.cpp
SOURCES	+= my_surf_objs.cpp
SOURCES      += ../sort_neuron_swc/openSWCDialog.cpp


TARGET	= $$qtLibraryTarget(inter_node_pruning)
DESTDIR	= $$VAA3DPATH/../bin/plugins/neuron_utilities/inter_node_pruning/
