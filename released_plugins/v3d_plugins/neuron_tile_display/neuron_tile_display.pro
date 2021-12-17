
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH =  ../../../../v3d_external/v3d_main
INCLUDEPATH	+= $$VAA3DPATH/v3d_main/basic_c_fun
INCLUDEPATH     += $$VAA3DPATH/v3d_main/common_lib/include
QT += widgets
HEADERS	+= neuron_tile_display_plugin.h \
    neuron_tile_display_dialog.h
SOURCES	+= neuron_tile_display_plugin.cpp \
    neuron_tile_display_dialog.cpp
SOURCES	+= $$VAA3DPATH/v3d_main/basic_c_fun/v3d_message.cpp
HEADERS += $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.h
SOURCES += $$VAA3DPATH/v3d_main/basic_c_fun/basic_surf_objs.cpp
HEADERS += $$VAA3DPATH/v3d_main/neuron_editing/neuron_xforms.h
SOURCES += $$VAA3DPATH/v3d_main/neuron_editing/neuron_xforms.cpp

TARGET	= $$qtLibraryTarget(neuron_tile_display)
DESTDIR	= $$VAA3DPATH/bin/plugins/neuron_utilities/tile_display_multiple_neurons/
