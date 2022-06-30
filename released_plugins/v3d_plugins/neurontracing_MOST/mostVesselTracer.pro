
TEMPLATE      = lib
CONFIG       += qt plugin warn_off
QT+=widgets
#CONFIG      += x86_64

V3DMAINPATH = ../../../../v3d_external/v3d_main
QT_DIR = $$[QT_INSTALL_PREFIX]

INCLUDEPATH  += $$V3DMAINPATH/basic_c_fun 
INCLUDEPATH  += $$V3DMAINPATH/v3d
INCLUDEPATH  += $$V3DMAINPATH/neuron_editing
INCLUDEPATH  += $$V3DMAINPATH/common_lib/include

USE_Qt5 {
  INCLUDEPATH += $$QT_DIR/lib/QtConcurrent.framework/Versions/5/Headers  # for QtConcurrent, by PHC 2015May
  #SHARED_FOLDER = $$QT_DIR/demos/shared # for arthurwidgets
  SHARED_FOLDER = ./painting/shared/ # for arthurwidgets
  include($$SHARED_FOLDER/shared.pri)
  INCLUDEPATH += $$SHARED_FOLDER
  LIBS += -L$$SHARED_FOLDER
} else {
  SHARED_FOLDER = D:/software/Qt/Examples/Qt-6.2.2/activeqt # for arthurwidgets
  include($$SHARED_FOLDER/shared.pri)
  INCLUDEPATH += $$SHARED_FOLDER
  LIBS += -L$$SHARED_FOLDER
}

HEADERS       = mostimage.h
HEADERS      += tools.h
HEADERS      += mostVesselTracer.h
HEADERS      += voxelcluster.h
HEADERS      += srb.h
HEADERS      += $$V3DMAINPATH/neuron_editing/neuron_format_converter.h
HEADERS      += $$V3DMAINPATH/neuron_editing/v_neuronswc.h

SOURCES       = mostimage.cpp
SOURCES      += tools.cpp
SOURCES      += $$V3DMAINPATH/basic_c_fun/v3d_message.cpp
SOURCES      += mostVesselTracer.cpp
SOURCES      += voxelcluster.cpp
SOURCES      += srb.c
SOURCES      += $$V3DMAINPATH/neuron_editing/neuron_format_converter.cpp
SOURCES      += $$V3DMAINPATH/neuron_editing/v_neuronswc.cpp



TARGET        = $$qtLibraryTarget(mostVesselTracer)

DESTDIR       = $$V3DMAINPATH/../bin/plugins/neuron_tracing/MOST_tracing/



