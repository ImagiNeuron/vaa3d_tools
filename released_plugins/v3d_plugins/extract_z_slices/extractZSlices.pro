
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH = ../../../../v3d_external/v3d_main
INCLUDEPATH += $$VAA3DPATH/basic_c_fun
INCLUDEPATH += $$VAA3DPATH/common_lib/include
QT += widgets
HEADERS	+= extractZSlices_plugin.h
SOURCES	+= extractZSlices_plugin.cpp

HEADERS += $$VAA3DPATH/basic_c_fun/basic_memory.h

SOURCES	+= $$VAA3DPATH/basic_c_fun/v3d_message.cpp
SOURCES += $$VAA3DPATH/basic_c_fun/basic_4dimage_create.cpp

SOURCES += $$VAA3DPATH/basic_c_fun/basic_memory.cpp

TARGET	= $$qtLibraryTarget(extractZSlices)
DESTDIR = $$VAA3DPATH/../bin/plugins/image_geometry/extract_Z_Slices/

