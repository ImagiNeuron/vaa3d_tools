
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
V3DMAINPATH = ../../../../v3d_external/v3d_main
INCLUDEPATH	+= $$V3DMAINPATH/basic_c_fun
INCLUDEPATH += $$V3DMAINPATH/common_lib/include
INCLUDEPATH += $$V3DMAINPATH/jba/newmat11

macx{
    LIBS += -L$$V3DMAINPATH/common_lib/lib_mac64 -lv3dtiff
#    CONFIG += x86_64
}
unix:!macx {
    #LIBS += -L$$V3DMAINPATH/common_lib/lib -lv3dtiff
    LIBS += -L$$V3DMAINPATH/common_lib/lib -ltiff
    LIBS += -L$$V3DMAINPATH/jba/c++ -lv3dnewmat
}
win32 {
    contains(QMAKE_HOST.arch, x86_64) {
    LIBS     += -L$$VAA3DPATH/common_lib/winlib64 -llibtiff
    } else {
    LIBS     += -L$$VAA3DPATH/common_lib/winlib -llibtiff
    }
}

INCLUDEPATH += main


HEADERS += $$V3DMAINPATH/basic_c_fun/basic_memory.h
HEADERS += $$V3DMAINPATH/basic_c_fun/mg_utilities.h
HEADERS += $$V3DMAINPATH/basic_c_fun/mg_image_lib.h
HEADERS += $$V3DMAINPATH/basic_c_fun/stackutil.h
HEADERS	+= $$V3DMAINPATH/../../vaa3d_tools/released_plugins/v3d_plugins/swc_to_maskimage/filter_dialog.h

HEADERS	+= image_blend_plugin.h
SOURCES	+= image_blend_plugin.cpp

SOURCES	+= $$V3DMAINPATH/basic_c_fun/v3d_message.cpp
SOURCES += $$V3DMAINPATH/basic_c_fun/basic_memory.cpp
SOURCES += $$V3DMAINPATH/basic_c_fun/mg_utilities.cpp
SOURCES += $$V3DMAINPATH/basic_c_fun/mg_image_lib.cpp
SOURCES += $$V3DMAINPATH/basic_c_fun/stackutil.cpp
SOURCES	+= $$V3DMAINPATH/../../vaa3d_tools/released_plugins/v3d_plugins/swc_to_maskimage/filter_dialog.cpp
SOURCES	+= $$V3DMAINPATH/basic_c_fun/basic_surf_objs.cpp
SOURCES	+= $$V3DMAINPATH/v3d/colormap.cpp


TARGET	= $$qtLibraryTarget(blend_two_images)
DESTDIR	=  ../../../../v3d_external/bin/plugins/image_blending/blend_two_images
