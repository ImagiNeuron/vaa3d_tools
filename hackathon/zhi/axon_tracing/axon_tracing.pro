
TEMPLATE	= lib
CONFIG	+= qt plugin warn_off
#CONFIG	+= x86_64
VAA3DPATH =  ../../../../v3d_external/v3d_main
INCLUDEPATH	+= $$VAA3DPATH/basic_c_fun
#INCLUDEPATH	+= $$VAA3DPATH/common_lib/include
INCLUDEPATH     += $$VAA3DPATH/common_lib/include/boost/graph

INCLUDEPATH     += $$VAA3DPATH/neuron_editing
INCLUDEPATH     += $$VAA3DPATH/worm_straighten_c
INCLUDEPATH     += $$VAA3DPATH/cellseg
INCLUDEPATH     += ../../../released_plugins/v3d_plugins/neurontracing_vn2/app2
INCLUDEPATH     += $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/include
HEADERS += $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/include/ANN/ANN.h
HEADERS += $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/bd_tree.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_fix_rad_search.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_pr_search.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_search.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_split.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_tree.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_util.h \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/pr_queue.h
SOURCES += $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/ANN.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/bd_fix_rad_search.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/bd_pr_search.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/bd_search.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/bd_tree.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/brute.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_dump.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_fix_rad_search.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_pr_search.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_search.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_split.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_tree.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/kd_util.cpp \
    $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/src/perf.cpp


macx{
    LIBS += -L$$VAA3DPATH/common_lib/lib_mac64 -lv3dtiff
    LIBS += -L$$VAA3DPATH/jba/c++ -lv3dnewmat

    PRE_TARGETDEPS = $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/lib/libANN.a
    LIBS            += $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/lib/libANN.a
    ANN.target = $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/lib/libANN.a
    ANN.commands = cd $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann && make macosx-g++
    ANN.depends = $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/Makefile
    QMAKE_EXTRA_TARGETS += ANN
}

win32 {
    contains(QMAKE_HOST.arch, x86_64) {
    LIBS     += -L$$VAA3DPATH/common_lib/winlib64 -llibtiff
    LIBS     += -L$$VAA3DPATH/common_lib/winlib64 -llibnewmat
    } else {
    LIBS     += -L$$VAA3DPATH/common_lib/winlib -llibtiff
    LIBS     += -L$$VAA3DPATH/common_lib/winlib -llibnewmat
    }
}

unix:!macx {
    LIBS += -L$$VAA3DPATH/common_lib/lib -ltiff
    LIBS += -L$$VAA3DPATH/jba/c++ -lv3dnewmat
    QMAKE_CXXFLAGS += -fopenmp
    LIBS += -fopenmp

    PRE_TARGETDEPS = $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/lib/libANN.a
    LIBS            += $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/lib/libANN.a
    ANN.target = $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/lib/libANN.a
    ANN.commands = cd $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann && make linux-g++
    ANN.depends = $$VAA3DPATH/../../vaa3d_tools/hackathon/xiaoxiaol/consensus_skeleton_2/ann/Makefile
    QMAKE_EXTRA_TARGETS += ANN
}


HEADERS	+= axon_tracing_plugin.h
HEADERS	+= $$VAA3DPATH/jba/c++/histeq.h
HEADERS += $$VAA3DPATH/basic_c_fun/mg_utilities.h \
    $$VAA3DPATH/basic_c_fun/basic_memory.h \
    $$VAA3DPATH/basic_c_fun/stackutil.h\
    $$VAA3DPATH/basic_c_fun/mg_image_lib.h\
    $$VAA3DPATH/basic_c_fun/img_definition.h \
    $$VAA3DPATH/basic_c_fun/volimg_proc_declare.h \
    $$VAA3DPATH/basic_c_fun/volimg_proc.h \
    $$VAA3DPATH/basic_c_fun/v3d_message.h \
    $$VAA3DPATH/basic_c_fun/color_xyz.h \
    $$VAA3DPATH/basic_c_fun/basic_surf_objs.h \
    $$VAA3DPATH/basic_c_fun/basic_landmark.h \
    $$VAA3DPATH/graph/graph.h \
    $$VAA3DPATH/graph/graph_basic.h \
    $$VAA3DPATH/graph/dijk.h \
    $$VAA3DPATH/worm_straighten_c/bdb_minus.h \
    $$VAA3DPATH/worm_straighten_c/mst_prim_c.h \
    $$VAA3DPATH/worm_straighten_c/bfs.h \
    $$VAA3DPATH/worm_straighten_c/spline_cubic.h \
    $$VAA3DPATH/neuron_editing/apo_xforms.h \
    $$VAA3DPATH/neuron_editing/neuron_xforms.h \
    $$VAA3DPATH/neuron_editing/neuron_sim_scores.h \
    $$VAA3DPATH/neuron_editing/v_neuronswc.h \
    $$VAA3DPATH/cellseg/template_matching_seg.h \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/FL_bwdist.h \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/gd.h \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/tip_detection.h \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/autoseed.h \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/v3dneuron_gd_tracing.h


SOURCES	+= axon_tracing_plugin.cpp
SOURCES	+= $$VAA3DPATH/basic_c_fun/v3d_message.cpp
SOURCES	+= $$VAA3DPATH/basic_c_fun/basic_surf_objs.cpp
SOURCES	+= $$VAA3DPATH/jba/c++/histeq.cpp

SOURCES += $$VAA3DPATH/basic_c_fun/mg_utilities.cpp \
    $$VAA3DPATH/basic_c_fun/basic_memory.cpp \
    $$VAA3DPATH/basic_c_fun/stackutil.cpp\
    $$VAA3DPATH/basic_c_fun/mg_image_lib.cpp\
    $$VAA3DPATH/graph/dijk.cpp \
    $$VAA3DPATH/worm_straighten_c/bdb_minus.cpp \
    $$VAA3DPATH/worm_straighten_c/mst_prim_c.cpp \
    $$VAA3DPATH/worm_straighten_c/bfs_1root.cpp \
    $$VAA3DPATH/worm_straighten_c/spline_cubic.cpp \
    $$VAA3DPATH/neuron_editing/apo_xforms.cpp \
    $$VAA3DPATH/neuron_editing/neuron_xforms.cpp \
    $$VAA3DPATH/neuron_editing/neuron_sim_scores.cpp \
    $$VAA3DPATH/neuron_editing/v_neuronswc.cpp \
    $$VAA3DPATH/cellseg/template_matching_seg.cpp \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/gd.cpp \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/tip_detection.cpp \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/autoseed.cpp \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/calculate_cover_scores.cpp \
    ../../../released_plugins/v3d_plugins/neurontracing_vn2/app1/v3dneuron_gd_tracing.cpp \



TARGET	= $$qtLibraryTarget(axon_tracing)
DESTDIR	= $$VAA3DPATH/../bin/plugins/axon_tracing/
