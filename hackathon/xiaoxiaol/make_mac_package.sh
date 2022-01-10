#!/bin/bash
cmake_build_dir="/Users/xiaoxiaoliu/work/src/v3d_external/build_macosx-x86_64"
qmake_build_dir="/Users/xiaoxiaoliu/work/v3d/v3d_external/bin"
package_dir="/Users/xiaoxiaoliu/Downloads/Vaa3d_V3.20_MacOSX10.9_64bit"

# use cmake to build the main v3d executable
cd $cmake_build_dir
make -j8

# copy qmake built plugins to the cmake build older
cd $qmake_build_dir
find . -type d -empty -delete
cd $cmake_build_dir/v3d/Mac_Fat/Vaa3d.app/Contents/MacOS
rm -rf plugins
cp -r  $qmake_build_dir/plugins ./


# mac deploy all plugins: the change the qt dependency path
cd $cmake_build_dir
make DeployPlugins

# copy the executables and plugins into the pre-existing package
cp -rf $cmake_build_dir/v3d/Mac_Fat/Vaa3d.app $package_dir/
rm -rf $package_dir/plugins
cp -r $package_dir/Vaa3d.app/Contents/MacOS/plugins  $package_dir/

#run testing scripts
#cp ~/work/v3d/v3d_external/testing ~/Downloads/Vaa3d_V3.100_MacOSX10.9_64bit/
#cd testing
#./test_plugin_installation.sh


