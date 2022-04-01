#!/bin/bash
# batch build all released V3D plugin projects


cd v3d_plugins
 
QMAKE_CONFIG=
MAKE_ARGS=
MYDIR=

for arg in $*; do
  #echo $arg		
  if [ $arg == "-m" ]; then
  	QMAKE_CONFIG="CONFIG+=x86_64"
  elif [ $arg == "-n" ]; then
  	QMAKE_CONFIG="CONFIG+=x86"
  elif [ ${arg:0:2} == "-d" ]; then
  	MYDIR="${arg:2}"
  else
  	MAKE_ARGS+=" $arg"
  fi
done

if [ ${#MYDIR} -gt 0 ]; then
  ALLDIRS=$MYDIR
else
  ALLDIRS=$( ls -d */ )
  #ALLDIRS=$( cat ../linux_plugins.txt )  
fi

# CMB 01 Dec, 2010
# Need to define QMAKESPEC on Mac
# because recent Qt installs default to creating xcode project files.
# We want Makefiles for this script.
if [[ `uname` == 'Darwin' ]]; then
   QMAKE_ARGS='-spec macx-g++'
else
   QMAKE_ARGS=''
fi

for mydir in $ALLDIRS; do
  echo 
  echo $mydir
  echo ===============================  
  cd $mydir

  if [ -f build.sh ]; then
    sh build.sh $MAKE_ARGS;
  else
    #if [ -f *.pro ]; then
    for mypro in $( ls *.pro ); do
  	  qmake $QMAKE_ARGS $mypro
    	  make $MAKE_ARGS 
    done;
    #fi
  fi;

  cd ..
done  

cd ..

cd ../../v3d_external/bin/plugins
echo $( find -depth -empty -type d -delete )

