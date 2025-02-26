:: build the MIND plugin in the v3d_plugins directory
:: (window version) by Athmane Benarous
:: 2024-11-16
:: Best add release argument to make release\v3d.exe can recognize plugin
:: revised from the compile_plugin bat file
@echo off
set PATH=%PATH%;
cd v3d_plugins

:: Automatically set the directory name to MIND
set DIR_NAME=MIND

:: Check if the directory exists
if not exist %DIR_NAME% (
    echo Error: Directory "%DIR_NAME%" does not exist.
    goto :eof
)

:: Change to the specified directory
cd %DIR_NAME%

:: Build the plugin
if exist *.pro (
    qmake
    mingw32-make clean
    mingw32-make -f Makefile.Release
) else (
    echo Error: No .pro file found in the directory "%DIR_NAME%".
)

cd ..
goto :eof
