:: build a particular plugin in the v3d_plugins directory
:: (window version) by Shidan Javaheri
:: 2024-11-16
:: Best add release argument to make release\v3d.exe can recognize plugin
:: revised from the original bat file
@echo off
set PATH=%PATH%;
cd v3d_plugins

:input_loop
:: Prompt the user for a directory name if not provided as an argument
if "%1"=="" (
    set /p DIR_NAME="Enter the directory name: "
) else (
    set DIR_NAME=%1
)

:: Check if the directory exists
if not exist %DIR_NAME% (
    echo Error: Directory "%DIR_NAME%" does not exist.
    set DIR_NAME=
    set "1="
    goto input_loop
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

:: Build the plugin
if exist *.pro (
    qmake
    mingw32-make clean
    mingw32-make -f Makefile.Release
)

cd ..
goto :eof
