
@echo off
cls

@SETLOCAL
@SET CMAKE_EXEC="c:\prog\cmake-gui\bin\cmake.exe"
@IF EXIST %CMAKE_EXEC% GOTO START
@GOTO ERROR

:START
mkdir build
cd build
%CMAKE_EXEC% ^
-DPYTHON_VERSION=37 ^
-DGUI=ON ^
-DWIN_QT_PATH=c:\Qt\online\5.10.0\msvc2017_64\ ^
-Wno-dev -G"Visual Studio 16" -A"x64" ..

rem -Wno-dev -G"Visual Studio 16" -A"x64" -T v141 ..

cd ..

@GOTO EOF

@:ERROR
@echo ERROR: CMake not found.

:EOF
pause
