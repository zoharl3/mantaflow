
@echo off
cls

@SETLOCAL

mkdir build
cd build

cmake ^
-DCMAKE_TOOLCHAIN_FILE=C:\prj-external-libs\vcpkg\scripts\buildsystems\vcpkg.cmake ^ -DVCPKG_TARGET_TRIPLET=x64-windows ^
-DPYTHON_VERSION=37 ^
-DGUI=ON ^
-DWIN_QT_PATH=c:\Qt\online\5.10.0\msvc2017_64\ ^
-DOPENVDB=1 ^
-DOPENVDB_ROOT=c:\Program Files\OpenVDB\ ^
-DTBB=1 ^
..

rem -A x64 -T v141

cd ..

pause
