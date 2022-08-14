
@echo off
cls

@SETLOCAL

mkdir build
cd build

cmake ^
-DCMAKE_TOOLCHAIN_FILE=C:/prj-external-libs/vcpkg/scripts/buildsystems/vcpkg.cmake ^ -DVCPKG_TARGET_TRIPLET=x64-windows ^
-DPYTHON_VERSION=37 ^
-DGUI=ON ^
-DDEBUG=OFF ^
-DDOUBLEPRECISION=ON ^
-DPREPDEBUG=ON ^
-DDEBUG_PYTHON_WITH_RELEASE=ON ^
-DPython_LIBRARY="c:/python37/libs/python37.lib" ^
-DPython_INCLUDE_DIR="c:/python37/include/" ^
-DOPENVDB=ON ^
-DOPENVDB_BLOSC=ON ^
-DOPENVDB_ROOT="c:/Program Files/OpenVDB/" ^
-DHOUDINI_ROOT="C:/Program Files/Side Effects Software/Houdini 19.0.657/" ^
-DTBB=ON ^
-DNUMPY=ON ^
-DMatlab_ROOT_DIR="%MATLAB_DIR%" ^
-Wno-dev -G "Visual Studio 17 2022" -T "v142" ..

cd ..

pause
