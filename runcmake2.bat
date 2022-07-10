
@echo off
cls

@SETLOCAL

mkdir build
cd build

cmake ^
-DCMAKE_TOOLCHAIN_FILE=C:/prj-external-libs/vcpkg/scripts/buildsystems/vcpkg.cmake ^ -DVCPKG_TARGET_TRIPLET=x64-windows ^
-DPYTHON_VERSION=37 ^
-DGUI=ON ^
-DDEBUG=ON ^
-DDOUBLEPRECISION=ON ^
-DPREPDEBUG=ON ^
-DDEBUG_PYTHON_WITH_RELEASE=ON ^
-DOPENVDB=ON ^
-DOPENVDB_BLOSC=ON ^
-DOPENVDB_ROOT="c:/Program Files/OpenVDB/" ^
-DHOUDINI_ROOT="C:/Program Files/Side Effects Software/Houdini 19.0.657/" ^
-DTBB=ON ^
-Wno-dev ..

rem -A x64 -T v141

cd ..

pause
