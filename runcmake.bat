
@echo off
cls

@SETLOCAL

mkdir build
cd build

cmake ^
-DCMAKE_TOOLCHAIN_FILE=d:/vcpkg/scripts/buildsystems/vcpkg.cmake ^ -DVCPKG_TARGET_TRIPLET=x64-windows ^
-DPYTHON_VERSION=311 ^
-DGUI=ON ^
-DDEBUG=OFF ^
-DDOUBLEPRECISION=ON ^
-DPREPDEBUG=ON ^
-DDEBUG_PYTHON_WITH_RELEASE=ON ^
-DOPENVDB=ON ^
-DOPENVDB_BLOSC=ON ^
-DHOUDINI_ROOT="C:/Program Files/Side Effects Software/Houdini 20.5.410/" ^
-DTBB=ON ^
-DNUMPY=ON ^
-Wno-dev -G "Visual Studio 17 2022" ..

cd ..

pause
