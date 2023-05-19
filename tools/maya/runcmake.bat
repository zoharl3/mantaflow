
@echo off
cls

@SETLOCAL

mkdir build
cd build

cmake ^
-DCMAKE_TOOLCHAIN_FILE=C:/prj-external-libs/vcpkg/scripts/buildsystems/vcpkg.cmake ^ -DVCPKG_TARGET_TRIPLET=x64-windows ^
-DMAYA_VERSION=2023 ^
-Wno-dev -G "Visual Studio 17 2022" -T "v142" ..

cd ..

pause
