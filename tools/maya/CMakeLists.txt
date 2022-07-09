
cmake_minimum_required(VERSION 3.5)

project(maya_menta)

# OSX specific flags
if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

# Add path to Chad Vernons FindMaya.cmake modules
# https://github.com/chadmv/cgcmake
SET(AUTODESK_DIR /Applications/Autodesk)
list(APPEND CMAKE_MODULE_PATH ${AUTODESK_DIR}/Resources/CMake "${CMAKE_SOURCE_DIR}/cmake/")

# Output folder
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/lib)

# Find Maya and declare variables:
#   MAYA_FOUND          Defined if a Maya installation has been detected
#   MAYA_EXECUTABLE     Path to Maya's executable
#   MAYA_<lib>_FOUND    Defined if <lib> has been found
#   MAYA_<lib>_LIBRARY  Path to <lib> library
#   MAYA_INCLUDE_DIR    Path to the devkit's include directories
#   MAYA_LIBRARIES      All the Maya libraries
find_package(Maya REQUIRED)

set(SILENT_SOURCES)
set(INCLUDE_PATHS)

# Z lib compression
if(1)
	# default: build from own sources
	set(ZLIB_SRC adler32.c compress.c crc32.c deflate.c gzclose.c gzlib.c gzread.c gzwrite.c
				 inflate.c infback.c inftrees.c inffast.c trees.c uncompr.c zutil.c)
	foreach(it ${ZLIB_SRC})
		list(APPEND SILENT_SOURCES ../../../dependencies/zlib-1.2.8/${it}) 
	endforeach(it)

	set(ZLIB_ADDFLAGS "-Dverbose=-1")
	if(NOT WIN32)
		# otherwise we get warnings that we could only fix by upgrading zlib to a version > 1.2.8
		set(ZLIB_ADDFLAGS "-Wno-implicit-function-declaration -Wno-shift-negative-value -Dverbose=-1")
	endif()
	set_source_files_properties(${SILENT_SOURCES} PROPERTIES COMPILE_FLAGS "${ZLIB_ADDFLAGS}")

	list(APPEND INCLUDE_PATHS ../../dependencies/zlib-1.2.8)
endif()

# Build options
include_directories( ${MAYA_INCLUDE_DIR} ${INCLUDE_PATHS} )
link_directories( ${MAYA_LIBRARY_DIR} )

add_subdirectory( densityloader )
add_subdirectory( bobjloader )

