
# Source files
set(SOURCE_FILES "../bobjFluidObject.cpp")

set( NAME1 "bobjloader" )

add_library( ${NAME1} SHARED ${SOURCE_FILES} ${SILENT_SOURCES} )
target_link_libraries( ${NAME1} ${MAYA_LIBRARIES} )
MAYA_PLUGIN( ${NAME1} )


