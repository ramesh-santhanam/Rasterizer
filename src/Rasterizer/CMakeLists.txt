message("Bin Path: ${SDK_VTK_BIN} ...")

set(CMAKE_PREFIX_PATH ${SDK_VTK_BIN} ${CMAKE_PREFIX_PATH})

find_package(VTK COMPONENTS
             vtkCommonDataModel  vtkFiltersGeometry NO_MODULE)

include(${VTK_USE_FILE})

include_directories(${VTK_INCLUDE_DIRS})

message("include dirs...: ${VTK_INCLUDE_DIRS} ...")

set( HEADER_FILES ${PROJECT_SOURCE_DIR}/src/Rasterizer/Rasterize2d.h )
set( HEADER_FILES ${HEADER_FILES} ${PROJECT_SOURCE_DIR}/src/Rasterizer/QuadTree.h ${PROJECT_SOURCE_DIR}/src/Rasterizer/Rasterer2d.h  ${PROJECT_SOURCE_DIR}/src/Rasterizer/Tree2ps.h )

set( HEADER_FILES ${HEADER_FILES} ${VTK_INCLUDE_DIRS} )

message("public include dirs...: ${HEADER_FILES} ...")

add_library(Rasterizer SHARED Rasterize2d.cpp QuadTree.cpp Rasterer2d.cpp Tree2PS.cpp ${HEADER_FILES})

target_include_directories(Rasterizer PUBLIC ${HEADER_FILES})
target_link_libraries("Rasterizer" ${VTK_LIBRARIES})

