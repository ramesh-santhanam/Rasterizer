include_directories("${PROJECT_SOURCE_DIR}/src/Rasterizer")
include_directories(${VTK_INCLUDE_DIRS})

add_executable(testRasterizer main.cpp)
add_executable(testEval testEval.cpp)

set( PROJECT_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src/Rasterizer/Rasterize2d.h")
message("... include: ${PROJECT_INCLUDE_DIRS}..")

set( PROJECT_LINK_LIBS Rasterizer )
target_include_directories( testRasterizer PUBLIC ${PROJECT_INCLUDE_DIRS} )
target_link_libraries( testRasterizer ${PROJECT_LINK_LIBS} )

set( PROJECT_LINK_LIBS Rasterizer gtest)
target_include_directories( testEval PUBLIC ${PROJECT_INCLUDE_DIRS} )
target_link_libraries( testEval ${PROJECT_LINK_LIBS} )
