cmake_minimum_required(VERSION 3.9)

# ===============================================
# Project name

# Set project name if not build within another project
if("${PROJECT_NAME}" STREQUAL "")
    project(HomologyInference)
endif()

# ===============================================
# Global settings

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Disable some compiler warnings
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-sign-compare")

# Turn some compiler warnings into errors
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=return-type")

# Export symbols of all functions so we get useful stack traces
if (CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_ENABLE_EXPORTS ON)
endif()

# Set glow bin dir
set(GLOW_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin/glow)

# Set output dir for meshes, screenshots, etc...
set(ISM_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/output)

# Disable glfw additionals
option(GLFW_BUILD_EXAMPLES "" OFF)
option(GLFW_BUILD_TESTS "" OFF)
option(GLFW_BUILD_DOCS "" OFF)
option(GLFW_INSTALL "" OFF)

# ===============================================
# Find packages

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

find_package(Gurobi REQUIRED)
find_package(Eigen3 REQUIRED)

# ===============================================
# Add submodules

add_subdirectory(libs/OpenMesh)

# Supress an OpenMesh warning
if (${QT_QMAKE_EXECUTABLE})
endif()

#add_subdirectory(libs/eigen-lean)
add_subdirectory(libs/typed-geometry)
add_subdirectory(libs/polymesh)
add_subdirectory(libs/glfw)
add_subdirectory(libs/glow)
add_subdirectory(libs/imgui-lean)
add_subdirectory(libs/glow-extras)

# ===============================================
# Target InterSurfaceMaps

file(GLOB_RECURSE ISM_HEADERS src/*.hh src/*.h)
file(GLOB_RECURSE ISM_SOURCES src/*.cc src/*.c)
add_library(HomologyInference ${ISM_SOURCES} ${ISM_HEADERS})
target_compile_features(HomologyInference PUBLIC cxx_std_17)
target_compile_definitions(HomologyInference PUBLIC SOURCE_PATH_STR="${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(HomologyInference PUBLIC DATA_PATH_STR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_compile_definitions(HomologyInference PUBLIC OUTPUT_PATH_STR="${ISM_OUTPUT_DIR}")
target_include_directories(HomologyInference PUBLIC
    ${OPENMESH_INCLUDE_DIRS}
    src
)
target_link_libraries(HomologyInference PUBLIC
    Eigen3::Eigen
    OpenMeshCore
    typed-geometry
    polymesh
    glfw
    glow
    imgui
    glow-extras
    stdc++fs
    gurobi_c++
)

# ===============================================
# Target apps (Executables)

file(GLOB_RECURSE ISM_APPS apps/*.cc)
foreach(main ${ISM_APPS})
    get_filename_component(name ${main} NAME_WE)
    get_filename_component(dir ${main} DIRECTORY)
    get_filename_component(dir_name ${dir} NAME)
    message(STATUS "Adding app ${name}")

    add_executable(${name} ${main})
    target_compile_definitions(${name} PRIVATE APP_EXECUTABLES)
    add_dependencies(${name} HomologyInference)
    target_include_directories(${name} PUBLIC
        "src"
    )
    target_link_libraries(${name} PUBLIC
        HomologyInference
    )
endforeach()
