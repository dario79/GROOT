#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GROOT)
set(CMAKE_CXX_FLAGS "-Wall -std=c++17")
set(CMAKE_CXX_VERSION 17)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# CADMesh
# -Dcadmesh_DIR=/home/lattuadad/geom_import/cadmesh/lib/cmake/cadmesh-./
#find_package(cadmesh)
#SET( cadmesh_INCLUDE_DIR "${cadmesh_DIR}/../../../include")
#SET( cadmesh_INCLUDE_DIR "/opt/cadmesh/cadmesh/include")
#
#INCLUDE_DIRECTORIES(${cadmesh_INCLUDE_DIR})
#FIND_LIBRARY(cadmesh_LIBRARY cadmesh)
#SET( cadmesh_libs ${cadmesh_libs} ${cadmesh_LIBRARY})


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
#message(${CMAKE_MODULE_PATH})
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /home/lattuadad/geant4.10.04.p02/cmake/Modules/)
find_package(ROOT REQUIRED)
SET( ROOT_INCLUDE_DIR "/usr/local/include/root")
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

message( STATUS ${ROOT_INCLUDE_DIR})
#message(${ROOT_LIBRARIES})


include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
#                    ${CADMESH_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(GROOT GROOT.cc ${sources} ${headers})
target_link_libraries(GROOT ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
	
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Hadr03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(GROOT_SCRIPTS
    vis.mac run.mac inputFileDetectors.txt
  )

foreach(_script ${GROOT_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS GROOT DESTINATION bin)

