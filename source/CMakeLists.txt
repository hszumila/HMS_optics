#----------------------------------------------------------------------------
# Setup the project.
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(HMS_optics)

#----------------------------------------------------------------------------
# Setup ROOT.
find_package(ROOT REQUIRED COMPONENTS Spectrum)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#----------------------------------------------------------------------------
# Setup include directories.
include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Locate sources and headers for this project.
set(sources
  ${PROJECT_SOURCE_DIR}/src/cmdOptions.cpp
  ${PROJECT_SOURCE_DIR}/src/myConfig.cpp
  ${PROJECT_SOURCE_DIR}/src/myEvent.cpp
  ${PROJECT_SOURCE_DIR}/src/myMath.cpp
  ${PROJECT_SOURCE_DIR}/src/myOther.cpp
  ${PROJECT_SOURCE_DIR}/src/myRecMatrix.cpp
)
set(headers
  ${PROJECT_SOURCE_DIR}/inc/cmdOptions.hpp
  ${PROJECT_SOURCE_DIR}/inc/myConfig.hpp
  ${PROJECT_SOURCE_DIR}/inc/myEvent.hpp
  ${PROJECT_SOURCE_DIR}/inc/myMath.hpp
  ${PROJECT_SOURCE_DIR}/inc/myOther.hpp
  ${PROJECT_SOURCE_DIR}/inc/myRecMatrix.hpp
)

#----------------------------------------------------------------------------
# Copy all data files to build directory.
set(DATA_FILES
  hms_recon_coeff__dep.dat
  hms_recon_coeff__indep.dat
  setup_optics_example.txt
  hmsoptics/hms_replay_production_1235_50000.root
)

foreach(_file ${DATA_FILES})
  configure_file(
    ${PROJECT_SOURCE_DIR}/data/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Set compile flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra -Weffc++ -Wconversion -Wsign-conversion -Wsign-promo")

#----------------------------------------------------------------------------
# Add the executable, and link it.
add_executable(reconstruct reconstruct.cpp ${sources})
target_link_libraries(reconstruct ${ROOT_LIBRARIES} ${ROOT_TSpectrum_LIBRARY})

add_executable(hms_optics hms_optics.cpp ${sources})
target_link_libraries(hms_optics ${ROOT_LIBRARIES} ${ROOT_TSpectrum_LIBRARY})
