cmake_minimum_required (VERSION 3.0)
project (MoopCorr)

### General setup
set(CMAKE_CXX_STANDARD 11)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

#Enable debug build 
# alternative from command-line:
# cmake -DCMAKE_BUILD_TYPE=debug
#set(CMAKE_BUILD_TYPE debug)
#set(ROOT_CONFIG_DEBUG 1)

### Setup dependencies
# ROOT
if(EXISTS $ENV{ROOTSYS}/cmake/ROOTConfig.cmake)
  list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
else()
  list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
endif()
find_package(ROOT REQUIRED COMPONENTS ASImage)

# Silly: simply image library
find_package(SILLY) #custom-writen

# libXML2
find_package(LibXml2 REQUIRED)

# Latex (optional)
#set(LATEX_OUTPUT_PATH doc)
#include(UseLATEX)

### Include folder 
include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(${ROOT_INCLUDE_DIRS})
include_directories(${SILLY_INCLUDE_DIRS})
include_directories(${LIBXML2_INCLUDE_DIR})

#Load additional cmake macros utilities
include(${ROOT_USE_FILE})

### Libararies
ROOT_GENERATE_DICTIONARY(G__MoopMap MoopMap.h LINKDEF MoopMapLinkDef.h)
ROOT_GENERATE_DICTIONARY(G__MoopMapPic MoopMapPic.h LINKDEF MoopMapPicLinkDef.h)
add_library(MoopMap SHARED src/MoopMap.cxx src/MoopMapPic.cxx G__MoopMap.cxx G__MoopMapPic.cxx)
target_link_libraries(MoopMap ${ROOT_LIBRARIES} ${SILLY_LIBRARIES} ${LIBXML2_LIBRARIES})
set_target_properties(MoopMap PROPERTIES COMPILE_FLAGS ${SILLY_CXX_FLAGS})

add_library(MoopAnalyzer SHARED src/DataDefs.cxx src/MoopDataAnalyzer.cxx)

install (TARGETS MoopMap MoopAnalyzer DESTINATION lib)

### Executables
add_executable(generateToys util/generateToys.cxx)
target_link_libraries(generateToys MoopMap MoopAnalyzer)

ROOT_GENERATE_DICTIONARY(G__checkMap checkMap.h LINKDEF checkMapLinkDef.h)
add_executable(checkMap util/checkMap.cxx G__checkMap.cxx)
target_link_libraries(checkMap MoopMap)

add_executable(analyzeToyData util/analyzeToyData.cxx)
target_link_libraries(analyzeToyData MoopMap MoopAnalyzer)

install (TARGETS generateToys checkMap DESTINATION bin)

### Documentation
#add_latex_document(
#    ${CMAKE_SOURCE_DIR}/doc/stat.tex
#    IMAGE_DIRS ${CMAKE_SOURCE_DIR}/doc/images
#    USE_INDEX
#    FILTER_OUTPUT
#    )


