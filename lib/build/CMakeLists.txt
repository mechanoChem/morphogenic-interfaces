##
#  CMake script
##
#Code structure and usage
##---------------------------------------------------------------------------------------------------------------------
# How to install
# in build folder where you found this CMakeList.txt, change it following guides below
# $cmake CMakeList.txt
# $make release install
# please check "libopenIGA.a" is installed in lib(or other) folder
# delete this build folder if necessary
##---------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------


CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

set(CMAKE_CXX_STANDARD 14)


#Include deal.II
##---------------------------------------------------------------------------------------------------------------------
#usually DEAL_II_DIR for explicitly deal.II path; deal.II_DIR and DEAL_II_DIR could be preloaded enviromental variable

#SET(ENV{DEAL_II_DIR} "/nfs/mcfs_comp/home/rudraa/software/dealii-8.4.1/installWithTrilinos")
FIND_PACKAGE(deal.II 8.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()
##---------------------------------------------------------------------------------------------------------------------




##Include all source files, add library dealMultiPhysics
##---------------------------------------------------------------------------------------------------------------------
##Usually you do NOT need to change anything here
# Set the name of the project and target

SET(TARGET "main")
#this lib is NOT installed as a shared lib
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)

FILE(GLOB mechanoChemFEM_SRC "../src/FEMdata/*.cc" "../src/initBoundValProbs/*.cc" "../src/hpFEMbase/*.cc" "../src/Residual/*.cc" "../src/solve/*.cc" "../src/supplementary/*.cc" )
ADD_LIBRARY(mechanoChemFEM ${mechanoChemFEM_SRC})
DEAL_II_SETUP_TARGET(mechanoChemFEM RELEASE)
FILE (GLOB TARGET_SRC "main.cc")

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
set(CMAKE_STATIC_LIBRARY_PREFIX "")




SET (sacado_INC $ENV{TRILINOS_DIR}/include)
SET (sacado_LIB $ENV{TRILINOS_DIR}/lib/libkokkoscore.so $ENV{TRILINOS_DIR}/lib/libsacado.so)
INCLUDE_DIRECTORIES(${sacado_INC})
#TARGET_LINK_LIBRARIES(main mechanoChemFEM)
##---------------------------------------------------------------------------------------------------------------------


##Set dirdctory for openIGA, please give a full path and make sure the directory is valid 
##---------------------------------------------------------------------------------------------------------------------
install(TARGETS mechanoChemFEM DESTINATION $ENV{mechanoChemFEM_DIR}/v0.71/lib)
FILE(COPY ../include DESTINATION $ENV{mechanoChemFEM_DIR}/v0.71/)

