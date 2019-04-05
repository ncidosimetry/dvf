cmake_minimum_required (VERSION 2.8.11)
project (dvf)

find_package (GDCM REQUIRED)
include (${GDCM_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -std=c++11 --stdlib=libc++" )

file(GLOB gdcmlibs "${GDCM_LIBRARY_DIRS}/*.so*" "${GDCM_LIBRARY_DIRS}/*.a" "${GDCM_LIBRARY_DIRS}/*.dylib")
file(GLOB excludegdcmlibs "${GDCM_LIBRARY_DIRS}/libgdcmMEXD.*" "${GDCM_LIBRARY_DIRS}/libsocketxx.*")

if(APPLE)
	find_library(COREFOUNDATION_LIBRARY CoreFoundation )
endif()


#include_directories (
#    ${PROJECT_SOURCE_DIR}/../../src/base
#    ${PROJECT_SOURCE_DIR}/../../src/4ddc
#)

add_executable(dvf dvf_smc.cpp)

install(PROGRAMS dvf DESTINATION .)

link_directories (${GDCM_DIR} )

target_link_libraries(dvf
    ${COREFOUNDATION_LIBRARY}
    gdcmCommon
    gdcmDSED     
    gdcmMEXD     
    gdcmcharls   
    gdcmjpeg12   
    gdcmjpeg8    
    gdcmuuid     
    socketxx
    gdcmDICT     
    gdcmIOD      
    gdcmMSFF     
    gdcmexpat    
    gdcmjpeg16   
    gdcmzlib   
#    gdcmopenjpeg
#    vtkgdcm
)


