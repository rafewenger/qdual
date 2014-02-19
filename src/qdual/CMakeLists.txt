PROJECT(QDUAL)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF (NOT DEFINED ${QDUAL_DIR})
  GET_FILENAME_COMPONENT(QDUAL_ABSOLUTE_PATH "../.." ABSOLUTE)
  SET(QDUAL_DIR ${QDUAL_ABSOLUTE_PATH} CACHE PATH "QDUAL directory")
ENDIF (NOT DEFINED ${QDUAL_DIR})

SET(CMAKE_INSTALL_PREFIX "${QDUAL_DIR}/")
SET(LIBRARY_OUTPUT_PATH ${QDUAL_DIR}/lib CACHE PATH "Library directory")
SET(SRC_QDUAL_DIR "src/qdual")
SET(NRRD_LIBDIR "${QDUAL_DIR}/lib")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

#Find Expat
include(FindEXPAT)
find_package(EXPAT REQUIRED)
include_directories(${EXPAT_INCLUDE_DIRS})

find_library (ITKZLIB_LIBRARY ITKZLIB PATHS "${QDUAL_DIR}/lib")
find_library (ZLIB_FOUND ZLIB PATHS "${QDUAL_DIR}/lib")

#Find itkzlib or z library
if ("${ITKZLIB_LIBRARY}" STREQUAL "ITKZLIB_LIBRARY-NOTFOUND")
  if ("${ZLIB_LIBRARY}" STREQUAL "ZLIB_LIBRARY-NOTFOUND")
    message ("ITKZLIB and ZLIB not found.  Exiting.")
    return(10)
  else()
  message ("ZLIB FOUND")
    set (LIB_ZLIB "z")
  endif()
else()
  message ("ITKZLIB FOUND")
  set (LIB_ZLIB "ITKZLIB")
endif()

INCLUDE_DIRECTORIES("${QDUAL_DIR}/include")
LINK_DIRECTORIES("${NRRD_LIBDIR}")
#LINK_LIBRARIES(expat NrrdIO z)

ADD_EXECUTABLE(qdual qdual_main.cxx qdualIO.cxx qdual.cxx 
                     qdual_datastruct.cxx qdual_extract.cxx 
                     qdual_position.cxx 
                     ijkdualtable.cxx ijkdualtable_ambig.cxx
					 qdual_remove_degen.cxx qdual_table.cxx
					 qdual_collapse.cxx qdual_restrictions.cxx)
target_link_libraries(qdual ${EXPAT_LIBRARIES} NrrdIO ${LIB_ZLIB})
SET(CMAKE_INSTALL_PREFIX ${QDUAL_DIR})
INSTALL(TARGETS qdual DESTINATION "bin")

ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY ../.. COMMAND tar cvfh ${SRC_QDUAL_DIR}/qdual.tar ${SRC_QDUAL_DIR}/*.cxx ${SRC_QDUAL_DIR}/*.h ${SRC_QDUAL_DIR}/CMakeLists.txt)
