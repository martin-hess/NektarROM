SET(BOOST_INCLUDE_SEARCH_PATH ${CMAKE_SOURCE_DIR}/../ThirdParty/include /usr/local/include /usr/include /sw/include )

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )
    SET(BOOST_INCLUDE_SEARCH_PATH ${BOOST_INCLUDE_SEARCH_PATH}
    "C:\\Program Files\\Microsoft Visual Studio .NET 2003\\Vc7\\include"
    c:\\Boost\\include )
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 7 .NET 2003" )

IF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )
    SET(BOOST_INCLUDE_SEARCH_PATH ${BOOST_INCLUDE_SEARCH_PATH}
        "C:\\Program Files (x86)\\Microsoft Visual Studio 8\\VC\\include"
        "C:\\Program Files\\Microsoft Visual Studio 8\\VC\\include"
        c:\\Boost\\include)
ENDIF( ${CMAKE_GENERATOR} STREQUAL "Visual Studio 8 2005" )

FIND_PATH(BOOST_BASE_DIR boost/weak_ptr.hpp ${BOOST_INCLUDE_SEARCH_PATH}
NO_DEFAULT_PATH )
FIND_PATH(BOOST_BASE_DIR boost/weak_ptr.hpp ${BOOST_INCLUDE_SEARCH_PATH} )
SET(BOOST_INCLUDE_DIR ${BOOST_BASE_DIR})
    
SET(BoostFileSystemName "")
SET(BoostFileSystemDebugName "")

SET(BoostFileSystemName "libboost_filesystem-vc71-mt" "libboost_filesystem-vc80-mt"
			boost_filesystem-gcc boost_filesystem-gcc41 boost_filesystem)

SET(BoostFileSystemDebugName "libboost_filesystem-vc71-mt-gd" "libboost_filesystem-vc80-mt-gd"
					boost_filesystem-gcc-d
					boost_filesystem-d
					boost_filesystem-gcc41-d
					${BoostFileSystemName} )

SET(BoostThreadName "boost_thread-vc71-mt" "boost_thread-vc80-mt"
			boost_thread-gcc-mt boost_thread-gcc41-mt
			boost_thread-gcc boost_thread-gcc41 boost_thread-mt boost_thread)

SET(BoostThreadDebugName "boost_thread-vc71-mt-gd" "boost_thread-vc80-mt-gd"
				boost_thread-gcc-mt-d  boost_thread-gcc41-mt-d 
				boost_thread-mt-d
				boost_thread-d ${BoostThreadName} )

FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES ${BoostFileSystemName}
          PATHS
          ${BOOST_BASE_DIR}/lib
	  ${BOOST_BASE_DIR}/../lib
	  NO_DEFAULT_PATH )
FIND_LIBRARY( BOOST_FILESYSTEM_LIB NAMES ${BoostFileSystemName}
          PATHS
          C:\\Boost\\lib)

FIND_LIBRARY( BOOST_FILESYSTEM_DEBUG_LIB NAMES ${BoostFileSystemDebugName}
          PATHS
          ${BOOST_BASE_DIR}/lib
          ${BOOST_BASE_DIR}/../lib
	  NO_DEFAULT_PATH)
FIND_LIBRARY( BOOST_FILESYSTEM_DEBUG_LIB NAMES ${BoostFileSystemDebugName}
          PATHS
          C:\\Boost\\lib)

FIND_LIBRARY( BOOST_THREAD_LIB NAMES ${BoostThreadName}
          PATHS
          ${BOOST_BASE_DIR}/lib
	    ${BOOST_BASE_DIR}/../lib
	  NO_DEFAULT_PATH)
FIND_LIBRARY( BOOST_THREAD_LIB NAMES ${BoostThreadName}
          PATHS
          C:\\Boost\\lib)
FIND_LIBRARY( BOOST_THREAD_DEBUG_LIB NAMES ${BoostThreadDebugName}
          PATHS
          ${BOOST_BASE_DIR}/lib
          ${BOOST_BASE_DIR}/../lib
	  NO_DEFAULT_PATH)
FIND_LIBRARY( BOOST_THREAD_DEBUG_LIB NAMES ${BoostThreadDebugName}
          PATHS
          C:\\Boost\\lib )
IF (BOOST_INCLUDE_DIR)
  SET(BOOST_FOUND TRUE)
ENDIF (BOOST_INCLUDE_DIR)

#SET (BOOST_LIB_DIR ${BOOST_BASE_DIR}/lib )
GET_FILENAME_COMPONENT(BOOST_LIB_DIR ${BOOST_THREAD_LIB} PATH CACHE)
LINK_DIRECTORIES(${BOOST_LIB_DIR})
    
IF (BOOST_FOUND)
  IF (NOT Boost_FIND_QUIETLY)
     MESSAGE(STATUS "Found Boost: ${BOOST_INCLUDE_DIR}")
  ENDIF (NOT Boost_FIND_QUIETLY)
ELSE(BOOST_FOUND)
  IF (Boost_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Boost")
  ENDIF (Boost_FIND_REQUIRED)
ENDIF (BOOST_FOUND)
