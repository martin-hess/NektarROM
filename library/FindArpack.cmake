
SET(ARPACK_SEARCH_PATHS ${CMAKE_SOURCE_DIR} /opt/local/lib/  /../ThirdParty/ARPACK)

FIND_LIBRARY(ARPACK_LIB NAMES arpack.1 arpack PATHS ${ARPACK_SEARCH_PATHS})

SET(ARPACK_FOUND FALSE)
IF (ARPACK_LIB)
  SET(ARPACK_FOUND TRUE)
    MARK_AS_ADVANCED(ARPACK_LIB)
ENDIF (ARPACK_LIB)

IF (ARPACK_FOUND)
  IF (NOT ARPACK_LIB_FIND_QUIETLY)
     MESSAGE(STATUS "Found Arpack : ${ARPACK_LIB}")
     GET_FILENAME_COMPONENT(ARPACK_PATH ${ARPACK_LIB} PATH CACHE)
     SET(ARPACK_INCLUDE_DIR ${ARPACK_PATH}/../include CACHE FILEPATH "ARPACK include directory.")
  ENDIF (NOT ARPACK_LIB_FIND_QUIETLY)
ELSE(ARPACK_FOUND)
  IF (ARPACK_LIB_FIND_REQUIRED)
     MESSAGE(FATAL_ERROR "Could not find Arpack")
  ENDIF (ARPACK_LIB_FIND_REQUIRED)
ENDIF (ARPACK_FOUND)



