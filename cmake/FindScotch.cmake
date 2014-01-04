
SET(SCOTCH_SEARCH_PATHS 
	${CMAKE_SOURCE_DIR}/ThirdParty/scotch_6.0.0/
	${CMAKE_SOURCE_DIR}/ThirdParty/scotch_6.0.0/build/
	${CMAKE_SOURCE_DIR}/../ThirdParty/scotch_6.0.0/
	${CMAKE_SOURCE_DIR}/../ThirdParty/scotch_6.0.0/build
    ${CMAKE_SOURCE_DIR}/ThirdParty/dist/lib 
    ${CMAKE_SOURCE_DIR}/../ThirdParty/dist/lib)

FIND_LIBRARY(SCOTCH_LIB NAMES scotch PATHS ${SCOTCH_SEARCH_PATHS})
FIND_LIBRARY(SCOTCHMETIS_LIB NAMES scotchmetis PATHS ${SCOTCH_SEARCH_PATHS})
FIND_LIBRARY(SCOTCHERR_LIB NAMES scotcherr PATHS ${SCOTCH_SEARCH_PATHS})

SET(SCOTCH_FOUND FALSE)
IF (SCOTCH_LIB AND SCOTCHMETIS_LIB AND SCOTCHERR_LIB)
    SET(SCOTCH_FOUND TRUE)
    GET_FILENAME_COMPONENT(SCOTCH_PATH ${SCOTCH_LIB} PATH)
    INCLUDE_DIRECTORIES(${SCOTCH_PATH}/../../include)
    MARK_AS_ADVANCED(SCOTCH_PATH)
    MARK_AS_ADVANCED(SCOTCH_LIB)
    MARK_AS_ADVANCED(SCOTCHERR_LIB)
    MARK_AS_ADVANCED(SCOTCHMETIS_LIB)
ENDIF (SCOTCH_LIB AND SCOTCHMETIS_LIB AND SCOTCHERR_LIB)

IF (SCOTCH_FOUND)
    IF (NOT SCOTCH_LIB_FIND_QUIETLY)
        MESSAGE(STATUS "Found Scotch")
    ENDIF (NOT SCOTCH_LIB_FIND_QUIETLY)
ELSE(SCOTCH_FOUND)
    IF (SCOTCH_LIB_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find Scotch")
    ENDIF (SCOTCH_LIB_FIND_REQUIRED)
ENDIF (SCOTCH_FOUND)
