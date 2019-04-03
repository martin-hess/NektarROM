# FindScotch
# ----------
# Finds the graph partitioning library Scotch.
#
# This module defines the following variables
#   SCOTCH_INCLUDE_DIRS - Location of scotch.h
#   SCOTCH_LIBRARY - Location of main scotch library
#   SCOTCHERR_LIBRARY - Location of scotcherr library
#   SCOTCH_FOUND - TRUE if Scotch found
#   SCOTCH_VERSION - Version of scotch library found

# Determine if we are looking for PT-scotch, or just scotch
SET(SCOTCH_COMPONENTS SERIAL PARALLEL)
SET(PARALLEL OFF)
LIST(FIND Scotch_FIND_COMPONENTS "scotch" FIND_SERIAL)
LIST(FIND Scotch_FIND_COMPONENTS "ptscotch" FIND_PARALLEL)
IF (FIND_PARALLEL GREATER -1)
    SET(PARALLEL ON)
ENDIF()

MESSAGE(STATUS "Searching for Scotch:")
SET(TEST_SCOTCH_DIR $ENV{SCOTCH_DIR})
SET(TEST_SCOTCH_INCLUDE_DIR $ENV{SCOTCH_INCDIR})

SET(SCOTCH_HEADER_DIRS "SCOTCH_HEADERS_DIR-NOTFOUND")
IF(TEST_SCOTCH_INCLUDE_DIR)
    FIND_PATH(SCOTCH_HEADER_DIRS NAMES scotch.h 
                HINTS ${TEST_SCOTCH_INCLUDE_DIR})
ELSE()
    IF(TEST_SCOTCH_DIR)
        FIND_PATH(SCOTCH_HEADERS_DIRS NAMES scotch.h
                    HINTS ${TEST_SCOTCH_DIR}
                    PATH_SUFFIXES "include" "include/scotch")
    ELSE()
        FIND_PATH(SCOTCH_HEADERS_DIRS NAMES scotch.h
                    HINTS ${MACPORTS_PREFIX}/include
                    PATH_SUFFIXES "scotch")
    ENDIF()
ENDIF()
MARK_AS_ADVANCED(SCOTCH_HEADERS_DIRS)

IF (SCOTCH_HEADERS_DIRS)
    SET(SCOTCH_INCLUDE_DIR ${SCOTCH_HEADERS_DIRS})

    TRY_RUN(
        RUN_RESULT COMPILE_RESULT
        ${CMAKE_CURRENT_BINARY_DIR}/
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts/get-scotch-version.c
        CMAKE_FLAGS -DINCLUDE_DIRECTORIES=${SCOTCH_INCLUDE_DIR}
        COMPILE_OUTPUT_VARIABLE COMPILER_OUTPUT
        RUN_OUTPUT_VARIABLE SCOTCH_VERSION)

    IF("${COMPILE_RESULT}" AND ("${RUN_RESULT}" EQUAL 0))
        STRING(STRIP "${SCOTCH_VERSION}" SCOTCH_VERSION)
    ELSE()
        MESSAGE(STATUS "${COMPILER_OUTPUT}")
        MESSAGE(STATUS "${SCOTCH_VERSION}")
        MESSAGE(ERROR "ERROR checking for Scotch version.")
    ENDIF()
ELSE ()
    SET(SCOTCH_INCLUDE_DIR "SCOTCH_INCLUDE_DIRS-NOTFOUND")
    IF (NOT SCOTCH_FIND_QUIETLY)
        MESSAGE(STATUS "Looking for Scotch -- scotch.h not found")
    ENDIF ()
ENDIF ()
LIST(REMOVE_DUPLICATES SCOTCH_INCLUDE_DIR)

FIND_LIBRARY(SCOTCH_LIBRARY    NAMES scotch PATHS ${MACPORTS_PREFIX}/lib)
FIND_LIBRARY(SCOTCHERR_LIBRARY NAMES scotcherr PATHS ${MACPORTS_PREFIX}/lib)
GET_FILENAME_COMPONENT(SCOTCH_LIBRARY_DIR ${SCOTCH_LIBRARY} PATH)

IF (SCOTCH_LIBRARY AND SCOTCHERR_LIBRARY AND SCOTCH_INCLUDE_DIR)
    SET(Scotch_scotch_FOUND TRUE)

    IF (PARALLEL)
        FIND_LIBRARY(PTSCOTCH_LIBRARY NAMES ptscotch PATHS ${SCOTCH_LIBRARY_DIR}
            NO_DEFAULT_PATH)
        FIND_LIBRARY(PTSCOTCHERR_LIBRARY NAMES ptscotcherr PATHS
            ${SCOTCH_LIBRARY_DIR}
            NO_DEFAULT_PATH)

        IF (PTSCOTCH_LIBRARY AND PTSCOTCHERR_LIBRARY)
            SET(Scotch_ptscotch_FOUND TRUE)
        ENDIF()
    ENDIF()
ENDIF()

SET(SCOTCH_FOUND TRUE)
IF (Scotch_scotch_FOUND)
    MESSAGE(STATUS "-- Found Scotch version ${SCOTCH_VERSION}: ${SCOTCH_LIBRARY}")
    IF (PARALLEL)
        IF (Scotch_ptscotch_FOUND)
            MESSAGE(STATUS "-- Found PT-Scotch version ${SCOTCH_VERSION}:${PTSCOTCH_LIBRARY}")
        ELSE()
            SET(SCOTCH_FOUND FALSE)
            MESSAGE(STATUS "-- Could not find PT-Scotch library")
        ENDIF()
    ENDIF()
    IF (Scotch_FIND_VERSION)
        IF (${Scotch_FIND_VERSION} VERSION_GREATER ${SCOTCH_VERSION})
            SET(SCOTCH_FOUND FALSE)
            MESSAGE(STATUS "-- This scotch is too old")
        ENDIF()
    ENDIF()
ELSE()
    SET(SCOTCH_FOUND FALSE)
    MESSAGE(STATUS "-- Could not find Scotch library")
ENDIF()
