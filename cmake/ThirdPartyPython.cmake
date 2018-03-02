########################################################################
#
# ThirdParty configuration for Nektar++
#
# Python interfaces
#
########################################################################

IF (NEKTAR_BUILD_PYTHON)
    # Find Python
    FIND_PACKAGE(PythonInterp  REQUIRED)
    FIND_PACKAGE(PythonLibsNew REQUIRED)
    INCLUDE_DIRECTORIES(SYSTEM ${PYTHON_INCLUDE_DIRS})

    # Include headers from root directory for config file.

    # Try to find Boost.NumPy
    FIND_LIBRARY(BOOST_NUMPY_LIB boost_numpy PATHS ${Boost_LIBRARY_DIRS})

    # If we can't find it, pull it from git and compile it
    IF (NOT BOOST_NUMPY_LIB)
        INCLUDE(ExternalProject)
        EXTERNALPROJECT_ADD(
            boost-numpy
            PREFIX ${TPSRC}
            GIT_REPOSITORY http://github.com/ndarray/Boost.NumPy.git
            STAMP_DIR ${TPBUILD}/stamp
            DOWNLOAD_DIR ${TPSRC}
            SOURCE_DIR ${TPSRC}/boost-numpy
            BINARY_DIR ${TPBUILD}/boost-numpy
            TMP_DIR ${TPBUILD}/boost-numpy-tmp
            INSTALL_DIR ${TPDIST}
            CONFIGURE_COMMAND ${CMAKE_COMMAND}
                -G ${CMAKE_GENERATOR} -DCMAKE_INSTALL_PREFIX:PATH=${TPDIST} -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF -DLIBRARY_TYPE=STATIC
                ${TPSRC}/boost-numpy
            )

        SET(BOOST_NUMPY_LIB ${TPDIST}/lib64/${CMAKE_STATIC_LIBRARY_PREFIX}boost_numpy${CMAKE_STATIC_LIBRARY_SUFFIX})
        INCLUDE_DIRECTORIES(SYSTEM ${TPDIST}/include)
        CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/python/setup.py.in ${CMAKE_BINARY_DIR}/setup.py)

        ADD_CUSTOM_TARGET(nekpy-install-user
            DEPENDS _MultiRegions
            COMMAND python setup.py install --user
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

        FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/NekPy)
        FILE(WRITE ${CMAKE_BINARY_DIR}/NekPy/__init__.py "# placeholder")
    ELSE()
        ADD_CUSTOM_TARGET(boost-numpy ALL)
        ADD_DEFINITIONS(-DBOOST_HAS_NUMPY)
    ENDIF()
ENDIF()
