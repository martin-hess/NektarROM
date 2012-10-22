SET(NEKTAR_USE_MPI    OFF CACHE BOOL 
    "Use MPICH2 for parallelisation.")
SET(THIRDPARTY_BUILD_GSMPI OFF CACHE BOOL
    "Build GSMPI if needed")

IF( NEKTAR_USE_MPI )
    INCLUDE (FindMPI)
    MARK_AS_ADVANCED(MPI_LIBRARY)
    MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY)
    MARK_AS_ADVANCED(file_cmd)
    ADD_DEFINITIONS(-DNEKTAR_USE_MPI)
    INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )
    
    IF (THIRDPARTY_BUILD_GSMPI)
        EXTERNALPROJECT_ADD(
            gsmpi-1.1
            PREFIX ${TPSRC}
            URL ${TPURL}/gsmpi-1.1.tar.bz2
            URL_MD5 "f2c1f7695f361c6d87365e2ea63aece1"
            DOWNLOAD_DIR ${TPSRC}
            CONFIGURE_COMMAND 
                ${CMAKE_COMMAND} 
                -DCMAKE_INSTALL_PREFIX:PATH=${TPSRC}/dist 
                ${TPSRC}/src/gsmpi-1.1
        )
        SET(GSMPI_LIBRARY gsmpi CACHE FILEPATH
            "GSMPI path" FORCE)
        SET(XXT_LIBRARY xxt CACHE FILEPATH
            "XXT path" FORCE)
    ELSE (THIRDPARTY_BUILD_GSMPI)
        INCLUDE (FindGSMPI)
        INCLUDE (FindXXT)
    ENDIF (THIRDPARTY_BUILD_GSMPI)

ENDIF( NEKTAR_USE_MPI )

