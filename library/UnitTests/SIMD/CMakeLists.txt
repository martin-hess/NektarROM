SET(VecDataTestSources
    main.cpp
    TestSimdLib.cpp
)

ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
ADD_NEKTAR_EXECUTABLE(VecDataUnitTests
    COMPONENT unit-test DEPENDS MatrixFreeOps SOURCES ${VecDataTestSources})
ADD_TEST(NAME VecDataUnitTests COMMAND VecDataUnitTests --detect_memory_leaks=0)
