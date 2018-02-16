SET(src main.cpp)

ADD_EXECUTABLE(parallel-test ${src})

TARGET_LINK_LIBRARIES(parallel-test eema ${MPI_LIBRARIES})


