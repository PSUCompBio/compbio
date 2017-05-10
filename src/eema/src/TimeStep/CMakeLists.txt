SET(TimeStep_src 
    fe_getTimeStep.cpp
    fe_calWaveSpeed.cpp
    fe_calTimeStep.cpp
     )

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) 
ADD_LIBRARY(TimeStep OBJECT ${TimeStep_src})
