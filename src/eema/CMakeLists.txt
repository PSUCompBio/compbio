project("EEMA")
#Bring the headers, such as *.h into the project, should be renamed to include
include_directories(include)
add_subdirectory(src)

#===============================================================================
# Add Parallel Capabilties
#===============================================================================
option(EEMA_PARALLEL "Activate parallelization" OFF)
if(EEMA_PARALLEL)
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
  include_directories(${MPI_INCLUDE_PATH})
endif()

#===============================================================================
# Examples and tests
#===============================================================================
option(EEMA_EXAMPLES "Activate examples" OFF)
if(EEMA_EXAMPLES)
  add_subdirectory(examples)
endif()

#===============================================================================
# Documentation
#===============================================================================
option(EEMA_DOCUMENTATION "Build source documentation using Doxygen." OFF)
if(EEMA_DOCUMENTATION)
  add_subdirectory(doc)
endif()

#===============================================================================
# Graphic User Interface
#===============================================================================
option(EEMA_GUI "Build the EEMA GUI using Qt." OFF)
if(EEMA_GUI)
  add_subdirectory(gui)
endif()
