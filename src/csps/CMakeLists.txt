find_program(PYTHON "python")


#Bring the headers, such as *.h into the project, should be renamed to include
#include_directories(include)


add_subdirectory(src)

#===============================================================================
# Examples and tests
#===============================================================================
#option(CSPS_EXAMPLES "Activate examples" OFF)
#if(CSPS_EXAMPLES)
#  add_subdirectory(examples)
#endif()

#===============================================================================
# Documentation
#===============================================================================
#option(CSPS_DOCUMENTATION "Build source documentation using Doxygen." OFF)
#if(CSPS_DOCUMENTATION)
#  add_subdirectory(doc)
#endif()

#===============================================================================
# Graphic User Interface
#===============================================================================
#option(CSPS_GUI "Build the CSPS GUI using Qt." ON)
#if(CSPS_GUI)
#  add_subdirectory(gui)
#endif()
