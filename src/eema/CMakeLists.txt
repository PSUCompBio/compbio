#Bring the headers, such as *.h into the project, should be renamed to include
include_directories(include)
add_subdirectory(src)

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
