cmake_minimum_required (VERSION 2.6)

# find python
find_package(PythonInterp REQUIRED)

find_package (PythonLibs REQUIRED)

if(NOT PYTHONINTERP_FOUND)
    message(FATAL_ERROR "ERROR: Python interpreter not found. Cannot continue - sorry.")
endif()

include_directories ( ${PYTHON_INCLUDE_DIRS} )

add_subdirectory(src)
