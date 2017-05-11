project("CSPS")

macro(add_python_target tgt)
  foreach(file ${ARGN})
    set(OUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.pyo)
    list(APPEND OUT_FILES ${OUT})
    add_custom_command(OUTPUT ${OUT}
        COMMAND <python command you use to byte-compile .py file>)
  endforeach()

  add_custom_target(${tgt} ALL DEPENDS ${OUT_FILES})
endmacro()

add_python_target(csps main.py)
#SET(csps_src 
#    main.py  
#   )
#ADD_LIBRARY(csps ${csps_src})
