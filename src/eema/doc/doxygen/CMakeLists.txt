set(DOXYGEN_INPUT ${CMAKE_CURRENT_BINARY_DIR}/eema.dox)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(DOXYGEN_OUTPUT ${DOXYGEN_OUTPUT_DIR}/index.html)

string(REGEX REPLACE ":" " " EEMA_DOXYGEN_DEFINTIONS "${EEMA_DEFINITIONS}")

make_directory(${DOXYGEN_OUTPUT_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/eema.dox.in
  ${CMAKE_CURRENT_BINARY_DIR}/eema.dox)

add_custom_command(
  OUTPUT ${DOXYGEN_OUTPUT} #${CMAKE_CURRENT_BINARY_DIR}/latex/refman.tex
  COMMAND ${CMAKE_COMMAND} -E echo "${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}"
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building eema Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${DOXYGEN_INPUT}
  )

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.pdf
  COMMAND ${CMAKE_COMMAND} -E echo "Building eema RefMan..."
  COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/latex
  #COMMAND ${CMAKE_COMMAND} -E rename ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.pdf ${CMAKE_CURRENT_BINARY_DIR}/refman.pdf
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/latex/refman.tex
)


add_custom_target(eema-doc ALL DEPENDS ${DOXYGEN_OUTPUT})

add_custom_target(eema-doc-forced
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  )

install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}
  DESTINATION share/eema-${EEMA_VERSION}/doc)
