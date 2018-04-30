add_subdirectory(TimeStep)
add_subdirectory(Math)
add_subdirectory(Materials)
add_subdirectory(Input-Output)
add_subdirectory(Elements)
add_subdirectory(FEM)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
add_library(eema $<TARGET_OBJECTS:TimeStep> $<TARGET_OBJECTS:Math> $<TARGET_OBJECTS:Materials> $<TARGET_OBJECTS:Input-Output> $<TARGET_OBJECTS:Elements> $<TARGET_OBJECTS:FEM>)

