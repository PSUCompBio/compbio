SET(2-elt-hex_src main.cpp)

ADD_EXECUTABLE(2-elt-hex ${2-elt-hex_src})

TARGET_LINK_LIBRARIES(2-elt-hex eema)

file(COPY input.inp DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)
