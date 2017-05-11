SET(example2_src main.cpp)

ADD_EXECUTABLE(eema_example2 ${example2_src})

TARGET_LINK_LIBRARIES(eema_example2 eema)

file(COPY input.inp DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)
