SET(example2_src main.cpp)

ADD_EXECUTABLE(eema_example2 ${example2_src})

TARGET_LINK_LIBRARIES(eema_example2 eema)

file(COPY input.inp DESTINATION .)
file(COPY input_nofiber_current.inp DESTINATION .)
file(COPY input_case-1-mp-16elements.inp DESTINATION .)
file(COPY case-1-mphysics.inp DESTINATION .)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/vtu)
