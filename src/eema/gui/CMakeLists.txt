# configure Qt

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Core Widgets REQUIRED)

# Testing configuration
enable_testing()
set(TEST_LINK_LIBRARIES Qt5::Test)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defines.h.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/defines.h)

file(GLOB_RECURSE UI_FILES *.ui)
file(GLOB_RECURSE CODE_FILES *.cpp)


qt5_wrap_ui(UI_HEADERS ${UI_FILES})
qt5_add_resources(RESOURCE_FILES ${CMAKE_SOURCE_DIR}/src/eema/gui/resources.qrc)

add_executable(eema_gui
  ${UI_HEADERS}
  ${CODE_FILES}
  ${RESOURCE_FILES}
)

target_link_libraries(eema_gui
  Qt5::Widgets
)
