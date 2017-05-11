option (EXAMPLE1  "Build library with Explicit Dynamic Example 1" OFF)
if (EXAMPLE1)
   add_subdirectory(example-1)
endif (EXAMPLE1)
option (EXAMPLE2  "Multiphysics Simulation" OFF)
if (EXAMPLE2)
   add_subdirectory(example-2)
endif (EXAMPLE2)


