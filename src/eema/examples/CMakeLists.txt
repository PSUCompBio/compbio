option (EXAMPLE1  "Build library with Explicit Dynamic Example 1" OFF)
if (EXAMPLE1)
   add_subdirectory(example-1)
endif (EXAMPLE1)

option (2-ELT-HEX-EXAMPLE  "Build library with Explicit Dynamic 2 Element Hex Example" OFF)
if (2-ELT-HEX-EXAMPLE)
   add_subdirectory(2-elt-hex)
endif (2-ELT-HEX-EXAMPLE)

option (EXAMPLE2  "Multiphysics Simulation" OFF)
if (EXAMPLE2)
   add_subdirectory(example-2)
endif (EXAMPLE2)


