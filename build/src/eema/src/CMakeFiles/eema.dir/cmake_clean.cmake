file(REMOVE_RECURSE
  "../../../lib/libeema.pdb"
  "../../../lib/libeema.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/eema.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
