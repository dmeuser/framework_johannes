root_generate_dictionary(TreeParticles_dict tree/TreeParticles.hpp LINKDEF tree/LinkDef.h)

add_library(tree SHARED TreeParticles.cpp TreeParticles_dict)
add_custom_command(TARGET tree POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/src/tree/TreeParticles_dict_rdict.pcm ${LIBRARY_OUTPUT_PATH}
  )

