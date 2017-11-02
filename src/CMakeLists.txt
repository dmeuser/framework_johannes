include_directories(${CMAKE_SOURCE_DIR}/src)

set(SRCS
  main.cpp
  Dataset.cpp
  Config.cpp
  )

set(MYLIBS
  tools
  tree
  style
  )
set(EXTLIBS
  ${ROOT_LIBRARIES}
  -lRooFit
  ${Boost_LIBRARIES}
  ${CMAKE_DL_LIBS}
  )

add_executable(run.x ${SRCS})
target_link_libraries(run.x ${MYLIBS} ${EXTLIBS})

add_subdirectory(tools)
add_subdirectory(modules)
add_subdirectory(tree)
add_subdirectory(style)
