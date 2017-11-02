set(SRCS
  io.cpp
  util.cpp
  physics.cpp
  gfx.cpp
  hist.cpp
  efficiency.cpp
  )

add_library(tools SHARED ${SRCS})

