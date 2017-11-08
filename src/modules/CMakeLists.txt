set(MODULES
  danilo_plot_combined
  danilo_datacards
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
