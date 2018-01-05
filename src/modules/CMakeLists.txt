set(MODULES
  danilo_plot_combined
  danilo_datacards
  danilo_distributions
  signal_scan
  distributions
  SFapplication
  test_2d
  danilo_check_scan
  danilo_SFapplication
  danilo_check_overlap
  danilo_plot2d_scan
  danilo_signal_scan
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
