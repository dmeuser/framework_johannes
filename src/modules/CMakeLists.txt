set(MODULES
  danilo_plot_combined
  danilo_datacards
  danilo_distributions
  danilo_check_scan
  danilo_SFapplication
  danilo_check_overlap
  danilo_plot2d_scan
  danilo_signal_scan
  danilo_cutflow_GGM
  danilo_acceptanceHist
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
