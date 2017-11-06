set(MODULES
  weighting
  fit_weighting
  SR_variance
  SR_JES
  JES_shift
  fit_JES
  JES_studies
  tests
  distributions
  trigger
  prescales
  matching
  genComp
  efake
  8TeVXcheck
  eventPicker
  zinvisible
  roofit
  SFapplication
  SRbinning
  signal_scan
  limits
  correlation_test
  pu
  acceptance
  mets_gain
  xsec_plots
  wz_ratio
  80XXcheck
  plot_limit
  plot_uncert
  cutflow
  )

foreach(MODULE ${MODULES})
  add_library(${MODULE} SHARED ${MODULE}.cpp)
endforeach()
