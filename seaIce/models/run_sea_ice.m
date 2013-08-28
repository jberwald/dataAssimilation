function zf = run_sea_ice( E0, delta_t, t0 )
  % function zf = run_sea_ice( E0, delta_t, t0 )
  %
  % (Opaque) Wrapper to run sea_ice_EW09.m.
  %
  % E0 - initial value for model. If 0, then sea_ice_EW09 will run in
  % spin_up mode for delta_t time 'years'. Else, E0 should be the last
  % value of E in z_forecast, where
  %
  % z_forecast = [t, E, T_srf, -E/Li, E/Cml*Hml, F_top, F0, F_T, F_sw]
  %
  % Thus, E0 = z_forecast(end,2);, which for efficiency should be
  % extracted from the full z_forcast matrix before passing it to
  % run_sea_ice.m
  %
  % delta_t - time (in years) between observational data points
  % 
  % t0 - end time of previous forecast. If t0==0, sea_ice_EW09.m
  % assumes a spin_up run.
  %
  % author: Jesse Berwald, August 2013
  z_forecast = sea_ice_EW09( E0, delta_t ,t0 );
  zf = z_forecast;
     
