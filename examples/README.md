## Examples
1. `ex1a_generate_wl_resample.m`: use the *boot_tide.m* function for generating MBB and SPB water level resamples given the water level reconstructions (from tidal harmonic analysis regressions) at two stations. 
2. `ex1b_compute_boot_stats_on_resamples.m`: use the *boot_tide.m* function for generating MBB water level resamples at two stations and compute bootstrap plug-in statistics using *boot_tide_param.m*.
3. `ex3a_compare_boot_ci.m`: compute and compare the bootstrap estimators and CIs of the tidal HA amplitudes and phases using the percentile, gaussian, and studentized bootstrap methods. 


* `auxf`: folder with auxiliary functions used to load data and create figures
* `load_one_station.m` and `load_two_stations.m`: load data (hourly water levels and the results of HA regressions applied on these series)

