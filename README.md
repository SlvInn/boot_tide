# boot_tide: MATLAB scripts for tide model resampling-based uncertainty analysis 
Set of functions implementing two residual boostrap algorithms for tidal analyses.  
1. a Moving Block Bootstrap (MBB) that resamples fixed or random length blocks of the tidal model resisuals, and 
2. a Semi-Parametric Bootstrap (SPB) that simulates noise realizations from the Fast Fourier Transform of the residuals. 

Both methods are intended to produce similar results with comparable accuracy. The MBB can be applied in virtually all sampling conditions but it is more computationally demanding than SPB, especially for long records. MBB can also handle the construction of the residual resamples at various spatial locations (e.g., stations) while accounting for (part of) the spatial autocorrelation beetween them. Similarly, it could be in principle used to work with bi-dimentional data (current series, treating 2D series as 2-location input matrices) but it hasn't been tested in this use. The SPB is faster than the MBB but constructs independent residual resamples at the various spatial locations. Also, some sort of interpolation is needed for resampling via SPB irregularly sampled records or series with missing data, since it involves a Fourier analysis of residuals.  

For a discussion of the MBB and SPB residual resampling, including their comparison with other existing tidal analysis packages: Innocenti et al. 2022, DOI:10.1175/JTECH-D-21-0060.1

The project use `fttnoise.m` by Aslak Grinsted (2009; https://www.mathworks.com/matlabcentral/fileexchange/32111-fftnoise-generate-noise-with-a-specified-power-spectrum) and other contributed functions.

---- 
silvia.innocenti@ec.gc.ca, August 2022.

Innocenti, S. and Matte, P. and Fortin, V. and Bernier, N. B., (2022), Residual bootstrap methods for parameter uncertainty assessment in tidal analysis with temporally correlated noise, JTECH, DOI:10.1175/JTECH-D-21-0060.1

#### Contents
1. [boot_tide.m:](#boot_tide.m) construct the residual bootstrap resamples/replicates.  
2. [boot_tide_param.m:](#boot_tide_param.m) compute the bootstrap plug-in estimators of the parameters and corresponiding confidence intervals.
3. examples and templates folder


### boot_tide.m
Construct residual bootstrap resamples given an or some observed series and 
a tidal reconstruction (e.g., tidal HA regression predictions) or the corresponding
residuals (i.e., observed series-reconstruction).

This function implements the MBB and SPB algorithms (`method` option) to generate the resamples of the water level series and, if requested, estimate the parameters of a user-specifed model (`tidal_model` option) for each of the generated resamples. 
   
```MATLAB
    % general syntax:
    out = boot_tide(time, y, {'option_name',option_value})
```

### boot_tide_param.m
Calculate the plug-in bootstrap estimators of the parameter replicates of a tidal model (not necessarily produced via `boot_tide.m`, it only needs a matrix of parameters replicates) and the corresponding Confidence Intervals (CIs). 
Via the specification of the `circular` option, the user can request to use circular statistics for the parameters expressed in radiant or degrees units [e.g., tidal phases]. 


```MATLAB
    % general syntax:
    [m_theta]  = boot_tide_param(theta, {'option_name',option_value}) 
    [m_theta,s_theta]  = boot_tide_param(theta, {'option_name',option_value})
    [m_theta,s_theta,ci_theta]  = boot_tide_param(theta, {'option_name',option_value})
    [m_theta,s_theta,ci_theta,opt]  = boot_tide_param(theta, {'option_name',option_value})
```

