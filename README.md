## boot_tide: tide model uncertainty analysis based on residual resampling
Source code implementing the residual bootstrap algorithms presented in 
Innocenti et al. [2022](https://journals.ametsoc.org/view/journals/atot/aop/JTECH-D-21-0060.1/JTECH-D-21-0060.1.xml):

1. a Moving Block Bootstrap (MBB) that resamples fixed or random length blocks of the tidal model resisuals, and 
2. a Semi-Parametric Bootstrap (SPB) that simulates noise realizations from the Fast Fourier Transform of the residuals. 

Both methods are intended to produce similar results with comparable accuracy. The MBB can be applied in virtually all sampling conditions, but it is more computationally demanding than SPB, especially for long records. MBB can also handle the construction of the residual resamples at various spatial locations (e.g., stations) while accounting for (part of) the spatial autocorrelation between them. Similarly, it could be in principle used to work with bi-dimensional data (current series, treating 2D series as 2-location input matrices) but it hasn't been tested in this use. The SPB is faster than some MBB methods, but constructs independent residual resamples at various spatial locations. Also, some interpolation is needed for resampling via SPB irregularly sampled records or series with missing data since it involves a Fourier analysis of residuals.  

For a discussion of the MBB and SPB residual resampling, including their comparison with other existing tidal analysis packages, see Innocenti et al. [2022](DOI:10.1175/JTECH-D-21-0060.1)

The project uses (adapted versions) of contributed MATLAB functions, including: 
- `fttnoise.m` by Aslak Grinsted ([2009](https://www.mathworks.com/matlabcentral/fileexchange/32111-fftnoise-generate-noise-with-a-specified-power-spectrum)), and 
- `CircStat` by Berens ([2009](https://www.jstatsoft.org/article/view/v031i10))

---- 
silvia.innocenti@ec.gc.ca

Innocenti, S. and Matte, P. and Fortin, V. and Bernier, N. B., (2022), Residual bootstrap methods for parameter uncertainty assessment in tidal analysis with temporally correlated noise, JTECH, [DOI:10.1175/JTECH-D-21-0060.1](https://journals.ametsoc.org/view/journals/atot/aop/JTECH-D-21-0060.1/JTECH-D-21-0060.1.xml)

# Contents
1. [boot_tide.m:](#boot_tide.m) construct the residual bootstrap resamples/replicates.  
2. [boot_tide_param.m:](#boot_tide_param.m) compute the bootstrap plug-in estimators of the parameters and corresponding confidence intervals.
3. examples and templates folder


### boot_tide.m
Construct residual bootstrap resamples given an or some observed series and 
a tidal reconstruction (e.g., tidal HA regression predictions) or the corresponding
residuals (i.e., observed series - reconstruction).

This function implements the MBB and SPB algorithms (`method` option) to generate the resamples of the water level series and, if requested, estimate the parameters of a user-specified model (`tidal_model` option) for each of the generated resamples. 
   
```MATLAB
    % compute resamples
    out = boot_tide(time, y, {'option_name',option_value})
```

general options: 
- `n_boot`: number of bootstrap resamples to generate
- `seed`: pseudo-number generator seed: generate an odd number as seed  
- `method`: bootstrap type: MBB or SPB

options for MBB only:
- `circular`: if true, use circular MBB 
- `lblock_dist`: distribution (MATLAB name) for simulating the block length
- `lblock_par`: block length distribution parameter
- `lblock`: MBB block length [in hours] (only for 'fix' MBB)
- `block_output`: if true, return the MBB matrices used for constructing the resamples (e.g., time indices)

options for SPB only (EXPERIMENTAL):
- `noise`: method to estimate the residual spectrum: fft, fftnoise, skfft, sknoise, cpfft, cpnoise, dllftt, or dllnoise

option for computing the parameter resamples: 
- `tide_model`: function handle for estimating the tidal parameters on each y resample 
- `theta`: vector of parameters calculated on the original sample

### boot_tide_param.m
Calculate the plug-in bootstrap estimators of the parameter replicates of a tidal model (not necessarily produced via `boot_tide.m`, it only needs a matrix of parameter replicates) and the corresponding Confidence Intervals (CIs). 
Via the specification of the `circular` option, the user can request to use circular statistics for the parameters expressed in radiant or degrees units [e.g., tidal phases]. 

```MATLAB
    % bootstrap estimates of the parameters (bootstrap distribution means)
    [m_theta]  = boot_tide_param(theta, {'option_name',option_value}) 

    % bootstrap estimates of the parameters and their standard error / circular variance
    [m_theta,s_theta]  = boot_tide_param(theta, {'option_name',option_value})

    % bootstrap estimates of the parameters and their CIs
    [m_theta,s_theta,ci_theta]  = boot_tide_param(theta, {'option_name',option_value})
```

options: 
- `circular`: vector indicating which component of theta are circular variables (e.g., phases)  
- `circ_units`: units of circular parameters: 'radians' or 'degrees'
- `ci`: method for computing the parameter CI: 'percentile', 'gaussian', or 'studentized'
- `alpha`: CI confidence level
- `theta0`: vector of parameters computed on the original samples (used for studentized boot CI)
