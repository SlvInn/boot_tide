# MATLAB scripts for applying UTide HA and a resampling-based uncertainty analysis. 
S. Innocenti - Last modification 30/06/2021
Largely adapted from the UTide MATLAB package, Codiga (2011)
[UTide v1p0 9/2011 d.codiga@gso.uri.edu; http://www.po.gso.uri.edu/~codiga/utide/utide.htm]

Users should be primarily interested in calling the following 4 functions:
1. `cut_solv.m`     - Adamptation of the *_ut_solv()_* UTide function. 
2. `cut_reconstr.m` - Adamptation of the *_ut_reconstr()_* UTide function.
3. `cut_boot.m`     - Contructs *_n_boot_* bootstrap resamples of the HA residuals from a cut_solv-HA estimation and re-estimate the HA regression on each resample. 
4. `cut_boot_reconstr.m` - Generates hindcast and/or forecast/prediction at user-specified times
using the HA output from *_cut_reconstr()_*.

the other scripts generally contains sub-functions for these 4 main programs. 
The *_"cut"_* prefix in file name is used to differentiate the adapted scripts from the original UTide ones. Through the scripts, comments are used to identify and possibly explain the modifications made 
with respect to the original UTide code [the #(SI) and #(PM) tags usually identify these modifications]. 



## Modified UTide package
Execute tidal HA using a modified version of UTide. 

Major modifications to UTide consist in:
* Allowing to solve the HA based on a non-complex formulation of the regression model 
in 1D analyses [i.e. as in Foreman and Henry, 1989, and Innocenti et al., submitted to JTECH].
*  Adding some useful output for 1D analyses [e.g., return some statistics for uncertainty assessments in the CoefDist and coef.aux structures]. **NOTE:** not all these quantities have been added (yet) to the 2D-analysis output.
* Allowing to use ridge, lasso, or elastic-net estimators for the HA regerssion. 
**NOTE:** these code options have not been thoroughly tested (use at your own risk) but the code includes many warnings related to the use of these estimators.
* Separating *_ut_solv.m_* in various scripts and add comments to improve readability of subfuctions. 

**NOTE:** `cut_solv.m` and  `cut_reconstr.m` are intended to work with bi-dimentional data (current series) but they haven't been tested in this use. 

---- 
Foreman and Henry, 1989 - DOI:10.1016/0309-1708(89)90017-1


### cut_solv.m 
Define the analysis setup and execute the leat-square optimization of the regression model. 
Various statistics intended to be used in model selection and residual analysis are also computed. 

Specific analysis steps consist in:
- choice of the estimated and inferred constituents;  
- regression model definition (e.g., trend/no-trend, use nodal corrections, calculation of the design matrix / basis of frequencies, etc);
- choice of the regression estimator (e.g., ols vs irls, irls parameters, etc) and least square optimization; 
- computation of basic statistics for model evaluation, further constituent selection, and residual analysis (e.g., CoefDist output structure).

<!-- #### Subfunctions:
... write docs .... -->

### cut_reconstr.m 
Generate hindcasts and/or forecasts/predictions at user-specified times using the HA coefficient structure computed by CUT_SOLV(). Specific steps of the reconstruction are:
* select a list (e.g., a subset) of constituents from those estimated by cut_solv.m and compute the corresponding model basis function 
* reconstruct the regression prediction based the selected list of constituents

<!-- #### Subfunctions:
... write docs .... -->

## Residual resampling functions 
Functions implementing two residual boostrap algorithms for the analysis of residuals in HA uncertainty assessment:  

1. a Moving Block Bootstrap (MBB) that directly resamples fixed or random length blocks of the HA regression resisuals, and 
2. a Semi-Parametric Bootstrap (SPB) that simulates noise realizations from the Fast Fourier Transform of the observed residuals. 

Both methods are intended to produce similar results with comparable accuracy. The MBB can be applied in virtually all sampling conditions but it is more computationally demanding than SPB, especially for long records. The SPB is faster than the MBB but some sort of interpolation is needed for analyzing  irregularly sampled records or series with missing data, since it involves a Fourier analysis of residuals.  

A detailed description of the two bootstraps can be found in the manuscript from Innocenti et al. (under review), while a teoretical summary is reported in https://www.youtube.com/watch?v=voEg6NEgT8M.

---- 
Innocenti, S. and Matte, P. and Fortin, V. and Bernier, N. B., (under review), Residual bootstrap methods for parameter uncertainty assessment in tidal analysis with temporally correlated noise, JTECH


### cut_boot.m 
Construct residual resamples accoring to the specified algorith and execute the leat-square optimization of the regression model for each bootstrap resample. HA regression parameters for each resample are then expressed as tidal amplitude and phases. 

**Note:** the SPB use fttnoise.m (included in the Boot_Tide repo) by Aslak Grinsted (2009; https://www.mathworks.com/matlabcentral/fileexchange/32111-fftnoise-generate-noise-with-a-specified-power-spectrum).

### cut_boot_reconstr.m 
Generate hindcasts and/or forecasts/predictions at  user-specified times for each bootstrap replicate of the HA parameters. This function is similar to cut_reconstr.m but avoids recomputing the model basis function for each set of parameters. Accordingly, the HA model basis function is computed once, or it is read from the optional input arguments to further increase the computational efficiency.