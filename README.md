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
[i.e. as in Foreman and Henry, 1989, and Innoncenti et al., submitted to JTECH].
*  Adding some useful output for 1D analyses [e.g., return some statistics for uncertainty assessments in the CoefDist structure]. NOTE: not all these quantities have been added (yet) to the 2D-analysis output.
* Separating *_ut_solv.m_* in various scripts to improve readability of subfuctions.

Note that `cut_solv.m` and  `cut_reconstr.m` are intended to work with bi-dimentional data (current series) but they haven't been tested in this use. 

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

#### Subfunctions:
... write docs ....

### cut_reconstr.m 
Generate hindcasts and/or forecasts/predictions at user-specified times using the HA coefficient structure computed by CUT_SOLV(). Specific steps of the reconstruction are:
* select a list (e.g., a subset) of constituents from those estimated by cut_solv.m and compute the corresponding model basis function 
* reconstruct the regression prediction based the selected list of constituents

#### Subfunctions:
... write docs ....

## Residual resampling functions 

### cut_boot.m 
### cut_boot_reconstr.m 