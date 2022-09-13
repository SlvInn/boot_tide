function out = boot_tide(time,y,varargin) 
% Construct residual bootstrap resamples given an or some observed series and 
% a tidal reconstruction (e.g., tidal HA regression predictions) or the corresponding
% residuals (i.e., observed series-reconstruction).
%
%     out = boot_tide(time,y,{'option_name',option_value})
%
% INPUT:
% time  - Tx1 vector of arbitrary times [datenum UCT/GMT],
%         It may include NaNs and time steps may be irregularly distributed. 
%  
% y     - TxS time series of observed water levels at the T time steps for S 
%         spatial locations. If S>1, all resamples will be constructed using
%         the same time index resampling to conserve the spatial structure of
%         data. y may include NaNs.
% 
% VARARGIN >> {'option_name',option_value}: 
% ** Must contain at least one of the following pairs: **
% 'yhat', yhat - TxS time series of reconstructed water levels at the T time steps for S 
%                spatial locations. yhat must be s.t. yhat = y - yres
% 'yres', yres - TxS time series of residuals (y-yhat) at the T time steps for S 
%                spatial locations. yres must be s.t. yres = y - yhat
%
%
% Other VARARGIN >> {'option_name',option_value}:
% 'n_boot', int - Number of bootstrap resamples to generate. Default: n_boot = 10^3
% 
% 'method', mtdname - string indicating which bootstrap algorithm must be used 
%           to construct the residual resamples. Accepted strings are:
%          'mbb' for implementing the Moving Block Boostrap (either with fixed or random block length')
%          'spb' for the Semi-Parametric Bootstrap based on the residual power spectrum.
%           Default: 'mbb'.
%
%           !! NOTE !! 'MBB' try to preserve the spatial consistency of the residuals at S locations,
%                   while 'SPB' ignores any spatial dependence. Also, when considering S locations
%                   with 'SPB', missing values in one series entail missing values at all other locations.
%
% 'circular', true/false  - boolean flag determining if a circular strategy 
%           [Politis and Romano (1994)] must be used in constructing the resamples.
%           Default: false, ignored when using SPB. 
%
% 'lblock_dist', distname - string defining the probability distribution 
%            to extract/simulate the block lengths for each MBB resample. 
%            Default: 'geom', ignored when using SPB. 
%
%            Allowed distname:
%            * 'fix'  - blocks with fixed length = nhours (see the 'blocklength' input below), 
%                       no distribution involved. 
%            * 'geom' - blocks with random lengths esimulated from a geometric distribution with 
%                       parameter 1/p = nhours (see 'lblock_par' definition below)
%            * 'unif' - blocks with random lengths esimulated from an uniform distribution with 
%                       parameters [0,b] = nhours (see 'lblock_par' definition below)
%
% 'lblock_par', lblkpar - block length distribution parameter(s) when using 
%           random blocks MBB. Specifically, theta has the following interpretations:
%            * inverse of the average block length for blockdist ='geom'. Default: 1/(30*24 +1) 
%            * max block length for blockdist = 'unif' (i.e. nhours is a 2-element vector). 
%            Default: 60*24, ignored when using SPB.  
%                      
% 'lblock', nhours - block length [in hours] for fixed length block in MBB 
%           Default: 30*24 + 1, ignored when using SPB.  
%
% 'block_output', bool - if true return the matrices used for construction the MBB resamples  
%                     (time indices, block start, and block length) in the OUTPUT structure. See below, 
%                     the OUTPUT section.
%
%
% 'seed', seeds - 2-element vector containing the random seeds to be used in the 
%           resampling. For MBB, the first seed is used in the simulation of the block 
%           random starts, while the second is used in the simulation of the block lengths.
%           For the SPB, only the first seed is used to initialize the fttnoise function.      
%          
% 'tide_model', @(x) tide_model(x) - function handle returning an estimate (vector) of tidal parameters for a given vector
%                      of water level observations. If provided, boot_tide.m would return an estimate of the tidal parameters 
%                      for each bootstrap resample. Default: [] (no parameter estimate is made on the bootstrap resamples).
%                      NOTE: if 'tide_model' is provided, 'theta' is also needed (see below) 
%  
%
% EXPERIMENTAL OPTIONS (presently deprecated, need more tests):
% 'noise', fname - name (string) of the ftt method (noise generating function) to be used to 
%            simulate semi-parametrically residual resamples with the same spectrum as the observed residuals. 
%            Allowed values are:
%               - 'fft'    | 'fftnoise': use the fftnoise.m function  >> Default: it is the only option that don't print a warning at screen
%               - 'skfft'  | 'sknoise'  | 'skewed': use the sknoise.m function 
%               - 'cpfft'  | 'cpnoise': use the cpnoise.m function 
%               - 'dllftt' | 'dllnoise' | 'daniell': use the dllnoise.m function 
%            >> See get_spb_noise.m for a detailed description of each function. 
%
%
% OUTPUT >> structure with fields:
% yhat       - original water level reconstructions: yhat = y- yres. 
% yres       - original water level residuals: yres = y - yhat.
% options    - structure with option values used for generating MBB or SPB resamples.
% y_boot     - (T x nboot x S) matrix of the nboot simulated resamples for each spatial location S. 
%               Returned if no 'tide_model' has been provided.
% theta_boot - (P x nboot x S)  matrix of the nboot replicates of tide model parameters estimated on the resamples.
%              Returned if a 'tide_model' has been provided.
%
% >>>> optional OUTPUT for MBB (options.method = 'mbb'): 
% If options.block_output=true the following fields are also returned: 
% i_boot  - (T x nboot) matrix of time indices that can be used for constructing the y bootstrap resamples: 
%           for each replicate, y_boot_b = yhay + yres(i_boot(:,b),:), b = 1,2, ..., n_boot; y_boot_b = (T x S) matrix
% blocks - (nboot x 2) cell with the series of block lengths and block start indices for each bootstrap resample
%
% >>>> OUTPUT for SPB (options.method = 'spb'): 
% res_ftt - (T x S) matrix of the estimated S residual spectra
%
% TODO:
% - rename tidal_model to something better
% - test for non-hourly data
% - test for hourly data in UTide and NS_Tide
% - test works if the estimation at multiple location works when prescribing a tidal_model: the model must be the same funtion handle at all locations
% - change the way of dealing with missing values > add the possibility of choosing valid time steps one location at time
% - change the SPB to account for the spatial dependece of S locations 
%
%   
% !! NOTE !! If working with irregular times, attention in setting the bootstrap
% parameters (check the block length, etc., since internal variables and 
% scripts have only been tested for hourly data)
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04 - 2022/09
%

    % define the data structure with fields t, y, yhat, and yres and check their consistency
    data = set_options('data',varargin{:},'t',time(:),'y',y);
    data = check_data(data);

    % set the bootstrap options based on default and user-defined values
    boot_opt = set_options('boot',varargin{:});
    boot_opt = check_boot_options(data,boot_opt);
    
    % intialize the output structure based on the input data and options
    out = set_options('output','options',boot_opt, 'yhat',data.yhat, 'yres',data.yres);

    % get the resamples (time indices and or simulated y)
    if strcmp(boot_opt.method,'mbb')
        out_boot = get_mbb(data,boot_opt);
    else
        out_boot = get_spb(data,boot_opt);
    end


    % add the bootstrap fields to the out structure
    outfields = fieldnames(out_boot);
    for f = 1 : length(outfields)
        out.(outfields{f}) = out_boot.(outfields{f});
    end

end

% OLD / UNUSED OPT:
% 'theta', theta - Px1 vector of parameter estimated on the observed sample [i.e., theta = tide_model(y)] used 
%                  to check the tide_model output and for computing the bootstrap estimates of the tide_model parameters.              
%
% 'theta_circ', bool - Px1 vector of boolean values indicating which theta components are circular variables (e.g., phases).
%                      Default: [true for each theta],  if 'tide_model' and 'theta' are provided 