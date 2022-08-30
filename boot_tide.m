function out = boot_tide(time,y,varargin) 
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04 - 2022/07
%
% Construct residual bootstrap resamples given an or some observed series and 
% a tidal reconstruction (e.g., tidal HA regression predictions) or the corresponding
% residuals (i.e., observed series-reconstruction).
%
% INPUT:
% time  - Tx1 vector of arbitrary times [datenum UCT/GMT],
%         It may include NaNs and time steps may be irregularly distributed. 
%  
% y     - TxS time series of observed water levels at the T time steps for S 
%         spatial locations. If S>1 all resamples will be constructed using
%         the same time index resampling to conserve the spatial structure of
%         data. y may include NaNs.
% 
% {OPTIONS/VARARGIN} must contain at least one of the following pairs:
% 'yhat', yhat - TxS time series of reconstructed water levels at the T time steps for S 
%                spatial locations. yhat must be s.t. yhat = y - yres
% 'yres', yres - TxS time series of resuduals (y-yhat) at the T time steps for S 
%                spatial locations. yres must be s.t. yres = y - yhat
%
%
% {OPTIONS/VARARGIN}:
% 'n_boot', int - Number of bootstrap resamples to generate. Default: n_boot = 10^3
% 
% 'method', mtdname - string indicating which bootstrap algorithm must be used 
%           to construct the residual resamples. Accepted strings are:
%          'mbb' for implementing the Moving Block Boostrap (either with fixed or random block length')
%          'spb' for the Semi-Parametric Bootstrap based on the residual power spectrum.
%           Default: 'mbb'.
%
%   !! NOTE !! only 'MBB' preserves the spatial consistency of the estimates at S locations,
%           while 'SPB' ignores the spatial dependence of the residals. 
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
%            * 'pois' - blocks with random lengths esimulated from a poisson distribution with 
%                       parameter lambda = nhours (see 'lblock_par' definition below)
%            * 'unif' - blocks with random lengths esimulated from an uniform distribution with 
%                       parameters [a,b] = nhours (see 'lblock_par' definition below)
%
% 'lblock_par', theta - block length distribution parameter(s) when using 
%           random blocks MBB. Specifically, theta has the following interpretations:
%            * inverse of the average block length for blockdist ='geom'. Default: 1/(30*24 +1) 
%            * average block length for blockdist = 'pois'. Default: 30*24 
%            * min and max block length for blockdist = 'unif' (i.e. nhours is a 2-element vector). 
%            Default: [0,60*24], ignored when using SPB.  
%                      
% 'lblock', nhours - block length [in hours] for fixed length block in MBB 
%           Default: 30*24 + 1, ignored when using SPB.  
%
%
% 'seed', seeds - 2-element vector containing the random seeds to be used in the 
%           resampling. For MBB, the first seed is used in the simulation of the block 
%           random starts, while the second is used in the simulation of the block lengths.
%           For the SPB, only the first seed is used to initialize the fttnoise function.      
%     
% 'spnoise', fname - name (string) of the ftt method or noise generating function to be used to 
%            simulate semi-parametrically residual resamples with the same spectrum as the observed residuals. 
%            Allowed values are:
%               - 'fft'    | 'fftnoise': use the fftnoise.m function adapred by 
%               - 'skfft'  | 'sknoise'  | 'skewed': use the sknoise.m function 
%               - 'cpfft'  | 'cpnoise': use the cpnoise.m function 
%               - 'dllftt' | 'dllnoise' | 'daniell': use the dllnoise.m function 
%            >> See get_spb_noise.m for a detailed description of each function. 
%          
% 
% OUTPUT: 
% structure with fields:
% yhat    - original water level reconstructions: yhat = y- yres. 
% yres    - original water level residuals: yres = y - yhat.
% options - structure with option values used for generating MBB or SPB resamples.
%
% For options.method = 'mbb', the following fields are also returned: 
% iboot  - (T x nboot) matrix of time indices that can be used for constructing the y bootstrap resamples: 
%           for each replicate, yboot_b = yhay + yres(iboot(:,b),:), b = 1,2, ..., n_boot; yboot_b = (T x S) matrix
% blocks - (nboot x 2) cell with the series of block lengths and block start indices for each bootstrap resample
%
%   !! NOTE !! iboot and blocks values apply to the S spatial locations to preserve the spatial dependence structure of residuals 
%
% For options.method = 'spb', the following fields are also returned: 
% eboot   - (T x nboot x S) matrix of the simulated residuals at each spatial location (nboot x S 
%            independent realizations of noise series of length T).
% spectra - (T x S) matrix of the estimated residual spectrum at each spatial location.
%
% TODO:
% - test for non-hourly data
% - test for hourly data in UTide and NS_Tide
% - add the possibility of estimating a HA regression or other tidal models 
% - add the computation of bootstrap statistics (estimators) and CI
% - extend the account for the spatial dependece of S locations to the SPB  
%
%   !! NOTE !!   
% If working with irregular times, attention in setting the bootstrap
% parameters (check the block length, etc., since internal variables and 
% scripts have only been tested for hourly data)
%

    % define the data structure with fields t, y, yhat, and yres
    % and check the consistency of these variables
    data = check_data(time,y,varargin);


    % set the bootstrap options based on default and user-defined values
    boot_opt = set_options('boot',varargin);
    boot_opt = check_boot_opt(data,boot_opt);

    % intialize the output structure based on the input data and options
    out = set_options('output',{'options',boot_opt, 'yhat',data.yhat, 'yres',data.yres});


    % get the resamples (time indices and or simulated y)
    if strcmpi(opt.method,'mbb')
        [out.iboot,out.blocks] = get_mbb(data,boot_opt);
    else
        [out.eboot,out.spectra] = get_spb(data,boot_opt);
    end


end