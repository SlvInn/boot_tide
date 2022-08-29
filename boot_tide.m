function out = boot_tide(time,y,varargin) 
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04 - 2022/07
%
% Construct residual bootstrap resamples given an or some observed series and 
% a tidal reconstruction (e.g., tidal HA regression predictions)  
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
%                spatial locations.
% 'yres', yres - TxS time series of resuduals (y-yhat) at the T time steps for S 
%                spatial locations
%
%
% {OPTIONS/VARARGIN}:
% 'n_boot', int - Number of bootstrap resamples to generate. Default: n_boot = 10^3
% 
% 'method', mtdname - string indicating which bootstrap algorithm must be used 
%           to construct the residual resamples. Accepted strings are:
%          'mbb' for implementing the Moving Block Boostrap (either with fixed or random block length')
%          'spb' for the Semi-Parametric Bootstrap based on the residual power spectrum 
%           Default: 'mbb'
%
% 'circular', true/false  - boolean flag determining if a circular strategy 
%           [Politis and Romano (1994)] must be used in constructing the resamples.
%           Default: false
%
% 'lblock_dist', distname - string defining the probability distribution 
%            to extract/simulate the block lengths for each MBB resample. 
%            Default: 'geom'
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
%            * min and max block length for blockdist = 'unif' (i.e. nhours is a 2-element vector). Default: [0,60*24] 
%                      
% 'lblock', nhours - block length [in hours] for fixed length block in MBB 
%           Default: 30*24 + 1 
%
% NOTE: The 'circular', 'lblock_dist', 'lblock_par', and 'lblock' options
%       are ignored when using SPB. 
%
%
% 'seed', seeds - 2-element vector containing the random seeds to be used in the 
%           resampling. For MBB, the first seed is used in the simulation of the block 
%           random starts, while the second is used in the simulation of the block lengths.
%           For the SPB, only the first seed is used to initialize the fttnoise function.      
%       
% 'yres', yres - Tx1 vector or Tx2 matrix  of known regression residuals 
%       to be used to construct the bootstrap resamples. Default: [] - Estimated by cut_reconstr1.m
%
% 'yhat' , yhat - Tx1 vector or Tx2 matrix  of known regression predictions  
%       to be used to construct the bootstrap resamples. Default: [] - Estimated by cut_reconstr1.m
%
%
% 
% OUTPUT: 
% structure with fields:
% options -
% yhat    - 
% yres    -
%
% For options.method = 'mbb', the following fields are also returned: 
% iboot  - 
% blocks -
%
% For options.method = 'spb', the following fields are also returned: 
% yboot  - 
% psd    - 
%
% TODO:
% - test for non-hourly data
% - test for hourly data in UTide and NS_Tide
%
% - allow y to be a TxSx2 matrix for handling 2D currents
% - ?? add the possibility of estimating a regression model ???
% 
% % % ====== NOTES ===== % % % 
% if working with irregular times, attention in setting the bootstrap
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
        [out.yboot,out.psd] = get_spb(data,boot_opt);
    end


end