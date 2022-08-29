function [yboot, psd] = get_spb(data,opt)
% Get the simulated y resamples with the Semi-Parametric Bootstrap (SPB) method.
%
% INPUT:
% data  - structure of processed data series   >> See help check_data
% opt   - structure with the bootstrap options >> See help boot_tide.m for a description of the opt structure fields 
% 
% OUTPUT: 
% yboot - (T x nboot x S) matrix of the nboot simulated resamples for each spatial location S
% psd   - ...
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2022/08

% time and residual vector without nan
t    = data.t(valid);
lt   = numel(data.t);
yres = data.yres(valid,:);
ns   = size(yres,2);

# n of resamples
nb = opt.nboot

% empty vector of seeds
seeds = nan(nb,1);

% empty matrix for the residual FTT 
yres_ftt = nan(lt,ns);

% SILVIA: sono at line 165 de cut_bootcore.m