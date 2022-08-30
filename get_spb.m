function [eboot, eftt] = get_spb(data,opt)
% Get the simulated y resamples with the Semi-Parametric Bootstrap (SPB) method.
%
% INPUT:
% data  - structure of processed data series   >> See help check_data
% opt   - structure with the bootstrap options >> See help boot_tide.m for a description of the opt structure fields 
% 
% OUTPUT: 
% eboot - (T x nboot x S) matrix of the nboot simulated resamples for each spatial location S
% eftt  - (T x S) matrix of the estimated residual spectra
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2022/08

% time and residual vector without nan
t    = data.t(valid);
yres = data.yres(valid,:);

% number of spatial locations
ns   = size(yres,2); 
if ns>1
    warning(['SPB residuals are simulated independently at the ' num2str(ns) ' locations: spatial dependency is not accounted for in SPB'])
end

% n of resamples
nb = opt.nboot

% empty vector of seeds
% seeds = nan(nb,1);

% geenerate and store random seeds for each boot 
rng(opt.seed(1))
seeds = 2*floor(rand(nb,ns)*(10^5)) +1;

% empty matrix for the residual FTT 
lt = numel(data.t);
eftt  = nan(lt,ns);
eboot = nan(lt,nb,ns);

% To handle irregularly spaced observations, residuals are interpolated over a time vector with fixed time step
teq = roundn(linspace(min(t), max(t), lt),-4); % rounded at the hour resolution

if ~isequal(teq,roundn(t,-4))
    disp("** warning ** interpolating non-regularly spaced obs and/or missind data")
end

for s = 1 : ns % for each spatial location

    % interpolate non-missing values over a regular time vector [to be changed]
    yreseq = interp1(t,yres(:,s),teq,'pchip'); %'pchip',0


    % generate nb random residuals based on the fft and store residual FFT
    [eboot(:,:,s),eftt(:,s)] = get_spb_noise(y_residual,method,n_series)
    
    
end

