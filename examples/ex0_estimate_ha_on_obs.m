% ex0_estimate_ha_on_obs.m
% matlab -nodisplay -batch ex0_estimate_ha_on_obs.m
% 
% Apply UTide harmonic analysis to get water level reconstructions 
% for next exxemples.
% Note: aux() is a function container used to store miscellaneaous functions in aux.m
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 

clc
clear

addpath('../.boot_UTide./')

% fix the period
t0 = datenum(2000,01,01,0,0,0);  % date start
te = datenum(2000,12,31,23,0,0); % date end

% load data at 2 stations 
stations = {'Halifax', 'St.Francois'};
      ns = length(stations);
   coord = [47.000252, -70.808219; 44.666667, -63.583333;];

[hobs(:,1),t] = aux('read_hourly_data','./data/hobs_HAL_1999_2009.mat',t0,te); % water level and time @ Halifax
[hobs(:,2),~] = aux('read_hourly_data','./data/hobs_STF_1999_2009.mat',t0,te); % water level @ St.Francois

% utide parameters 
optut = {'NoTrend',... % remove secular trend (detrend linearly)
        'Rmin',0.1,... % Rmin for Rayleigh criterion (automated constituent selection)
        'Nrlzn',10^4,.... % use Nrlzn repetitions in the Monte Carlo approximation of error
        'DiagnMinSNR',2,... % SRNmin value for  the signal-to-noise
        'OrderCnstit','frq',...  % order by increasing frequency the coef output structure
        'RunTimeDisp','nnn', ... % suppress all output printet at screen
        'sincosf'};


% empties
cf = cell(2,1);
cfstats = cell(2,1);



% solve at station 1 
[cf,cfstats]   = cut_solv(t,hobs(:,1),[],coord(1,1),'auto',optut{:});

% reconstruction
cf.boot_opt.aux   = cf.aux;
[h,~,~,B] = cut_boot_reconstr(t,cf); 

% save data
var2save = {'t','h','B','cf','optut'};
save('./data/hhat_HAL_2000.mat',var2save{:})

% solve at station 2 
[cf,cfstats]   = cut_solv(t,hobs(:,2),[],coord(2,1),'auto',optut{:});

% reconstruction
cf.boot_opt.aux   = cf.aux;
[h,~,~,B] = cut_boot_reconstr(t,cf); 

% save data
save('./data/hhat_STF_2000.mat',var2save{:})

