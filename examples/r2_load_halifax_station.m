% Load water level data for one station and the corresponding harmonic analysis 
% Note: aux() is a function container used to store miscellaneaous (sub)-functions in aux.m
addpath('../')

% fix the period
t0 = datenum(2000,01,01,0,0,0);  % date start
te = datenum(2000,12,31,23,0,0); % date end

% load data at 2 stations 
station = 'Halifax';


% load the observed water levels
[h,t] = aux('read_hourly_data','./data/hobs_HAL_1999_2009.mat',t0,te); % water level and time @ Halifax


% load some reconstructions at these two stations (from HA)
[hhat,that] = aux('read_hourly_data','./data/hhat_HAL_2000.mat',t0,te); % estimated water level and time @ Halifax


% load the classical HA 
load('./data/hhat_HAL_2000.mat','cf','B')


