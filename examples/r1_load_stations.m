% Load water level data for two stations and the corresponding harmonic analysis 
% Note: aux() is a function container used to store miscellaneaous (sub)-functions in aux.m

addpath('../')

% fix the period
t0 = datenum(2000,01,01,0,0,0);  % date start
te = datenum(2000,12,31,23,0,0); % date end

% load data at 2 stations 
stations = {'Halifax', 'St.Francois'};
nS = length(stations);

% load the observed water levels
[h(:,1),t] = aux('read_hourly_data','./data/hobs_HAL_1999_2009.mat',t0,te); % water level and time @ Halifax
[h(:,2),~] = aux('read_hourly_data','./data/hobs_STF_1999_2009.mat',t0,te); % water level @ St.Francois

% load some reconstructions at these two stations (from HA)
[hhat(:,1),that] = aux('read_hourly_data','./data/hhat_HAL_2000.mat',t0,te); % estimated water level and time @ Halifax
[hhat(:,2),~] = aux('read_hourly_data','./data/hhat_STF_2000.mat',t0,te);  % estimated water level @ St.Francois


% load the classical HA 
load('./data/hhat_HAL_2000.mat','cf','B')
coef.hal = cf;
coef.hal.B = B;
clear B cf

load('./data/hhat_STF_2000.mat','cf','B')
coef.stf = cf;
coef.stf.B = B;
clear B cf 
