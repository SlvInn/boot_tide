% load_two_stations.m
% Load one year of water level data at two stations and the corresponding 
% harmonic analysis (UTide estimates)
% 

addpath('../')
addpath('./auxf/')

% fix the period
t0 = datenum(2000,01,01,0,0,0);  % date start
te = datenum(2000,12,31,23,0,0); % date end

% load data at 2 stations 
stations = {'Halifax', 'St.Francois'};
ns = length(stations);


% load the observed water levels
[h(:,1),t] = read_hourly_data('./data/hobs_HAL_1999_2009.mat',t0,te); % water level and time @ Halifax
[h(:,2),~] = read_hourly_data('./data/hobs_STF_1999_2009.mat',t0,te); % water level @ St.Francois

% load some reconstructions at these two stations (from HA)
[hhat(:,1),that] = read_hourly_data('./data/hhat_HAL_2000.mat',t0,te); % estimated water level and time @ Halifax
[hhat(:,2),~] = read_hourly_data('./data/hhat_STF_2000.mat',t0,te);  % estimated water level @ St.Francois


% load the output from a classical HA 
load('./data/hhat_HAL_2000.mat','cf','B')
coef.hal = cf;
coef.hal.B = B;
clear B cf

load('./data/hhat_STF_2000.mat','cf','B')
coef.stf = cf;
coef.stf.B = B;
clear B cf 
