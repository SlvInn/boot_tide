% ex1_generate_wl_resample.m
% 
% Use the boot_tide.m fonction for generating MBB and SPB water level resamples
% given the water level reconstructions (from tidal HA) at two stations. 
%
% silvia.innocenti@ec.gc.ca, 2022/08/30. 
clc; clear

% load one year of hourly data at two stations:
load_two_stations
% h: water level and time @ Halifax and St.Francois de l'Ile d'Orleans
% hhat: estimated/reconstructed water level and time @ Halifax and St.Francois de l'Ile d'Orleans

% MBB method
disp('>> apply boot_tide MBB without options (default values)')
tic; mbb_out = boot_tide(t(:,1),h,'yhat',hhat); tc = toc;
% this is equivalent to using the following options (see mbb_out.options): 
% {'nboot', 10^4, ...
%  'seed', [89045777, 6442509], ... % randomly generated by default
%  'method', 'mbb', ...
%  'circular', true, ...,  
%  'lblock_dist', 'geom', ...
%  'lblock_par', 0.0014, ...
%  'block_output', false }

disp('used options (defaults):')
disp(mmb_out.options)

% print the computing time
disp(['constructing ' num2str(spb_out.options.nboot) ' MBB resamples took ' num2str(tc) ' sec.'])

% constructed resamples in mbb_out.y_boot:
disp('constructed resamples of warer levels in mbb_out.y_boot:')
disp(mbb_out)

% SPB method
disp('>> apply boot_tide SPB with some options')
tic; spb_out = boot_tide(t(:,1),h,'yhat',hhat,'nboot', 10001,'method','spb','noise','fftnoise');tc = toc;

disp('used options:')
disp(spb_out.options)

% print the computing time
disp(['constructing ' num2str(spb_out.options.nboot) ' SPB resamples took ' num2str(tc) ' sec.'])





