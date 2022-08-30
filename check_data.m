function data = check_data(data)
% Apply basic check on input data
%
% INPUT: data - struct with fields:
% time - Tx1 vector of arbitrary times [datenum UCT/GMT],
%        It may include NaNs and time steps may be irregularly distributed 
%  
% y    - TxS time series of observed water levels at the T time steps for S 
%        spatial locations. If S>1 all resamples will be constructed using
%        the same time index resampling to conserve the spatial structure of
%        data. y may include NaNs.
%
% varargin must contain at least one of the following pairs:
% 'yhat', yhat   - TxS time series of reconstructed water levels at the T time steps for S 
%                  spatial locations.
% 'yres', yres   - TxS time series of resuduals (y-yhat) at the T time steps for S 
%                  spatial locations.
% 'tresol', tr   - data (avg) temporal resolution in hours. Default: 1 [h]
%
% OUTPUT:
% data - structure with fields       
%        t    - Tx1 vector of times [datenum UCT/GMT]
%        y    - TxS time series of observed water levels (T time steps, S locations)
%        yhat - TxS time series of reconstructed water levels 
%        yres - TxS time series of residuals 
%       valid - TxS time series of valid data index [non NaN time, observation, and residual values]
%      tresol - 1x1 temporal resolution of data [h], may be approximate or avg resolution
%
% S.Innocenti, silvia.innocenti@ec.gc.ca, 2022/07
        
% make time to be a col vector
t  = data.t;
lt = length(t);
y  = data.y;

% get the observation matrix size
[ry,cy] = size(y);

% chekc time/y dim
assert(lt==ry, 'inconsistent dimensions for time and y') 


% check that at least a tidal reconstruction or the residuals are in input
[ryh,cyh] = size(data.yhat);
[rr,cr]   = size(data.yres);

assert(ryh>0 || rr>0, 'at least one of yhat or yres (reconstructions or residuals) must be given in input')


% if yhat in input
if ryh>0
   assert(ryh==lt, 'inconsistent dimensions for time and yhat ') 
   assert(cyh==cy, 'inconsistent dimensions for y and yhat ') 
   
   % compute the resuals based on yhat and y
   yres = y-data.yhat; 
   
   if rr==0
       data.yres = yres;
   else
       
       % check that the computed and input residual are the same
       for s = 1 : cyh
           if any((roundn(data.yres(:,s),-3)==roundn(yres(:,s),-3))~=0)
               warning('check_data: yres may be different than y-yhat at location ' + num2str(s))
           end
       end    
   end 
end


% if residuals in input
if rr>0
   assert(rr==lt, 'inconsistent dimensions for time and yres ') 
   assert(cr==cy, 'inconsistent dimensions for y and yres ') 
   
   % compute yhat based on the resuals 
   yhat = y-data.yres; 
   
   if ryh==0
       data.yhat = yhat;
   else
       
       % check that the computed and input yhat are the same
       for s = 1 : cr
           if any((roundn(data.yhat(:,s),-3)==roundn(yhat(:,s),-3))~=0)
               warning('check_data: yhat may be different than y-yres at location ' + num2str(s))
           end
       end     
   end     
   
end


% compute valid data indices:
data.valid = isfinite(data.t) & (data.t>=0) & all(isfinite(data.y),2) & all(isfinite(data.yhat),2);
% TODO: change the way of dealing with missing values > add the possibility of choosing valid time steps at 
% all the locations or one location at time


  
end