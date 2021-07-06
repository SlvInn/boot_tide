%%--------------------------------------------------------- 
function [E,F,U,V] = cut_E(t,tref,frq,lind,lat,ngflgs,opt)
% Silvia Innocenti, adapted from ut_E [UTide v1p0 9/2011 d.codiga@gso.uri.edu]. 
% 2019/10/17
%
% Compute complex exponential basis function
% INPUT
%      t = times [datenum UTC] (nt x 1)
%   tref = reference time [datenum UTC] (1 x 1)
%    frq = frequencies [cph] (nc x 1)
%   lind = list indices of constituents in ut_constants.mat (nc x 1)
%    lat = latitude [deg N] (1 x 1)
%   ngflgs = [NodsatLint NodsatNone GwchLint GwchNone] each 0/1
%       ([0 1 0 1] case not allowed, and not needed, in ut_E)
%   prefilt = 'prefilt' input to ut_solv (correction to account for pre-filtering that was applied to the raw inputs before the analysis)
%
% OUTPUT
%   E = complex exponential basis function [unitless] (nt x nc)
%  

prefilt = opt.prefilt;

nt = length(t);
nc = length(lind);
if ngflgs(2) && ngflgs(4)
    F = ones(nt,nc);
    U = zeros(nt,nc);
    V = 24*(t-tref)*frq';
else
    [F,U,V] = ut_FUV(t,tref,lind,lat,ngflgs);
end
E = F.*exp(1i*(U+V)*2*pi);
if ~isempty(prefilt)
    P = interp1(prefilt.frq,prefilt.P,frq)';
    P( P>max(prefilt.rng) | P<min(prefilt.rng) | isnan(P) )=1;
    E = E.*P(ones(nt,1),:);
    
    if ~opt.complexf
    error 'not done yet: prefiltering treatement must be adapted to the sin/cos definition of the model >> try with complexf option = true'
    % see the other part of the IF statement (when prefilt is ~empty)
    end
end
