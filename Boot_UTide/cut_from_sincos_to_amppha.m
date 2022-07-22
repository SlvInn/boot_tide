function [coef,Xu,Yu] = cut_from_sincos_to_amppha(coef,m,nNR,nR,Q,beta,cnstit,lor,opt)
% S. Innocenti adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% 2019/09
% Pass from ssin-cos regression parameters to tidal amplitude and phases 

ap = m([1:nNR       2*nNR+(1:nR)]);
am = m([nNR+(1:nNR) 2*nNR+nR+(1:nR)]);
if ~isempty(opt.infer)
    error 'this must be done: see cut_fromComplex2AmpPha.m ~lines 6-13'
    % Q,beta, and cnstit would eb used here
end

% this must be checked:
Xu = ap;
Yu = am;


coef.A = sqrt(ap.^2 + am.^2);
   pha = atan2d(am,ap);
%    pha = angle(ap/am);
   pha = ut_cluster(pha,360);
   pha(pha<0) = 360+pha(pha<0);
coef.g = pha;



%% mean and trend
if opt.notrend
    coef.mean = real(m(end));
else
    coef.mean  = real(m(end-1));
    coef.slope = real(m(end))/lor;
end


end