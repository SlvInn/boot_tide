function [coef,Xu,Yu,Xv,Yv] = cut_from_compx_to_amppha(coef,m,nNR,nR,Q,beta,cnstit,lor,opt)
% S. Innocenti adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% 2019/09
% Pass from complex regression parameters to tidal amplitude and phases 

% nc = nNR+nR;
ap = m([1:nNR       2*nNR+(1:nR)]);
am = m([nNR+(1:nNR) 2*nNR+nR+(1:nR)]);
if ~isempty(opt.infer)
    if opt.inferaprx
        for k = 1:nR
              ap(nNR+k) = ap(nNR+k)./(1+beta(k)*cnstit.R{k}.I.Rp.*Q(k));
              am(nNR+k) = am(nNR+k)./(1+beta(k)*cnstit.R{k}.I.Rm.*conj(Q(k)));
        end
    end
end
Xu = real(ap+am);
Yu = -imag(ap-am);
if ~opt.twodim
    [coef.A,~,~,coef.g] = ut_cs2cep([Xu Yu]);
    Xv = [];
    Yv = [];
else
    Xv = imag(ap+am);
    Yv = real(ap-am);
    [coef.Lsmaj,coef.Lsmin,coef.theta,coef.g] = ut_cs2cep([Xu Yu Xv Yv]);
end
    
%% mean and trend
if opt.twodim
    if opt.notrend
        coef.umean = real(m(end));
        coef.vmean = imag(m(end));
    else
        coef.umean = real(m(end-1));
        coef.vmean = imag(m(end-1));
        coef.uslope = real(m(end))/lor;
        coef.vslope = imag(m(end))/lor;
    end
else
    if opt.notrend
        coef.mean = real(m(end));
    else
        coef.mean  = real(m(end-1));
        coef.slope = real(m(end))/lor;
    end
end


end