%%--------------------------------------------------------- 
function [u,v] = cut_reconstr1(tin,coef,varargin)
% S. Innocenti adapted from ut_reconstr1.m
% 2019/10
%
% UT_RECONSTR1()
% Reconstruction for a single record. See comments for UT_RECONSTR().
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

% fprintf('ut_reconstr: '); % (SI): Commented

% parse inputs and options
[t,opt] = ut_rcninit(tin,varargin);

% determine constituents to include
if ~isempty(opt.cnstit)
    [~,ind] = ismember(cellstr(opt.cnstit),coef.name);
    if ~isequal(length(ind),length(cellstr(opt.cnstit)))
        error('ut_reconstr: one or more of input constituents Cnstit not found in coef.name');
    end
    
else % include constituents with high SNR and high Percent Energy
% (SI): NOTE  THAT for evaluating the component significance, the SNR or equivalent statistics 
% should be computed on complex coefficients, not on amplitudes 
    ind = 1:length(coef.aux.frq);
    if coef.aux.opt.twodim
        SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
         PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
         PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
    else
        SNR = (coef.A.^2)./((coef.A_ci/1.96).^2); % signal-to-noise
         PE = 100*coef.A.^2/sum(coef.A.^2); % Percent Energy
    end
    ind = ind(SNR(ind)>=opt.minsnr & PE(ind)>=opt.minpe);
end


% get complex coefficients
rpd = pi/180;
if coef.aux.opt.twodim
    ap = 0.5*(coef.Lsmaj(ind) + coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) - coef.g(ind))*rpd);
    am = 0.5*(coef.Lsmaj(ind) - coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) + coef.g(ind))*rpd);
else
    ap = 0.5*coef.A(ind).*exp(-1i*coef.g(ind)*rpd);
    am = conj(ap);
end

% options to compute the basis function
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone coef.aux.opt.gwchlint coef.aux.opt.gwchnone];
% fprintf('prep/calcs ... '); % (SI) Commented

% compute basis function 
E = cut_E(t,coef.aux.reftime,coef.aux.frq(ind),coef.aux.lind(ind),coef.aux.lat,ngflgs,coef.aux.opt);

% fit: Xhat
fit = E*ap + conj(E)*am;

% mean (& trend)
  u = nan*ones(size(tin));
whr = ~isnan(tin);
if coef.aux.opt.twodim
    v = u;
    if coef.aux.opt.notrend
        u(whr) = real(fit) + coef.umean;
        v(whr) = imag(fit) + coef.vmean;
    else
        u(whr) = real(fit) + coef.umean + coef.uslope*(t-coef.aux.reftime);
        v(whr) = imag(fit) + coef.vmean + coef.vslope*(t-coef.aux.reftime);
    end
else
    if coef.aux.opt.notrend
        u(whr) = real(fit) + coef.mean;
    else
        u(whr) = real(fit) + coef.mean + coef.slope*(t-coef.aux.reftime);
    end
    v = [];
end
% fprintf('done.\n');% Commented by SI

%%--------------------------------------------------------- 
function [t,opt] = ut_rcninit(tin,args)
% UT_RCNINIT()
% parse inputs and options for UT_RECONSTR1()
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

          t = tin(:);
t(isnan(t)) = [];
 opt.cnstit = [];
 opt.minsnr = 2;
  opt.minpe = 0;
while ~isempty(args)
    switch(lower(args{1}))
        case 'cnstit'
            opt.cnstit = args{2};    
            args(1:2) = [];
        case 'minsnr'
            opt.minsnr = args{2};
            args(1:2) = [];
        case 'minpe'
            opt.minpe = args{2};
            args(1:2) = [];
        otherwise 
            error(['ut_reconstr: unrecognized input: ' args{1}]);
    end
end
