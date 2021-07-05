function [u,v,opt,B] =  cut_boot_reconstr(tin,coef_boot ,varargin)
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04
% Based on ut_reconstr() and ut_solv() [UTide v1p0 2011/09, 
% http://www.po.gso.uri.edu/~codiga/utide/utide.htm], and their 
% adapted versions cut_reconstr() and cut_solv()
%
% TO DO:
% - ...  
%
% OVERVIEW: 
% Generate hindcasts and/or forecasts/predictions at  
% user-specified times for each bootstrap replicate of the HA parameters. 
% This function is similar to cut_reconstr.m but avoids recomputing the model 
% basis function for each set of parameters. Accordingly, the HA model basis 
% function is computed once, or it is read from the optional input arguments
% to further increase the computational efficiency.
%
% INPUT: 
% t_fit = Tx1 vector of arbitrary times [datenum UCT/GMT] distributed uniformly 
%         or irregularly. Hindcasts are produced at these time steps. If t_fit include 
%         NaNs, the outputs (u_fit/v_fit or sl_fit) will put NaNs at the corresponding time steps. 
% coef_boot = results output structure from cut_boot() containing various bootstrap replicates
%         of the HA model parameters.
%
%
% {OPTIONS}:
% 'cnstit', cnstit = (cell list of 4-character strings) list of 
%       constituents, which must be used for the bootstrap hindcasts. The list of 
%       names must include all or a subset of constituents in the coef.name 
%       field of the  cut_solv() output structure.  
%
% 'basis'/'B'/'expbasis', Bmat = Tx(2K+1) matrix specifying the model basis function 
%        (dimensions are specified assuming no trend) to be used for reconstructions 
%        at the specified time steps. 
%        Default: Bmat=[], B is computed via the UTide cut_E() function.
%       
%
%
% OUTPUT:
% u & v =  Txn_boot matrix/matrices of series/hindcasts reconstructed for each 
%           bootstrap parameter replicate. 
% B     = Tx(2K+1) matrix specifying the model basis function
% opt   = structure containing the 'cnstit' and 'basis' input parameters
%

% remove nan
      t = tin(:);
   ntin = numel(tin(:));
    whr = ~isnan(tin);
      t = t(whr);
     

 % checks
 assert(isfield(coef_boot,'M'),'the input coef strcuture must have a field ''M'' with parameter valuess')
 nM    = size(coef_boot.M,1);
 nBoot = size(coef_boot.M,2);
 nC    = size(coef_boot.g,1);
 
% defaults
       aux = coef_boot.boot_opt.aux;
opt.cnstit = coef_boot.name;
opt.B      = [];


% Optional arguments
optargin = size(varargin,2);

i = 1;
while i <= optargin
    switch lower(varargin{i})
        case 'cnstit'
            opt.cnstit = args{2};    
            args(1:2) = [];
            
       case {'basis' 'b' 'expbasis'}  % complex exponential basis function 
            opt.B = args{2};    
            args(1:2) = [];
            
    otherwise 
            error(['cut_boot_reconstr: unrecognized input: ' args{1}]);        
    end
end

% select constituents
[~,ind] = intersect(coef_boot.name,cellstr(opt.cnstit),'stable');
assert(~isempty(ind),'None of the input constituents is in coef_boot.name ')    
    
% model basis function 
if ~isempty(opt.B)
    
    assert(size(opt.B,2)==nM, 'dimension of the input basis function and coef_boot.M matrices are not consistent')
    assert(size(opt.B,1)==nt, 'dimension of the input basis function and time vector are not consistent')
    B = opt.B;
   
else

     ngflgs = [aux.opt.nodsatlint aux.opt.nodsatnone aux.opt.gwchlint aux.opt.gwchnone];
  [E,F,U,V] = cut_E(t,aux.reftime,aux.frq(ind),aux.lind(ind),aux.lat,ngflgs,aux.opt);
         
  if ~ aux.opt.complexf 
        B = [F.*cos(2*pi*(U+V)), F.*sin(2*pi*(U+V))];  
  else % complex formulation 
        B = [E conj(E)];       
  end
 
end

% reconstructions
   u = nan*ones(ntin,nBoot);
if aux.opt.twodim
   v = u;
else
   v = []; 
end

% for each bootstrap
  iap = ind;
  iam = nC+ind; 
for b = 1 : nBoot

        
     ap = coef_boot.M(iap,b);
     am = coef_boot.M(iam,b);
    fit = B(:,iap)*ap +  B(:,iam)*am;

    % mean (& trend)
    if aux.opt.twodim
        if aux.opt.notrend
            u(whr,b) = real(fit) + coef_boot.umean(b);
            v(whr,b) = imag(fit) + coef_boot.vmean(b);
        else
            u(whr,b) = real(fit) + coef_boot.umean(b) + coef_boot.uslope(b)*(t-aux.reftime);
            v(whr,b) = imag(fit) + coef_boot.vmean(b) + coef_boot.vslope(b)*(t-aux.reftime);
        end
    else
        if aux.opt.notrend
            u(whr,b) = real(fit) + coef_boot.mean(b);
        else
            u(whr,b) = real(fit) + coef_boot.mean(b) + coef_boot.slope(b)*(t- aux.reftime);
        end
    end

end

end