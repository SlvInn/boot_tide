%%--------------------------------------------------------- 
function [nt,t,u,v,tref,lor,elor,opt,tgd,uvgd] = cut_slvinit(tin,uin,vin,cnstit,args)
% S. Innocenti. Adapted from UT_SLVINIT()
% 2019/09
%
% Initial input parsing and conditioning for UT_SOLV1()
% inputs
%        tin = input times [datenum UTC] passed in to ut_solv (ntin x 1)
%   uin, vin = raw inputs [units arb.] passed in to ut_solv (ntin x 1)
%     cnstit = cell array of 4-char strings passed in to ut_solv (nc x 1)
%       args = varargin from call to ut_solv
% outputs
%     nt = number of t/u/v values ( b/c of nans, can be < length(tin) )
%      t = tin values having non-nan uin, vin [datenum UTC] (nt x 1)
%    u,v = uin,vin elements [units arb.] where both are non-nan (nt x 1)
%   tref = central reference time [datenum UTC] (1 x 1)
%    lor = length of record, (max(t) - min(t)) [days] (1 x 1)
%   elor = effective length of record, lor*nt/(nt-1) [days] (1 x 1)
%    opt = option flags information
%    tgd = result of ~isnan(tin) (ntin x 1)
%   uvgd = result of ~isnan(uin) & tgd [ & ~isnan(vin) ] (ntin x 1)
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

%% input checking/reshaping
assert(isequal(size(tin),size(uin)), 'ct_solv: vectors of input times and input values must be same size.')

    tgd = ~isnan(tin);
    uin = uin(tgd);

if ~isempty(vin)
    assert(isequal(size(tin),size(vin)), 'cut_solv: vectors of input times and input values must be same size.')
    vin = vin(tgd);
    vin = vin(:);
    opt.twodim = 1;
    
else
    opt.twodim = 0;
    v = [];
end

    tin = tin(tgd);
    tin = tin(:);
   
% some basic checks    
assert(isequal(sort(tin),tin) || ~isequal(length(unique(tin)),length(tin)), 'cut_solv: input times must be monotonically increasing.')    
assert(~(min(tin)<datenum(500,0,0) || max(tin)>datenum(3500,0,0)),'cut_solv: input times should be from the present millenium.')

uin = uin(:);
assert(isequal(size(uin,2),1) || size(uin,1)<1,'cut_solv: multiple input values required.')


% only keep t,u,v values with non-nan u(&v) 
if opt.twodim
    uvgd = ~isnan(uin) & ~isnan(vin);
       v = vin(uvgd);
else
    uvgd = ~isnan(uin);
end

 t = tin(uvgd);
nt = length(t);
 u = uin(uvgd);
 
if var(unique(diff(tin)))<eps
    
    opt.equi = 1; % based on times; u/v can still have nans ("gappy")
         lor = (max(tin)-min(tin));
        elor = lor*length(tin)/(length(tin)-1);
        tref = 0.5*(tin(1)+tin(end));
else
    
    opt.equi = 0;
         lor = (max(t) - min(t));
        elor = lor*nt/(nt-1);
        tref = 0.5*(t(1)+t(end));
end

%% options - set defaults
opt.complexf   = 1; % Added (SI)

opt.notrend    = 0;
opt.prefilt    = [];
opt.nodsatlint = 0;
opt.nodsatnone = 0;

opt.gwchlint = 0;
opt.gwchnone = 0;

opt.infer = [];
opt.inferaprx = 0;
opt.rmin = 1;


opt.method  = 'cauchy';
methnotset  = 1;
opt.tunrdn  = 1;
allmethods = {'ols','andrews','bisquare','fair','huber','logistic','talwar','welsch',...
              'cauchy','ridge','ridgesvd','lasso','enet'}; % Modified (SI): added the 'ridge','ridgesvd','lasso','enet' regression but only used in experimental mode
          
% opt.regalgo = ''; % Added (SI)
opt.alpha   = []; % Added (SI): alpha parameters of the lasso and elastic net regressions

opt.linci = 0;
opt.white = 0;
opt.nrlzn = 200; 

opt.lsfrqosmp   = 1;
opt.nodiagn     = 0;
opt.diagnplots  = 0;
opt.diagnminsnr = 2;

opt.ordercnstit = [];
opt.runtimedisp = 'yyy';


while ~isempty(args)
    switch(lower(args{1}))
        case 'sincosf'        % Added (SI)
            opt.complexf = 0; % Added (SI)
                 args(1) = [];% Added (SI)
   
        case 'complexf'      % Added (SI)
            opt.complexf = 1;% Added (SI)
            args(1) = [];    % Added (SI)    
            
        case 'notrend'
            opt.notrend = 1;
            args(1) = [];
            
        case 'prefilt'
            opt.prefilt = args{2};
            args(1:2) = [];
            
        case 'nodsatlint'
            opt.nodsatlint = 1;
            args(1) = [];
            
        case 'nodsatnone'
            opt.nodsatnone = 1;
            args(1) = [];
            
        case 'gwchlint'
            opt.gwchlint = 1;
            args(1) = [];
            
        case 'gwchnone'
            opt.gwchnone = 1;
            args(1) = [];
            
        case 'infer'
            opt.infer = args{2};
            args(1:2) = [];
            
        case 'inferaprx'
            opt.inferaprx = 1;
            args(1) = []; 
            
        case 'rmin'
            opt.rmin = args{2};
            args(1:2) = [];
            
        case allmethods
            if methnotset
                opt.method = lower(args{1});
                   args(1) = [];
                methnotset = 0;
            else
                error('cut_solv: only one ''method'' option allowed.');
            end
            
        case 'alpha' % Added (SI) but experimental only - Elastic net parameter
             assert(numel(args{2})<=1,'alpha must be empty or scalar')
             opt.alpha = args{2}; % Added (SI)
             args(1:2) = [];      % Added (SI)
             warning('''alpha'' parameter found in input >> The elastic-net estimator has not been thoroughly tested')
             
        case 'regalgo' % Added (SI) but experimental only
             assert(any(strcmpi(args{2},{'' 'svd' 'bic'})),'regalgo must be empty or {svd | bic}')
             opt.regalgo = args{2}; % Added (SI)
             args(1:2) = [];        % Added (SI)
             warning('''regalgo'' parameter found in input >> The ridge, lasso, and elastic-net estimators have not been thoroughly tested')
            
             
        case 'tunrdn'
            opt.tunrdn = args{2};
            args(1:2) = [];
            
        case 'linci'
            opt.linci = 1;
            args(1) = [];
            
        case 'white'
            opt.white = 1;
            args(1) = [];
            
            
        case 'nrlzn'
            opt.nrlzn = args{2};
            args(1:2) = [];
            
        case 'lsfrqosmp'
            opt.lsfrqosmp = round(args{2});
            args(1:2) = [];
            
        case 'nodiagn'
            opt.nodiagn = 1;
            args(1) = [];
            
        case 'diagnplots'
            opt.diagnplots = 1;
            args(1) = [];
            
        case 'diagnminsnr'
            opt.diagnminsnr = args{2};
            args(1:2) = [];
            
        case 'ordercnstit'
            opt.ordercnstit = args{2};
            args(1:2) = [];
            
        case 'runtimedisp'
            opt.runtimedisp = lower(args{2});
            args(1:2) = [];
        otherwise
            error(['cut_solv: unrecognized input: ' args{1}]);
    end
end

% check input consistency
if  opt.twodim == 1
    assert(opt.complexf == 1, 'the real sin/cos formulation cannot be used for the twodim case')   
end

% checks on nodal corrections
assert((opt.nodsatlint + opt.nodsatnone)<2,'cut_solv: only one of ''NodSatNone'' or ''NodSatLinT'' is allowed.')
assert((opt.gwchlint   + opt.gwchnone<2), 'cut_solv: inconsistent inputs (only one gwch option allowed).')
assert(opt.lsfrqosmp>=1,'cut_solv: Lomb-Scargle frequency oversampling factor must be >=1.' )

% checks on the constituent selection parameters
if ~isempty(opt.ordercnstit) 
    if ~isequal(opt.ordercnstit,'frq') && ~isequal(opt.ordercnstit,'snr')        
        assert(~isequal(lower(cnstit),'auto'),'cut_solv: OrderCnstit argument must be ''frq'' or ''snr'' since cnstit=''auto''.' )
        assert(~isequal(sort(opt.ordercnstit),sort(cnstit)),'cut_solv: OrderCnstit argument (if not ''SNR'' nor ''Frq'') must include same elements as ''cnstit''.')         
    end
end
assert(~(opt.nodiagn && isequal(opt.ordercnstit,'snr')),['cut_solv: OrderCnstit cannot be ''snr'' if NoDiagn is selected.']);

% check user-define input for inference parameters
if ~isempty(opt.infer)
    ninf = length(opt.infer.infnam);
    
    if opt.twodim 
        assert(~(isequal(size(opt.infer.amprat,1),2*ninf) || ~isequal(size(opt.infer.amprat,1),2*ninf)),...
            'cut_solv: opt.infer.amprat/phsoff must each have 2*length(opt.infer.infnam) elements.')
    else
       assert(~(~isequal(size(opt.infer.amprat,1),ninf) || ~isequal(size(opt.infer.amprat,1),ninf)), ...
           'cut_solv: opt.infer.amprat/phsoff must each have length(opt.infer.infnam) elements.')
    end
end
if ~isempty(opt.infer) && opt.inferaprx
    assert(  size(opt.infer.refnam,1) <= length(unique(cellstr(opt.infer.refnam))), ...
        'cut_solv: cannot infer multiple constituents from one reference constituent with ''InfrAprx''.')
    assert((opt.nodsatlint||opt.nodsatnone) && (opt.gwchlint||opt.gwchnone), ...
        'cut_solv: ''InfrAprx'' requires ''nodsatlint'' or ''nodsatnone'', as well as ''gwchlint'' or ''gwchnone''.') 
end
       
% check estimator input parameters % Modified (SI)
if ~any(strcmp(opt.method,{'cauchy' 'ridge' 'lasso' 'enet' 'conha'})) % Modified (SI)
         [~,ind] = ismember(opt.method,allmethods);
        allconst = [NaN 1.339 4.685 1.400 1.345 1.205 2.795 2.985]; 
    opt.tunconst = allconst(ind);
else
    opt.tunconst = 2.385; % for cauchy 
    % opt.tunconst = 0.1; % (SI) test 2020.04.17
    
    % check alpha (SI)
    switch opt.method
        case 'ridge'
             opt.alpha = 0;
                 if ~(isempty(opt.alpha) || opt.alpha==0)
                      opt.alpha = 0;
                      disp('warning: Ridge SVD parameters: alpha set = 0')
                 end
            
        case 'lasso'  
                if ~(isempty(opt.alpha) || opt.alpha==1)
                  disp('warning: Lasso and Lasso SVD parameters: alpha set = 1')
                end
             opt.alpha = 1;
        case 'enet'   
             
      if isempty(opt.regalgo)
          % assert(~isempty(opt.alpha),'Elastic net parameters: alpha can not be empty')  
          if isempty(opt.alpha)
              disp('warning: alpha ENET set to alpha = 0.5') 
              opt.alpha = 0.5;
          end
      else
          if ~isempty(opt.alpha)
            assert(opt.alpha>0 && opt.alpha<1,'Elastic net parameters: alpha must be 0<alpha<1') 
          end
      end   

    end           
end


% check output preferences parameters
assert((opt.nodiagn + opt.diagnplots)<2,'cut_solv: cannot use both ''NoDiagn'' and ''DiagnPlots''  options.')

nf  = length(fieldnames(opt));
opt = orderfields(opt,[1:12 nf 13:(nf-1)]);
assert(isequal(length(opt.runtimedisp),3),'cut_solv: ''RunTimeDisp'' must be a 3-character string.')

if ~ismember(opt.runtimedisp(1),{'y','n'}) || ...
   ~ismember(opt.runtimedisp(2),{'y','n'}) || ...
   ~ismember(opt.runtimedisp(3),{'y','n'})
    error('cut_solv: ''RunTimeDisp'', while not case-sensitive, must include only y or n characters.');
end



% (SI): add warnings for some specific configurations
if opt.white == 1
   if ~strcmpi(opt.method, 'ols')
      warning(['''white'' and ' opt.method ' options found in input >>' ... 
          'using a WLS estimator (var-cov matrix ~= withe noise)'])
   end
end 

% if opt.nrlzn<100
%      warning([num2str(opt.nrlzn) ' MC replicate may not suffice to approximate parameter distribution ' ...
%          ' >> consider to increase the nrlzn parameter'])
% end
end
