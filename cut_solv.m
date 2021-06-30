function [coef,CoefDist,NodCorr] = cut_solv(tin,uin,vin,lat,cnstit,varargin)
% S. Innocenti adapted from UT_SOLV(). UTide v1p0 9/2011 d.codiga@gso.uri.edu 
% 2019/09
% 
% Execute harmonic tidal analysis using a modified version of UTide. 
% Major modifications conisist in:
%
% - Allowing to solve the HA based on a non complex formulation of the regression model 
% [i.e. as in Foreman and Henry, 1989, and Innoncenti et al., submitted to JTECH]
%
% - Adding some useful output for 1D-analyses [e.g., return some statistics
% for uncertainty assessments in the CoefDist structure]. 
% NOT ALL these quantities have been added (yet) to the 2D-analysis output.
%
% - Separating ut_solv.m in various scripts to improve readability of subfuctions.
%
%
% OVERVIEW 
% 
% Syntax for two-dimensional raw input, such as velocities:
% 	coef = cut_solv ( tin, uin, vin, lat, cnstit , {options} ); 
% 
% Syntax for one-dimensional raw input, such as sea level:
% 	coef = cut_solv ( tin, uin,  [], lat, cnstit , {options} ); 
% 
% Analysis of groups of records
%   The descriptions that follow next for INPUTS, DEFAULTS, OPTIONS, and 
%   OUTPUT are for treatment of a single record. Following that, 
%   explanations are given for modifications that enable treating a group 
%   of records with a single execution.
%
% INPUT : 
%   * tin      = (ntin x 1) vector of raw times [datenum UTC] passed in to ut_solv 
%   * uin, vin = (ntin x 1) vecotrs of raw input values [units arb.] passed in to ut_solv 
%   * lat      = latitude (required for use in default nod./sat. corrections)
%   * cnstit   = specification of constituents to include in model: 
%                (nc x 1) cell array of 4-char strings passed in to ut_solv 
%
% OUTPUT:
%     nt = number of t/u/v values ( b/c of nans, can be < length(tin) )
%      t = tin values having non-nan uin, vin [datenum UTC] (nt x 1)
%    u,v = uin,vin elements [units arb.] where both are non-nan (nt x 1)
%   tref = central reference time [datenum UTC] (1 x 1)
%    lor = length of record, (max(t) - min(t)) [days] (1 x 1)
%   elor = effective length of record, lor*nt/(nt-1) [days] (1 x 1)
%    opt = option flags information
%    tgd = result of ~isnan(tin) (ntin x 1)
%   uvgd = result of ~isnan(uin) & tgd [ & ~isnan(vin) ] (ntin x 1)
%
% OPTIONS:
% Option flags are not case-sensitive but cannot be abbreviated. The order
% of the option flags is not important but they must be passed in 
% after all other arguments. See report for more complete explanations.
% See ut_solv.m and Codiga (2011) for UTide options. 
%
% New optional input (not in UTide)are:
%
%   * 'sincosf' - string that impose the use of the sin/cos formulation of
%               the HA regression model. Default: non-specified (ie. use
%               the complex formulation)
%   * few experimantal additional optios are explained directly in the code
%      but not officially presented here
%    
%
% OUTPUT
%
%   * coef = results structure from cut_solv()
%
% NEW OUTPUT (compared to ut_solv.m):
%   * coef.Std  - structure in coef containing the vectors of the standard error 
%             estimated for each parameter >> See cut_ci
%
%   * CoefDist  - structure containing the var-cov estimates and MC simulations
%       and the Monte-Carlo simulations for complex parameters, as well as for
%       the amplitudes and phases >> See cut_ci
%
%
% Foreman and Henry, 1989 - DOI:10.1016/0309-1708(89)90017-1
%
%

if isequal(sum(size(uin)>2),1)  % single record
    [coef,CoefDist,NodCorr] = cut_solv1(tin,uin,vin,lat,cnstit,varargin{:});   
    
else % group of records

    % n_t = number of times in each sequence in the group
    n_t = size(uin,1); 
	
    % alln = [n_1 n_2 n_3 ... n_n_d]
    alln = size(uin);
    alln = alln(2:end);
	
    % n_d = number of dimensions of group
    n_d = length(alln);
	
    % n_s = number of time sequences in group
    n_s = prod(alln); 
	
    % size of group
    if isequal(n_d,1)
        gsz = [1 alln];
    else
        gsz = alln;
    end
	
    % check and condition inputs
    if (~isequal(size(tin,1),size(uin,1)) && ...
            ~isequal(size(tin),size(uin))) || (~isempty(vin) && ...
            ~isequal(size(vin),size(uin)))
        error(['ut_solv(group): sizes of t & u, or of their first '...
            'dimensions, must be equal; u & v sizes must also be equal.']);
    end
	
    if ~isequal(size(lat),[1 1]) && ~isequal(size(lat),gsz)
        error('ut_solv(group): lat has neither acceptable size.');
    end
	
    if isequal(lower(cnstit),'auto')
        error(['ut_solv(group): ''cnstit'' cannot be ''auto'' when '... 
            'analyzing groups of records.']);
    end
	
    ind = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'ordercnstit');
                ind = i;
            end
        end
    end
    if ~isempty(ind)
        if ~isequal(varargin{ind+1},'frq')
            error(['ut_solv(group): if ''OrderCnstit'' is passed in '...
                'it must be ''frq''. ']);
        end
        varargin{ind:ind+1} = [];
    end
    ind = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'runtimedisp');
                ind = i;
            end
        end
    end
    if ~isempty(ind)
        if ~isequal(varargin{ind+1},'nnn')
            error(['ut_solv(group): if ''RunTimeDisp'' is passed in '...
                'it must be ''nnn''. ']);
        end
        varargin{ind:ind+1} = [];
    end
	
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'diagnplots');
                 error('ut_solv(group): ''DiagnPlots'' is not allowed.');
            end
        end
    end
	
    infind = [];
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'infer');
                 infind = i;
            end
        end
    end
    if ~isempty(infind)
        infer0 = varargin{infind+1};
        ninf = length(infer0.infnam);
        if ~isempty(vin)
            if ( ~isequal(size(infer0.amprat),[2*ninf 1]) && ...
                    ~isequal(size(infer0.amprat),[2*ninf alln]) ) || ...
                    ( ~isequal(size(infer0.phsoff),[2*ninf 1]) && ...
                    ~isequal(size(infer0.phsoff),[2*ninf alln]) )
                error(['ut_solv(group): infer.amprat and/or '...
                    'infer.phsoff have neither of the acceptable sizes.']);
            end
        else
            if ( ~isequal(size(infer0.amprat),[ninf 1]) && ...
                    ~isequal(size(infer0.amprat),[ninf alln]) ) || ...
                    ( ~isequal(size(infer0.phsoff),[ninf 1]) && ...
                    ~isequal(size(infer0.phsoff),[ninf alln]) )
                error(['ut_solv(group): infer.amprat and/or '...
                    'infer.phsoff have neither of the acceptable sizes.']);
            end
        end
    end
	
    nttst = 0;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'notrend')
                nttst = 1;
            end
        end
    end
    ndtst = 0;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            if isequal(lower(varargin{i}),'nodiagn')
                ndtst = 1;
            end
        end
    end
	
    % reshape group to vectors (from multi-dim arrays) if needed
    if n_d > 1
        if isequal(size(tin),size(uin))
            tin = reshape(tin,n_t,n_s);
        end
        uin = reshape(uin,n_t,n_s);
        if ~isempty(vin)
            vin = reshape(vin,n_t,n_s);
        end
        if ~isempty(infind)
            if ~isempty(vin)
                if ~isequal(size(infer0.amprat),[2*ninf 1])
                    infer0.amprat = reshape(infer0.amprat,2*ninf,n_s);
                    infer0.phsoff = reshape(infer0.phsoff,2*ninf,n_s);
                end
            else
                if ~isequal(size(infer0.amprat),[ninf 1])
                    infer0.amprat = reshape(infer0.amprat,ninf,n_s);
                    infer0.phsoff = reshape(infer0.phsoff,ninf,n_s);
                end
            end
        end
        if ~isequal(size(lat),[1 1])
            lat = reshape(lat,1,n_s);
        end
    end
	
    % initialize storage (as vectors)
    if ~isempty(infind)
        nallc = length(unique([cnstit; infer0.infnam; infer0.refnam;]));
    else
        nallc = length(cnstit);
    end
	
       coef.g = nan*ones(nallc,n_s);
    coef.g_ci = coef.g;
	
    if ~isempty(vin)
        coef.Lsmaj = coef.g;
        coef.Lsmaj_ci = coef.g;
        coef.Lsmin = coef.g;
        coef.Lsmin_ci = coef.g;
        coef.theta = coef.g;
        coef.theta_ci = coef.g;
    else
        coef.A = coef.g;
        coef.A_ci = coef.g;
    end
	
    if ~isempty(vin)
        coef.umean = nan*ones(1,n_s);
        coef.vmean = coef.umean;
        if ~nttst
            coef.uslope = coef.umean;
            coef.vslope = coef.vmean;
        end
    else
        coef.mean = nan*ones(1,n_s);
        if ~nttst
            coef.slope = coef.mean;
        end
    end
	
    coef.results = cell(1,n_s);
    coef.aux.rundescr = cell(1,n_s);
    coef.aux.opt.equi = nan*ones(1,n_s);
    coef.aux.reftime = nan*ones(1,n_s);
    if ~ndtst
       coef.diagn.PE = coef.g;
       coef.diagn.SNR = coef.g;
       coef.diagn.TVraw = nan*ones(1,n_s);
       coef.diagn.TVallc = nan*ones(1,n_s);
       coef.diagn.TVsnrc = nan*ones(1,n_s);
       coef.diagn.PTVallc = nan*ones(1,n_s);
       coef.diagn.PTVsnrc = nan*ones(1,n_s);
       coef.diagn.SNRallc = nan*ones(1,n_s);
       coef.diagn.K = nan*ones(1,n_s);
       coef.diagn.lo.CorMx = coef.g;
       coef.diagn.hi = coef.diagn.lo;
       coef.diagn.table = cell(1,n_s); 
    end
	
    % main loop, over sequences in group
    for is = 1:n_s
	
        % prep inputs to cut_solv1
        uin1 = uin(:,is);
        if ~isempty(vin)
            vin1 = vin(:,is);
        else
            vin1 = [];            
        end
        if isequal(size(uin),size(tin))
            tin1 = tin(:,is);
        else
            tin1 = tin;
        end
        if ~isequal(size(lat),[1 1])
            lat1 = lat(is);
        else
            lat1 = lat;
        end
        if ~isempty(infind)
            infer1 = infer0;
            if ~isempty(vin) 
                if ~isequal(size(infer0.amprat),[2*ninf 1])
                    infer1.amprat = infer0.amprat(:,is);
                    infer1.phsoff = infer0.phsoff(:,is);
                end
            else
                if ~isequal(size(infer0.amprat),[ninf 1])
                    infer1.amprat = infer0.amprat(:,is);
                    infer1.phsoff = infer0.phsoff(:,is);
                end
            end
            varargin{infind+1} = infer1;
        end
		
        % execute cut_solv1
        fprintf(sprintf('[%d/%d]',is,n_s));
        try
            coef1 = cut_solv1(tin1,uin1,vin1,lat1,cnstit,...
                'OrderCnstit','frq','RunTimeDisp','nnn',varargin{:});
        catch err
            disp(err.message);
        end
        % store results from 1st run that apply to all runs
        if is == 1
            coef.name = coef1.name;
            coef.aux.frq = coef1.aux.frq;
            coef.aux.lind = coef1.aux.lind;
            coef.aux.opt.twodim = coef1.aux.opt.twodim;
            coef.aux.opt.notrend = coef1.aux.opt.notrend;
            coef.aux.opt.prefilt = coef1.aux.opt.prefilt;
            coef.aux.opt.nodsatlint = coef1.aux.opt.nodsatlint;
            coef.aux.opt.nodsatnone = coef1.aux.opt.nodsatnone;
            coef.aux.opt.gwchlint = coef1.aux.opt.gwchlint;
            coef.aux.opt.gwchnone = coef1.aux.opt.gwchnone;
            coef.aux.opt.inferaprx = coef1.aux.opt.inferaprx;
            coef.aux.opt.rmin = coef1.aux.opt.rmin;
            coef.aux.opt.method = coef1.aux.opt.method;
            coef.aux.opt.tunconst = coef1.aux.opt.tunconst;
            coef.aux.opt.tunrdn = coef1.aux.opt.tunrdn;
            coef.aux.opt.linci = coef1.aux.opt.linci;
            coef.aux.opt.white = coef1.aux.opt.white;
            coef.aux.opt.nrlzn = coef1.aux.opt.nrlzn;
            coef.aux.opt.lsfrqosmp = coef1.aux.opt.lsfrqosmp;
            coef.aux.opt.nodiagn = coef1.aux.opt.nodiagn;
            coef.aux.opt.diagnplots = coef1.aux.opt.diagnplots;
            coef.aux.opt.diagnminsnr = coef1.aux.opt.diagnminsnr;
            coef.aux.opt.ordercnstit = coef1.aux.opt.ordercnstit;
            coef.aux.opt.runtimedisp = coef1.aux.opt.runtimedisp;
            coef.aux.opt.cnstit = coef1.aux.opt.cnstit;
            if ~ndtst
                coef.diagn.name = coef1.diagn.name;
                coef.diagn.lo.name = coef1.diagn.lo.name;
                coef.diagn.lo.RR = coef1.diagn.lo.RR;
                coef.diagn.lo.RNM = coef1.diagn.lo.RNM;
                coef.diagn.hi.name = coef1.diagn.hi.name;
                coef.diagn.hi.RR = coef1.diagn.hi.RR;
                coef.diagn.hi.RNM = coef1.diagn.hi.RNM;
            end
        end
		
        % store results from this particular (is'th) run 
        if ~isempty(vin)            
            coef.Lsmaj(:,is) = coef1.Lsmaj;
            coef.Lsmaj_ci(:,is) = coef1.Lsmaj_ci;
            coef.Lsmin(:,is) = coef1.Lsmin;
            coef.Lsmin_ci(:,is) = coef1.Lsmin_ci;
            coef.theta(:,is) = coef1.theta;
            coef.theta_ci(:,is) = coef1.theta_ci;
            coef.umean(is) = coef1.umean;
            coef.vmean(is) = coef1.vmean;
            if ~nttst
                coef.uslope(is) = coef1.uslope;
                coef.vslope(is) = coef1.vslope;
            end 
        else
            coef.A(:,is) = coef1.A;
            coef.A_ci(:,is) = coef1.A_ci;
            coef.mean(is) = coef1.mean;
            if ~nttst
                coef.slope(is) = coef1.slope;
            end
        end
        coef.g(:,is) = coef1.g;
        coef.g_ci(:,is) = coef1.g_ci;
        coef.results{is} = coef1.results;
        coef.aux.rundescr{is} = coef1.aux.rundescr;
        coef.aux.reftime(is) = coef1.aux.reftime;
        coef.aux.opt.equi(is) = coef1.aux.opt.equi;
        if ~ndtst
            coef.diagn.PE(:,is) = coef1.diagn.PE;
            coef.diagn.SNR(:,is) = coef1.diagn.SNR;
            coef.diagn.TVraw(is) = coef1.diagn.TVraw;
            coef.diagn.TVallc(is) = coef1.diagn.TVallc;
            coef.diagn.TVsnrc(is) = coef1.diagn.TVsnrc;
            coef.diagn.PTVallc(is) = coef1.diagn.PTVallc;
            coef.diagn.PTVsnrc(is) = coef1.diagn.PTVsnrc;
            coef.diagn.SNRallc(is) = coef1.diagn.SNRallc;
            coef.diagn.K(is) = coef1.diagn.K;
            coef.diagn.lo.CorMx(:,is) = coef1.diagn.lo.CorMx;
            coef.diagn.hi.CorMx(:,is) = coef1.diagn.hi.CorMx;
            coef.diagn.table{is} = coef1.diagn.table;
        end
    end
    if ~isempty(infind)
        coef.aux.opt.infer = infer0;
    else
        coef.aux.opt.infer = [];
    end
    coef.aux.lat = lat;
	
    % reshape back to original shapes (from vectors) if needed
    if n_d > 1
        coef.g = reshape(coef.g,[nallc alln]);
        coef.g_ci = reshape(coef.g_ci,[nallc alln]);
        if ~isempty(vin)
            coef.Lsmaj = reshape(coef.Lsmaj,[nallc alln]);
            coef.Lsmaj_ci = reshape(coef.Lsmaj_ci,[nallc alln]);
            coef.Lsmin = reshape(coef.Lsmin,[nallc alln]);
            coef.Lsmin_ci = reshape(coef.Lsmin_ci,[nallc alln]);
            coef.theta = reshape(coef.theta,[nallc alln]);
            coef.theta_ci = reshape(coef.theta_ci,[nallc alln]);
            coef.umean = reshape(coef.umean,alln);
            coef.vmean = reshape(coef.vmean,alln);
            if ~nttst
                coef.uslope = reshape(coef.uslope,alln);
                coef.vslope = reshape(coef.vslope,alln);
            end
        else
            coef.A = reshape(coef.A,[nallc alln]);
            coef.A_ci = reshape(coef.A_ci,[nallc alln]);
            coef.mean = reshape(coef.mean,alln);
            if ~nttst
                coef.slope = reshape(coef.slope,alln);
            end
        end
        coef.results = reshape(coef.results,alln);
        if ~isempty(infind)
            if ~isempty(vin)
                if ~isequal(size(coef.aux.opt.infer.amprat),[2*ninf 1])
                    coef.aux.opt.infer.amprat = ...
                        reshape(coef.aux.opt.infer.amprat,[2*ninf alln]);
                    coef.aux.opt.infer.phsoff = ...
                        reshape(coef.aux.opt.infer.phsoff,[2*ninf alln]);
                end
            else
                if ~isequal(size(coef.aux.opt.infer.amprat),[ninf 1])
                    coef.aux.opt.infer.amprat = ...
                        reshape(coef.aux.opt.infer.amprat,[ninf alln]);
                    coef.aux.opt.infer.phsoff = ...
                        reshape(coef.aux.opt.infer.phsoff,[ninf alln]);
                end
            end
        end
        if ~isequal(size(lat),[1 1])
            coef.aux.lat = reshape(lat,alln);
        end
        coef.aux.rundescr = reshape(coef.aux.rundescr,alln);
        coef.aux.reftime = reshape(coef.aux.reftime,alln);
        coef.aux.opt.equi = reshape(coef.aux.opt.equi,alln);
        if ~ndtst
            coef.diagn.PE = reshape(coef.diagn.PE,[nallc alln]);
            coef.diagn.SNR = reshape(coef.diagn.SNR,[nallc alln]);
            coef.diagn.TVraw = reshape(coef.diagn.TVraw,alln);
            coef.diagn.TVallc = reshape(coef.diagn.TVallc,alln);
            coef.diagn.TVsnrc = reshape(coef.diagn.TVsnrc,alln);
            coef.diagn.PTVallc = reshape(coef.diagn.PTVallc,alln);
            coef.diagn.PTVsnrc = reshape(coef.diagn.PTVsnrc,alln);
            coef.diagn.SNRallc = reshape(coef.diagn.SNRallc,alln);
            coef.diagn.K = reshape(coef.diagn.K,alln);
            coef.diagn.lo.CorMx = reshape(coef.diagn.lo.CorMx,...
                [nallc alln]);
            coef.diagn.hi.CorMx = reshape(coef.diagn.hi.CorMx,...
                [nallc alln]);
            coef.diagn.table = reshape(coef.diagn.table,alln);
        end
    end
	
    % reorder fields
    if coef.aux.opt.twodim
        fldord = {'name'; 'Lsmaj'; 'Lsmaj_ci';...
            'Lsmin'; 'Lsmin_ci'; 'theta'; 'theta_ci'; 'g'; 'g_ci';...
            'umean'; 'vmean';};
    else
        fldord = {'name'; 'A'; 'A_ci'; 'g'; 'g_ci';'mean';};
    end
	
    if ~coef.aux.opt.notrend
        if coef.aux.opt.twodim
            fldord{end+1} = 'uslope';
            fldord{end+1} = 'vslope';
        else
            fldord{end+1} = 'slope';
        end
    end
	
    fldord{end+1} = 'results';
    fldord{end+1} = 'aux';
    if ~coef.aux.opt.nodiagn
        fldord{end+1} = 'diagn';
    end
            coef = orderfields(coef,fldord);
        coef.aux = orderfields(coef.aux,{'rundescr';'opt';'frq';'lind';'lat';'reftime';});    
    coef.aux.opt = orderfields(coef.aux.opt,{'twodim';'equi';'notrend';...
        'prefilt';'nodsatlint';'nodsatnone';'gwchlint';'gwchnone';...
        'infer';'inferaprx';'rmin';'method';'tunconst';'tunrdn';'linci';...
        'white';'nrlzn';'lsfrqosmp';'nodiagn';'diagnplots';...
        'diagnminsnr';'ordercnstit';'runtimedisp';'cnstit';});
		
    if ~ndtst
        coef.diagn = orderfields(coef.diagn,{'name';'PE';'SNR';'TVraw';...
            'TVallc';'TVsnrc';'PTVallc';'PTVsnrc';'SNRallc';'K';'lo';...
            'hi';'table';});
        coef.diagn.lo = orderfields(coef.diagn.lo,{'name';'RR';'RNM';...
            'CorMx';});
        coef.diagn.hi = orderfields(coef.diagn.hi,{'name';'RR';'RNM';...
            'CorMx';});
    end
end

