function [u,v] = cut_reconstr(tin,coef,varargin)
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2019/09 - 2020/12
% Adapted from ut_reconstr() - UTide v1p0 2011/09, http://www.po.gso.uri.edu/~codiga/utide/utide.htm. 
% 
% TO DO (warning messages have been set up for most of these issues):
% - complete the help with a complete description of the options;
% - complete the comments in the sub-function scripts.
%
% OVERVIEW:
% Generate hindcasts and/or forecasts/predictions at 
% user-specified times using the HA coefficient structure computed by CUT_SOLV(). 
% Specific steps of the reconstruction are:
% - select a list (e.g., a subset) of constituents from those estimated by 
%   cut_solv.m and compute the corresponding model basis function 
% - reconstruct the regression prediction based the selected list of constituents
%
% Syntax for two-dimensional raw input, such as velocities:
%     [ u_fit, v_fit ] = cut_reconstr ( t_fit, coef , {options} );
%
% Syntax for one-dimensional raw input, such as sea level:
%     [ sl_fit, ~ ] = cut_reconstr ( t_fit, coef , {options} ); 
%
%     Analysis of groups of records (TO BE TESTED)
% The description below considers INPUTS, DEFAULTS, OPTIONS, and OUTPUT for 
% treating a single record analysis. Modifications for treating a group of records 
% with a single execution are given in Codiga (2011) and UTide help. 
%
% INPUT: 
%     * t_fit = Column vector of arbitrary times [datenum UCT/GMT] distributed uniformly 
%               or irregularly. Hindcasts are produced at these time steps. If t_fit include 
%               NaNs, the outputs (u_fit/v_fit or sl_fit) will put NaNs at the corresponding time steps. 
%     * coef  = results output structure from CUT_SOLV()
%
% {OPTIONS}:
% Option flags are not case-sensitive but cannot be abbreviated. The option order is not relevant. Please refer to ut_solv.m help and Codiga (2011) 
% for complete explanations of UTide options. 
%     * ‘MinSNR’, MinSNR (num) to select constituents with SNR >= MinSNR. 
%        Default MinSNR=2.
%     * ‘MinPE’, MinPE (num) to select constituents with PE > MinPE. 
%        Default MinPE=0. 
%
%        If both ‘MinSNR’ and ‘MinPE’ are specified, no constituent with either 
%        SNR or PE values lower than the defined thresholds will be used.
%
%     * ‘Cnstit’, Cnstit (cell list of 4-character strings) list of 
%       constituents, which must be used for regerssion hindcasts. The list of 
%       names must include all or a subset of constituents in the coef.name 
%       field of the  cut_solv() output structure. If ‘Cnstit’ is specified, the 
%       MinSNR and MinPE options are ignored.
%
% DEFAULTS: 
%     * Constituents with SNR > 2 are included in the superposition. 
%     * All options used for the HA regression are read from the coef structure, output of cut_solv().
%
% OUTPUT:
%     * t_fit    = Tx1 vector of times [datenum UTC] defining the dates at which hindcasts must be produced 
%     * u_fit & v_fit, or sl_fit =  Tx1 vector(s) of reconstructed series/hindcasts. 
%
% GROUPS OF RECORDS (TEMPORARILY ELIMINATED)
%
% REFERENCES:
% Codiga, 2011 - http://www.po.gso.uri.edu/~codiga/utide/utide.htm



if isequal(size(coef.g,2),1) % single record 
    [u,v] = cut_reconstr1(tin,coef,varargin{:});
    
else % group of records
    n_t = size(tin,1); 
    % alln = [n_1 n_2 n_3 ... n_n_d]
    alln = size(coef.g);
    alln = alln(2:end);
    % check inputs
    if (~isequal(size(tin),[n_t 1]) && ~isequal(size(tin),[n_t alln]))
        error(['cut_reconstr(group): input times have neither acceptable size.']);
    end
    % n_d = number of dimensions of group
    n_d = length(alln);
    % n_s = number of time sequences in group
    n_s = prod(alln); 
    % nallc = number of all constituents
    nallc = size(coef.g,1);
    % reshape t to a vector
    if size(tin,2)>1
        tin = reshape(tin,n_t,n_s);
    end
    % reshape fields of coef to vectors
    if n_d > 1
        coef.g = reshape(coef.g,nallc,n_s);
        coef.g_ci = reshape(coef.g_ci,nallc,n_s);
        if coef.aux.opt.twodim 
            coef.Lsmaj = reshape(coef.Lsmaj,nallc,n_s);
            coef.Lsmaj_ci = reshape(coef.Lsmaj_ci,nallc,n_s);
            coef.Lsmin = reshape(coef.Lsmin,nallc,n_s);
            coef.Lsmin_ci = reshape(coef.Lsmin_ci,nallc,n_s);
            coef.theta = reshape(coef.theta,nallc,n_s);
            coef.theta_ci = reshape(coef.theta_ci,nallc,n_s);
            coef.umean = reshape(coef.umean,1,n_s);
            coef.vmean = reshape(coef.vmean,1,n_s);
            if ~coef.aux.opt.notrend
                coef.uslope = reshape(coef.uslope,1,n_s);
                coef.vslope = reshape(coef.vslope,1,n_s);
            end
        else
            coef.A    = reshape(coef.A,nallc,n_s);
            coef.A_ci = reshape(coef.A_ci,nallc,n_s);
            coef.mean = reshape(coef.mean,1,n_s);
            if ~coef.aux.opt.notrend
                coef.slope = reshape(coef.slope,1,n_s);
            end
        end
        coef.aux.reftime = reshape(coef.aux.reftime,1,n_s);
        if size(coef.aux.lat,2)>1
            coef.aux.lat = reshape(coef.aux.lat,1,n_s);
        end
    end
    % initialize storage
    u = nan*ones(n_t,n_s);
    if coef.aux.opt.twodim
        v = u;
    else
        v = [];
    end
    % main loop
    coef1.aux.frq = coef.aux.frq;
    coef1.aux.lind = coef.aux.lind;
    coef1.aux.opt = coef.aux.opt;
    for is = 1:n_s % for each time sequence
        if size(tin,2) > 1
            tin1 = tin(:,is);
        else
            tin1 = tin;
        end
        % select coef results for the is'th record
        coef1.name = coef.name;
        coef1.g    = coef.g(:,is);
        coef1.g_ci = coef.g_ci(:,is);
        if coef.aux.opt.twodim
            coef1.Lsmaj = coef.Lsmaj(:,is);
            coef1.Lsmaj_ci = coef.Lsmaj_ci(:,is);
            coef1.Lsmin = coef.Lsmin(:,is);
            coef1.Lsmin_ci = coef.Lsmin_ci(:,is);
            coef1.theta = coef.theta(:,is);
            coef1.theta_ci = coef.theta_ci(:,is);
            coef1.umean = coef.umean(is);
            coef1.vmean = coef.vmean(is);
            if ~coef.aux.opt.notrend
                coef1.uslope = coef.uslope(is);
                coef1.vslope = coef.vslope(is);
            end
        else
            coef1.A = coef.A(:,is);
            coef1.A_ci = coef.A_ci(:,is);
            coef1.mean = coef.mean(is);
            if ~coef.aux.opt.notrend
                coef1.slope = coef.slope(is);
            end
        end
        if size(coef.aux.lat,2)>1
            coef1.aux.lat = coef.aux.lat(is);
        else
            coef1.aux.lat = coef.aux.lat;
        end
        coef1.aux.reftime = coef.aux.reftime(is);
        % execute reconstruct for one record
        fprintf('[%d/%d]',is,n_s);
        [u1,v1] = cut_reconstr1(tin1,coef1,varargin{:});
        % store results
        u(:,is) = u1;
        if coef.aux.opt.twodim
            v(:,is) = v1;
        end
    end
    % reshape back to original array sizes
    if n_d > 1
        u = reshape(u,[n_t alln]);
        if coef.aux.opt.twodim
            v = reshape(v,[n_t alln]);
        end
    end
end

