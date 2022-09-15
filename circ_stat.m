function out  = circ_stat(fname,varargin)
% Container circular statistics functions, mostly adapted from the  
% Circular Statistics Toolbox for Matlab - P. Berens, 2009
% berens@tuebingen.mpg.de - https://www.jstatsoft.org/article/view/v031i10
%
% Major changes: 
% - select only some functions
% - use operations handling/ignoring nan (e.g., nansum, nanmean, etc)
%
% ---- 
% Berens (2009),  CircStat: A MATLAB Toolbox for Circular Statistics, DOI:0.18637/jss.v031.i10
%
% Based on:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
%
% INPUT:
%    fname: function name, see each specific subfunction
% varargin: function args, see each specific subfunction
%
% silvia.innocenti@ec.gc.ca, 08/2022

    switch fname
        case 'mean'
            out = circ_mean(varargin{:});

        case 'mean_ci'
            out = circ_mean_ci_width(varargin{:}); 

        case 'r'
            out = circ_r(varargin{:});  

        case {'var', 'variance'}  
            out = circ_var(varargin{:});     

        case {'med', 'median'}
            out = circ_median(varargin{:});     

        case {'std', 'std0'}
            [out,out0] = circ_std(varargin{:}); 
            if strcmpi(fname,'std0')
                out = out0;
            end
            
        case {'skewness','skewness0'}
            [out,out0] = circ_skewness(varargin{:});   
            if strcmpi(fname,'skewness0')
                out = out0;
            end
            
        case {'kurtosis','kurtosis0'} 
            [out,out0] = circ_kurtosis(varargin{:}); 
            if strcmpi(fname,'kurtosis0')
                out = out0;
            end

        case {'quartiles'}
            out = circ_quartiles(varargin{:});     
           
        case {'corr'}
            out = circ_corr(varargin{:}); 
            
        case {'lincorr','circlincorr', 'clcorr', 'corrcl'}      
            out = circlin_corr(varargin{:}); 
             
        case {'ang2rad', 'ang_to_rad'}
            out = varargin{1} * pi ./180;     

        case {'rad2ang', 'rad_to_ang'}
            out = varargin{1} ./ pi * 180;
    end
end


%% Circ Stat funtions 
function mu = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%
% Computes the mean direction for circular data.
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%



    if nargin < 3
       dim = 1;
    end
    
   

    if nargin < 2 || isempty(w)
        % if no specific weighting has been specified
        % assume no binning has taken place
          w = ones(size(alpha));
    else
        if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
            error('Input dimensions do not match');
        end 
    end

    % compute weighted sum of cos and sin of angles
    r = nansum(w.*exp(1i*alpha),dim);

    % obtain mean by
    mu = angle(r);


end

function t = circ_mean_ci_width(alpha, xi, w, d, dim)
%
% t = circ_mean_ci(alpha, xi, w, d, dim)
%   Computes the confidence limits on the mean for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [xi   (1-xi)-confidence limits are computed, default 0.05]
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%   Output:
%     t     mean +- t yields upper/lower (1-xi)% confidence limit
%
%


    if nargin < 5
        dim = 1;
    end

    if nargin < 4 || isempty(d)
        % per default do not apply correct for binned data
        d = 0;
    end

    if nargin < 3 || isempty(w)
        % if no specific weighting has been specified
        % assume no binning has taken place
        w = ones(size(alpha));
    else
        if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
            error('Input dimensions do not match');
        end 
    end

    % set confidence limit size to default
    if nargin < 2 || isempty(xi)
        xi = 0.05;
    end

    % compute ingredients for conf. lim.
    r  = circ_r(alpha,w,d,dim);
    n  = nansum(w,dim);
    R  = n.*r;
    c2 = chi2inv((1-xi),1);

    % check for resultant vector length and select appropriate formula
    t = zeros(size(r));

    for i = 1:numel(r)
        if r(i) < .9 && r(i) > sqrt(c2/2/n(i))
            t(i) = sqrt((2*n(i)*(2*R(i)^2-n(i)*c2))/(4*n(i)-c2));  % equ. 26.24
        elseif r(i) >= .9
            t(i) = sqrt(n(i)^2-(n(i)^2-R(i)^2)*exp(c2/n(i)));      % equ. 26.25
        else 
            t(i) = NaN;
            warning('Requirements for confidence levels not met.');
        end
    end

    % apply final transform
    t = acos(t./R);
end

function S = circ_var(alpha, w, d, dim)
%
% s = circ_var(alpha, w, d, dim)
%   
% Computes circular variance for circular data (equ. 26.17/18, Zar).   
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_var(alpha, [], [], dim)
%
%   Output:
%     S     circular variance 1-r


    if nargin < 4
        dim = 1;
    end

    if nargin < 3 || isempty(d)
        % per default do not apply correct for binned data
        d = 0;
    end

    if nargin < 2 || isempty(w)
        % if no specific weighting has been specified
        % assume no binning has taken place
        w = ones(size(alpha));
    else
        if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
            error('Input dimensions do not match');
        end 
    end

    % compute mean resultant vector length
    r = circ_r(alpha,w,d,dim);

    % apply transformation to var
    S = 1 - r;
end

function [S, s0] = circ_std(alpha,varargin)
%
% s = circ_std(alpha, w, d, dim)
%   
% Computes circular standard dev for circular data (equ. 26.17/18, Zar).   
%
%   Input:
%     alpha	sample of angles in radians
%
%   {varargin}:
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_var(alpha, [], [], dim)
%
%   Output:
%     S     circular std 2*(1-r)

    r = 1 - circ_var(alpha, varargin{:});
    
    S  = sqrt(2*(1-r));      % 26.20
    s0 = sqrt(-2*log(r));    % 26.21

end

function med = circ_median(alpha,dim)
%
% med = circ_median(alpha)
%
% Computes the median direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [dim  compute along this dimension, default is 1, must 
%           be either 1 or 2 for circ_median]
%
%   Output:
%     mu		mean direction
%
%   circ_median can be slow for large datasets
%

    if nargin < 2
        dim = 1;
    end

    Beta = alpha;  
      
    if dim == 2
        Beta = alpha';
    elseif all(dim~=[1,2])
        error('circ_med only works along first two dimensions')
    end




      M = size(alpha);
    med = NaN(M(2),1);
    for i = 1 : M(2)
        
        beta = mod(Beta(:,i),2*pi);
           n = size(beta,1);
          dd = circ_dist3(beta,beta); %circ_dist2(beta,beta);
        
        
        m1 = nansum(dd>=0,1);
        m2 = nansum(dd<=0,1);
        dm = abs(m1-m2);
    
    
        if mod(n,2)==1
            [m, idx] = nanmin(dm);
        else
              m = nanmin(dm);
            idx = find(dm==m,2);
        end  

        if m > 1
            warning('circ_med:ties','Ties detected.') %#ok<WNTAG>
        end

        % compute circ mean
        md = angle(nansum(exp(1i*beta(idx))));
    
        % compute circ dist between md and other angles
         mb = angle(nansum(exp(1i*beta))); 
         cd = angle(exp(1i*mb)./exp(1i*md));
        cd2 = angle(exp(1i*mb)./exp(1i*(md+pi)));
    

        if abs(cd) > abs(cd2)
            md = mod(md+pi,2*pi);
        end
    
        med(i) = md;
    end

    if dim == 2
        med = med';
    end
end

function [q25,q50,q75] = circ_quartiles(x)
%
% [q25,q50,q75] = circ_quartiles(alpha)
%
% Compute pseudoquartiles as in Innocenti et al. (2019)
% DOI:10.1029/2019JD031210
%
%   Input:
%     alpha	sample of angles in radians
%
%   Output:
%   [q25,q50,q75]   1st pseudoquartile, circumar median, 2nd pseudoquartile
%
%   circ_quartiles can be slow for large datasets
%

% compute the circ median
    x = mod(x(:),2*pi);
   me = circ_med(x);
  
   % divide the points in two grousp according to the median defintiion
    me1 = mod(me,2*pi);
    me2 = me1 + pi;
   me12 = sort([me1 me2]);

   up =  me12(1) <= x & x < me12(2);
   do =  x <= me12(1) |  me12(2) < x;
   
    % compute the medians of the two grousp
    dy1 = x(up);
    if ~isempty(dy1)
        q25 = mod(circ_med(dy1),2.*pi);
    else
        q25 = NaN;
    end

    dy2 = x(do);
    if ~isempty(dy2)
        q75 = mod(circ_med(dy2),2.*pi);
    else
        q75 = NaN;
    end
    
    q50 = mod(me,2.*pi);
        
end

function [b, b0] = circ_skewness(alpha, w, dim)
%
% [b b0] = circ_skewness(alpha, w, dim)
%
% Computes two measures of angular skewness.
%
%
%   Input:
%     alpha     sample of angles
%     [w        weightings in case of binned angle data]
%     [dim      statistic computed along this dimension, 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_skewness(alpha, [], dim)
%
%   Output:
%     b         skewness (from Pewsey)
%     b0        alternative skewness measure (from Fisher, formula 2.29 in Berens 2009)
%

    
    if nargin < 3
      dim = 1;
    end
    
    if nargin < 2 || isempty(w)
       % if no specific weighting has been specified
       % assume no binning has taken place
         w = ones(size(alpha));
    else
       if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
          error('Input dimensions do not match');
       end 
    end
    
    
    % compute neccessary values
               R = circ_r(alpha,w,[],dim);
           theta = circ_mean(alpha,w,dim);
    [~,rho2,mu2] = circ_moment(alpha,w,2,true,dim);
    
    % compute skewness 
    theta2 = repmat(theta, size(alpha)./size(theta));
        b  = nansum(w.*(sin(2*(circ_dist(alpha,theta2)))),dim)./nansum(w,dim);
        b0 = rho2.*sin(circ_dist(mu2,2*theta))./(1-R).^(3/2);  % (formula 2.29)
end    
    
function [k, k0] = circ_kurtosis(alpha, w, dim)
%
% [k k0] = circ_kurtosis(alpha, w, dim)
%
% Computes two measures of angular kurtosis.
%
%
%   Input:
%     alpha     sample of angles
%     [w        weightings in case of binned angle data]
%     [dim      statistic computed along this dimension, 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_kurtosis(alpha, [], dim)
%
%   Output:
%     k         kurtosis (from Pewsey)
%     k0        kurtosis (from Fisher, formula 2.30 of Berens 2009)
    
    if nargin < 3
        dim = 1;
    end
    
    if nargin < 2 || isempty(w)
       % if no specific weighting has been specified
       % assume no binning has taken place
         w = ones(size(alpha));
    else
         if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
            error('Input dimensions do not match');
         end 
    end
    
    % compute mean direction
              R = circ_r(alpha,w,[],dim);
          theta = circ_mean(alpha,w,dim);
      [~, rho2] = circ_moment(alpha,w,2,true,dim);
    [~, ~, mu2] = circ_moment(alpha,w,2,false,dim);
    
    % compute skewness 
    theta2 = repmat(theta, size(alpha)./size(theta));
         k = nansum(w.*(cos(2*(circ_dist(alpha,theta2)))),dim)./nansum(w,dim);
        k0 = (rho2.*cos(circ_dist(mu2,2*theta))-R.^4)./(1-R).^2;    % (formula 2.30)

end    
   
function rho = circ_corr(alpha1, alpha2)
%
% rho = circ_corr(alpha1, alpha2)
%
% Compute circ correlation coefficients for two circular random variables.
%
%   Input:
%     alpha1	sample of angles in radians
%     alpha2	sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     

    if size(alpha1,2) > size(alpha1,1)
        alpha1 = alpha1';
    end

    if size(alpha2,2) > size(alpha2,1)
        alpha2 = alpha2';
    end

    if length(alpha1)~=length(alpha2)
      error('Input dimensions do not match.')
    end

    % compute mean directions
    n = length(alpha1);
    alpha1_bar = circ_mean(alpha1);
    alpha2_bar = circ_mean(alpha2);

    % compute correlation coeffcient from p. 176
    num = nansum(sin(alpha1 - alpha1_bar) .* sin(alpha2 - alpha2_bar));
    den = sqrt(nansum(sin(alpha1 - alpha1_bar).^2) .* nansum(sin(alpha2 - alpha2_bar).^2));
    rho = num / den;	

    % compute pvalue
    l20 = nanmean(sin(alpha1 - alpha1_bar).^2);
    l02 = nanmean(sin(alpha2 - alpha2_bar).^2);
    l22 = nanmean((sin(alpha1 - alpha1_bar).^2) .* (sin(alpha2 - alpha2_bar).^2));

        ts = sqrt((n * l20 * l02)/l22) * rho;
    % pval = 2 * (1 - normcdf(abs(ts)));

end

function rho = circlin_corr(alpha, x)
%
% rho  = circlin_corr(alpha, x)
%
%   Correlation coefficient between one circular and one linear random variable.
%
%   Input:
%     alpha   sample of angles in radians
%     x       sample of linear random variable
%
%   Output:
%     rho     correlation coefficient


    if size(alpha,2) > size(alpha,1)
        alpha = alpha';
    end

    if size(x,2) > size(x,1)
        x = x';
    end

    if length(alpha)~=length(x)
        error('Input dimensions do not match.')
    end

    n = length(alpha);

    % compute correlation coefficent for sin and cos independently
    rxs = corr(x,sin(alpha),'rows','complete');
    rxc = corr(x,cos(alpha),'rows','complete');
    rcs = corr(sin(alpha),cos(alpha),'rows','complete');

    % compute angular-linear correlation (equ. 27.47)
    rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));

end



%% INTERNAL CIRC FUNCTIONS
function r = circ_r(alpha, w, d, dim)
%
% r = circ_r(alpha, w, d)
%
% Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied 
%           correction factor is used to correct for bias in 
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
    if nargin < 4
        dim = 1;
    end

    if nargin < 2 || isempty(w) 
        % if no specific weighting has been specified
        % assume no binning has taken place
        w = ones(size(alpha));
    else
        if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
            error('Input dimensions do not match');
        end 
    end

    if nargin < 3 || isempty(d)
        % per default do not apply correct for binned data
        d = 0;
    end

    % compute weighted sum of cos and sin of angles
    r = nansum(w.*exp(1i*alpha),dim);

    % obtain length 
    r = abs(r)./nansum(w,dim);

    % for data with known spacing, apply correction factor to correct for bias
    % in the estimation of r (see Zar, p. 601, equ. 26.16)
    if d ~= 0
        c = d/2/sin(d/2);
        r = c*r;
    end
end    

function r = circ_dist(x,y)
%
% r = circ_dist3(alpha, beta)
%
% Computes all the x_i-y_i pairwise differences x_i-y_j around the circle.
%
%   Input:
%     alpha      sample of linear random variable
%     beta       sample of linear random variable or one single angle
%
%   Output:
%     r       matrix with differences

    if size(x,1)~=size(y,1) && size(x,2)~=size(y,2) && length(y)~=1
       error('Input dimensions do not match.')
    end

    r = angle(exp(1i*x)./exp(1i*y));
end

function r = circ_dist3(x,y)
%
% r = circ_dist3(alpha, beta)
%
% Computes all the pairwise differences x_i-y_j around the circle.
%
%   Input:
%     alpha      sample of linear random variable
%     beta       sample of linear random variable
%
%   Output:
%     r       matrix with pairwise differences
%


  x = x(:);
  y = y(:);
  y = y';

  nx = numel(x);
  ny = numel(y);
  r  = nan(nx,ny);
  

    for i = 1 : nx
        ri = exp(1i*x(i))./exp(1i*y);
        r(i,:) = angle(ri);
    end
end 

function [mp, rho_p, mu_p] = circ_moment(alpha, w, p, cent, dim)
%
% [mp cbar sbar] = circ_moment(alpha, w, p, cent, dim)
%
% Computes the complex p-th centred or non-centred moment 
%   of the angular data in angle.
%
%   Input:
%     alpha     sample of angles
%     [w        weightings in case of binned angle data]
%     [p        p-th moment to be computed, default is p=1]
%     [cent     if true, central moments are computed, default = false]
%     [dim      compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_moment(alpha, [], [], [], dim)
%
%   Output:
%     mp        complex p-th moment
%     rho_p     magnitude of the p-th moment
%     mu_p      angle of th p-th moment


    if nargin < 5
      dim = 1;
    end

    if nargin < 4
      cent = false;
    end

    if nargin < 3 || isempty(p)
        p = 1;
    end

    if nargin < 2 || isempty(w)
      % if no specific weighting has been specified
      % assume no binning has taken place
        w = ones(size(alpha));
    else
      if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
        error('Input dimensions do not match');
      end 
    end


    if cent
      theta = circ_mean(alpha,w,dim);
          v = size(alpha)./size(theta);
      alpha = circ_dist(alpha,repmat(theta,v));
    end
  

   n = size(alpha,dim);
cbar = nansum(cos(p*alpha).*w,dim)/n;
sbar = nansum(sin(p*alpha).*w,dim)/n;
  mp = cbar + 1i*sbar;

rho_p = abs(mp);
 mu_p = angle(mp);
end


