function out  = circ_stat(fname,varargin)
% Container circular statistics functions. Most functions are from the  
% Circular Statistics Toolbox for Matlab - P. Berens, 2009
% berens@tuebingen.mpg.de - https://www.jstatsoft.org/article/view/v031i10
%
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

        case {'std'}
            out = circ_std(varargin{:});     
            
        case {'skewness'}
            [out,~] = circ_skewness(varargin{:});      
            
        case {'kurtosis'} 
            [out,~] = kurtosis(varargin{:});      

        case {'circ_quartiles'}
            out = circ_quartiles(varargin{:});     

        case {'ang2rad', 'ang_to_rad'}
            out = varargin{1} * pi ./180;     

        case {'rad2ang', 'rad_to_ang'}
            out = varargin{1} ./ pi * 180;
    end
end



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
    r = sum(w.*exp(1i*alpha),dim);

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
    n  = sum(w,dim);
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


function S = circ_std(alpha, w, d, dim)
%
% s = circ_std(alpha, w, d, dim)
%   
% Computes circular standard dev for circular data (equ. 26.17/18, Zar).   
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
%     S     circular std 2*(1-r)
    S = 2 * circ_var(alpha, w, d, dim);
end


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
    r = sum(w.*exp(1i*alpha),dim);

    % obtain length 
    r = abs(r)./sum(w,dim);

    % for data with known spacing, apply correction factor to correct for bias
    % in the estimation of r (see Zar, p. 601, equ. 26.16)
    if d ~= 0
        c = d/2/sin(d/2);
        r = c*r;
    end
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
        
        
        m1 = sum(dd>=0,1);
        m2 = sum(dd<=0,1);

        dm = abs(m1-m2);
    
    
        if mod(n,2)==1
            [m, idx] = min(dm);
        else
            m = min(dm);
            idx = find(dm==m,2);
        end  

        if m > 1
            warning('circ_med:ties','Ties detected.') %#ok<WNTAG>
        end

        % compute circ mean
        md = angle(sum(exp(1i*beta(idx))));
    
        % compute circ dist between md and other angles
        mb = angle(sum(exp(1i*beta))); 
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


function [b b0] = circ_skewness(alpha, w, dim)
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
    [~, rho2 mu2] = circ_moment(alpha,w,2,true,dim);
    
    % compute skewness 
    theta2 = repmat(theta, size(alpha)./size(theta));
    b  = sum(w.*(sin(2*(circ_dist(alpha,theta2)))),dim)./sum(w,dim);
    b0 = rho2.*sin(circ_dist(mu2,2*theta))./(1-R).^(3/2);    % (formula 2.29)
end    
    


function [k k0] = circ_kurtosis(alpha, w, dim)
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
    k = sum(w.*(cos(2*(circ_dist(alpha,theta2)))),dim)./sum(w,dim);
    k0 = (rho2.*cos(circ_dist(mu2,2*theta))-R.^4)./(1-R).^2;    % (formula 2.30)

end    
    
    