function [m_theta, s_theta,  ci_theta,  opt] = boot_tide_param(theta,varargin) 
% Compute the plug-in bootstrap estimators of the tidal model parameter replicates, 
% allowing for the use of circular statistics for (some) parameters (i.e. theta rows). 
%
%     [m_theta]  = boot_tide_param(theta, {'option_name',option_value}) 
%     [m_theta,s_theta]  = boot_tide_param(theta, {'option_name',option_value})
%     [m_theta,s_theta,ci_theta]  = boot_tide_param(theta, {'option_name',option_value})
%     [m_theta,s_theta,ci_theta,opt]  = boot_tide_param(theta, {'option_name',option_value})
%
% INPUT:
% theta - (P x nboot x S) matrix of the nboot replicates of tide model parameters.
%
%
% VARARGIN >> {'option_name',option_value}: 
% circular   - (P x 1) vector indicating which theta components are angles (e.g., phases)
%               Default: [false for each theta].
%
% circ_units - 'radians' or 'degrees' indicating the units of the circular components of theta. 
%               Default: 'radians.
%
% ci         - string indicating the method for compunting the Confidence Interval (CI) of non-circular theta.
%              Possible values: 
%                   * 'percentile'  - percentile bootstrap CI: [theta_alpha/2, theta_1-alpha/2]
%                   * 'gaussian'    - gaussian approximation of CIs: [mean_theta +- z_alpha/2 * sigma_theta]
%                   * 'studentized' - Studentized bootstrap CI: ...
%                   * 'accelerated' - Accelerated bootstrap CI (not ready yet)
%                   * 'bca'         - Bias-Correct Accelerated bootstrap CI (not ready yet)   
%              Default: 'percentile'.  
%              
% alpha      - confidence level for constructing CI. Default: 0.05
% 
% theta0     - (P x S) matrix of the original tide model parameters (i.e. parameters computed on 
%              the original data sample, not the bootstrap replicates).  theta0 is mandatory for  
%              computing studentized bootstrap CIs, and ignored by all other CI methods.
% 
% OUTPUT: 
% m_theta - (P x S) matrix of the parameter mean overt the nboot replicates.
%           circular mean calculated from circ_mean.m are returned for circular parameters.
%
% {OPTIONAL OUTPUT/ VARARGOUT}:
% >> [m_theta,s_theta]  = boot_tide_param(theta)
% s_theta - (P x S) matrix of the parameter standard errors over the nboot replicates.
%            circular variance calculated from circ_var.m are returned for circular parameters.
%
% >> [m_theta,s_theta,ci_theta]  = boot_tide_param(theta)
% ci_theta - (P x 2 x S) matrix of the parameter confidence interval bounds.
%
% >> [m_theta,s_theta,ci_theta,opt]  = boot_tide_param(theta)
% opt - Structures containing the option values used in the boot_tide_param.m call
%
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2020/04 - 2022/09
%

    % make theta being a 3D matrix
    thsh = size(theta);
    assert(length(thsh)>=2, 'theta must be a matrix with bootstrap parameter replicates')
    if length(thsh)<3
        theta =  reshape(theta,[thsh, 1]);
    end
    


    % set the options based on default and user-defined values
    opt = set_options('boot_param',varargin{:});
    opt = check_boot_param_options(theta,opt);

    % find the index of the linear- and circular-scale components of theta
    lin = find(opt.circular==0);
    cir = find(opt.circular==1);
    

    % express the circular theta components in radians
    if strcmp(opt.circ_units,'degrees')
        theta(cir,:) = circ_stat('ang2rad',theta(cir,:));
    end

    % % % compute plug-in estimators % % %

    % mean 
    m_theta = mean_theta(theta,lin,cir); 


    % variance
    if nargout > 1
        s_theta = std_theta(theta,lin,cir);
    end

    % % % compute confidence intervals for the mean % % %
    if nargout > 2
        ci_theta = mean_theta_ci(theta,m_theta,s_theta,lin,cir,opt);
    end

    % put the computed circular quantities back in degrees
    if strcmp(opt.circ_units,'degrees')
        m_theta(cir,:) = circ_stat('rad2ang',m_theta(cir,:));
        if nargout>2
            ci_theta(cir,:) = circ_stat('rad2ang',ci_theta(cir,:));
        end
    end

end



function m_theta =  mean_theta(theta,lin,cir)

    % empty
    [np,~,ns] = size(theta);
    m_theta = nan(np,ns);

    
    
    % compute the mean theta for linear-scale parameters
    m_theta(lin,:) = nanmean(theta(lin,:,:) ,2);

    % compute the mean theta for circular-scale parameters
    for cc = 1 : length(cir)
        c = cir(cc);
        for s = 1 : ns
            thcs = theta(c,:,s);
            m_theta(c,s) = circ_stat('mean',thcs(:)); 
        end
    end

end



function s_theta = std_theta(theta,lin,cir)

    % empty
    [np,~,ns] = size(theta);
    s_theta = nan(np,ns);


    % compute the mean theta for linear-scale parameters
    s_theta(lin,:,:) = sqrt(nanvar(theta(lin,:,:),[],2));

    % compute the mean theta for circular-scale parameters
    for cc = 1 : length(cir)
        c = cir(cc);
    
        for s = 1 : ns
            thcs = theta(c,:,s);
            s_theta(c,s) = circ_stat('var',thcs(:)); 
        end
    end

   
end


function ci_theta =  mean_theta_ci(theta,m_theta,s_theta,lin,cir,opt)

    % empty CI matrix
    [np,nb,ns] = size( theta);
    ci_theta  = nan(np,2,ns);

    % compute the mean theta for linear-scale parameters
    a2 = opt.alpha/2;
    switch opt.ci
        case 'gaussian' 

            z = norminv(1-a2);
            ci_theta(lin,1,:) = m_theta(lin,:,:) - z.*sqrt(s_theta(lin,:,:)); % lower CI bound
            ci_theta(lin,2,:) = m_theta(lin,:,:) + z.*sqrt(s_theta(lin,:,:)); % upper CI bound

        case 'percentile'

            ci_theta(lin,1,:) = quantile(theta(lin,:,:),a2,2);   % lower CI bound
            ci_theta(lin,2,:) = quantile(theta(lin,:,:),1-a2,2); % upper CI bound


        case 'studentized'
            
            std_thetab = s_theta(lin,:,:)/sqrt(nb);
            stud_theta = (theta(lin,:,:)-m_theta(lin,:,:))./std_thetab ;
             
            ci_theta(lin,1,:) = opt.theta0(lin) - std_thetab.*quantile(stud_theta,1-a2,2);   % lower CI bound
            ci_theta(lin,2,:) = opt.theta0(lin) - std_thetab.*quantile(stud_theta,a2,2); % upper CI bound

        otherwise 
            disp('accelerated and bca methods not implemented yet')
            
    end    

    
    % compute the mean theta for circular-scale parameters
    for cc = 1 : length(cir)
        c = cir(cc);
     
        for s = 1 : ns
            thcs = theta(c,:,s);
            w = circ_stat('mean_ci',thcs(:),a2); 
            ci_theta(c,1,s) = m_theta(c,s) - w; % lower CI bound
            ci_theta(c,2,s) = m_theta(c,s) + w; % upper CI bound
        end
    end


end



function opt = check_boot_param_options(theta,opt)

    % check the length of the circular vector
    [np,~,ns] = size(theta);
    if ~isempty(opt.circular)
        opt.circular = opt.circular(:); %make circular  being a col vector
        assert(np==numel(opt.circular),'inconsistent theta and theta_circ dimensions')
    else
        opt.circular = zeros(np,1); % set false for each param in theta
    end

    % string indicating the method for computing CI
    legal_ci = {'percentile', 'gaussian', 'studentized', 'accelerated', 'bca'};
    opt.ci   = lower(opt.ci);
    assert(any(strcmp(opt.ci,legal_ci)), ['Not a valid ci, possible types are ' legal_ci{:}])

    % tmp error
    if ~any(strcmp(opt.ci,{'percentile', 'gaussian','studentized'}))
        error('tmp: accelerated and bca methods not implemented yet') 
    end

    % if studentized bootstrap, check for the availability of 
    % the original theta value (theta computed on the original series, not resamples)
    if strcmp(opt.ci,'studentized')
        [np0,ns0] = size(opt.theta0);
        assert(np0*ns0>0,'theta0 must be provided for studentized bootstrap CIs')
        assert(np==np0 && ns==ns0 , 'inconsistent theta and theta0 dimensions')
        ncirc = sum( opt.circular);
        if ncirc >0
            warning(['studentized bootstrap not applied for the ' num2str(ncirc) 'circular parameters'])
        end
    end

    % confidence level for constructing CI
    assert(opt.alpha>0 && opt.alpha<1,'alpha must be in (0,1)')

    % check the circular units
    legal_circ_units = {'rad', 'radians', 'deg', 'degrees'};
    opt.circ_units   = lower(opt.circ_units );
    assert(any(strcmp(opt.circ_units,legal_circ_units)), ['Not a valid circ_units, possible units are ' legal_circ_units{:}])

    if strcmp(opt.circ_units,'deg')
        opt.circ_units = 'degrees';
    end
    
end