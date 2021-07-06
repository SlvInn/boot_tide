%%--------------------------------------------------------- 
function [coef,m,W,Conv,solnstats] = cut_solve_reg(xraw,B,opt,uin)
% S. Innocenti adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% Apply the least square estimation on the model defined by the B basis
% function and with the given opt

        nt = size(xraw,1);
        nu = size(uin,1);
        
% initialize the weight at W = I and lambda = 0
          W = sparse(1:nt,1:nt,1);  % assign unit weights
     coef.W = ones(nu,1);     % moved here (SI)
coef.constr = [0 NaN];        % [Lambda Alpha]; % Added (SI)


% estimate parameters
lastwarn(''); % initialise 
switch opt.method
    
    case 'ols'
              m = B\xraw;
            % m = regress(xraw,B);

            solnstats.w     = [];
            solnstats.resid = [];

    case {'cauchy','andrews','bisquare','fair','huber', 'logistic','talwar','welsch'}
        
        [m,solnstats] = robustfit(B,ctranspose(xraw),opt.method,opt.tunconst,'off');
      
      
        % save weights in the output struct,W~=I:
                W = sparse(1:nt,1:nt,solnstats.w);
           coef.W(~isnan(uin))= solnstats.w; % moved here (SI)
           solnstats.resid = []; % (SI): remove residuals since are computed otside the function ((for all estimators)



    case {'ridge' 'lasso' 'enet' } % added SI
       
        warning(['using ' opt.method ' estimator BUT the method has not been thoroughly tested']) % Added (SI)
        [m, coef.constr] = cut_ridlasenet(xraw,B,opt);
        
        solnstats.w     = [];
        solnstats.resid = [];
        % warning(['solnstats.w set to NaN for ' opt.method ' estimation >> check if right'])
        
      
end

coef.M = m;    
Conv = 1;
      

end
