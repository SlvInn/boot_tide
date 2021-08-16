function [coef, CoefDist,NodCorr] = cut_solv1(tin,uin,vin,lat,cnstit,varargin)
% S. Innocenti adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% 2019/09
%
% Tidal analysis of a single record. See comments for cut_solv().
%
% fprintf('ut_solv: '); % (SI): Commented
%
%% input checking/conditioning, options parsing
[nt,t,u,v,tref,lor,elor,opt,tgd,uvgd] = cut_slvinit(tin,uin,vin,cnstit,varargin);

%% constituent selection, initial storage
             opt.cnstit = cnstit;
[nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref,opt.rmin/(24*lor),opt.cnstit,opt.infer);
      coef.aux.rundescr = ut_rundescr(opt,nNR,nR,nI,t,tgd,uvgd,lat);
          coef.aux.opt  = opt;
          coef.aux.lat  = lat;
		  coef.aux.rpds = []; % initialize Residual Power Spectral Density (band averaged). Added (SI)

%% basis function matrix
[B,nm,Q,beta,NodCorr] = cut_EB(t,tref,cnstit,opt,lat,nt,nNR,nR,lor);


%% HA regression solution
xraw = u;
if opt.twodim
   xraw = complex(u,v);
end

% estimate regression parameters:
[coefreg,m,W,Conv,solnstats] = cut_solve_reg(xraw,B,opt,uin); 
solnstats.resid = [];
solnstats.M     = coefreg.M;
solnstats.B     = B;

         % add new fields to the structure coef
            fcreg = fieldnames(coefreg);
          [~,icr] = setxor(fcreg,fieldnames(coef));
         for ficr = 1 : length(icr)
             coef.(fcreg{ficr}) = coefreg.(fcreg{ficr});
         end

         % if the optimization did not converge
         if Conv==0
            % nan-fill, create coef.results, reorder coef fields, do runtime display
             coef = cut_finish(coef,nNR,nR,nI,elor,cnstit); % abort remainder of calcs
            return;        
         end
     
% Xhat: 
xmod = B*m;
if ~opt.twodim
   xmod = real(xmod);
end

%% Residuals
solnstats.res = (xraw-xmod); % Added (SI): unweighted residuals, to be used!
e = W*(xraw-xmod); % weighted Residuals used by UTide


if opt.complexf
    
    % For NR & R cnstits: from complex coeffs to current ellipse &  Ampl/pha params
    [coef,Xu,Yu,Xv,Yv] = cut_from_compx_to_amppha(coef,m,nNR,nR,Q,beta,cnstit,lor,opt);
    
else
  
   % For NR & R cnstits: from cos/sin coeffs to current ellipse &  Ampl/pha params
   [coef,Xu,Yu] = cut_from_sincos_to_amppha(coef,m,nNR,nR,Q,beta,cnstit,lor,opt);
             Xv = [];
             Yv = []; 
end

  % spectral power (Puu, Pvv, Puv) of residual as in UTide
    if ~opt.white
       [coef,Puu,Puv,Pvv] = cut_spectpow(coef,opt,e,t,tin,tgd,uvgd,elor);
    else
        Puu = [];
        Puv = [];
        Pvv = [];
    end

    % NR & R cnstits: confidence intervals
    [coef,varMSM,varcov_mCw,varcov_mCc,Gall,Hall,coefdiststr] = ... % Modified (SI) by adding ciefdiststr output
          cut_ci(xraw,m,B,W,coef,opt,Xu,Yu,Xv,Yv,Puu,Puv,Pvv,nm,nNR,nR);
      
      % create an output structure for the estimated var-cov and param distrib
       CoefDist = coefdiststr; % Added (SI)
       CoefDist.IRLSstats = solnstats;
     

    % Infer cnstits: complex & cos/sin coeffs, c.e.p., & conf ints
    if ~isempty(opt.infer)
        coef = cut_infer_and_ci(coef,varcov_mCw,varcov_mCc,nNR,nI);   
        warning(['CoefDist does not include the var-cov information on inferred components. ' ... % Added (SI)
                 'It must be modified accordingly']) % Added (SI)
    end
    
    % diagnostics
    if ~opt.nodiagn
         coef = cut_disgnostic(coef,opt,e,cnstit,t,u,v,xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,tin,uin,vin);
    end

%% re-order constituents and prepare output
[coef,CoefDist] = cut_prepareout(coef,opt,CoefDist);

%% create coef.results, reorder coef fields, do runtime display
coef = cut_finish(coef,nNR,nR,nI,elor,cnstit);

end
