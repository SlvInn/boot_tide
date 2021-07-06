function [coef,varMSM,varcov_mCw,varcov_mCc,Gall,Hall,CoefDist] = cut_ci(xraw,m,B,W,coef,opt,Xu,Yu,Xv,Yv,Puu,Puv,Pvv,nm,nNR,nR)
% S. Innocenti adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% 2019/09
% NEW OUTPUT: 
% coef.Std  - structure in coef containing the vectors of the standard error 
%             estimated for each parameter
%
% CoefDist  - structure containing the var-cov estimates and MC simulations
%             and the Monte-Carlo simulations for complex parameters, as well 
%             as for the amplitudes and phases.
% 
%    CoefDist.Sig.Se -  UTide estimation of the scalar residual standard error (sigma epsilon)
%    CoefDist.Sig.W  -  white    component of the UTide var-cov matrix 
%    CoefDist.Sig.C  -  coloured component of the UTide var-cov matrix
%    
% 
%    CoefDist.M_mc         -  Monte-Carlo replicates of M parameters 
%    CoefDist.A_mc / g_mc  -  Monte-Carlo replicates of amplitude and phases
%    OR 
%    CoefDist.Lsmaj_mc /Lsmin_mc / theta_mc / g_mc  -  Monte-Carlo replicates of 2D amplitude and phases 

    nt = size(xraw,1);
    nc = nNR+nR;
    
%% NR & R cnstits: confidence intervals
varMSM = real((ctranspose(xraw)*W*xraw - ctranspose(m)*ctranspose(B)*W*xraw)/(nt-nm)); % sigma-e 

  gamC = inv(ctranspose(B)*W*B)*varMSM; %#ok
  gamP = 0.*gamC; % Added (SI): if we use the sin-cos model formulation the pseudo-cov must be 0
  if opt.complexf~=0 % if complex formulation - Added (SI)
     gamP = inv(transpose(B)*W*B)*((transpose(xraw)*W*xraw - transpose(m)*transpose(B)*W*xraw)/(nt-nm)); %#ok
  end
  
  Gall = gamC+gamP; % sum of covariance and pseudo-covariance
  Hall = gamC-gamP; % diff
   
            
% inits
  coef.g_ci   = nan*ones(size(coef.g));
  coef.g_rlzn = nan(size(coef.g,1),opt.nrlzn); % Added (SI)

if opt.twodim
    coef.Lsmaj_ci = coef.g_ci;
    coef.Lsmin_ci = coef.g_ci;
    coef.theta_ci = coef.g_ci;
       varcov_mCw = nan*ones(nc,4,4);
else
      coef.A_ci   = coef.g_ci;
      coef.A_rlzn = coef.g_rlzn; % Added (SI)
       varcov_mCw = nan*ones(nc,2,2);
end

if ~opt.white
      varcov_mCc = varcov_mCw;
else
      varcov_mCc = [];
end

  % % Added (SI):
  % define the  CoefDist structure to contain the output for uncertainty
  % analysis (var-cov estimates and MC simulations)
   CoefDist.Sig.Se = varMSM; % scalar residual standard error (sigma epsilon)
   CoefDist.Sig.W  = [];     % white component of the var-cov
   CoefDist.Sig.C  = [];     % coloured component of the var-cov
   CoefDist.M_mc   = nan(size(coef.M,1),opt.nrlzn); % produced MC simulations
       coef.M_rlzn = [];     % same as CoefDist.M_mc but in the coef output

 
% main loop on constituents
for c = 1:nc % for each constituent
    
        G = [Gall(c,c) Gall(c,c+nc); Gall(c+nc,c) Gall(c+nc,c+nc);];
        H = [Hall(c,c) Hall(c,c+nc); Hall(c+nc,c) Hall(c+nc,c+nc);];
    varXu = real(G(1,1)+G(2,2)+2*G(1,2))/2;
    varYu = real(H(1,1)+H(2,2)-2*H(1,2))/2;
    
    if opt.twodim
        varXv = real(H(1,1)+H(2,2)+2*H(1,2))/2;
        varYv = real(G(1,1)+G(2,2)-2*G(1,2))/2;
    end
    
   
    
    if opt.linci % linearized
        
        
        [coef,varcov_mCw, varcov_mCc] = cut_linear_ci(opt,varcov_mCw, varcov_mCc,varXu,varYu, varXv, varYv);
        
    else % monte carlo
        
                      covXuYu = imag(H(1,1)-H(1,2)+H(2,1)-H(2,2))/2;
                          Duu = [varXu covXuYu; covXuYu varYu;];
        varcov_mCw(c,1:2,1:2) = Duu;
        
        
        
        % if colored noise
        if ~opt.white

                              Duu = Puu(c)*Duu/trace(Duu);
            varcov_mCc(c,1:2,1:2) = Duu;
 
        end
        
        % if one dimension
        if ~opt.twodim
           
             muC = [Xu(c) Yu(c)];
             
            if ~opt.white % if colored noise
                 varcov_mCc(c,:,:) = ut_nearposdef(squeeze(varcov_mCc(c,:,:)));
                              sigC = squeeze(varcov_mCc(c,:,:));
                             mCall = mvnrnd(muC,sigC,opt.nrlzn);
                             
                            
            else % if withe noise
               % varcov_mCw(c,:,:) = ut_nearposdef(varcov_mCw(c,:,:));
                             sigC = squeeze(varcov_mCw(c,:,:));
                             sigC = ut_nearposdef(sigC); 
                            mCall = mvnrnd(muC,sigC,opt.nrlzn);    
            end
            

             coef.Std.M([c c+nc],1)  = diag(sigC); % Added (SI) 2020/04/17
             
             CoefDist.M_mc(c,:)    = mCall(:,1); % Added (SI)
             CoefDist.M_mc(nc+c,:) = mCall(:,2); % Added (SI)

            % Compute coeff in MC and save values
            % then compute the half-length of the CI  
                   [A,~,~,g] = ut_cs2cep(mCall); % compute current ellipse parameters from cosine-sine coefficients
            coef.A_rlzn(c,:) = A; % added (PM)
                         mad = median(abs(A-median(A))); 
             coef.Std.A(c,1) = mad/0.6745;   % Added (SI)
               coef.A_ci(c)  = 1.96*coef.Std.A(c);
             
            
                        g(1) = coef.g(c);
                           g = ut_cluster(g,360);
            coef.g_rlzn(c,:) = g; %added (PM)
                         mad = median(abs(g-median(g)));  % LINEAR meadian! no sense (SI)
             coef.Std.g(c,1) = mad/0.6745; % Added (SI)   
               coef.g_ci(c)  = 1.96*coef.Std.g(c);
            
               % (SI): The way A_ci and g_ci are computed should be changed: 
               % - variable significance level
               % - asymmetric interval 
               % - no linear inteval for g (angular variable) 
               % => significance (or SNR) must be computed for the M coefficients  
               % - (The factor 1/0.6745 = 1.48 makes this an approximately
               % unbiased estimate of scale when the error model is
               % Gaussian) [Holland and welsch 1997] => evident from the definition of the estimator!
               
              
              
        else  % if two dimensions 
                          covXvYv = imag(G(1,1)-G(1,2)+G(2,1)-G(2,2))/2;
                              Dvv = [varXv covXvYv; covXvYv varYv;];
            varcov_mCw(c,3:4,3:4) = Dvv;
            
            if ~opt.white
                Dvv = Pvv(c)*Dvv/trace(Dvv);
                varcov_mCc(c,3:4,3:4) = Dvv;
            end
            
                          covXuXv = imag(-H(1,1)-H(1,2)-H(2,1)-H(2,2))/2;
                          covXuYv = real(G(1,1)-G(2,2))/2;
                          covYuXv = real(-H(1,1)+H(2,2))/2;
                          covYuYv = imag(-G(1,1)+G(1,2)+G(2,1)-G(2,2))/2;
                              Duv = [covXuXv covXuYv; covYuXv covYuYv;];
            varcov_mCw(c,1:2,3:4) = Duv;
            varcov_mCw(c,3:4,1:2) = transpose(Duv);
            
            if ~opt.white
                if sum(abs(Duv(:)))>0
                    Duv = Puv(c)*Duv/sum(abs(Duv(:)));
                    varcov_mCc(c,1:2,3:4) = Duv;
                    varcov_mCc(c,3:4,1:2) = transpose(Duv);
                else
                    varcov_mCc(c,1:2,3:4) = 0;
                    varcov_mCc(c,3:4,1:2) = 0;
                end
                varcov_mCc(c,:,:) = ut_nearposdef(squeeze(varcov_mCc(c,:,:)));
                sigMc = squeeze(varcov_mCc(c,:,:)); % Added (SI)
                mCall = mvnrnd([Xu(c) Yu(c) Xv(c) Yv(c)],sigMc,opt.nrlzn);
            else
                sigMc = squeeze(varcov_mCw(c,:,:)); % Added (SI)
                mCall = mvnrnd([Xu(c) Yu(c) Xv(c) Yv(c)],sigMc,opt.nrlzn);
            end
            
            % Added (SI)
          
            sigMc = diag(sigMc); % Added (SI) 2020/04/17 >>  to be tested for 2D series
            for jj = 1 : 4 
                CoefDist.M_mc((jj-1)*nc + c,:) = mCall(:,jj ); % Added (SI)
                coef.Std.M((jj-1)*nc + c,1)  = sigMc(jj); % Added (SI) 2020/04/17 >>  to be tested for 2D series
            end
            [Lsmaj,Lsmin,theta,g] = ut_cs2cep(mCall);
            
            
            
                              mad = median(abs(Lsmaj-median(Lsmaj))); 
              coef.Std.Lsmaj(c,1) = mad/0.6745; % Added (SI)   
                 coef.Lsmaj_ci(c) = 1.96*coef.Std.Lsmaj(c);
             coef.Lsmaj_rlzn(c,:) = Lsmaj; % Added (SI)
                 
                              mad = median(abs(Lsmin-median(Lsmin)));
              coef.Std.Lsmin(c,1) = mad/0.6745; % Added (SI)   
                 coef.Lsmin_ci(c) = 1.96*coef.Std.Lsmin(c);
             coef.Lsmin_rlzn(c,:) = Lsmin; % Added (SI)
        
                         theta(1) = coef.theta(c);
                            theta = ut_cluster(theta,360);
                            
                              mad = median(abs(theta-median(theta)));
              coef.Std.theta(c,1) = mad/0.6745;  % Added (SI)   
                coef.theta_ci(c)  = 1.96*coef.Std.theta(c);
             coef.theta_rlzn(c,:) = theta; % Added (SI)
                 
                             g(1) = coef.g(c);
                                g = ut_cluster(g,360);
                              mad = median(abs(g-median(g)));  
                  coef.Std.g(c,1) = mad/0.6745;  % Added (SI)   
                  coef.g_ci(c)    = 1.96*coef.Std.g(c);
                 coef.g_rlzn(c,:) = g; % Added (SI)
        end
    end
    
end

 % create an output structure CoefDist for estimated variances and param distributions
 CoefDist.Sig.W = varcov_mCw;
 CoefDist.Sig.C = varcov_mCc;
    coef.M_rlzn = CoefDist.M_mc;

 % If monte-carlo add the simulated values of parameters in the output
 % structure CoefDist: % % Added (SI)
 if ~opt.linci % linearized
         CoefDist.g_mc = coef.g_rlzn;
     if ~opt.twodim
         CoefDist.A_mc = coef.A_rlzn;
     else
         CoefDist.Lsmaj_mc = coef.Lsmaj_rlzn;
         CoefDist.Lsmin_mc = coef.Lsmin_rlzn;
         CoefDist.theta_mc = coef.theta_rlzn;
     end
 end

end

