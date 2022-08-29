function [m, constr] = cut_ridlasenet(xraw,B,opt)
% resolve the tidal HA regression using a ridge, lasso, or elastic net
% estimators. WARNING: the code has not been thoroughly tested and must 
% be further commented.

% Add a/some column/s if a trenb is included in the model
if opt.notrend
    Bin = B(:,1:end-1);
else
    % error 'check the code for models with trends'
    Bin = [B(:,1:end-2) B(:,end)];
    warning('check the code for models with trends')
end
 
% set parameters for the numerical optimization 
optfms = optimset('MaxFunEvals',800,'MaxIter',400,'display','off');
         
% define a set of lambas-0 to find the best lambda in the optimization
  el0 = 0.1; %-0.5;
   l0 = 10^(el0);
   l1 = 10^(el0+1);
 Lam0 = unique([l0:l0:l1 l1*2:l1:l1*10 l1*100:l1*100:1 1:10 11:4:137]);

if isempty(opt.regalgo)% % use LC-Curve criterion

       % estimate the models for several lambdas   
         if isempty(opt.alpha)
            error('warning: this part of code is incomple (alpha set = 0.5)') 
            alp = 0.5;

         elseif opt.alpha==0 % if ridge

          m0RRtmp = ridge(xraw,Bin,Lam0,0);
             m0RR = [m0RRtmp(:,2:end), m0RRtmp(:,1)];
             alp  = opt.alpha;

         else % if lasso or elastic net
              alp = opt.alpha;
 [m0RRtmp,interc] = CLasso(Bin,xraw,Lam0,'alpha',alp,'standardize',false);%
             m0RR = [m0RRtmp;interc];
         end

        % apply L-Curve criterion      
        [LC,Lam1] = LCurvePM(B,xraw,m0RR,Lam0,alp);
       [~,ilamRR] = max(LC);

            lamCT = Lam1(ilamRR);
            alpCT = opt.alpha;

         if roundn(lamCT,-3)== roundn(nanmin(Lam1),-3)
            warning('creg:solve:minlambda',[opt.method ':penality: lambda = min allowed smooth parameter'])
         end

 elseif strcmp(opt.regalgo,'svd') % % use LC-Curve criterion with SVD

     % optimize lambda and m at the same time
     % L-Curve criterion with SVD approach 
       [U,S] = svd(B,'econ'); % apply singular value decomp

       if isempty(opt.alpha)  % if elastic net and we want to optimize lambda and alpha 

              Clax = @(LA) CLasso(Bin,xraw,LA(1),'alpha',LA(2),'standardize',false);
            BestBL = @(LA) -LCurveSVD(Bin,U,S,xraw,Clax(LA),LA(1),LA(2));
            LamAlp = fminsearch(BestBL, [rand 0.01+(1-0.01).*rand],optfms);
             lamCT = LamAlp(1);
             alpCT = LamAlp(2);

       else  % if we know alpha   
           if opt.alpha==0 % if ridge optimize lambda and m
               BestBL = @(L) -LCurveSVD(Bin,U,S,xraw,ridge(xraw,Bin,L,0),L,0);

           else % if lasso or elastic net
                 Clax = @(L) CLasso(Bin,xraw,L,'alpha',opt.alpha,'standardize',false);
               BestBL = @(L) -LCurveSVD(Bin,U,S,xraw,Clax(L),L,opt.alpha); 
           end   
           lamCT = fminsearch(BestBL,rand,optfms);
           alpCT = opt.alpha;
       end

 elseif strcmp(opt.regalgo,'bic')  % % use BIC

       % estimate the models for several lambdas  
        if isempty(opt.alpha)
             alpha = 0.01:0.01:0.99;
            nalpha = numel(alpha);
            for aa = 1 : nalpha
                [m0RRtmp,interc] = CLasso(Bin,xraw,Lam0,'alpha',alpha(aa),'standardize',false);%
                            m0RR = [m0RRtmp;interc];
                              LC = LassoBIC(B,xraw,m0RR,Lam0);
          [LBicAlpha(aa),ilamRR] = nanmin(LC);
                      lamCTA(aa) = Lam0(ilamRR);
            end

             [~,ialphRR] = nanmin(LBicAlpha);
                   lamCT = lamCTA(ialphRR);
                   alpCT = alpha(ialphRR);

       elseif opt.alpha==0 % if ridge

          m0RRtmp = ridge(xraw,Bin,Lam0,0);
             m0RR = [m0RRtmp(:,2:end), m0RRtmp(:,1)];   

       else % if lasso ro elastic net but alpha is known 

            [m0RRtmp,interc] = CLasso(Bin,xraw,Lam0,'alpha',opt.alpha,'standardize',false);%
             m0RR = [m0RRtmp;interc];
        end

             % compute BIC
                 LC = LassoBIC(B,xraw,m0RR,Lam0);

         [~,ilamRR] = nanmin(LC);
              lamCT = Lam0(ilamRR);
              alpCT = opt.alpha;

 end
         
% save the value of best lambda in the output struct:
constr  = [lamCT alpCT];

% re-estimate the regression with the best lambda
if alpCT < eps
       m = ridge(xraw,Bin,lamCT,0);
       m = [m(2:end); m(1)];
else
    [m,interc] = CLasso(Bin,xraw,lamCT,'alpha',alpCT,'standardize',false);
             m = [m; interc];
end

end