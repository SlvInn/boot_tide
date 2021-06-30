function [coefr,covb] = NS_Tide_Cornoise(X,coef,xres,n_rep)
% Correlated noise model. Matte et al. (2013), 10.1175/JTECH-D-12-00016.1, Section 2.e
% Function isolated from the ns_tide package (Version 1.03, experimental, P. Matte,  06/2012)
% Modified by S. Innocenti (2020) to use unweighted residuals and regressor matrix. 
% 
% INPUT:
% X     - Txp matrix: basis functions used as HA regressor. 
%         For instance, for conventional stationary HA models, 
%         p=2K with K= num. of tidal constituents (no intercept column)
% coef  - px1 vector of HA regression coefficients WITHOUT the intercept
% xres  - Tx1 vector of regression unweighted residuals. 
%         Note: IRLS weight should not be multiplied to the residuals         
% n_rep - Number of Monte-Carlo repetitions to be generated
% 
% OUTPUT:
% coefr - n_rep x p  matrix containing the random MC reaplicates of coef
% covb  - p x p regression variance-covariance matrix (i.e. coef var matrix) 

 % get matrix dimensions
  [T,p] =  size(X);

% check for intercept 
if p <length(coef) 
    X = [ones(length(xres),1), X];
    p = p+1;
end

% locate missings/NaNs
  dum = sum(X,2) + xres; 
  Xgd = ~isnan(dum); % non-NaNs for both X and xres
  
 
  
  % Modify ns_tide to use unweighted residuals, S. Innocenti, 2020
   Pf = fft(X(Xgd,:));  % (SI): in ns_tide.m Pf is computed on W.*X but here there are no weights 
 Perr = fft(xres(Xgd)); % (SI): in ns_tide.m Perr is computed on W.*xres but here there are no weights 
 Terr = length(Perr);


% Daniell avg the spectrum (boxcar avg)
  B   = ones(16,1) ;
  B   = B./sum(B) ;
  msk = ones(Terr ,1) ;

  
% apply a running mean of 16 values
  Perr = conv2(abs(Perr).^2,B,'same')./conv2(msk,B,'same'); %for older Matlab versions
 %Perr = conv(abs(Perr).^2,B,'same')./conv(msk,B,'same');   %for newer Matlab versions

% clone the Perr vector p times 
PerrRep = Perr(:,ones(1,p));  
invPerr = (1./PerrRep); % residual var-cov

% inv of the regression var-cov
covmat = Pf'*(invPerr.*Pf);
covmat = (covmat + covmat')/2; % symmetrizing the matrix


% SIGMA = inv(Pf'*(diag(1./Perr)*Pf));
  SIGMA = real(inv(covmat));
     MU = coef(:);
 
  
% (SI) : the following lines should be eliminated since 
% MU and SIGMA must have consistent sizes without this correction.
  covb = SIGMA; % get a SIGMA copy to be sure to return the var-cov uncorrected for size
if length(SIGMA)~=length(MU)
    warning('Length of SIGMA is different than length of MU');
    for i=1:length(MU)-length(SIGMA)
        SIGMA(end+1,end+1) = 0;
    end
end

% Simulate n_rep random realizations of M:
% coefr = mvnrnd(MU,SIGMA,n_rep-1);
coefr = mvnrnd(MU,SIGMA,n_rep);

end