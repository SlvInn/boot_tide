%%--------------------------------------------------------- 
function [coef_boot,B] = cut_boot1(tin,uin,vin,coef,n_boot,varargin)
% S. Innocenti 
% 2020/04/25
%
% Bootstrap for a single record. 
% See comments and help for cut_boot(), cut_renconstr(), and cut_solv().
%
% TO DO (warning messages have been set up for most of these issues):
% - finish controls on user-defined inputs
% - test the functions for 2D analyses; 
% - 

% heritate estimation options
caux = coef.aux;
opt  = coef.aux.opt;

% Bootstrap options >> Set defaults and read the ser-defined inputs
bopt = cut_bootinit(coef,varargin)
print(bopt)

% get the HA reconstructions and residuals
if isempty(bopt.Yhat) % if we have no Yhat calculated
    % reconstruct the original series based on the estiamted parameters (i.e. regression fit)
    [uhat,vhat] = cut_reconstr1(tin,coef,'cnstit',coef.name);

    % check series dimension
    if opt.twodim   
        Yhat = [uhat,vhat];
        Yres = [uhat-uin,vhat-vin];
        warning('the code must be checked and tested for the two dimentional case')   
    else 
        Yhat = uhat;
        Yres = uhat-uin;   
    end

    % remore NaNs
    if opt.twodim
        uvgd = ~isnan(uin) & ~isnan(vin);
    else
        uvgd = ~isnan(uin);
    end


else % if they are provided by the user

    Yhat = bopt.Yhat;
    Yres = bopt.Yres;  

    % remore NaNs
    uvdg = ~any(isnan(Yhat),2)
    
end


% time vector without nan:
t  = tin(uvgd);    
nt = numel(t);

% define the LOR for equispaced or irregularly sampled series        
if opt.equi == 1 % based on times; u/v can still have nans ("gappy")
    lor = (max(tin)-min(tin)); 
else
    lor = (max(t) - min(t));   
end
    

%% ------ Prepare matrix for the bootstrap regressions ------- %%
         
% options to compute the basis function
ngflgs = [opt.nodsatlint opt.nodsatnone opt.gwchlint opt.gwchnone];

% compute the model basis function 
[E,F,U,V] = cut_E(t,caux.reftime,caux.frq,caux.lind,caux.lat,ngflgs,opt);
     
if ~opt.complexf % if sin/cos formulation 
 B = [F.*cos(2*pi*(U+V)), F.*sin(2*pi*(U+V))];
 
 if ~isempty(opt.infer)
     error 'not done yet: inference must be adapted to the sin/cos definition of the model >> try with complexf option = true'
     % see the other part of the IF statement (when opt.infer is ~empty)
 else
     Q = [];
  beta = [];
    nR = 0;
   nNR = numel(coef.name);
 end
else % if complex formulation 
 B = [E conj(E)];
 
 if ~isempty(opt.infer)
     error 'not done yet: inference must be adapted to the sin/cos definition of the model >> try with complexf option = true'
 else
     Q = [];
 beta = [];
     nR = 0;
     nNR = numel(coef.name);
 end
 
end

% include (a) column(s) for temporal trends
if opt.notrend
 B = [B ones(nt,1)];
else
 B = [B ones(nt,1) (t-caux.reftime)/lor];   
end   
% ------------------------------------- %

% construct resample and re-estimate the parameters on each one
coef_boot = cut_bootcore(coef,tin,uvgd,B,Yres,Yhat,bopt,n_boot)

% from M to A,g, mean, and trend:
coef_boot = cut_from_mboot_to_param(coef_boot,n_boot,opt,nNR,nR);


% output
coef_boot.boot_opt = bopt;     % bootstrap options
coef_boot.boot_opt.aux = caux; % HA model setup
end


            



