function B = cut_get_B_boot(coef,tin,t)
% return the basis function B according to the option specified 
% in coeff and the time vector t (time vector without nan)
opt = coef.aux.opt;


% t is the time vector without nan
nt = numel(t);
     

% define the LOR for equispaced or irregularly sampled series        
if opt.equi == 1 % based on times; u/v can still have nans ("gappy")
    lor = (max(tin)-min(tin)); 
else
    lor = (max(t) - min(t));   
end
    

% options to compute the basis function
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone coef.aux.opt.gwchlint coef.aux.opt.gwchnone];

% compute the model basis function 
[E,F,U,V] = cut_E(t,coef.aux.reftime,coef.aux.frq,coef.aux.lind,coef.aux.lat,ngflgs,coef.aux.opt);
        
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
    B = [B ones(nt,1) (t-coef.aux.reftime)/lor];   
end
    

