function [u,v,opt,B] =  cut_boot_reconstr(tin,coef_boot ,varargin)

% remove nan
      t = tin(:);
   ntin = numel(tin(:));
    whr = ~isnan(tin);
      t = t(whr);
     

 % checks
 assert(isfield(coef_boot,'M'),'the input coef strcuture must have a field ''M'' with parameter valuess')
 nM    = size(coef_boot.M,1);
 nBoot = size(coef_boot.M,2);
 nC    = size(coef_boot.g,1);
 
% defaults
       aux = coef_boot.boot_opt.aux;
opt.cnstit = coef_boot.name;
opt.B      = [];


% Optional arguments
optargin = size(varargin,2);

i = 1;
while i <= optargin
    switch lower(varargin{i})
        case 'cnstit'
            opt.cnstit = args{2};    
            args(1:2) = [];
            
       case {'basis' 'b' 'expbasis'}  % complex exponential basis function 
            opt.B = args{2};    
            args(1:2) = [];
            
    otherwise 
            error(['cut_boot_reconstr: unrecognized input: ' args{1}]);        
    end
end

% select constituents
[~,ind] = intersect(coef_boot.name,cellstr(opt.cnstit),'stable');
assert(~isempty(ind),'None of the input constituents is in coef_boot.name ')    
    
% complex exponential basis function 
if ~isempty(opt.B)
    
    assert(size(opt.B,2)==nM, 'dimension of the input basis function and coef_boot.M matrices are not consistent')
    assert(size(opt.B,1)==nt, 'dimension of the input basis function and time vector are not consistent')
                      B = opt.B;
   
else
    
    
     ngflgs = [aux.opt.nodsatlint aux.opt.nodsatnone aux.opt.gwchlint aux.opt.gwchnone];
  [E,F,U,V] = cut_E(t,aux.reftime,aux.frq(ind),aux.lind(ind),aux.lat,ngflgs,aux.opt);
         
 
  
  if ~ aux.opt.complexf 
      
        B = [F.*cos(2*pi*(U+V)), F.*sin(2*pi*(U+V))];
         
         
  else % complex formulation 
        B = [E conj(E)];
         
  end

    
end

% reconstructions
   u = nan*ones(ntin,nBoot);
if aux.opt.twodim
   v = u;
else
   v = []; 
end

% for each bootstrap
  iap = ind;
  iam = nC+ind; 
for b = 1 : nBoot

        
         ap = coef_boot.M(iap,b);
         am = coef_boot.M(iam,b);
         
         
% %          disp(size(B))
% %          disp(size(ap))
        fit = B(:,iap)*ap +  B(:,iam)*am;

    % mean (& trend)
    if aux.opt.twodim

        if aux.opt.notrend
            u(whr,b) = real(fit) + coef_boot.umean(b);
            v(whr,b) = imag(fit) + coef_boot.vmean(b);
        else
            u(whr,b) = real(fit) + coef_boot.umean(b) + coef_boot.uslope(b)*(t-aux.reftime);
            v(whr,b) = imag(fit) + coef_boot.vmean(b) + coef_boot.vslope(b)*(t-aux.reftime);
        end
    else
        if aux.opt.notrend
% %             disp(size(u(whr,b)))
% %             disp(size(real(fit)))
            u(whr,b) = real(fit) + coef_boot.mean(b);
        else
            u(whr,b) = real(fit) + coef_boot.mean(b) + coef_boot.slope(b)*(t- aux.reftime);
        end


    end

end


end