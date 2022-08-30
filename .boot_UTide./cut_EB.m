function [B,nm,Q,beta,NodCorr] = cut_EB(t,tref,cnstit,opt,lat,nt,nNR,nR,lor)
% S. Innocenti, adapted from ut_solv1() [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% preparation of the basis function matrix


% fprintf('matrix prep ... '); % (SI): commented

     ngflgs = [opt.nodsatlint opt.nodsatnone opt.gwchlint opt.gwchnone];
  [E,F,U,V] = cut_E(t,tref,cnstit.NR.frq,cnstit.NR.lind,lat,ngflgs,opt); % Modified (SI)
   NodCorr.F = F;
   NodCorr.U = U;
   NodCorr.V = V;



  if ~opt.complexf 
        B = [F.*cos(2*pi*(U+V)), F.*sin(2*pi*(U+V))];
        if ~isempty(opt.infer)
            error 'not done yet: inference must be adapted to the sin/cos definition of the model >> try with complexf option = true'
            % see the other part of the IF statement (when opt.infer is ~empty)
        else
              Q = [];
           beta = [];
        end
  else % complex formulation 
        B = [E conj(E)];
        
        if ~isempty(opt.infer)
            Etilp = nan*ones(nt,nR);
            Etilm = Etilp;

            if ~opt.inferaprx
                for k = 1:nR
                             E = ut_E(t,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,ngflgs,opt.prefilt);
                             Q = ut_E(t,tref,cnstit.R{k}.I.frq,cnstit.R{k}.I.lind,lat,ngflgs,opt.prefilt) ./ repmat(E,1,cnstit.R{k}.nI);
                    Etilp(:,k) = E.*(1+sum(Q.*cnstit.R{k}.I.Rp(ones(nt,1),:),2));
                    Etilm(:,k) = E.*(1+sum(Q.*conj(cnstit.R{k}.I.Rm(ones(nt,1),:)),2));
                end
            else
                   Q = nan*ones(1,nR);
                beta = Q;
                for k = 1:nR
                             E = ut_E(t,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,ngflgs,opt.prefilt);
                    Etilp(:,k) = E;
                    Etilm(:,k) = E;
                          Q(k) = ut_E(tref,tref,cnstit.R{k}.I.frq,cnstit.R{k}.I.lind,lat,ngflgs,opt.prefilt)  ./ ...
                                 ut_E(tref,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,ngflgs,opt.prefilt);
                           arg = pi*lor*24*(cnstit.R{k}.I.frq-cnstit.R{k}.frq)*(nt+1)/nt;
                       beta(k) = sin(arg)/arg;
                end    
            end

               B = [B Etilp conj(Etilm)];

        else
            Q = [];
            beta = [];
        end
        
  end
   
% (SI): moved inside the IF  
%{
% if ~isempty(opt.infer)
%     Etilp = nan*ones(nt,nR);
%     Etilm = Etilp;
%     
%     if ~opt.inferaprx
%         for k = 1:nR
%                      E = ut_E(t,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,ngflgs,opt.prefilt);
%                      Q = ut_E(t,tref,cnstit.R{k}.I.frq,cnstit.R{k}.I.lind,lat,ngflgs,opt.prefilt) ./ repmat(E,1,cnstit.R{k}.nI);
%             Etilp(:,k) = E.*(1+sum(Q.*cnstit.R{k}.I.Rp(ones(nt,1),:),2));
%             Etilm(:,k) = E.*(1+sum(Q.*conj(cnstit.R{k}.I.Rm(ones(nt,1),:)),2));
%         end
%     else
%            Q = nan*ones(1,nR);
%         beta = Q;
%         for k = 1:nR
%                      E = ut_E(t,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,ngflgs,opt.prefilt);
%             Etilp(:,k) = E;
%             Etilm(:,k) = E;
%                   Q(k) = ut_E(tref,tref,cnstit.R{k}.I.frq,cnstit.R{k}.I.lind,lat,ngflgs,opt.prefilt)  ./ ...
%                          ut_E(tref,tref,cnstit.R{k}.frq,cnstit.R{k}.lind,lat,ngflgs,opt.prefilt);
%                    arg = pi*lor*24*(cnstit.R{k}.I.frq-cnstit.R{k}.frq)*(nt+1)/nt;
%                beta(k) = sin(arg)/arg;
%         end    
%     end
%    
%        B = [B Etilp conj(Etilm)];
%     
% else
%     Q = [];
%     beta = [];
% end
%}

if opt.notrend
    B = [B ones(nt,1)];
   nm = 2*(nNR + nR) + 1;
else
    B = [B ones(nt,1) (t-tref)/lor];
   nm = 2*(nNR + nR) + 2;
end

end