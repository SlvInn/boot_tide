function [coef,CoefDist] = cut_prepareout(coef,opt,CoefDist)
%% re-order constituents and add the following fileds to CoefDist:
%    CoefDist.M_mc         -  Monte-Carlo replicates of M parameters 
%    CoefDist.A_mc / g_mc  -  Monte-Carlo replicates of amplitude and phases
%    OR 
%    CoefDist.Lsmaj_mc /Lsmin_mc / theta_mc / g_mc  -  Monte-Carlo replicates of 2D amplitude and phases 


% get the constituent order
ind = get_out_cnstit_order(coef,opt);

% prepare output
 coef.name      = coef.name(ind);
 coef.g         = coef.g(ind);
 coef.g_rlzn    = coef.g_rlzn(ind,:); %added by P.Matte
 coef.g_lin_ci  = coef.g_lin_ci(ind);  % modified (SI)
 coef.Std.g_lin = coef.Std.g_lin(ind); % Added (SI)
 coef.Std.g     = coef.Std.g(ind);     % Added (SI)
 
 % Added (SI) (to check if opt.twodim: I don't know if the dimensions are consistent in following lines, ~ll 76 to 88)
           nc = numel(coef.g);
           MS = coef.M(1: nc);
           MC = coef.M(nc+1: 2*nc);
           MT = coef.M(((2*nc)+1):end);
       coef.M = [MS(ind); MC(ind); MT];  
           MS = coef.M_rlzn(1: nc,:);
           MC = coef.M_rlzn(nc+1: 2*nc,:);
           MT = coef.M_rlzn(((2*nc)+1):end,:);
  coef.M_rlzn = [MS(ind,:); MC(ind,:); MT];
           MS = coef.Std.M(1: nc);
           MC = coef.Std.M(nc+1: 2*nc);
           MT = coef.Std.M(((2*nc)+1):end);
   coef.Std.M = [MS(ind,:); MC(ind,:); MT]; 

    % reorder the varcov matrices (to be checked):
                nb = size(CoefDist.Sig.B,1);
    CoefDist.Sig.B = CoefDist.Sig.B([ind;nc+ind;2*nc+1:nb],:);
    CoefDist.Sig.B = CoefDist.Sig.B(:,[ind;nc+ind;2*nc+1:nb]);
        varcov_mCw = CoefDist.Sig.W;
    CoefDist.Sig.W = varcov_mCw(ind,:,:);
    if ~opt.white
        varcov_mCc = CoefDist.Sig.C;
    CoefDist.Sig.C = varcov_mCc(ind,:,:);
    end
<<<<<<< HEAD

=======
>>>>>>> c4e9d54e9c63fa833110e157c9359e3978304faf
% (SI):  check for opt.twodim
% --- %

if opt.twodim
       coef.Lsmaj      = coef.Lsmaj(ind);
       coef.Lsmaj_ci   = coef.Lsmaj_ci(ind);
       coef.Lsmaj_rlzn = coef.Lsmaj_rlzn(ind,:); % Added (SI)
       coef.Std.Lsmaj  = coef.Std.Lsmaj(ind);    % Added (SI)
       
       coef.Lsmin      = coef.Lsmin(ind);
       coef.Lsmin_ci   = coef.Lsmin_ci(ind);
       coef.Lsmin_rlzn = coef.Lsmin_rlzn(ind,:); % Added (SI)
       coef.Std.Lsmin  = coef.Std.Lsmin(ind);    % Added (SI)
       
       coef.theta         = coef.theta(ind);
       coef.theta_lin_ci  = coef.theta_lin_ci(ind);
       coef.theta_rlzn    = coef.theta_rlzn(ind,:); % Added (SI)
       coef.Std.theta_lin = coef.Std.theta_lin(ind); % Added (SI)
       coef.Std.theta     = coef.Std.theta(ind);    % Added (SI)
       
       coef.aux.rpds{1} = coef.aux.rpds{1}(ind,:); % Added (SI)
       coef.aux.rpds{2} = coef.aux.rpds{2}(ind,:); % Added (SI)
else
      coef.A        = coef.A(ind);
      coef.A_rlzn   = coef.A_rlzn(ind,:); %added by P.Matte
      coef.A_ci     = coef.A_ci(ind);
      coef.Std.A    = coef.Std.A(ind); % Added (SI)
      if ~opt.white
      coef.aux.rpds = coef.aux.rpds(ind,:); % Added (SI)
      end
end
 coef.aux.frq  = coef.aux.frq(ind);
 coef.aux.lind = coef.aux.lind(ind);

% fprintf('done.\n');


% If monte-carlo add the simulated values of parameters in the output
% structure CoefDist: 
% 
 if ~opt.linci % linearized
    CoefDist.M_mc = coef.M_rlzn;
    CoefDist.g_mc = coef.g_rlzn;
     if ~opt.twodim
         CoefDist.A_mc = coef.A_rlzn;
     else
         CoefDist.Lsmaj_mc = coef.Lsmaj_rlzn;
         CoefDist.Lsmin_mc = coef.Lsmin_rlzn;
         CoefDist.theta_mc = coef.theta_rlzn;
     end
     coef.Std.S = CoefDist.Sig;
      
 end

end


function ind = get_out_cnstit_order(coef,opt)
% get the constituent order based on the user choice 

    if ~isempty(opt.ordercnstit)
    
        if isequal(opt.ordercnstit,'frq')
           [~,ind] = sort(coef.aux.frq);
           
        elseif isequal(opt.ordercnstit,'snr')
            
                if opt.twodim
                    SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
                else
                    SNR = (coef.A.^2)./((coef.A_ci/1.96).^2);
                end
                [~,ind] = sort(SNR,'descend');
            
    %{        
    %         if ~opt.nodiagn
    %             [~,ind] = sort(coef.diagn.SNR,'descend');
    %              disp(log(coef.diagn.SNR(ind))')
    %         else
    %             if opt.twodim
    %                 SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
    %             else
    %                 SNR = (coef.A.^2)./((coef.A_ci/1.96).^2);
    %             end
    %             [~,ind] = sort(SNR,'descend');
    %            % disp(log(SNR(ind))')
    %         end
    %}    
             
        else        
                [~,ind] = ismember(cellstr(opt.ordercnstit),cellstr(opt.cnstit));
        end
    else
      
           
         if opt.nodiagn
            if opt.twodim
                PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
                PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
            else
                PE = 100*coef.A.^2/sum(coef.A.^2);
            end
         else
             PE = coef.diagn.PE;
         end
        [~,ind] = sort(PE,'descend');
        
           
        % 
        %{
    %     if ~opt.nodiagn
    %           [~,indPE] = sort(coef.diagn.PE,'descend');
    %                 ind = indPE;
    %     else
    %         if opt.twodim
    %             PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
    %             PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
    %         else
    %             PE = 100*coef.A.^2/sum(coef.A.^2);
    %         end
    %        [~,ind] = sort(PE,'descend');
    %     end
        %}
    end
end    