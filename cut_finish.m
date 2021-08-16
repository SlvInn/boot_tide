
%%--------------------------------------------------------- 
function coef = cut_finish(coef,nNR,nR,nI,elor,cnstit)
% S. Innocenti adapted from ut_finish [UTide v1p0 9/2011 d.codiga@gso.uri.edu]
% Finalize the output of the cut_solv() function: 
% - replace parameter values and teh corresponding confidence intervals by Nans 
%   if IRLS did not converge (only for the portion that uses input other than coef)
% - create the coef.results char array
% - order the coef structure fields 
% - display coef.results, coef.aux.rundescr, and/or coef.diagn.table as specified by opt.RunTimeDisp
% 

%% Optimization did not converge, so nan-fill coefficients and conf-ints
if ~isfield(coef,'g')
    if coef.aux.opt.twodim
        
               % coef.M = NaN;  % Added (SI)
           coef.Lsmaj = nan*ones(size(coef.name));
        coef.Lsmaj_ci = coef.Lsmaj;
           coef.Lsmin = coef.Lsmaj;
        coef.Lsmin_ci = coef.Lsmaj;
           coef.theta = coef.Lsmaj;
        coef.theta_lin_ci = coef.Lsmaj;
               coef.g = coef.Lsmaj;
            coef.g_lin_ci = coef.Lsmaj;
           coef.umean = NaN;
           coef.vmean = NaN;
           
          coef.Std.Lsmaj  = NaN; % Added (SI)
          coef.Std.Lsmin  = NaN; % Added (SI)
          coef.Std.theta  = NaN; % Added (SI)
          coef.Std.g      = NaN; % Added (SI)
           
          
          coef.Lsmaj_rlzn = NaN;  % Added (SI)
          coef.Lsmin_rlzn = NaN;  % Added (SI)
          coef.theta_rlzn = NaN;  % Added (SI)
          coef.g_rlzn     = NaN;  % Added (SI)
          coef.M_rlzn     = NaN;  % Added (SI)
          coef.W          = NaN;  % coef.W*0+nan; % Added (SI)
           
        if ~coef.aux.opt.notrend
            coef.uslope = NaN;
            coef.vslope = NaN;
        end
    else
           %coef.M = nan*ones(size(coef.M));  % Added (SI)
           coef.A = nan*ones(size(coef.name));
        coef.A_ci = coef.A;
           coef.g = coef.A;
        coef.g_lin_ci = coef.A;
        coef.mean = NaN;
        
          coef.Std.A  = NaN;
          coef.Std.g  = NaN;

          coef.A_rlzn = NaN; % Added (SI)
          coef.g_rlzn = NaN; % Added (SI)
          coef.M_rlzn = NaN;  % Added (SI)
          coef.W      = NaN; %coef.W*0+nan; % Added (SI)
          
        if ~coef.aux.opt.notrend
            coef.slope = NaN;
        end
    end
    
    if ~coef.aux.opt.nodiagn
                     nallc = nNR+nR+nI;
                        nc = nNR+nR;
           coef.diagn.name = coef.name;
             coef.diagn.PE = nan*ones(nallc,1);
            coef.diagn.SNR = coef.diagn.PE;
          coef.diagn.TVraw = nan;
         coef.diagn.TVallc = nan;
         coef.diagn.TVsnrc = nan;
        coef.diagn.PTVallc = nan;
        coef.diagn.PTVsnrc = nan;
        coef.diagn.SNRallc = nan;
              coef.diagn.K = nan;
        coef.diagn.lo.name = cellstr(repmat('none',nallc,1));
        coef.diagn.hi.name = coef.diagn.lo.name;
          coef.diagn.lo.RR = nan*ones(nallc,1);
          coef.diagn.hi.RR = coef.diagn.lo.RR;
                   [~,ind] = sort(coef.aux.frq(1:nc),'ascend');
                   
        for i = 1 : nc-1
            c1 = ind(i);
            c2 = ind(i+1);
            coef.diagn.lo.name(i+1) = coef.name(c1);
              coef.diagn.hi.name(i) = coef.name(c2);
              coef.diagn.lo.RR(i+1) = 24*elor*(coef.aux.frq(c2) - coef.aux.frq(c1))/coef.aux.opt.rmin;
        end
        coef.diagn.hi.RR(1:nc-1) = coef.diagn.lo.RR(2:nc);
               coef.diagn.lo.RNM = nan*ones(nallc,1);
               coef.diagn.hi.RNM = nan*ones(nallc,1);
             coef.diagn.lo.CorMx = nan*ones(nallc,1);
             coef.diagn.hi.CorMx = nan*ones(nallc,1);
    end
    
    minsnr = coef.aux.opt.diagnminsnr;
       tbl = cell(nallc+6+sign(nI)*(nI+1),1);
    tbl{1} = 'UTide Summary Diagnostics Table:';
    tbl{2} = sprintf('   Rmin= %#-9.3g   MinSNR= %#-9.3g (* SNR >= MinSNR)',coef.aux.opt.rmin,minsnr);
    tbl{3} = sprintf('   K= %#-9.3g   SNRallc= %#-9.3g', coef.diagn.K,coef.diagn.SNRallc);
    tbl{4} = sprintf('   TVallc= %#-9.3g    TVsnrc= %#-9.3g  TVraw= %#-9.3g',coef.diagn.TVallc,coef.diagn.TVsnrc,coef.diagn.TVraw);
    tbl{5} = sprintf('   PTVallc= %5.1f%%   PTVsnrc= %5.1f%%',coef.diagn.PTVallc,coef.diagn.PTVsnrc);
    tbl{6} = ' NAME    PE       SNR loNAME   loRR    loRNM  loCorMx hiNAME   hiRR    hiRNM  hiCorMx';
      main = [blanks(nallc)' char(coef.diagn.name)];
      main = [main reshape(sprintf('%6.2f%%',coef.diagn.PE),7,nallc)'];
      main = [main reshape(sprintf(' %#8.2g',coef.diagn.SNR),9,nallc)'];
    
    if nI
        coef.diagn.lo.name(end-nI+1:end) = cellstr(repmat('(I) ',nI,1));
        coef.diagn.hi.name(end-nI+1:end) = cellstr(repmat('(I) ',nI,1));
    end
    
    lo = char(coef.diagn.lo.name);
    lo = [lo reshape(sprintf(' %#8.2g',coef.diagn.lo.RR),9,nallc)'];
    lo = [lo reshape(sprintf(' %#8.2g',coef.diagn.lo.RNM),9,nallc)'];
    lo = [lo reshape(sprintf(' %#8.2g',coef.diagn.lo.CorMx),9,nallc)'];
    hi = char(coef.diagn.hi.name);
    hi = [hi reshape(sprintf(' %#8.2g',coef.diagn.hi.RR),9,nallc)'];
    hi = [hi reshape(sprintf(' %#8.2g',coef.diagn.hi.RNM),9,nallc)'];
    hi = [hi reshape(sprintf(' %#8.2g',coef.diagn.hi.CorMx),9,nallc)'];
    
    tbl(7:7+nallc-1) = cellstr([main repmat(' ',nallc,1) lo repmat(' ',nallc,1) hi]);
    if nI
        tbl{end-nI} = '  Inferred constituents (I):';
         cntr = 0;
        for i = 1:nR
            for j = 1:cnstit.R{i}.nI
             cntr = cntr+1;
                if coef.aux.opt.twodim
                    tbl{end-nI+cntr} = sprintf(['    %s inferred from reference %s (r+,r- = %#.3g,%#.3g; ', 'zta+,zta- = %#.3g,%#.3g deg)'],...
                        deblank(cnstit.R{i}.I.name{j}),...
                        deblank(cnstit.R{i}.name),...
                        abs(cnstit.R{i}.I.Rp(j)),...
                        abs(cnstit.R{i}.I.Rm(j)),...
                        (180/pi)*angle(cnstit.R{i}.I.Rp(j)),...
                        -(180/pi)*angle(cnstit.R{i}.I.Rm(j)));
                else
                    tbl{end-nI+cntr} = sprintf(['    %s inferred from reference %s (r = %#.3g zta = %#.3g deg)'],...
                        deblank(cnstit.R{i}.I.name{j}),...
                        deblank(cnstit.R{i}.name),...
                        abs(cnstit.R{i}.I.Rp(j)),...
                        (180/pi)*angle(cnstit.R{i}.I.Rp(j)));
                end
            end
        end
    end
    coef.diagn.table = char(reshape(strrep(tbl(:)','NaN','---')',size(tbl,1),size(tbl,2)));
end

%% create coef.results
coef.results = sprintf('CUTide Results:');
if ~isnan(coef.g(1))
    nallc = length(coef.g);
    if coef.aux.opt.twodim
       coef.results = char(coef.results,' Cnstit  Lsmaj  Lsmaj_ci     Lsmin  Lsmin_ci     Theta  theta_lin_ci         g      g_lin_ci');
         main = [coef.Lsmaj coef.Lsmaj_ci coef.Lsmin coef.Lsmin_ci coef.theta coef.theta_lin_ci coef.g coef.g_lin_ci];
        for i = 1:nallc
            coef.results = char(coef.results,sprintf(['%4s %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g %#9.3g'],coef.name{i},main(i,:)));
        end
            coef.results = char(coef.results,sprintf([' Mean u,v (u/v units) = %#9.3g, %#9.3g'],coef.umean,coef.vmean));
        if ~coef.aux.opt.notrend
            coef.results = char(coef.results,sprintf([' Trend u,v slope (u/v units per day) = %#9.3g, %#9.3g'], coef.uslope,coef.vslope));            
        end
    else
        coef.results = char(coef.results,[' Cnstit      A      A_ci         g      g_lin_ci']);
         main = [coef.A coef.A_ci coef.g coef.g_lin_ci];
        for i = 1:nallc
            coef.results = char(coef.results,sprintf(['%4s %#9.3g %#9.3g %#9.3g %#9.3g'],coef.name{i},main(i,:)));
        end
            coef.results = char(coef.results,sprintf([' Mean (input units) = %#9.3g' ],coef.mean));
        if ~coef.aux.opt.notrend
            coef.results = char(coef.results,sprintf([' Trend slope (input units per day) = %#9.3g '],coef.slope));
        end
    end
else
    coef.results = char(coef.results,['Regression solution reached iteration limit without converging; '],...
                   '    coefficients and confidence intervals set to NaN.');
end

%% order fields
if coef.aux.opt.twodim
    fldord = {'name'; 'Lsmaj'; 'Lsmaj_ci'; 'W';'M'; 'M_rlzn'; ... % Added 'W', 'M', and 'M_rlzn' (SI)
              'Lsmin'; 'Lsmin_ci'; 'theta'; 'theta_lin_ci'; 'g'; 'g_lin_ci';...
              'Lsmaj_rlzn'; 'Lsmin_rlzn'; 'theta_rlzn'; 'g_rlzn'; ... % Added (SI)
              'umean'; 'vmean'; 'constr', 'Std'}; % Added 'constr' and 'Std' (SI)
else
    fldord = {'name'; 'A'; 'A_rlzn'; 'A_ci'; 'g'; 'g_rlzn'; 'g_lin_ci';'W';'M'; 'M_rlzn';...  % Added 'W' (PM) 'W','M', and 'M_rlzn' (SI)
              'mean';  'constr'; 'Std'}; % Added 'constr' and 'Std' (SI)
end
if ~coef.aux.opt.notrend
    if coef.aux.opt.twodim
        fldord{end+1} = 'uslope';
        fldord{end+1} = 'vslope';
    else
        fldord{end+1} = 'slope';
    end
end
   fldord{end+1} = 'results';
   fldord{end+1} = 'aux';
if ~coef.aux.opt.nodiagn
   fldord{end+1} = 'diagn';
end
%      disp(coef)
%      disp(fldord)
    coef = orderfields(coef,fldord);
coef.aux = orderfields(coef.aux,{'rundescr';'opt';'frq';'lind';'lat';'reftime';'rpds';}); % Added 'rpds' (SI)

%% runtime display
if ~isequal(coef.aux.opt.runtimedisp,'nnn')
    fprintf('\n');
end
if isequal(coef.aux.opt.runtimedisp(1),'y')
    disp(char(coef.results,' '));
end
if isequal(coef.aux.opt.runtimedisp(2),'y')
    disp(char(coef.aux.rundescr,' '));
end
if isequal(coef.aux.opt.runtimedisp(3),'y')
    if ~coef.aux.opt.nodiagn && ~isnan(coef.g(1))
        disp(char(coef.diagn.table,' '));
    end
end