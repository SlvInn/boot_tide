function coef = cut_disgnostic(coef,opt,e,cnstit,t,u,v,xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,tin,uin,vin)
%     fprintf('diagnostics ... ');

    if opt.twodim
         PE = sum(coef.Lsmaj.^2 + coef.Lsmin.^2);
         PE = 100*(coef.Lsmaj.^2 + coef.Lsmin.^2)/PE;
        SNR = (coef.Lsmaj.^2 +coef.Lsmin.^2)./((coef.Lsmaj_ci/1.96).^2 + (coef.Lsmin_ci/1.96).^2);
    else
         PE = 100*coef.A.^2/sum(coef.A.^2);
        SNR = (coef.A.^2)./((coef.A_ci/1.96).^2);
    end
    

    % Attention (SI): I think we should reorder the diagnostic as for the
    % output but the original UTide code does not do that:
            [~,indPE] = sort(PE,'descend');
      coef.diagn.name = coef.name(indPE);
      coef.diagn.PE   = PE(indPE);
      coef.diagn.SNR  = SNR(indPE); % used in ut_diagntable; ordered by PE there
    
    if opt.twodim
        [coef.diagn,usnrc,vsnrc] = ut_diagntable(coef,cnstit,t,u,v,xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,indPE);
    else
            [coef.diagn,usnrc,~] = ut_diagntable(coef,cnstit, t,u,[],xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mCw,indPE);
    end
    
    if opt.diagnplots
              tmp = nan*ones(size(uin));
        tmp(uvgd) = usnrc;
            usnrc = tmp;
              tmp = nan*ones(size(uin));
        tmp(uvgd) = e;
        
        e = tmp;
        if opt.twodim           
                  tmp = nan*ones(size(uin));
            tmp(uvgd) = vsnrc;
                vsnrc = tmp;
                ut_diagnfigs(coef,indPE,tin,uin,vin,usnrc,vsnrc,e);
        else
                ut_diagnfigs(coef,indPE,tin,uin,[],usnrc,[],e);
        end
    end
 
end
 
 
 
%%--------------------------------------------------------- 
function [diagn,usnrc,vsnrc] = ut_diagntable(coef,cnstit,t,u,v,xmod,m,B,W,varMSM,Gall,Hall,elor,varcov_mC,indPE)
% UT_DIAGNTABLE()
% compute and store diagnostics, then create character array table
% inputs: see ut_solv1()
% outputs:
%   diagn = diagnostic quantities results (see report)
%   usnrc, vsnrc = reconstructed superposed harmonics using constituents
%       with SNR above the threshold
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

nNR = length(cnstit.NR.frq);
 nR = 0;
 nI = 0;
if ~isempty(coef.aux.opt.infer)
    nR = length(cnstit.R);
    nI = length(coef.aux.opt.infer.infnam);
end    
nallc = nNR+nR+nI;
   nc = nNR+nR;
diagn = coef.diagn;
 tref = coef.aux.reftime;
[usnrc,vsnrc] = ut_diagnrcn(t,coef,diagn.SNR,coef.aux.opt.diagnminsnr);

% TIDAL VARIANCE DIAGNOSTICS 
% TVallc, TVraw, PTVallc
if coef.aux.opt.twodim
    urawtid  = u - coef.umean;
    vrawtid  = v - coef.vmean;
    umodtid  = real(xmod) - coef.umean;
    vmodtid  = imag(xmod) - coef.vmean;
    usnrctid = usnrc - coef.umean;
    vsnrctid = vsnrc - coef.vmean;
    if ~coef.aux.opt.notrend
        umodtid  =  umodtid - coef.uslope*(t-tref);
        vmodtid  =  vmodtid - coef.vslope*(t-tref);
        urawtid  =  urawtid - coef.uslope*(t-tref);
        vrawtid  =  vrawtid - coef.vslope*(t-tref);
        usnrctid = usnrctid - coef.uslope*(t-tref);
        vsnrctid = vsnrctid - coef.vslope*(t-tref);
    end
    diagn.TVraw  = mean(urawtid.*urawtid + vrawtid.*vrawtid);
    diagn.TVallc = mean(umodtid.*umodtid + vmodtid.*vmodtid);
    diagn.TVsnrc = mean(usnrctid.*usnrctid + vsnrctid.*vsnrctid);
else
    urawtid  = u - coef.mean;
    umodtid  = real(xmod) - coef.mean;
    usnrctid = usnrc - coef.mean;
    if ~coef.aux.opt.notrend
        urawtid  = urawtid - coef.slope*(t-tref);
        umodtid  = umodtid - coef.slope*(t-tref);
        usnrctid = usnrctid - coef.slope*(t-tref);
    end
    diagn.TVraw = mean(urawtid.*urawtid);
    diagn.TVallc = mean(umodtid.*umodtid);
    diagn.TVsnrc = mean(usnrctid.*usnrctid);
end
diagn.PTVallc = 100*diagn.TVallc/diagn.TVraw;
diagn.PTVsnrc = 100*diagn.TVsnrc/diagn.TVraw;

% SNRallc, K
diagn.SNRallc = real(((m'*B')*W*(B*m))/varMSM);
diagn.K = cond(B);

% RR, RNM: Rayleight criterion ratio and Noise modified Rayleight criterion
% CorMx: Max corr with the nearest lower-higher frequency 
diagn.lo.name  = cellstr(repmat('none',nallc,1));
diagn.lo.RR    = nan*ones(nallc,1);
diagn.lo.RNM   = diagn.lo.RR;
diagn.lo.CorMx = diagn.lo.RR;
diagn.hi       = diagn.lo;

    [~,ind] = sort(coef.aux.frq(1:nc),'ascend');
for i = 1:nc-1
   c1 = ind(i);
   c2 = ind(i+1);
   
    diagn.lo.name(i+1) = coef.name(c1);
    diagn.hi.name(i)   = coef.name(c2);
    diagn.lo.RR(i+1)   = 24*elor*(coef.aux.frq(c2) - coef.aux.frq(c1))/coef.aux.opt.rmin;
    diagn.lo.RNM(i+1)  = diagn.lo.RR(i+1)*  sqrt((diagn.SNR(c1)+diagn.SNR(c2))/2);
    
       den = sqrt(diag(squeeze(varcov_mC(c1,:,:)))*diag(squeeze(varcov_mC(c2,:,:)))');
         G = [Gall(c1,c2) Gall(c1,c2+nc); Gall(c1+nc,c2) Gall(c1+nc,c2+nc);];
         H = [Hall(c1,c2) Hall(c1,c2+nc); Hall(c1+nc,c2) Hall(c1+nc,c2+nc);];
      cDuu = nan*ones(2,2);
      
    if ~coef.aux.opt.twodim
        
        cDuu(1,1) = real(sum(sum(G)))/2;
        cDuu(1,2) = imag(H(1,1)-H(1,2)+H(2,1)-H(2,2))/2;
        cDuu(2,1) = imag(-G(1,1)-G(1,2)+G(2,1)+G(2,2))/2;
        cDuu(2,2) = real(H(1,1)-H(1,2)-H(2,1)+H(2,2))/2;
        
        diagn.lo.CorMx(i+1) = max(max(abs(cDuu)./den));
    else
        cDuv(1,1) = imag(sum(sum(-H)))/2;
        cDuv(1,2) = real(G(1,1)-G(1,2)+G(2,1)-G(2,2))/2;
        cDuv(2,1) = real(-H(1,1)-H(1,2)+H(2,1)+H(2,2))/2;
        cDuv(2,2) = imag(-G(1,1)+G(1,2)+G(2,1)-G(2,2))/2;
        cDvu(1,1) = imag(sum(sum(G)))/2;
        cDvu(1,2) = real(-H(1,1)+H(1,2)-H(2,1)+H(2,2))/2;
        cDvu(2,1) = real(-G(1,1)-G(1,2)+G(2,1)+G(2,2))/2;
        cDvu(2,2) = imag(H(1,1)-H(1,2)-H(2,1)+H(2,2))/2;
        cDvv(1,1) = real(sum(sum(H)))/2;
        cDvv(1,2) = imag(G(1,1)-G(1,2)+G(2,1)-G(2,2))/2;
        cDvv(2,1) = imag(-H(1,1)-H(1,2)+H(2,1)+H(2,2))/2;
        cDvv(2,2) = real(G(1,1)-G(1,2)-G(2,1)+G(2,2))/2;
        
               cv = [cDuu cDuv; cDvu cDvv;];
        diagn.lo.CorMx(i+1) = max(max(abs(cv)./den));
    end
end

diagn.hi.RR(1:nc-1)    = diagn.lo.RR(2:nc);
diagn.hi.RNM(1:nc-1)   = diagn.lo.RNM(2:nc);
diagn.hi.CorMx(1:nc-1) = diagn.lo.CorMx(2:nc);

% sort
ind(ind) = 1:length(ind);

diagn.lo.name(1:nc)  = diagn.lo.name(ind);
diagn.lo.RR(1:nc)    = diagn.lo.RR(ind);
diagn.lo.RNM(1:nc)   = diagn.lo.RNM(ind);
diagn.lo.CorMx(1:nc) = diagn.lo.CorMx(ind);

diagn.hi.name(1:nc)  = diagn.hi.name(ind);
diagn.hi.RR(1:nc)    = diagn.hi.RR(ind);
diagn.hi.RNM(1:nc)   = diagn.hi.RNM(ind);
diagn.hi.CorMx(1:nc) = diagn.hi.CorMx(ind);

diagn.lo.name  = diagn.lo.name(indPE);
diagn.lo.RR    = diagn.lo.RR(indPE);
diagn.lo.RNM   = diagn.lo.RNM(indPE);
diagn.lo.CorMx = diagn.lo.CorMx(indPE);
diagn.hi.name  = diagn.hi.name(indPE);
diagn.hi.RR    = diagn.hi.RR(indPE);
diagn.hi.RNM   = diagn.hi.RNM(indPE);
diagn.hi.CorMx = diagn.hi.CorMx(indPE);

diagn.SNR = diagn.SNR(indPE);

% summary table char array
minsnr = coef.aux.opt.diagnminsnr;
tbl    = cell(nallc+6+sign(nI)*(nI+1),1);
tbl{1} = 'UTide Summary Diagnostics Table:';
tbl{2} = sprintf('   Rmin= %#-9.3g   MinSNR= %#-9.3g  (* SNR >= MinSNR)',coef.aux.opt.rmin,minsnr);
tbl{3} = sprintf('   K= %#-9.3g   SNRallc= %#-9.3g', diagn.K, diagn.SNRallc);
tbl{4} = sprintf(['   TVallc= %#-9.3g    TVsnrc= %#-9.3g    TVraw= %#-9.3g'],diagn.TVallc,diagn.TVsnrc,diagn.TVraw);
tbl{5} = sprintf('   PTVallc= %5.1f%%   PTVsnrc= %5.1f%%', diagn.PTVallc,diagn.PTVsnrc);
tbl{6} = [' NAME    PE       SNR loNAME   loRR    loRNM  loCorMx hiNAME   hiRR    hiRNM  hiCorMx'];

  main = [blanks(nallc)' char(diagn.name)];
  main(diagn.SNR>=minsnr,1) = repmat('*', sum(diagn.SNR>=minsnr),1);
  main = [main reshape(sprintf('%6.2f%%',diagn.PE),7,nallc)'];
  main = [main reshape(sprintf(' %#8.2g',diagn.SNR),9,nallc)'];
  
if nI
    indPE(indPE) = 1:nallc;
    diagn.lo.name(indPE(end-nI+1:end)) = cellstr(repmat('(I) ',nI,1));
    diagn.hi.name(indPE(end-nI+1:end)) = cellstr(repmat('(I) ',nI,1));
end

lo = char(diagn.lo.name);
lo = [lo reshape(sprintf(' %#8.2g',diagn.lo.RR),9,nallc)'];
lo = [lo reshape(sprintf(' %#8.2g',diagn.lo.RNM),9,nallc)'];
lo = [lo reshape(sprintf(' %#8.2g',diagn.lo.CorMx),9,nallc)'];
hi = char(diagn.hi.name);
hi = [hi reshape(sprintf(' %#8.2g',diagn.hi.RR),9,nallc)'];
hi = [hi reshape(sprintf(' %#8.2g',diagn.hi.RNM),9,nallc)'];
hi = [hi reshape(sprintf(' %#8.2g',diagn.hi.CorMx),9,nallc)'];

tbl(7:7+nallc-1) = cellstr([main repmat(' ',nallc,1) lo repmat(' ',nallc,1) hi]);

if nI
    tbl{end-nI} = '  Inferred constituents (I):';
    cntr =0;
    for i=1:nR
        for j = 1:cnstit.R{i}.nI
            cntr = cntr+1;
            if coef.aux.opt.twodim
                tbl{end-nI+cntr} = sprintf(['    %s inferred from'...
                    ' reference %s (r+,r- = %#.3g,%#.3g; ',...
                    'zta+,zta- = %#.3g,%#.3g deg)'],...
                    deblank(cnstit.R{i}.I.name{j}),...
                    deblank(cnstit.R{i}.name),...
                    abs(cnstit.R{i}.I.Rp(j)),...
                    abs(cnstit.R{i}.I.Rm(j)),...
                    (180/pi)*angle(cnstit.R{i}.I.Rp(j)),...
                    -(180/pi)*angle(cnstit.R{i}.I.Rm(j)));
            else
                tbl{end-nI+cntr} = sprintf(['    %s inferred from'...
                    ' reference %s (r = %#.3g zta = %#.3g deg)'],...
                    deblank(cnstit.R{i}.I.name{j}),...
                    deblank(cnstit.R{i}.name),...
                    abs(cnstit.R{i}.I.Rp(j)),...
                    (180/pi)*angle(cnstit.R{i}.I.Rp(j)));
            end
        end
    end
end
diagn.table = char(reshape(strrep(tbl(:)','NaN','---')',size(tbl,1),size(tbl,2)));
end

%%--------------------------------------------------------- 
function [u,v] = ut_diagnrcn(t,coef,SNR,MinSNR)
% UT_DIAGNRCN()
% abbreviated version of ut_reconstr1 solely for use by ut_diagntable:
%       accepts only non-nan times
%       instead of varargin, requires that you pass in
%           SNR (in NR-R-I order, same as coef.Lsmaj etc, or coef.A etc) 
%           and MinSNR
% does not allow for subsetting by PE or specified constits
%   (as ut_reconstr1 does)
% see comments for ut_reconstr()
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

ind = SNR>=MinSNR;
rpd = pi/180;
% complex coefficients
if coef.aux.opt.twodim
    ap = 0.5*(coef.Lsmaj(ind) + coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) - coef.g(ind))*rpd);
    am = 0.5*(coef.Lsmaj(ind) - coef.Lsmin(ind)) .* ...
        exp(1i*(coef.theta(ind) + coef.g(ind))*rpd);
else
    ap = 0.5*coef.A(ind).*exp(-1i*coef.g(ind)*rpd);
    am = conj(ap);
end
% exponentials
ngflgs = [coef.aux.opt.nodsatlint coef.aux.opt.nodsatnone ...
    coef.aux.opt.gwchlint coef.aux.opt.gwchnone];
E = cut_E(t,coef.aux.reftime,coef.aux.frq(ind),coef.aux.lind(ind),...
    coef.aux.lat,ngflgs,coef.aux.opt);
% fit calc
fit = E*ap + conj(E)*am;
% mean (& trend)
if coef.aux.opt.twodim
    if coef.aux.opt.notrend
        u = real(fit) + coef.umean;
        v = imag(fit) + coef.vmean;
    else
        u = real(fit) + coef.umean + ...
            coef.uslope*(t-coef.aux.reftime);
        v = imag(fit) + coef.vmean + ...
            coef.vslope*(t-coef.aux.reftime);
    end
else
    if coef.aux.opt.notrend
        u = real(fit) + coef.mean;
    else
        u = real(fit) + coef.mean + ...
            coef.slope*(t-coef.aux.reftime);
    end
    v = [];
end
end

%%--------------------------------------------------------- 
function ut_diagnfigs(coef,indPE,t,u,v,usnrc,vsnrc,e)
% UT_DIAGNFIGS()
% create two diagnostic figures
% see report for a description of their contents
% UTide v1p0 9/2011 d.codiga@gso.uri.edu

if coef.aux.opt.twodim
    % first figure
    figure;
    % description
    subplot(4,1,1);
    text(0.2,0.6,coef.aux.rundescr,'fontsize',8);
    set(gca,'visible','off');
    % u time sequences
    subplot(4,1,2);
    plot(t,u,'r.:');
    hold;
    plot(t,usnrc,'b.:');
    plot(t,real(e),'g.:');
    set(gca,'xlim',[min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #1:   u    Red= raw input;'...
        ' Blue= reconstructed fit (SNR>=%.1f); Green= residual'],...
        coef.aux.opt.diagnminsnr));
    grid on;
    % v time sequences
    subplot(4,1,3);
    plot(t,v,'r.:');
    hold;
    plot(t,vsnrc,'b.:');
    plot(t,imag(e),'g.:');
    set(gca,'xlim',[min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    ylabel('[raw input units]');
    title(sprintf(['v    Red= raw input; Blue= reconstructed fit'...
        ' (SNR>=%.1f); Green= residual'],coef.aux.opt.diagnminsnr));
    grid on;
    % all constituents, increasing frequency, rel SNR thresh
    subplot(4,1,4);
    [~,ind] = sort(coef.aux.frq,'ascend');
    SNRnum = coef.Lsmaj(ind).^2 + coef.Lsmin(ind).^2;
    SNRden = (coef.Lsmaj_ci(ind)/1.96).^2 + (coef.Lsmin_ci(ind)/1.96).^2;
    thresh = coef.aux.opt.diagnminsnr*SNRden;
    semilogy(thresh,'g');
    hold;
    ymin = 10.^(floor(log10(min([SNRnum; thresh])))-1);
    name = coef.name(ind);
    ind = find(SNRnum>=thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','r');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','r','verticalalign','bottom');
    end
    ind = find(SNRnum<thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','b');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','b','verticalalign','bottom');
    end
    set(gca,'ylim',[ymin max(get(gca,'ylim'))],'xlim',[0 length(SNRnum)+1],...
        'xticklabel',repmat(blanks(size(get(gca,'xticklabel'),2)),...
        size(get(gca,'xticklabel'),1),1));
    ylabel('L_{smaj}^{2} + L_{smin}^{2}');
    title(sprintf(['All constits, by increasing frequency; '...
        'Red= SNR>=%.1f, Blue= SNR<%.1f; Green= SNR thresh = %.1f('...
        '\\sigma_{L_{smaj}}^{2}+\\sigma_{L_{smin}}^{2})'],...
        coef.aux.opt.diagnminsnr,coef.aux.opt.diagnminsnr,...
        coef.aux.opt.diagnminsnr));
    % second figure
    figure;
    indPE = indPE(coef.diagn.SNR>=coef.aux.opt.diagnminsnr);
    % Lsmaj
    subplot(4,1,1);
    errorbar(coef.Lsmaj(indPE),coef.Lsmaj_ci(indPE),'.');
    grid on;
    set(gca,'xtick',1:length(indPE),'xticklabel',coef.name(indPE),...
        'xlim',[0 length(indPE)+1]);
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #2: Constits with SNR>%.1f, in '...
        'decreasing-PE order // Major axis L_{smaj} with 95%%'...
        ' confidence intervals'],coef.aux.opt.diagnminsnr));
    % Lsmin
    subplot(4,1,2);
    errorbar(coef.Lsmin(indPE),coef.Lsmin_ci(indPE),'.');
    grid on;
    set(gca,'xtick',1:length(indPE),'xticklabel',coef.name(indPE),...
        'xlim',[0 length(indPE)+1]);
    ylabel('[raw input units]');
    title('Minor axis L_{smin} with 95% confidence intervals');
    % orientation angle theta
    subplot(4,1,3);
    errorbar(coef.theta(indPE),coef.theta_ci(indPE),'.');
    grid on;
    set(gca,'ylim',[-45 225],'ytick',-45:45:225,'xlim',[0 length(indPE)+1],...
        'xtick',1:length(indPE),'xticklabel',coef.name(indPE));
    ylabel('[degrees]');
    title('Orientation angle theta with 95% confidence intervals');
    % phase lag g
    subplot(4,1,4);
    errorbar(coef.g(indPE),coef.g_ci(indPE),'.');
    grid on;
    set(gca,'ylim',[-90 450],'ytick',-90:90:450,'xlim',[0 length(indPE)+1],...
        'xtick',1:length(indPE),'xticklabel',coef.name(indPE));
    ylabel('[degrees]');
    title('Phase lag g with 95% confidence intervals');
    figure(1);
else % one-dimensional case
    % first figure
    figure;
    % description
    subplot(3,1,1);
    text(0.2,0.6,coef.aux.rundescr,'fontsize',8);
    set(gca,'visible','off');
    % time sequences
    subplot(3,1,2);
    plot(t,u,'r.:');
    hold;
    plot(t,usnrc,'b.:');
    plot(t,real(e),'g.:');
    set(gca,'xlim',[min(t) max(t)]);
    datetick('x','keepticks','keeplimits');
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #1:  Red= raw input;'...
        ' Blue= reconstructed fit (SNR>=%.1f); Green= residual'],...
        coef.aux.opt.diagnminsnr));
    grid on;
    % all constituents, increasing frequency, rel SNR thresh
    subplot(3,1,3);
    [~,ind] = sort(coef.aux.frq,'ascend');
    SNRnum = coef.A(ind).^2;
    SNRden = (coef.A_ci(ind)/1.96).^2;
    thresh = coef.aux.opt.diagnminsnr*SNRden;
    semilogy(thresh,'g');
    hold;
    ymin = 10.^(floor(log10(min([SNRnum; thresh])))-1);
    name = coef.name(ind);
    ind = find(SNRnum>=thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','r');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','r','verticalalign','bottom');
    end
    ind = find(SNRnum<thresh);
    for i = 1:length(ind)
        line([ind(i) ind(i)],[ymin SNRnum(ind(i))],...
            'marker','*','color','b');
        text(ind(i),SNRnum(ind(i)),['  ' name(ind(i))],'rotation',60,...
            'color','b','verticalalign','bottom');
    end
    set(gca,'ylim',[ymin max(get(gca,'ylim'))],'xlim',[0 length(SNRnum)+1],...
        'xticklabel',repmat(blanks(size(get(gca,'xticklabel'),2)),...
        size(get(gca,'xticklabel'),1),1));
    ylabel('A^{2}');
    title(sprintf(['All constits, by increasing frequency; '...
        'Red= SNR>=%.1f, Blue= SNR<%.1f; '...
        'Green= SNR thresh = %.1f \\sigma_{A}^{2} '],...
        coef.aux.opt.diagnminsnr,coef.aux.opt.diagnminsnr,...
        coef.aux.opt.diagnminsnr));
    % second figure
    figure;
    indPE = indPE(coef.diagn.SNR>=coef.aux.opt.diagnminsnr);
    % amplitude
    subplot(2,1,1);
    errorbar(coef.A(indPE),coef.A_ci(indPE),'.');
    grid on;
    set(gca,'xtick',1:length(indPE),'xticklabel',coef.name(indPE),...
        'xlim',[0 length(indPE)+1]);
    ylabel('[raw input units]');
    title(sprintf(['UTide Diagnostic Plot #2: Constits w/ SNR>%.1f, decreasing-PE order '...
        ' // Amplitude with 95%% confidence intervals'],...
        coef.aux.opt.diagnminsnr));
    % phase lag g
    subplot(2,1,2);
    errorbar(coef.g(indPE),coef.g_ci(indPE),'.');
    grid on;
    set(gca,'ylim',[-90 450],'ytick',-90:90:450,'xlim',[0 length(indPE)+1],...
        'xtick',1:length(indPE),'xticklabel',coef.name(indPE));
    ylabel('[degrees]');
    title('Phase lag g with 95% confidence intervals');
    figure(1);
end
end


