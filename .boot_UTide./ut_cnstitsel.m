%%--------------------------------------------------------- 
function [nNR,nR,nI,cnstit,coef] = ut_cnstitsel(tref,minres,incnstit,infer)
% UT_CNSTITSEL()
% carry out constituent selection
% inputs
%       tref = reference time (datenum UTC)
%     minres = freq separation (cph) used in decision tree
%   incnstit = 'cnstit' input to ut_solv
%      infer = 'opt.infer' input to ut_solv
% outputs
%        nNR,nR,nI = number non-reference, reference, inferred constituents
%   cnstit.NR.name = cellstr of 4-char names of NR constits
%   cnstit.NR.frq  = frequencies (cph) of NR constits
%   cnstit.NR.lind = list indices (in ut_constants.mat) of NR constits
%   cnstit.R       = empty if no inference; otherwise, for each (i'th) R constit:
%                    cnstit.R{i}.name, .frq, .lind = as above, but for R constits
%                    cnstit.R{i}.I{j}.name, .frq, .lind = as above for j'th I constit
%   coef.name        = cellstr of names of all constituents (NR, R, and I)
%   coef.aux.frq     = frequencies (cph) of all constituents
%   coef.aux.lind    = list indices of all constituents
%   coef.aux.reftime = tref
% UTide v1p0 9/2011 d.codiga@gso.uri.edu 

    %% freq vals, all avail cnstits
    load('ut_constants.mat','const','shallow');
           [~,ader] = ut_astron(tref);
                 ii = isfinite(const.ishallow); % #ok
    const.freq(~ii) = const.doodson(~ii,:)*ader/24;
    for  k = find(ii)'
        ik = const.ishallow(k) + (0:const.nshallow(k)-1);
             const.freq(k) = sum(const.freq(shallow.iname(ik)).*shallow.coef(ik));
    end

    %% cnstit.NR
    if isequal(lower(incnstit),'auto')
        cnstit.NR.lind = find(const.df>=minres); 
    else
        cnstit.NR.lind = nan*ones(size(incnstit,1),1);
        for j = 1:length(incnstit)
            lind1 = strmatch(incnstit{j},const.name);
            if isempty(lind1)
                error(['ut_solv: unrecognized non-reference constituent: ' incnstit(j,:) '.']);
            elseif lind1==1
                error(['ut_solv: do not pass in constituent Z0, mean is removed by default.']);
            else
                cnstit.NR.lind(j) = lind1;
            end
        end
               [~,seq] = sort(const.freq(cnstit.NR.lind));
        cnstit.NR.lind = cnstit.NR.lind(seq);
    end

    if isempty(cnstit.NR.lind)
        error('ut_solv: no constituents specified nor auto-selected.');
    end

    if ~isempty(infer)
        % remove from NR any R or I constits
              nam = cellstr(char(unique([infer.infnam; infer.refnam;])));
        [tst,ind] = ismember(nam,cellstr(const.name(cnstit.NR.lind,:)));
              tst = find(tst);
        if ~isempty(tst)
            cnstit.NR.lind(ind(tst)) = [];
        end
    end

    % output of NR const
    cnstit.NR.frq  = const.freq(cnstit.NR.lind);
    cnstit.NR.name = cellstr(const.name(cnstit.NR.lind,:));
               nNR = length(cnstit.NR.frq); 

    %% cnstit.R 
          nR = 0;
          nI = 0;
    cnstit.R = [];
    if ~isempty(infer)
                 nI = length(infer.infnam);
        allrefnames = unique(infer.refnam);
                 nR = length(allrefnames);
        for k = 1 : nR
            
            % reference
            cnstit.R{k}.name = allrefnames{k};
                       lind1 = strmatch(cnstit.R{k}.name,const.name);
                       
            if isempty(lind1)
                error(['ut_solv: unrecognized reference constituent: ' cnstit.R{k}.name '.']);
            elseif lind1==1
                error(['ut_solv: do not pass in constituent Z0, mean removal is handled automatically.']);
            else
                cnstit.R{k}.lind = lind1;
            end
            cnstit.R{k}.frq = const.freq(cnstit.R{k}.lind);

            % inferred
                       ind = strmatch(cnstit.R{k}.name,infer.refnam);
            cnstit.R{k}.nI = length(ind);
            for lk = 1:cnstit.R{k}.nI
                   cnstit.R{k}.I.Rp(lk) = infer.amprat(ind(lk)) .* ...
                                          exp( 1i * infer.phsoff(ind(lk))*pi/180 );
                if length(infer.amprat)>size(infer.infnam,1)
                   cnstit.R{k}.I.Rm(lk) = infer.amprat(ind(lk)+nI).*...
                                          exp( -1i * infer.phsoff(ind(lk)+nI)*pi/180 );
                else
                   cnstit.R{k}.I.Rm(lk) = conj(cnstit.R{k}.I.Rp(lk));
                end
                cnstit.R{k}.I.name{lk} = infer.infnam{ind(lk)};
                                 lind1 = strmatch(cnstit.R{k}.I.name{lk},const.name);
                if isempty(lind1)
                    error(['ut_solv: unrecognized inference constituent: ' cnstit.R{k}.I.name{lk} '.']);
                elseif lind1==1
                    error(['ut_solv: do not pass in constituent Z0, mean removal is handled automatically.']);
                else
                    cnstit.R{k}.I.lind(lk) = lind1;
                end
                cnstit.R{k}.I.frq(lk) = const.freq(cnstit.R{k}.I.lind(lk));
            end
        end
    end

        nallc     = nNR+nR+nI;
    coef.name     = cnstit.NR.name; 
    coef.aux.frq  = nan*ones(nallc,1);
    coef.aux.lind = coef.aux.frq;
    coef.aux.frq(1:nNR)  = cnstit.NR.frq;
    coef.aux.lind(1:nNR) = cnstit.NR.lind;

    if ~isempty(infer)
        for k = 1:nR
            coef.aux.frq(nNR+k)  = cnstit.R{k}.frq; 
            coef.aux.lind(nNR+k) = cnstit.R{k}.lind;
            coef.name(nNR+k)     = cellstr(cnstit.R{k}.name); 
        end
        ind = nNR+nR+1;
        for k = 1:nR
            coef.aux.frq(ind:ind+cnstit.R{k}.nI-1)  = cnstit.R{k}.I.frq; 
            coef.aux.lind(ind:ind+cnstit.R{k}.nI-1) = cnstit.R{k}.I.lind; 
            coef.name(ind:ind+cnstit.R{k}.nI-1)     = cellstr(char(cnstit.R{k}.I.name));
            ind = ind + cnstit.R{k}.nI;
        end
    end
    coef.aux.reftime = tref;
end
