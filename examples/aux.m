function varargout  = aux(fname,varargin)
% create a container for calling misc sub-functions defined in the aux.m file
%    fname: sub-function name
% varargin: sub-function args

    switch fname
        case 'read_hourly_data'
                [h,t] = read_hourly_data(varargin{:});
            varargout = {h,t};
            
        case {'bsc2ag', 'from_bsincos_to_amppha'}
            
            [a,g] = from_bsincos_to_amppha(varargin{:});
            varargout = {a,g}; 
            
        case {'plot_amplitudes'}
            varargout{1} = plot_amplitudes(varargin{:});

    end
end



function [h,T] = read_hourly_data(filename,date_start,date_end)
% read hourly observation within the date_start-date_end period from the specified file

         raw_data = load(filename);
            raw_t = raw_data.t;

               T  = ((date_start*24:date_end*24)/24)'; 
               lT = length(T); 
        [~,it,iT] = intersect(roundn(raw_t,-3),roundn(T,-3));

                h = nan(lT,1); % hourly water level for the ref period
            h(iT) = raw_data.h(it);
end


function [a,pha] = from_bsincos_to_amppha(m)
% convert sin-cos amplitudes to tidal amplitudes and phases

% n of tidal constituents
 m = m(:);
nc = (length(m)-1)/2; % assume to have the intercept in m

ap = m(1:nc);
am = m(nc+(1:nc));

a   = sqrt(ap.^2 + am.^2);
pha = atan2d(am,ap);


% ut_cluster
      I = ones(1,size(pha,2));
    ii  = (pha-pha(:,I ))> 360/2;
pha(ii) = pha(ii)-360;
     ii = (pha-pha(:,I))< - 360/2;
pha(ii) = pha(ii)+360;

pha(pha<0) = 360+pha(pha<0);

end



function pl = plot_amplitudes(amp,comp,freq, varargin)


% check if the plot must set the tick style
if any(strcmpi(varargin,'ticks'))
    id = find(strcmpi(varargin,'ticks'));
    ticks = varargin{id+1};
    varargin(id:id+1) = [];
else
    ticks = true;
end

% check if the plot must be staicase
if any(strcmpi(varargin,'staircase'))
    id = find(strcmpi(varargin,'staircase'));
    staircase = varargin{id+1};
    varargin(id:id+1) = [];
else
    staircase = true;
end

% make the amp beaing a column vect
  amp  = amp(:);

% No of constituents
 nc  = length(amp);
 
 
  % use a stair function for showing values at each frequency
      xtck = 1:nc;  
    box_tk = [xtck'-0.5 xtck' xtck'+0.5]';
    box_tk = box_tk(:);
    
   % order the constit by freq 
[frq,ifrq] = sort(freq);
       per = 1./frq;

       
   % create the ticks and tick labels for the const ordered by freq    
       ftk = [0.0027 0.0357  .5    1   2 4 6  ];
       Per = {'1yr' '28d' '48h' '24h' '12h' '6h' '4h' }; 
    
      clear xtk       
for f = 1 : numel(ftk) 
[~,xtk(f)] = min(abs(24*frq-ftk(f)));
end
 [xtk,iftk] = unique(xtk);
        Per = Per(iftk);
      
        Clist =  {'M2' 'S2'  'N2' 'K2' 'K1' 'O1' 'P1'  'MM' 'SA'  'M4' 'S4' 'M6'}; %'SSA' 'S1'
   [Clist,ic] = intersect(comp(ifrq),Clist,'stable');
   
   
 % use a stair function for showing values at each frequency 
 x = repmat(amp',3,1);
 if staircase==false
        x(1,:) = NaN;
        x(3,:) = NaN;
 end
 
  x = x(:);
pl = plot(box_tk,x,varargin{:});


    % handle the x-axis ticks and labels, if requested
    if ticks
        ylim = get(gca,'ylim');
        for t = 1 : numel(Clist)
           text(ic(t)-0.5,ylim(2)+0.125,[Clist{t}(1) '_{' Clist{t}(2:end) '}'],'fontsize',5)
        end

      grid on     
      set(gca,'xtick',xtk+0.5,'xticklabel',Per) 
    end  
end