function pl = plot_amplitudes(amp,comp,freq, varargin)
% Create tidal amplitude or amplitude std graphics as in 
% Fig. 1 of Innocenti et al, DOI:10.1175/JTECH-D-21-0060.1
  
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
 
  pl = plot(box_tk,x(:),varargin{:});
    
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