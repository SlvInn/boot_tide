   
function options = set_options(function_name,varargin)
% Set the optional parameter of the funtion specified in function_name
% 
% S. Innocenti, slv.innoc@gmail.com 
% adapted from SLMtools by John D'Errico (2009)
% 2014-09-14, modified 2022-07-22
% 
    % if we d'ont have optional inputs or the options are in a structure
   if (nargin==1) || ~isstruct(varargin{1}) 
      defaults = set_default(function_name);
        
    % if we already have a structure of options    
   elseif isstruct(varargin{1})                              
      defaults    = varargin{1};
      varargin(1) = [];
   end   

    
   if isempty(varargin)
      options = defaults; % use the defaults already present in defaults
   else    
      options = paired_arguments(defaults, varargin); % update the default values with those in varargin
   end  

end

function options = paired_arguments(options, paired_param) %%% function [options] = paired_arguments(options, paired_param)
% Update the strcucture with default parameters (options) 
% with the user entered values (paired_param)
% 
% S. Innocenti, slv.innoc@gmail.com 
% adapted from SLMtools by John D'Errico (2009)
% 2014-09-14, modified 2022-07-22
% 
                                                                      
                                
if ~isstruct(options)
    error 'No structure supplied to set default parameters'
end

% check the num of input args
arg  = paired_param;
narg = length(arg);
n    = narg/2;

if n~=floor(n)
   error 'Property/value pairs must come in pairs.'
end

% get teh default param names                             
 propnames = fieldnames(options);
lpropnames = lower(propnames);

% match input values and original properies                                
for i = 1:n
  p_i = lower(arg{2*i-1});
  v_i = arg{2*i};

  ind = strmatch(p_i,lpropnames,'exact');
  if isempty(ind)
     ind = find(strncmp(p_i,lpropnames,length(p_i)));
     
   % OLD:   
   %   if isempty(ind)
   %      error(['No matching property found for: ', arg{2*i-1}])
   %   elseif length(ind)>1
   %      error(['Ambiguous property name: ', arg{2*i-1}])
   %   end

   % NEW (we don't want to draw an error if a property not found):
   if length(ind)>1
      error(['Ambiguous property name: ', arg{2*i-1}])
   end
     
  end

  if ~isempty(ind) % NEW (we don't want to draw an error if a property not found)
      p_i = propnames{ind};

      % override the corresponding default in params
      options = setfield(options,p_i,v_i);
   end   
end

end
