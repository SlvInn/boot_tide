function yboot = from_mbb_to_resamples(out_mbb,locations)
% Return a matrix of the bootstrap resamples for y based on the time index 
% bootstrappend via MBB in boot_tide.m.
%
% INPUT:
% out_mbb   -  output structure from boot_tide.m
% locations - (optional) vector of y columns id indicating the spatial location to consider  
% 
% OUTPUT:
%  ....
%
% S. Innocenti, silvia.innocenti@ec.gc.ca, 2022/08

    assert(strcmp(out_mbb.options.method,'mbb'),'only MBB boot_tide output can be processed')

    % get reconstructions and original residuals
    yhat = out_mbb.yhat;
    yres = out_mbb.yres;
    assert(all(size(yhat)==size(yres)),'the output MBB structure contains inconsistent yhat and yres')

    # tot n of spatial locations
    ns = size(yres,s);

    % get the index of the spatial location to sample
    if nargin<2
       locations = 1:ns;
    else
           
       assert(length(locations)<=ns,'too many locations requested')
       assert(all(locations>1) && all(locations<ns),['location indices must be in [1,' num2str(ns) ']'])
       ns = length(locations); 
    end

    % time vector length
    lt = size(yres,1);

    # n of boot resamples
    nb = out_mbb.options.nboot;

    % construct the y resample for each spatial location
    yboot = nan(lt,nb,ns)
    for s = 1 : ns
        l  = locations(s);
        
        for b = 1 : nb
            ib = out_mbb.iboot(:,b);
            yboot(:,b,s) = yhat(:,l) + yres(ib,l);
        end
    end

end