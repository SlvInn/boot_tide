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
    
    