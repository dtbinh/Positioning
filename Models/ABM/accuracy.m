function U = accuracy( xy_i, active_cf, cf, a_a)
%ACCURACY
%   Updates the accuracy of all the active condition/forecast rules.
%   
%   Jonas K. Sekamane. 
%   Version 0.01
   
    M = size(cf, 1);
    N = size(xy_i, 1);
    
    % Form active_cf find the indecies that match cf.
    idx = [];
    for firm=1:size(active_cf,2)
        for target=1:size(active_cf,1)
            if(target ~= firm) % Non diagonal elements
                %[target firm]
                idx_firm_target = [repmat([firm target], length(active_cf{target, firm}), 1) active_cf{target, firm}];
                idx = [idx; idx_firm_target];
            end
        end
    end
    
    % Split index such that the variables and remaining code is legible.
    idx_firm = idx(:,1);
    idx_target = idx(:,2);
    idx_cf = idx(:,3);
    
    % Expand the coordinates set so it fits the index.
    xy_i_target = xy_i( idx_target, : );
    
    % forecasts
    xy_i_forecast = NaN(size(xy_i_target));
    for a=1:length(idx_firm)
        xy_i_forecast(a,:) = cf(idx_cf(a), 14:15, idx_firm(a)) + xy_i_target(a,:) * reshape( cf(idx_cf(a), 16:19, idx_firm(a)), 2, 2);
    end
    
    % Distance from forecast to actual position.
    %distance = sqrt(sum(( xy_i_target-xy_i_forecast ).^2,2)); % distance
    distance = sum(( xy_i_target-xy_i_forecast ).^2,2); % squared distance
    
    % Update accuracy
    cf_updated = cf;
    for a=1:length(idx_firm)
        cf_updated(idx_cf(a), 24, idx_firm(a)) = a_a * cf(idx_cf(a), 24, idx_firm(a)) + (1-a_a) * distance(a);
        % Count the number of times a rule has been activated.
        cf_updated(idx_cf(a), 26, idx_firm(a)) = cf_updated(idx_cf(a), 26, idx_firm(a))+1;
    end
        
    % Output variables
    U = cf_updated;
end