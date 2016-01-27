function U = accuracy( xy_i, active_cf, cf, a_a)
%ACCURACY
%   Updates the accuracy of all the active condition/forecast rules.
%   
%   Jonas K. Sekamane. 
%   Version 0.01
    
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
    
    %% Forecasts
    
    % Total number of active rules
    s = size(idx_firm,1);
    % Intercepts in column 14-15 and coeffecients in column 16-19. We need
    % to find values for each active rule so we repeat s times.
    v = repmat(14:19,s,1);
    % Finding the index of each value in our 3D matrix, by using using the 
    % cf rule index, coloumn number and firm index.
    idx_value = sub2ind(size(cf), repmat(idx_cf,6,1), v(:), repmat(idx_firm,6,1));
    % Extract values and reformat into 2D matric, with a coloumn for each
    % intercept/coeffecient and a row for each active rule.
    values = reshape( cf(idx_value), s, 6);
    intercepts = values(:,1:2);
    coeffecients = values(:,3:6);
    % f(x, y) = (a, b) + (x, y)*(e, f; g, h) = (a, b) + (x*e+y*g, x*f+y*h)
    xy_i_forecast = intercepts + [sum(xy_i_target.*coeffecients(:,[1 3]), 2) sum(xy_i_target.*coeffecients(:,[2 4]), 2)];
    
    
    %% Distance from forecast to actual position.
    %distance = sqrt(sum(( xy_i_target-xy_i_forecast ).^2,2)); % distance
    distance = sum(( xy_i_target-xy_i_forecast ).^2,2); % squared distance
    
    
    %% Update accuracy and active count
    cf_updated = cf;
    
    % Accuracy in column 24 and active count in column 26.
    aa = repmat([24 26],s,1);
    % Finding the index of both accuracy and active count in our 3D matrix, 
    % by using using the cf rule index, coloumn number and firm index.
    idx_aa = sub2ind(size(cf), repmat(idx_cf,2,1), aa(:), repmat(idx_firm,2,1));
    % Extract accucary and active count. Reformat into 2D matric, with a 
    % coloumn for accucary and one for active count, and row for each rule.
    accuracy_active = reshape( cf(idx_aa), s, 2);
    
    % Update accuracy
    updated_accuracy = a_a * accuracy_active(:,1) + (1-a_a) * distance;
    % Count the number of times a rule has been activated.
    active_count = accuracy_active(:,2) + 1;
    
    cf_updated(idx_aa) = [updated_accuracy; active_count];
    
    
    %% Output variables
    U = cf_updated;
end