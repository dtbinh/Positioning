function [F, I, A] = forecast( firm, xy_i, J, cf)
%FORECAST
%   Finds the active forecast rules (A). Selects the forecast rule with best
%   accuracy (I). And creates forecast of other firms position (F).
%   
%   Jonas K. Sekamane. 
%   Version 0.01
   
    M = size(cf, 1);
    N = size(xy_i, 1);
    others = 1:N; others(firm) = [];
    
    f_xy = NaN(size(xy_i));
    f_rule = NaN(1, N);
    f_active = cell(N,1);
    
    for n = others
        
        %%  Find the conditions that satisfy the current state J.
        % Subtracting the conditions from the state J gives 0 if a particular
        % condition is statisfied, and 1 or -1 if it is not satisfied. Summing over
        % all the particular conditions reveals if any of the conditions are 
        % unfulfilled. By ignoring NaN when summing makes NaN a wildcard character, 
        % that will match both 0s and 1s in state J.
        count_unfulfilled = sum( abs( repmat(J(n,:), M, 1) - cf(:,1:13) ), 2, 'omitnan' );
        % Indecies of the conditions that satisfy the current state J.
        c_idx = find( count_unfulfilled == 0);
        
        if(length(c_idx) == 1)
            c_idx_accuracy = c_idx;
        else
            % Select the condition/forecast rule with best accuracy.
            c_idx_accuracy = find( cf(c_idx,24) == min(cf(c_idx,24)) );
            % If there are several condition/forecast rules with the same
            % accuracy then randomly pick one.
            c_idx_accuracy_l = length(c_idx_accuracy);
            if (c_idx_accuracy_l > 1)
                c_idx_accuracy = c_idx_accuracy( randi(c_idx_accuracy_l) );
            end
        end
        
        % Use condition/forecast rule idx.
        idx = c_idx(c_idx_accuracy);
        % 2D linear forecast model
        f_xy(n,:) = cf(idx,14:15) + xy_i(n,:) * reshape(cf(idx,16:19), 2, 2);
        
        f_rule(n) = idx;
        f_active{n} = c_idx;
    end
    
    % Do not forecast own position
    f_xy(firm,:) = [];
    %f_xy(firm,:) = xy_i(firm,:)
    
    % Output variables
    F = f_xy;
    I = f_rule;
    A = f_active;
end