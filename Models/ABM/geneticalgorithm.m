function cf_new = geneticalgorithm(cf_i, C, prob_crossover, cf_range)
%GENETICALGORITHM
%   Replaces the 20% worst performing condition/forecast rules with new 
%   rules using crossover and mutation.
%   
%   Jonas K. Sekamane. 
%   Version 0.01

    %% 1. Setup

    % Total number of condition/forecast rules
    M = size(cf_i, 1);
    % Number of firms
    N = size(cf_i, 3);
    cf_n = cf_i;
    
    %% 2.  Find rules for replacement (worst), and the parent rules (elite).
    % Calculate specificity of each condition (count the number of 0s and 
    % 1s in the condition).
    specificity = sum( ~isnan( cf_i(:, 1:13, : )), 2);
    
    % Calculate fitness of remaining condition/forecast rules as
    % function of the forecast variance/accurracy and with a cost on
    % the specificity of the condition. The latter creates a slow drift
    % towards conditions with many wildcard characters (#/NaN).
    fitness = -cf_i(:, 24, :) - C*specificity;
        
    % Sort condition/forecast rules based on their accuracy.
    % We exclude the first rule (that matches any state) for the possible 
    % worst case candidates. We add it to the elite (so it may be parent).
    % The special case / first rule is handled at the end of the loop. 
    [~, order] = sort(fitness(2:end,:,:));
    order = order+1;
    order = squeeze(order);
    % Identify the 20% worst performing condition/forecast rules.
    idx_worst = order(1:M/5, :); % The worst performing rules.
    idx_elite = [ones(1, N); order(M/5+1:end, :)]; % The remaining elite rules.
    
    %% 3. Loop through firms
    for n = 1:N
        
        
        %% 3.1 Fitness
        fitness_n = fitness(idx_elite(:,n), :, n);
        
        %% 3.2 Replacement procedure
        
        % Randomly select two parents among remaining elite rules, for each 
        % offspring (each of the rules that needs a replacement).
        parents = NaN(M/5, 2);
        for replacement = 1:M/5
            parents(replacement,:) = randperm(M*4/5, 2);
        end
        %idx_parents_n = reshape(idx_elite(parents(:),n), M/5, 2);
        
        % Determining if offspring is created using the crossover operation
        % or the mutation operation. Using the crossover probability.
        num_crossover = sum( ( rand(M/5,1) < repmat(prob_crossover, M/5, 1) ) );
        
        %% 3.2.1 Mutation operation
            mutation = (num_crossover+1) : M/5;

            % Selecting the most fit parrent.
            [~, fit_parent] = max( fitness_n(parents(mutation,:)), [], 2);
            parent = parents( mutation' + M/5*(fit_parent - 1) );
            parent_bits = cf_i(parent, 1:13, n);
            
            % The offspring inherits the mutated condition of the parent.
            offspring_bits = parent_bits;
            % We mutate the condition bits of the parent by randomly 
            % flipping the bits. There is a 0.03 probability that an 
            % individual bits is mutated (on average 0.39 of 13 bits are 
            % mutated). 
            bits_flipped = find(rand(length(mutation), 13) < 0.03);
            % If the individual bit is mutated then the flips 
            % happens with the following probabilities:
            % * Flips 0 to # with prob. 2/3. Flips 0 to 1 with prob. 1/3
            % * Flips 1 to # with prob. 2/3. Flips 1 to 0 with prob. 1/3
            % * Flips # to 0 with prob. 1/3. Flips # to 1 with prob. 1/3.
            % and with 1/3 prob that # stays at #. 
            % The average fraction of #'s / specificity is unchanged.
            flip = [NaN     NaN     1; ...
                    NaN     NaN     0; ...
                    0       1       NaN];
            % Draw random number and determine if number is inbetween 
            % 0 - 1/3, 1/3 - 2/3 or inbetween 2/3 - 1. The three bins are
            % represented with the numbers 1, 2 and 3.
            bits_draw = rand(length(bits_flipped),1);
            bits_drawbin = NaN(length(bits_flipped),1);
            for bit = 1 : length(bits_flipped)
                bits_drawbin(bit) = find(bits_draw(bit) <= [1/3 2/3 1], 1, 'first');
            end
            % Transform the bits from [0, 1, NaN] into [1, 2, 3]
            bits_transform = parent_bits(bits_flipped);
            bits_transform(isnan(bits_transform)) = 2;
            bits_transform = bits_transform + 1;
            % Match the randomly drawn number with the respective flip.
            flip_idx = sub2ind(size(flip), bits_transform, bits_drawbin);
            % Add the mutated values to offspring.
            offspring_bits(bits_flipped) = flip(flip_idx);
            
        cf_n(idx_worst(mutation,n), 1:13, n) = offspring_bits;
            
        
            % The offspring inherits the mutated values of the parent.
            parent_values = cf_i(parent, 14:19, n);
            offspring_values = parent_values;
            % We mutate the values by randomly choosing new values from
            % range, or by changing the values with a random ammount. Each
            % case happens with probability 0.2.
            values_draw = rand(size(parent_values));
            
            % Randomly choose new value from range.
            values_replace_idx = find(values_draw <= 0.2);
            [~, values_replace_value] = ind2sub(size(parent_values), values_replace_idx);
            offspring_values(values_replace_idx) = cf_range(values_replace_value, 1) + rand(length(values_replace_value),1) .* diff(cf_range(values_replace_value, :), 1, 2);
            
            % Change value with small uniformly randomly drawn amount from -0.5% to +0.5% of range.
            values_change_idx = find(values_draw > 0.2 & values_draw <= 0.4);
            [~, values_change_value] = ind2sub(size(parent_values), values_change_idx);
            values_change_range = diff(cf_range(values_change_value,:), 1, 2);
            drange = 0.005; % changed with up to +/- 0.5% of the range.
            dchange = (-drange*values_change_range) + rand(length(values_change_value), 1) .*  (drange*values_change_range*2);
            offspring_values(values_change_idx) = parent_values(values_change_idx) + dchange;
          
        cf_n(idx_worst(mutation,n), 14:19, n) = offspring_values;
        
            
        
        %% 3.2.2 Crossover operation
            crossover = 1 : num_crossover;
            
            % Uniform crossover of bits. The offspring inherits individual 
            % bits from one of its two parents. Where each parent is  
            % equally likely to be the one that passes the individual bit 
            % on to the offspring.
            % For each of the 13 bits randomly draw 0 or 1. Where 1 is the 
            % first parent, and 0 is the second parent.
            parent_draw = logical(randi([0 1], num_crossover, 13));
            % The 13 bits from each parent.
            parent1_bits = cf_i(parents(crossover,1), 1:13, n);
            parent2_bits = cf_i(parents(crossover,2), 1:13, n);
            % Create empty matrix and fill in bits from respective parent, 
            % based on the draw.
            offspring_bits = NaN(size(parent1_bits));
            offspring_bits(parent_draw) = parent1_bits(parent_draw);
            offspring_bits(~parent_draw) = parent2_bits(~parent_draw);
            
        cf_n(idx_worst(crossover,n), 1:13, n) = offspring_bits;
            
        
            % Crossover of parents' values. 
            parent1_values = cf_i(parents(crossover,1), 14:19, n);
            parent2_values = cf_i(parents(crossover,2), 14:19, n);
            
            % Calculating weights based on accuracy (better accuracy gives 
            % higher weight. Transformed so weights sum to one).
            parent1_weight = 1./(1+cf_i(parents(crossover,1), 24, n)')';
            parent2_weight = 1./(1+cf_i(parents(crossover,2), 24, n)')';
            parents_weightsum = parent1_weight + parent2_weight;
            parent1_w = parent1_weight./parents_weightsum;
            parent2_w = parent2_weight./parents_weightsum;
            
            % Three methods selected with equal probability.
            crossover_method = randi(3, num_crossover, 1);
            
            offspring_values = NaN(size(parent1_values));
            for offspring = 1 : num_crossover
               
                if(crossover_method(offspring) == 1)
                    % Method 1: Component-wise crossover of each parent value.
                    parent_draw = logical(randi([0 1], 1, 6));
                    offspring_values(offspring, parent_draw) = parent1_values(offspring, parent_draw);
                    offspring_values(offspring, ~parent_draw) = parent2_values(offspring, ~parent_draw);
                    
                elseif(crossover_method(offspring) == 2)
                    % Method 2: Weighted average of the values of parents.
                    offspring_values(offspring,:) = parent1_values(offspring,:).*parent1_w(offspring) + parent2_values(offspring,:).*parent2_w(offspring);
                    
                else 
                    % Method 3: Randomly pick parent to pass on all values
                    parent_draw = randi([0 1], 1, 1);
                    offspring_values(offspring,:) = parent1_values(offspring,:)*parent_draw + parent2_values(offspring,:)*(1-parent_draw);
                end
                
            end    
            
        cf_n(idx_worst(crossover,n), 14:19, n) = offspring_values;
            
        
        
        %% 3.2.3 Forecast variance of offspring.
            
            % Offspring inherits the average forecast variance of its
            % parents, unless parents have never been activated, in which
            % case the offspring gets the median forecast variance of all
            % elite forecast rules.
            parents_active = sum(reshape( cf_i(parents(:), 26, n), M/5, 2 ), 2) > 0;
            
            parents_accuracy = mean( [cf_i(parents(parents_active, 1), 24, n) cf_i(parents(parents_active, 2), 24, n)], 2);
            median_accuracy = median( cf_i(idx_elite(:,n), 24, n) );
            
            offspring_accuracy = NaN(M/5, 1);
            offspring_accuracy(parents_active) = parents_accuracy;
            offspring_accuracy(~parents_active) = median_accuracy;
            
        cf_n(idx_worst(:,n), 24, n) = offspring_accuracy;
        
        % Setting the number of times the rule has been used to zero.
        cf_n(idx_worst(:,n), 25, n) = 0;
        % Setting the number of times the rule has been active to zero.
        cf_n(idx_worst(:,n), 26, n) = 0;
        
        
        %% 3.3 Special case: Default rule that matches any state
        % The values of the default rule is set to be the weighted average
        % values of all other rules. The weight is based on the rules
        % accuracy.
        other_values = cf_n(2:end, 14:19, n);
        other_weight = 1./(1 + cf_n(2:end, 24, n));
        other_weightsum = sum( other_weight );
        other_w = other_weight/other_weightsum;
        
        cf_f(1, 14:19, n) = sum( other_values .* repmat(other_w, 1, 6) );
        
    end
    
    %% 4. Output variables
    cf_new = cf_n;
end