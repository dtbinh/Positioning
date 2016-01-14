%% ABM Endogenous number of firms
%   Jonas K. Sekamane. 
%   Version 0.01
%
%   Endogenous number of firms and decision rules (each rule selected with 
%   equal probability).
%
%   Inspired in part by: 
%       Lever and Sergenti (2011)

clearvars;

%% 1. PREFERENCES
pref.seed = rng('default');
pref.boundary = 10; % Number of standard deviations
pref.resolution = 70; % Length of the square (even number to include (0,0))
pref.N = 1; % Number of initial firms
pref.mu = 0; % Mean of subpopulation
pref.n_ratio = 1; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
pref.iterations = 100; % Number of system iterations/ticks
pref.ruleset = cellstr([{'STICKER'}, {'HUNTER'}, {'AGGREGATOR'}]); % The set of possible decision rules
pref.beta = 1; % The entry probability sensitivity to dissatisfaction.
pref.a_b = 0.5; % Customer memory parameter (weight on past dissatisfaction).
pref.a_f = 0.5; % Firm memory parameter (weight on past fitness).
pref.tau = 0.2; % The de facto survival threshold (market share/percentage).
pref.psi = 2; % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
pref.export_data = 0;
pref.export_fig = 0;


%% 2. SETUP
sd = 0.5; % Standard deviation of each subpopulation
b = pref.boundary/2;


    %%% 2.1 Population
    
    % Grid / Square
    [x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
    [X,Y] = meshgrid(x,y);
    % Create mask of the market that lies with in 3 std. dev. of (0,0).
    mask = sqrt(X.^2 + Y.^2) < 3; 
    X_mask = X(mask==1);
    Y_mask = Y(mask==1);
    
    % Subpopulation right
    mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
    sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
    F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
    F_r = reshape(F_r,length(y),length(x)); % Formated as grid

    % Subpopulation left
    mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
    sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
    F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
    F_l = reshape(F_l,length(y),length(x)); % Formated as grid

    % Total population 
    % The left subpopulation share of total population
    weight = pref.n_ratio/(1+pref.n_ratio);
    % Population probability density function (PDF). Left subpopulation is 
    % left unchanged, while the right subpopulation will scale according to 
    % n_ratio. Deviding by 2^2 has no effect on results, but scales PDF so 
    % its comparable to the baseline distribution with mean (0,0) and 
    % std. dev. (1,1). The subpopulations halves the std dev., thus their 
    % variance will be one-fourth.
    F = (F_l + F_r*1/pref.n_ratio)/4;
    %F = F_l*weight + F_r*(1-weight); % Population pdf

    
    %%% 2.2 Population descriptive statistics
    
    % Mean and covariance -- Johnson (1987)
    pop_mu = mu_l*weight + mu_r*(1-weight); % Population mean along each dimension; [mu_x mu_y]
    pop_sigma = sigma_l*weight + sigma_r*(1-weight) + weight*(1-weight)*(mu_l-mu_r)*(mu_l-mu_r)'; % Population covariance matrix
    %pop_sd_x = sqrt(pop_sigma(1,1)); % Population std. dev. along x-axis; sigma_x

   
    %%% 2.3 Prelocating / Preparing matrices for iteration
    
    % Many matrices have a row or coloumn for each firm. However, ex ante 
    % the total number of firms (both survining and dead) is unknown. 
    % So insure an excess amount of rows/coloumns that can accommodate all
    % possible firms. Using the initial number of firms + number of system
    % ticks (there is a max of one new firm per system tick). Ex post the 
    % excess rows and coloumns are removed.
    N_exante = pref.N + pref.iterations;
    
    % Creating a 3D matrix of firms. 
    % 1st dimension is firm, 2nd is (x,y) coordinates, and 3rd is iteration.
    xy                      = NaN(N_exante, 2, pref.iterations*pref.psi); 

    % Creating a 3D matrix of market and utility.
    % 1st dimension is x, 2nd is y, and 3rd is iteration.
    market                  = NaN(pref.resolution, pref.resolution, pref.iterations*pref.psi); % The cell/value is the closest firm
    %utility                 = NaN(pref.resolution, pref.resolution, pref.iterations); % The cell/value is utility

    shares                  = NaN(pref.iterations*pref.psi, N_exante);
    rank                    = NaN(pref.iterations*pref.psi, N_exante);
    eccentricity            = NaN(pref.iterations, N_exante);
    %mean_eccentricity       = NaN(pref.iterations, 1);
    ENP                     = NaN(pref.iterations, 1);
    %N                       = NaN(pref.iterations, 1);
    %mean_share              = NaN(pref.iterations, 1);
    age_death               = NaN(pref.iterations, N_exante);
    %perimeter               = NaN(pref.iterations, pref.N);
    %perimeter_extrema       = NaN(pref.iterations, pref.N);
    centroid                = NaN(N_exante, 2, pref.iterations*pref.psi);
    centroid_distance       = NaN(pref.iterations*pref.psi, N_exante);
    death                   = NaN(N_exante, 1);
    age                     = NaN(N_exante, 1);
    dissatisfaction         = NaN(pref.resolution, pref.resolution, pref.iterations);

    fitness                 = NaN(pref.iterations, N_exante);
    %rules                   = cell(1, N_exante);
    
    
    %%% 2.4 Firms
    
    % Active firms
    firms(1:pref.N) = 1:pref.N; % Array with active firms IDs/names.
    
    % Randomy draw initial position of firms
    % start at (0,0) move uniformly up to 3 std. dev. in random direction
    [x_0, y_0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*3 );
    xy_0 = [x_0 y_0];
    
    % Randomly draw the decision rule used by initial firms.
    % Draw decision rules with equal probability.
    rules(firms) = pref.ruleset( randi(length(pref.ruleset), pref.N,1) );
    
    % Set fitness of first firm equal 1. Set fitness of first firms equal 
    % one over the number for initial firms.
    fitness(1, firms) = 1/pref.N;
   
    % Set age of firms
    age(firms) = 0;
    
    newestfirm = pref.N; % Count the total number of firms added so far.
   
    
    %%% 2.5 Preparing for loop
    % TO-DO: There may be some preformance gains if some of the conditions 
    % are moved outside the for-loop.
    
    % For first iteration set initial firm position
    xy(firms, : ,1) = xy_0;
    
    
%% 3. EVOLUTION

for i = 1:pref.iterations*pref.psi
    
    %%% 3.1 Decision rules
    
    % Firm position determined by firm behaviour / decsion rule
    if(i>1) 
        for n = firms
            switch rules{n}
                case 'STICKER'
                    % Stay at current position
                    xy(n,:,i) = xy(n,:,i-1);

                    
                case 'AGGREGATOR'
                    % Move towards the center of its current market area.
                    % Maintain position is no customers
                    if(shares(i-1,n) > 0)
                        xy(n,:,i) = centroid(n,:,i-1);
                    else
                        xy(n,:,i) = xy(n,:,i-1);
                    end
                
                    
                case 'HUNTER'
                    % Continue same direction, if previous move proved fruitful,
                    % otherwise head the oppersite direction.
                     if(i==2)% | age(n) == 0)
                        % 0.1 std. dev. move in same direction, away from from (0,0)
                        [direction_in, ~] = cart2pol(xy(n,1,i-1), xy(n,2,i-1));
                        [dx_in, dy_in] = pol2cart(direction_in, 0.1);
                     else
                         [direction_in, ~] = cart2pol(xy(n,1,i-1)-xy(n,1,i-2), xy(n,2,i-1)-xy(n,2,i-2));
                         if(age(n) == 0)
                            direction_in = rand*2*pi; % If newly entered firm then heading is random.
                         end
                         if(shares(i-1,n) > shares(i-2,n))
                             % 0.1 std. dev. move in same direction as previous move
                             [dx_in, dy_in] = pol2cart(direction_in, 0.1);
                         else
                             % 0.1 std. dev. move in oppersite direction (90 degree + random 180 degree)
                             [dx_in, dy_in] = pol2cart(direction_in + pi/2 + rand * pi, 0.1);
                         end
                     end
                     %heading(i-1,n) = direction_in;
                     xy(n,:,i) = xy(n,:,i-1) + [dx_in dy_in];
  
                     
                case 'PREDATOR'
                    % Move towards its most profitable competitor.
                    if(rank(i-1,1) ~= n)
                        % Direction of most profitable competitor
                        [direction, rho] = cart2pol(xy(rank(i-1,1),1,i-1)-xy(n,1,i-1), xy(rank(i-1,1),2,i-1)-xy(n,2,i-1));
                        % 0.1 std. dev. move in direction of most profitable competitor, 
                        % if distance is less than 1 std. dev.
                        [dx_in, dy_in] = pol2cart(direction, 0.1);
                        if(rho>1)    
                            xy(n,:,i) = xy(n,:,i-1) + [dx_in dy_in];
                        else
                            xy(n,:,i) = xy(rank(i-1,1),:,i-1); % don't overshoot
                        end
                    else
                        xy(n,:,i) = xy(n,:,i-1);
                    end
                    
                
                %case 'EXPLORER'
                    % Survey several directions and move along the most profitable.

                
                %case 'INDUCTOR'
                    % Locate to maximise market share but subject to the predicted movements of other firms.
 					% Predictions are generated by its most reliable hypothesis. 
					% Gradually discards poorly performing hypotheses and forms new hypotheses.
                
            end
        end    
    end
    
    
    
    %%% 3.2 Propeties
    
    % Market and utility.
    [market_i, utility_i] = marketshare3(xy(:,:,i), [X(:) Y(:)]);
    market(:,:,i) = market_i;
    
    % Number of customers for each firm
    F_firm_i = arrayfun(@(firm) ...
            sum(F(find(market_i==firm))), ...
            firms);

    % Market shares
    shares(i,firms) = F_firm_i/sum(F(:)); % Calculate market share for each firm
    [~,rank(i,firms)] = sort(shares(i,firms), 'descend'); % Ranking firms according to market share
    
    % Market centroid coordinates
    for j =1:length(firms)
        firm = firms(j);
        % Index for each customer of the firm.
        idx = find(market_i==firm);
        % Each xy-coordinat within the firm's market weighted with the probability density
        centroid(firm,:,i) = ([X(idx) Y(idx)]' * F(idx))' ./ F_firm_i(j);
    end
    % Distance from each firm to the respective centroid of its market.
    centroid_distance(i,:) = diag(pdist2(xy(:,:,i), centroid(:,:,i), 'Euclidean'))';
    
    
    
    %%% 3.3 System tick
    % On system ticks determine entry and exit in market. Update memory.
    
    si = i/pref.psi; % System tick
    if ~mod(si,1) % if si integer
    
        %%% 3.3.1 Update fitness and dissatisfaction
        
        % Fitness: recusive updating algorithm where current fitness is
        % an weighted average of the fitness at the last system tick, and 
        % the current market share. Weight is the memory parameter a_f.
        
        % Dissatisfaction: recusive updating algorithm where current
        % dissatisfaction is an weighted average of the dissatisfaction at
        % the last system tick, and the current utility/distance to closest
        % firm. Weight is memory parameter a_b.
        
        if(si==1) 
            fitness(si,:) = pref.a_f + (1-pref.a_f) * shares(i,:); % fitness starts at 1
            dissatisfaction(:,:,si) = (1-pref.a_b) * utility_i; % dissatisfaction starts at 0
        else
            fitness(si,:) = pref.a_f * fitness(si-1,:) + (1-pref.a_f) * shares(i,:);
            dissatisfaction(:,:,si) = pref.a_b * dissatisfaction(:,:,si-1) + (1-pref.a_b) * utility_i;
        end
        dissatisfaction_i = dissatisfaction(:,:,si);
        
        % Set the age of firms
        age(firms) = age(firms)+1;
        
        
        
        %%% 3.3.2 Summary variables (in loop)
        
        % Utility
        %utility(:,:,si) = utility_i;
        
        % Eccentricity 
        % Firm euclidean distance from center/mean of voter distribution).
         eccentricity(si,:) = pdist2(xy(:,:,i), pop_mu, 'euclidean'); % in std. dev.
         %mean_eccentricity(si,:) = mean(eccentricity(si,:), 'omitnan');
        
        % Effective number of firms (ENP) -- Lever and Sergenti (2011), Laakso and Taagepera (1979)
         ENP(si,:) = sum(F(:))^2 / sum(F_firm_i(:).^2);
        % Herfindahl-Hirschman Index (HHI) -- Eiselt and Marianov (2011)
        % HHI(i,:) = sum( ( F_firm/sum(F(:)) ).^2 ); % HHI = 1/ENP
        
        % Number of surviving parties
        %N(si,:) = length(firms);
        
        % Mean market share
        %mean_share(si,:) = mean(shares(i,:), 'omitnan');
        
        % firm age at death
        %mean_age(si) = mean(age( ~isnan(death) ));
        age_death(si, ~isnan(death) ) = age( ~isnan(death) );
        
        % Misery
        % Quadratic loss function / representativeness
        %misery(si,:) = sum(utility_i(:).*F(:))/sum(F(:));
        
        % Market perimeter
        % market_props_i = regionpropsext(market(:,:,i), 'Extrema', 'Perimeter');
        % perimeter(i,:) = cat(1,market_props_i.Perimeter)';
        % perimeter_extrema(i,:) = cat(1,market_props_i.ExtremaPerimeter)';
        
        
        
        %%% 3.3.3 Exit
            
        % Firms with fitness below the de facto threshold exits market.
        deadfirms = find(fitness(si,:) < pref.tau);
        
        if any(deadfirms)
            % Need at least two remaning active firms.
            if( length(firms)-length(deadfirms) < 2 ) 
                % If fewer than two firms have fitness above threshold, then
                % randomly draw firms with fitness below threshold so exactly
                % two firms remain.
                deadfirms = deadfirms( randperm(length(deadfirms), length(firms)-2) );
            end
        end
        
        % Remove exiting firms from the list of active firms.
        firms( ismember(firms,deadfirms) ) = [];
        % Record time of death.
        death(deadfirms) = si;
        
        
        %%% 3.3.4 Entry
        
        % Randomly pick a patch within 3 std. dev. of (0,0) / in the mask. 
        % Where the probability of selecting the patch is proportinal to
        % the dissatisfaction and is weighted by the market share. Beta is 
        % the birth parameter that meassures sensitivity.
        P_mask = pref.beta * dissatisfaction_i(mask==1) .* F(mask==1)/sum(F(:));
        % P_mask = P_mask/sum(P_mask); % The probability should sum to 1.
        % Using a binary search algorithm. Create an arrray of cummulative 
        % probabilities, then draw a random number Uniform(0,1) and select 
        % the first cell in array/bin in cummulative probabilities that 
        % contains the drawn number. Returns the index of the patch inside 
        % the mask.
        i_mask = ( find(rand<cumsum(P_mask), 1) );
        
        if any(i_mask) % If new firm drawn
            % The ID of newest firm
            newestfirm = newestfirm+1;
            % Add to list of active firms
            firms = [firms newestfirm];
            % Coordinates for new firm
            xy(newestfirm, :, i) = [X_mask(i_mask) Y_mask(i_mask)];
            % Randomly draw decision rule of new firm (equal probability).
            rules(newestfirm) = pref.ruleset( randi(length(pref.ruleset)) );
            % Fitness of new firm
            fitness(si, newestfirm) = 1/length(firms);
            % Age of new firm
            age(newestfirm) = 0; 
        end
        
        
    end
    
end


%%% 3.4 Ex post adjustments

% Remove redundant rows and coloums from matrix.
xy(newestfirm+1:end, :, :)              = [];
shares(:, newestfirm+1:end)             = [];
fitness(:, newestfirm+1:end)            = [];
death(newestfirm+1:end, :)              = [];
age(newestfirm+1:end, :)                = [];
eccentricity(:, newestfirm+1:end)       = [];
age_death(:, newestfirm+1:end)          = [];
centroid(newestfirm+1:end, :, :)        = [];
centroid_distance(:, newestfirm+1:end)  = [];

% Mapping firms and rules
% Find the index corresponding to a decision rule for each firm. Length of 
% ruleidx is equal to the total number of firms. The value in each cell 
% corresponds to a decision rule in the ruleset, ie. 1 = first rule, 
% 2 = second rule, etc.
[~,ruleidx] = ismember(rules, pref.ruleset);

% Subset with the shares at system ticks only.
shares_sys = shares(pref.psi:pref.psi:pref.iterations*pref.psi,:);

% Colors
% Unique color for each firm.
color = linspecer(newestfirm); % Generates a color for each of the firms that is evenly space out / distinguishable.
% Unique color for each decision rule.
colorrule = linspecer(length(pref.ruleset)); % Generates a color for each of the firms that is evenly space out / distinguishable.
% Maps the decision rule colors to each firm, such that the firm is colored 
% according to its decision rule.
colorrulefirm = colorrule(ruleidx,:);


%%% 3.5 Summary variables (outside loop)
% If possible calculate the summary variables outside the for loop as
% opposed to inside the for loop, for the sake of preformance/speed.

% Mean eccentricity
mean_eccentricity = mean(eccentricity, 2, 'omitnan');
% Mean share
mean_share = mean(shares_sys, 2, 'omitnan');
% Number of surviving firms
N = sum(~isnan(shares_sys), 2, 'omitnan');
% Mean age of firm at exit
mean_age_death = mean(age_death, 2, 'omitnan');

% Calculate the summary variables for each rule.
mean_eccentricity_rule      = NaN(pref.iterations, length(pref.ruleset));
mean_share_rule             = NaN(pref.iterations, length(pref.ruleset));
N_rule                      = NaN(pref.iterations, length(pref.ruleset));
mean_age_death_rule         = NaN(pref.iterations, length(pref.ruleset));
for r=1:length(pref.ruleset)
    % Select all the firms that use the same rule
    rfirms = find(ruleidx==r);
    
    % Calculate the mean eccentricity for each rule
    mean_eccentricity_rule(:,r) = mean(eccentricity(:,rfirms), 2, 'omitnan');
    % Calculate the mean market share for each rule
    mean_share_rule(:,r) = mean(shares_sys(:,rfirms), 2, 'omitnan');
    % Calculate the number of surviving firms for each rule
    N_rule(:,r) = sum(~isnan(shares_sys(:,rfirms)), 2, 'omitnan');
    % Calculate the mean age of firm at death
    mean_age_death_rule(:,r) = mean(age_death(:,rfirms), 2, 'omitnan');
end
    

