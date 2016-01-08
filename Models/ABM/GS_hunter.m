%% Grid Sweep -- All-hunter
% version 0.01
% Jonas K. Sekamane
%
% Free paramter: N
%
% Exogenous number of firms.
% No polarization of subpopulation.
% Equal subpopulation size.
%
% Inspired in part by: 
%   Lever and Sergenti (2011)



%% 1. Introduction
% The All-hunter model can be representent as a stationary Markov chain / 
% time-homogenous Markov chain since the transition probabilities are
% independent of time / iteration.
%
% The All-hunter model contains random components. More specifically, if
% the previous move of the hunter-firm did not increase it's market share, 
% then it will head in a random oppersite direction (randomly drawn from a 
% 180 dregree arc oppersite of the current direction).
%
% So the Markov chain is stationary and stocastic, which implies that the
% Markov chain is ergodic (the process tends to a unquie limiting 
% distribution that is independent of the initial conditions).
%
% And since the state space is not too large and since the (off-diagonal) 
% transition probabilities are fairly high, we can use the time-average to 
% calculate representative estimates.
%

clearvars;

%% 2. Determining Burn-In
% Summary of procedure:
% # Specify Markov Chain representation; Specify vector of state space and summary variables.
% # Identify those runs that require most iterations to burn in.
% # Run several test repetitions of these; Use "second halves" to calculate R-hat (should be below 1.05) and calculate the estimate.
% # From estimate determine maximum burn-in of test repetition: within 1. std. dev. of estimate. 


    %%% 2.1 Markov chain representation
    % TO-DO: Specify Markov Chain
    %
    % Summary variables: mean eccentricity, effective number of firms (ENP)


    
    %%% 2.2 Idenfying runs that require most iteration
    % Run through all firm sizes. Run 5 repetitions of each size with 200 iterations.
    test.N = 2:12;
    
    pref.seed = rng('default'); % Seed such that the randomly generated results are repeatable
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 200; % Number of iterations
    pref.repetitions = 5; % Number of repetitions of run
    
    % Mean of subpopulation
    pref.mu = 0; % No polarization
    % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.n_ratio = 1; % Equal subpopulation size

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = zeros(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = zeros(pref.repetitions, pref.iterations, pref.runs);
    export_temp = zeros(pref.repetitions, 4, pref.runs); % The four extra coloumns are for: repetition number, N, mu, n_ratio
    
    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Number of firms
        pref.N = test.N(run);
        
        % Decision rules: All-hunter
        pref.rules = repmat( {'HUNTER'}, 1, pref.N);
        
        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP] = ABM(pref);

            % Store summary variables from each run
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
            
            export_temp(rep,:,run) = [pref.N pref.mu pref.n_ratio rep ];
        end

    end
    close(h);

    % Save summary variables
    % Reshape the data to required format before exporting
    export_temp_fmt = reshape(permute(export_temp,[1 3 2]), [pref.repetitions*pref.runs, 4]);
    export_mean_eccentricity    = [export_temp_fmt reshape( permute(data_mean_eccentricity,[1 3 2]), [pref.repetitions*pref.runs, pref.iterations] ) ];
    export_ENP                  = [export_temp_fmt reshape( permute(data_ENP,[1 3 2]), [pref.repetitions*pref.runs, pref.iterations] ) ];
    % Export
    csvwrite(strcat('data/data_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), export_mean_eccentricity);
    csvwrite(strcat('data/data_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), export_ENP);

    % Identify runs that require most iterations visually using filters:
    % <https://github.com/jsekamane/filter-time>
    
    % Preliminary conclusion
    % Parameter settings that require most iterations
    % mean eccentricity: Few firms (N: 2-3).
    % ENP: Many firms (N 10-12).
    
    
    
    %%% 2.3 Test repetitions
    test.N = [2 3 10 11 12];

    pref.seed = rng(88887852, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 1000; % Number of iterations
    pref.repetitions = 5; % Number of repetitions of run
    
    % Mean of subpopulation
    pref.mu = 0; % No polarization
    % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.n_ratio = 1; % Equal subpopulation size

    % Creating empty matrixes for summary variable
    data_mean_eccentricity = zeros(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = zeros(pref.repetitions, pref.iterations, pref.runs);
    
    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Number of firms
        pref.N = test.N(run);
        
        % Decision rules: All-hunter
        pref.rules = repmat( {'HUNTER'}, 1, pref.N);
        
        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP] = ABM(pref);

            % Store summary variables from each run
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
        end

    end
    close(h);
    
    % Second half
    data_2nd_mean_eccentricity = data_mean_eccentricity( :, size(data_mean_eccentricity,2)/2+1:end, :);
    data_2nd_ENP = data_ENP( :, size(data_ENP,2)/2+1:end, :);
    
    % Time average estimate (of second half)
    est_mean_eccentricity = mean(data_2nd_mean_eccentricity, 2);
    est_ENP = mean(data_2nd_ENP, 2);
    
    % Time average estimate standard deviation (of second half)
    est_std_mean_eccentricity = std(data_2nd_mean_eccentricity, 0, 2);
    est_std_ENP = std(data_2nd_ENP, 0, 2);
    
    % R-hat statistics / potential scale reduction factor (PSRF) -- Brooks and Gelman (1998)    
    rhat_mean_eccentricity = rhat( data_2nd_mean_eccentricity );
    rhat_ENP = rhat( data_2nd_ENP );
    % Should be below 1.05
    
    % Preliminary conclusion
    % R-hat statistics are below 1.05. The maximum value is 1.0340 and 1.0399.
    
    
    
    %%% 2.4 Maximum burn-in of test repetitions
    % within 1. std. dev. of estimate.
    
    burnin_ENP = burnin(data_ENP, est_ENP, est_std_ENP);
    burnin_mean_eccentricity = burnin(data_mean_eccentricity, est_mean_eccentricity, est_std_mean_eccentricity);
    
    % The maximum burn-in number of all test repetitions.
    max( burnin_ENP(:) )
    max( burnin_mean_eccentricity(:) )
    
    % Preliminary conclusion
    % Maximum burn-in for ENP is 77 [Test #2 (N=3) in rep #1].
    % Maximum burn-in for Mean Eccentricity is 120 [Test #2 (N=3) in rep #1].
    % So to be on the save side I will use a burn-in value of 150 for the all-hunter model.
    
    
    
%% 3. Final model    
% Using time-average to calculate estimates, so 1 single repetition. 
% Starting with 100 post-burnin iterations. 
% If this is not enough to pass our (five) sample size checks, then increase 
% the number of post-burnin iterations to 1000, etc.

    pref.burnin = 150; % Number of iterations before we have burnin 
    
    % Run through all firm sizes.
    test.N = 2:12;

    pref.seed = rng(80243125, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = pref.burnin + 1000; % Number of post-burnin iterations
    pref.repetitions = 1; % Number of repetitions of run

    % Mean of subpopulation
    pref.mu = 0; % No polarization
    % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.n_ratio = 1; % Equal subpopulation size

    % Creating empty matrixes for summary variable
    data_mean_eccentricity = zeros(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = zeros(pref.repetitions, pref.iterations, pref.runs);

    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Number of firms
        pref.N = test.N(run);

        % Decision rules: All-hunter
        pref.rules = repmat( {'HUNTER'}, 1, pref.N);

        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;

            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP] = ABM(pref);

            % Store summary variables from each run
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
        end

    end
    close(h);

    % Post-burnin iterations
    data_burnin_mean_eccentricity = data_mean_eccentricity( :, pref.burnin+1:end, :);
    data_burnin_ENP = data_ENP( :, pref.burnin+1:end, :);
    
    % Time average estimate (of post-burnin iterations)
    est_mean_eccentricity = mean(data_burnin_mean_eccentricity, 2);
    est_ENP = mean(data_burnin_ENP, 2);
    
    % Time average estimate standard deviation (of post-burnin iterations)
    est_std_mean_eccentricity = std(data_burnin_mean_eccentricity, 0, 2);
    est_std_ENP = std(data_burnin_ENP, 0, 2);
    
    % Time average estimate standard error (of post-burnin iterations)
    est_se_mean_eccentricity = est_std_mean_eccentricity ./ sqrt(pref.iterations - pref.burnin);
    est_se_ENP = est_std_ENP ./ sqrt(pref.iterations - pref.burnin);
    
    
    
    %%% 3.1 Check 1 -- Mapping / R-Hat statistics
    % TO-DO: ??? There is only 1 repetition here ???
    % R-hat statistics / potential scale reduction factor (PSRF) -- Brooks and Gelman (1998)    
    %rhat_mean_eccentricity = rhat( data_burnin_mean_eccentricity );
    %rhat_ENP = rhat( data_burnin_ENP );
    % Should be below 1.05
    
    
    %%% 3.2 Check 2 -- Convergence / F-test p-value
    % Not applicable. There is no prior expectation for ENP or mean eccentricty for the all-hunter model. 
    
    
    %%% 3.3 Check 3 -- Power zero
    % Calculating the power of the t-test where H_0: estimate is equal zero, H_A: different from zero
    power_zero_mean_eccentricity = powerzero(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.iterations-pref.burnin);
    power_zero_ENP = powerzero(squeeze(est_ENP), squeeze(est_std_ENP), pref.iterations-pref.burnin);
    % The power of the t-test should be at least 0.8
    
    
    %%% 3.4 Check 4 -- Power difference (grid sweeps)
    % Calculating the power of the two-sample t-test where H_0: estimate is equal estimate from adjacent grid, H_A: different
    power_diff_mean_eccentricity = powerdiff(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.iterations-pref.burnin);
    power_diff_ENP = powerdiff(squeeze(est_ENP), squeeze(est_std_ENP), pref.iterations-pref.burnin);
    % The power of the t-test should be at least 0.8
   
    
    %%% 3.5 Check 5 -- SE/SD ratio
    % 
    SESD_ratio_mean_eccentricity = squeeze(est_se_mean_eccentricity) ./ squeeze(est_std_mean_eccentricity);
    SESD_ratio_ENP = squeeze(est_se_ENP) ./ squeeze(est_std_ENP);
    % Have all summary variables been estimated with the same level of
    % precisions? This is (trivially) satisfied when the summary variables
    % have been estimated using the same number of post-burnin iterations,
    % since the SE/SD ratio simply returns one over the squareroot of the 
    % number of post-burnin iterations. 
   
   
   
    %%% 3.6 Export results
    % Format table before saving file
    export_mean_eccentricity = table(test.N', squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), squeeze(est_se_mean_eccentricity), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_mean_eccentricity, power_diff_mean_eccentricity, SESD_ratio_mean_eccentricity, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_ENP = table(test.N', squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_ENP, power_diff_ENP, SESD_ratio_ENP, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    %export_ENP.Properties.Description = ['burnin ' num2str(pref.burnin) ' iterations ' num2str(pref.iterations)];
    % Save file
    writetable(export_mean_eccentricity, strcat('data/GS_hunter_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/GS_hunter_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');

    % Preliminary conclusion
    %
    % With 100 post-burnin iterations the model satisfingly passes check 3 
    % and check 5. The model only passes check 4 for ENP, but not for mean 
    % eccentricity. Specifically for N = 6,8,12 the powers are below 0.8 
    % (and as low as 0.34 for N=8). Check 2 is not applicable. 
    % TO-DO: How to calculate R-hat when there is only 1 repetition of
    % model ???
    %
    % With 1000 post-burnin iterations the model satisfingly passes check 3 
    % and check 5. The model only passes check 4 for ENP, but not for mean 
    % eccentricity. Specifically for N = 6,9 the powers are below 0.8 
    % (and as low as 0.05 for N=6). Check 2 is not applicable. 
    % TO-DO: How to calculate R-hat when there is only 1 repetition of
    % model ???
    
    
    
