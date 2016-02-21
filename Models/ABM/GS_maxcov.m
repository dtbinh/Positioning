%% Grid Sweep -- All-maxcov
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
% The All-maxcov model can be representent as a stationary Markov chain / 
% time-homogenous Markov chain since the transition probabilities are
% independent of time / iteration.
%
% The All-maxcov model does not contain random components (besides the 
% initial position of firms). 
%
% So the Markov chain is stationary and deterministic. Thus process 
% converges to a single state, but the process is nonergodic (so the 
% limiting distribution depends on the initial conditions).
%
% Thus we use the ensemble-average to calculate representative estimates.
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
    % Run through all firm sizes. Run 10 repetitions of each size with 200 iterations.
    test.N = 2:12;
    
    pref.seed = rng(22589026, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 200; % Number of iterations
    pref.repetitions = 10; % Number of repetitions of run
    pref.psi = 1; % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
    pref.M = 1; % Number of condition/forecast rules that each firm holds.
    pref.a_a = 1-1/75; % Accuracy mememory paramenter.
    pref.C = 0.005; % Cost of specificity.
    pref.crossover = 0.3; % Probability that the offspring condition/forecast rule is created by crossover operations (rather than mutation). 
    
    % Mean of subpopulation
    pref.mu = 0; % No polarization
    % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.n_ratio = 1; % Equal subpopulation size

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    export_param = NaN(pref.repetitions, 4, pref.runs); % The four extra coloumns are for: repetition number, N, mu, n_ratio
    
    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Number of firms
        pref.N = test.N(run);
        
        % Decision rules: All-aggregator
        pref.rules = repmat( {'MAXCOV'}, 1, pref.N);
        
        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, ~] = ABM_ind(pref);

            % Store summary variables from each run
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
            
            export_param(rep,:,run) = [pref.N pref.mu pref.n_ratio rep ];
        end

    end
    close(h);
    
    % Save summary variables
    % Reshape the data to required format before exporting
    export_param_fmt = reshape(permute(export_param,[1 3 2]), [pref.repetitions*pref.runs, 4]);
    export_mean_eccentricity    = [export_param_fmt reshape( permute(data_mean_eccentricity,[1 3 2]), [pref.repetitions*pref.runs, pref.iterations] ) ];
    export_ENP                  = [export_param_fmt reshape( permute(data_ENP,[1 3 2]), [pref.repetitions*pref.runs, pref.iterations] ) ];
    % Export
    csvwrite(strcat('data/data_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), export_mean_eccentricity);
    csvwrite(strcat('data/data_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), export_ENP);

    % Identify runs that require most iterations visually using filters:
    % <https://github.com/jsekamane/filter-time>
    
    % Preliminary conclusion
    % Parameter settings that require most iterations
    % mean eccentricity: Converge within 25 iterations. Most slowly for many firms (N: 11-12).
    % ENP: Almost all converge within 25 iterations, however some runs
    % diverage again (especially for N: 9-12).
    
    
    
%% 3. Final model    
% Using ensemble-average to calculate estimates. Using 99 burnin iterations.
% Starting with 1000 repetitions. If this is not enough to pass our (five) 
% sample size checks, then increase the number of repetition to 2000, etc.

    pref.burnin = 99; % Number of iterations before we have burnin 
    
    % Run through all firm sizes.
    test.N = 2:12;

    pref.seed = rng(93433173, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = pref.burnin + 1; % Number of post-burnin iterations
    pref.repetitions = 1000; % Number of repetitions of run
    pref.psi = 1; % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
    pref.M = 1; % Number of condition/forecast rules that each firm holds.
    pref.a_a = 1-1/75; % Accuracy mememory paramenter.
    pref.C = 0.005; % Cost of specificity.
    pref.crossover = 0.3; % Probability that the offspring condition/forecast rule is created by crossover operations (rather than mutation). 


    % Mean of subpopulation
    pref.mu = 0; % No polarization
    % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.n_ratio = 1; % Equal subpopulation size

    % Creating empty matrixes for summary variable
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_mean_representation = NaN(pref.repetitions, pref.iterations, pref.runs);

    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Number of firms
        pref.N = test.N(run);

        % Decision rules: All-aggregator
        pref.rules = repmat( {'MAXCOV'}, 1, pref.N);

        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;

            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, o_mean_representation] = ABM_ind(pref);

            % Store summary variables from each run
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
            data_mean_representation(rep,:,run) = o_mean_representation';
        end

    end
    close(h);

    % Post-burnin iterations
    data_burnin_mean_eccentricity = data_mean_eccentricity( :, pref.burnin+1:end, :);
    data_burnin_ENP = data_ENP( :, pref.burnin+1:end, :);
    data_burnin_mean_representation = data_mean_representation( :, pref.burnin+1:end, :);
    
    % Ensemble average estimate
    est_mean_eccentricity = mean(data_burnin_mean_eccentricity, 1);
    est_ENP = mean(data_burnin_ENP, 1);
    est_mean_representation = mean(data_burnin_mean_representation, 1);
    
    % Ensemble average estimate standard deviation
    est_std_mean_eccentricity = std(data_burnin_mean_eccentricity, 0, 1);
    est_std_ENP = std(data_burnin_ENP, 0, 1);
    est_std_mean_representation = std(data_burnin_mean_representation, 0, 1);
    
    % Ensemble average estimate standard error
    est_se_mean_eccentricity = est_std_mean_eccentricity ./ sqrt(pref.repetitions);
    est_se_ENP = est_std_ENP ./ sqrt(pref.repetitions);
    est_se_mean_representation = est_std_mean_representation ./ sqrt(pref.repetitions);
    
    
    
    %%% 3.1 Check 1 -- Mapping / R-Hat statistics
    % Not applicable. Because we don't use time averages (only 1 post-burnin iteration).
    
    
    %%% 3.2 Check 2 -- Convergence / F-test p-value
    % Not applicable. There is no prior expectation for ENP or mean eccentricty for the all-maxcov model. 
    
    
    %%% 3.3 Check 3 -- Power zero
    % Calculating the power of the t-test where H_0: estimate is equal zero, H_A: different from zero
    power_zero_mean_eccentricity = powerzero(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.repetitions);
    power_zero_ENP = powerzero(squeeze(est_ENP), squeeze(est_std_ENP), pref.repetitions);
    power_zero_mean_representation = powerzero(squeeze(est_mean_representation), squeeze(est_std_mean_representation), pref.repetitions);
    % The power of the t-test should be at least 0.8
    
    
    %%% 3.4 Check 4 -- Power difference (grid sweeps)
    % Calculating the power of the two-sample t-test where H_0: estimate is equal estimate from adjacent grid, H_A: different
    power_diff_mean_eccentricity = powerdiff(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.repetitions);
    power_diff_ENP = powerdiff(squeeze(est_ENP), squeeze(est_std_ENP), pref.repetitions);
    power_diff_mean_representation = powerdiff(squeeze(est_mean_representation), squeeze(est_std_mean_representation), pref.repetitions);
    % The power of the t-test should be at least 0.8
   
    
    %%% 3.5 Check 5 -- SE/SD ratio
    % 
    SESD_ratio_mean_eccentricity = squeeze(est_se_mean_eccentricity) ./ squeeze(est_std_mean_eccentricity);
    SESD_ratio_ENP = squeeze(est_se_ENP) ./ squeeze(est_std_ENP);
    SESD_ratio_mean_representation = squeeze(est_se_mean_representation) ./ squeeze(est_std_mean_representation);
    % Have all summary variables been estimated with the same level of
    % precisions? This is (trivially) satisfied when the summary variables
    % have been estimated using the same number of repetitions,
    % since the SE/SD ratio simply returns one over the squareroot of the 
    % number of repetitions. 
   
   
   
    %%% 3.6 Export results
    % Format table before saving file
    export_mean_eccentricity = table(test.N', squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), squeeze(est_se_mean_eccentricity), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_mean_eccentricity, power_diff_mean_eccentricity, SESD_ratio_mean_eccentricity, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_ENP = table(test.N', squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_ENP, power_diff_ENP, SESD_ratio_ENP, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_mean_representation = table(test.N', squeeze(est_mean_representation), squeeze(est_std_mean_representation), squeeze(est_se_mean_representation), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_mean_representation, power_diff_mean_representation, SESD_ratio_mean_representation, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    %export_ENP.Properties.Description = ['burnin ' num2str(pref.burnin) ' iterations ' num2str(pref.iterations)];
    % Save file
    writetable(export_mean_eccentricity, strcat('data/GS_maxcov_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_r', num2str(pref.repetitions), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/GS_maxcov_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_r', num2str(pref.repetitions), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');
    writetable(export_mean_representation, strcat('data/GS_maxcov_mean_representation_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_r', num2str(pref.repetitions), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');

    % Preliminary conclusion
    %
    % ...
    
    
