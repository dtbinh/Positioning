%% Monte Carlo Parameterization -- All-maxcov
% version 0.01
% Jonas K. Sekamane
%
% Free paramter: N, mu, n_ratio
%
% Exogenous number of firms.
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
% Graphically.
    
    %%% 2.1 Idenfying runs that require most iteration
    % Idenfying from 500 runs with 200 iterations and 1 repetition.
    pref.seed = rng(34983315, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 200; % Number of iterations
    pref.repetitions = 1; % Number of repetitions of run
    pref.psi = 1; % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
    pref.M = 100; % Number of condition/forecast rules that each firm holds.
    pref.a_a = 1-1/75; % Accuracy mememory paramenter.
    pref.C = 0.005; % Cost of specificity.
    pref.crossover = 0.3; % Probability that the offspring condition/forecast rule is created by crossover operations (rather than mutation).    

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    export_param = NaN(pref.repetitions, 4, pref.runs); % The four extra coloumns are for: repetition number, N, mu, n_ratio

    poolobj = parpool('local',4)
    parfor_progress(pref.runs);
    parfor run=1:pref.runs
        pref2 = pref;
        pref2.run = run;

        % Randomly sample parameters from uniform distribution of the range of the parameters
        % Number of firms
        pref2.N = randi([2,12]); % N in [2,12]
        % Mean of subpopulation
        pref2.mu = rand * 1.5; % mu in [0,1.5]
        % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
        pref2.n_ratio = 1 + rand; % n_ratio in [1,2]
        
        % Decision rules: All-hunter
        pref2.rules = repmat( {'MAXCOV'}, 1, pref2.N);
        
        % Repetitions
        data_mean_eccentricity_run  = NaN(pref.repetitions, pref.iterations);
        data_ENP_run                = NaN(pref.repetitions, pref.iterations);
        export_param_run            = NaN(pref.repetitions, 4);
        for rep=1:pref.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, ~] = ABM_ind(pref3);

            % Store summary variables from each run
            data_mean_eccentricity_run(rep,:) = o_mean_eccentricity';
            data_ENP_run(rep,:) = o_ENP';
            
            export_param_run(rep,:) = [pref3.N pref3.mu pref3.n_ratio rep ];
        end
        data_mean_eccentricity(:,:,run) = data_mean_eccentricity_run;
        data_ENP(:,:,run)               = data_ENP_run;
        export_param(:,:,run)           = export_param_run;

        parfor_progress;
    end
    parfor_progress(0);
    delete(poolobj)

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
    % mean eccentricity: Converge within 50-75 iterations. Most slowly for many firms (N: 11-12).
    % ENP: Almost all converge within 25-75 iterations, however some runs
    % diverage again (especially for N: 9-12).
    
    
%% 3. Final model
% Using ensemble-average to calculate estimates. Using 150 burnin iterations.
% Starting with 50 repetitions. If this is not enough to pass our (five) 
% sample size checks, then increase the number of repetition to 100, etc.
    
    pref.burnin = 150; % Number of iterations before we have burnin 
    pref.seed = rng(28740504, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = pref.burnin + 1; % Number of post-burnin iterations
    pref.repetitions = 50; % Number of repetitions of run
    pref.psi = 1; % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
    pref.M = 100; % Number of condition/forecast rules that each firm holds.
    pref.a_a = 1-1/75; % Accuracy mememory paramenter.
    pref.C = 0.005; % Cost of specificity.
    pref.crossover = 0.3; % Probability that the offspring condition/forecast rule is created by crossover operations (rather than mutation).    
    
    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_mean_representation = NaN(pref.repetitions, pref.iterations, pref.runs);
    export_param = NaN(pref.repetitions, 4, pref.runs); % The four extra coloumns are for: repetition number, N, mu, n_ratio

    poolobj = parpool('local',4)
    parfor_progress(pref.runs);
    parfor run=1:pref.runs
        pref2 = pref;
        pref2.run = run;

        % Randomly sample parameters from uniform distribution of the range of the parameters
        % Number of firms
        pref2.N = randi([2,12]); % N in [2,12]
        % Mean of subpopulation
        pref2.mu = rand * 1.5; % mu in [0,1.5]
        % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
        pref2.n_ratio = 1 + rand; % n_ratio in [1,2]
        
        % Decision rules: All-aggregator
        pref2.rules = repmat( {'MAXCOV'}, 1, pref2.N);
        
        % Repetitions
        data_mean_eccentricity_run   = NaN(pref.repetitions, pref.iterations);
        data_ENP_run                 = NaN(pref.repetitions, pref.iterations);
        data_mean_representation_run = NaN(pref.repetitions, pref.iterations);
        export_param_run             = NaN(pref.repetitions, 4);
        for rep=1:pref2.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, o_mean_representation] = ABM_ind(pref3);
            
            % Store summary variables from each run and each repetition
            data_mean_eccentricity_run(rep,:) = o_mean_eccentricity';
            data_ENP_run(rep,:) = o_ENP';
            data_mean_representation_run(rep,:) = o_mean_representation';
            
            export_param_run(rep,:) = [pref3.N pref3.mu pref3.n_ratio rep ];
        end
        data_mean_eccentricity(:,:,run)   = data_mean_eccentricity_run;
        data_ENP(:,:,run)                 = data_ENP_run;
        data_mean_representation(:,:,run) = data_mean_representation_run;
        export_param(:,:,run)             = export_param_run;

        parfor_progress;
    end
    parfor_progress(0);
    delete(poolobj)
    
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

    
    
    %%# 3.1 Check 1 -- Mapping / R-Hat statistics
    % Not applicable. Because we don't use time averages (only 1 post-burnin iteration).
    
    
    %%# 3.2 Check 2 -- Convergence / F-test p-value
    % Does not apply. There is no prior expectation for ENP or mean eccentricty for the all-aggregator model. 
    
    
    %%# 3.3 Check 3 -- Power zero
    % Calculating the power of the t-test where H_0: estimate is equal zero, H_A: different from zero
    %power_zero_mean_eccentricity = powerzero(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.repetitions);
    %power_zero_ENP = powerzero(squeeze(est_ENP), squeeze(est_std_ENP), pref.repetitions);
    % The power of the t-test should be at least 0.8
    
    
    %%# 3.4 Check 4 -- Power difference (grid sweep)
    % Not applicable. Using Monte Carlo parameterisation and not the grid 
    % sweep method, so there is no adjacent grid/cell to compare with.
    
    %%# 3.5 Check 5 -- SE/SD ratio
    % 
    SESD_ratio_mean_eccentricity = squeeze(est_se_mean_eccentricity) ./ squeeze(est_std_mean_eccentricity);
    SESD_ratio_ENP = squeeze(est_se_ENP) ./ squeeze(est_std_ENP);
    SESD_ratio_mean_representation = squeeze(est_se_mean_representation) ./ squeeze(est_std_mean_representation);
    % Have all summary variables been estimated with the same level of
    % precisions? This is (trivially) satisfied when the summary variables
    % have been estimated using the same number of post-burnin iterations,
    % since the SE/SD ratio simply returns one over the squareroot of the 
    % number of post-burnin iterations. 

    % Reshape the data to required format before exporting
    export_param_fmt = reshape(permute(export_param(1,:,:),[1 3 2]), [pref.runs, 4]);
    
    %%% 3.6 Export results
    % Format table before saving file
    export_mean_eccentricity = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), squeeze(est_se_mean_eccentricity), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), SESD_ratio_mean_eccentricity, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_ENP = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), SESD_ratio_ENP, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_mean_representation = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_mean_representation), squeeze(est_std_mean_representation), squeeze(est_se_mean_representation), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), SESD_ratio_mean_representation, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    writetable(export_mean_eccentricity, strcat('data/MCP_maxcov_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/MCP_maxcov_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_mean_representation, strcat('data/MCP_maxcov_mean_representation_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');


