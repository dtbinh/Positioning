%% Monte Carlo Parameterization -- All-hunter
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


    %%% 2.1 Markov Chain representation
    % TO-DO: Specify Markov Chain
    %
    % Summary variables: mean eccentricity and ENP

    
    
    %%% 2.2 Idenfying runs that require most iteration
    % Idenfying from 500 runs with 200 iterations and 1 repetition.
    pref.seed = rng('default'); % Seed such that the randomly generated results are repeatable
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 200; % Number of iterations
    pref.repetitions = 1; % Number of repetitions of run

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
        pref2.rules = repmat( {'HUNTER'}, 1, pref2.N);
        
        % Repetitions
        data_mean_eccentricity_run  = NaN(pref.repetitions, pref.iterations);
        data_ENP_run                = NaN(pref.repetitions, pref.iterations);
        export_param_run            = NaN(pref.repetitions, 4);
        for rep=1:pref.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, ~] = ABM(pref3);

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
    % mean eccentricity: Few firms (N: 2-3), and large difference in ideal  point of subpopulation (mu close to 1.5). No visible effect from relative size of population (n_ratio).
    % ENP: Many firms (N 10-12), and large difference in ideal point (mu close to 1.5). No visible effect from relative size of population (n_ratio).
    
 
    
    %%% 2.3 Test repetitions
    % Making 10 runs with follow parameters setings
    test.N       = [2 2 3 3 10 10 11 11 12 12];
    test.mu      = [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5];
    test.n_ratio = [1 2 1 2 1 2 1 2 1 2];
    %test.N       = [2 3 10 11 12];
    %test.mu      = [0 0 0 0 0];
    %test.n_ratio = [1 1 1 1 1];
    
    pref.seed = rng(85387475, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now'); % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 1000; % Number of iterations
    pref.repetitions = 5; % Number of repetitions of run

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);

    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Number of firms
        pref.N = test.N(run);
        % Mean of subpopulation
        pref.mu = test.mu(run);
        % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
        pref.n_ratio = test.n_ratio(run);
        
        % Decision rules: All-hunter
        pref.rules = repmat( {'HUNTER'}, 1, pref.N);
        
        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, ~] = ABM(pref);
            
            % Store summary variables from each run and each repetition
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
        end

    end
    close(h);
    
    % Second half
    data_2nd_mean_eccentricity = data_mean_eccentricity( :, size(data_mean_eccentricity,2)/2+1:end, :);
    data_2nd_ENP = data_ENP( :, size(data_ENP,2)/2+1:end, :);
    
    % Time average estimate (second half)
    est_mean_eccentricity = mean(data_2nd_mean_eccentricity, 2);
    est_ENP = mean(data_2nd_ENP, 2);
    
    % Time average estimate std (second half)
    est_std_mean_eccentricity = std(data_2nd_mean_eccentricity, 0, 2);
    est_std_ENP = std(data_2nd_ENP, 0, 2);
    
    % R-hat statistics / potential scale reduction factor (PSRF) -- Brooks and Gelman (1998)    
    rhat_mean_eccentricity = rhat( data_2nd_mean_eccentricity );
    rhat_ENP = rhat( data_2nd_ENP );
    % Should be below 1.05

    
    %%% 2.4 Maximum burn-in of test repetitions
    % within 1. std. dev. of estimate.
    
    burnin_ENP = burnin(data_ENP, est_ENP, est_std_ENP);
    burnin_mean_eccentricity = burnin(data_mean_eccentricity, est_mean_eccentricity, est_std_mean_eccentricity);
    
    % The maximum burn-in number of all test repetitions.
    max( burnin_ENP(:) )
    max( burnin_mean_eccentricity(:) )
    
    % Preliminary conclusion
    % Maximum burn-in for ENP is 71 [Test #1 (N=2, mu=0, N=1) in rep #2].
    % Maximum burn-in for Mean Eccentricity is 99 [Test #2 (N=3, mu=0, N=1) in rep #5].
    % So to be on the save side I will use a burn-in value of 150 for the all-hunter model.

    
    
%% 3. Final model -- OUTDATED. SEE NEXT SECTION WITH ENSEMBLE AVERAGE
% Using time-average so a single repetition. Use 1000 runs of the model.
% Start with 100 post-burnin iterations. If this is not enough to pass our
% five sample size checks, then increase the number of post-burnin
% iterations to 1000, etc.
    
    pref.burnin = 150; % Number of iterations before we have burnin 
    pref.seed = rng(52591025, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 1000; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = pref.burnin + 100; % Number of post-burnin iterations
    pref.repetitions = 1; % Number of repetitions of run
    
    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_mean_representation = NaN(pref.repetitions, pref.iterations, pref.runs);
    export_param = NaN(pref.repetitions, 3, pref.runs); % The three extra coloumns are for: N, mu, n_ratio

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
        pref2.rules = repmat( {'HUNTER'}, 1, pref2.N);
        
        % Repetitions
        data_mean_eccentricity_run   = NaN(pref.repetitions, pref.iterations);
        data_ENP_run                 = NaN(pref.repetitions, pref.iterations);
        data_mean_representation_run = NaN(pref.repetitions, pref.iterations);
        export_param_run             = NaN(pref.repetitions, 3);
        for rep=1:pref.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, o_mean_representation] = ABM(pref3);
            
            % Store summary variables from each run and each repetition
            data_mean_eccentricity_run(rep,:) = o_mean_eccentricity';
            data_ENP_run(rep,:) = o_ENP';
            data_mean_representation_run(rep,:) = o_mean_representation';
            
            export_param_run(rep,:) = [pref3.N pref3.mu pref3.n_ratio];
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
    
    % Time average estimate (of post-burnin iterations)
    est_mean_eccentricity = mean(data_burnin_mean_eccentricity, 2);
    est_ENP = mean(data_burnin_ENP, 2);
    est_mean_representation = mean(data_burnin_mean_representation, 2);
    
    % Time average estimate standard deviation (of post-burnin iterations)
    est_std_mean_eccentricity = std(data_burnin_mean_eccentricity, 0, 2);
    est_std_ENP = std(data_burnin_ENP, 0, 2);
    est_std_mean_representation = std(data_burnin_mean_representation, 0, 2);
    
    % Time average estimate standard error (of post-burnin iterations)
    est_se_mean_eccentricity = est_std_mean_eccentricity ./ sqrt(pref.iterations - pref.burnin);
    est_se_ENP = est_std_ENP ./ sqrt(pref.iterations - pref.burnin);
    est_se_mean_representation = est_std_mean_representation ./ sqrt(pref.iterations - pref.burnin);

    
    
    %%# 3.1 Check 1 -- Mapping / R-Hat statistics
    % TO-DO: ??? There is only 1 repetition here ???
    
    
    %%# 3.2 Check 2 -- Convergence / F-test p-value
    % Does not apply. There is no prior expectation for ENP or mean eccentricty for the all-hunter model. 
    
    
    %%# 3.3 Check 3 -- Power zero
    % Calculating the power of the t-test where H_0: estimate is equal zero, H_A: different from zero
    power_zero_mean_eccentricity = powerzero(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.iterations-pref.burnin);
    power_zero_ENP = powerzero(squeeze(est_ENP), squeeze(est_std_ENP), pref.iterations-pref.burnin);
    power_zero_mean_representation = powerzero(squeeze(est_mean_representation), squeeze(est_std_mean_representation), pref.iterations-pref.burnin);
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
    export_param_fmt = reshape(permute(export_param,[1 3 2]), [pref.runs, 3]);
    
    %%% 3.6 Export results
    % Format table before saving file
    export_mean_eccentricity = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), squeeze(est_se_mean_eccentricity), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_mean_eccentricity, NaN(pref.runs,1), SESD_ratio_mean_eccentricity, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_ENP = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_ENP, NaN(pref.runs,1), SESD_ratio_ENP, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_mean_representation = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_mean_representation), squeeze(est_std_mean_representation), squeeze(est_se_mean_representation), NaN(pref.runs,1), NaN(pref.runs,1), power_zero_mean_representation, NaN(pref.runs,1), SESD_ratio_mean_representation, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    %export_ENP.Properties.Description = ['burnin ' num2str(pref.burnin) ' iterations ' num2str(pref.iterations)];
    % Save file
    writetable(export_mean_eccentricity, strcat('data/MCP_hunter_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/MCP_hunter_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');
    writetable(export_mean_representation, strcat('data/MCP_hunter_mean_representation_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '.csv'),'Delimiter',',');



%% 3. Final model
% Using ensemble-average to calculate estimates. Using 150 burnin iterations.
% Starting with 1000 repetitions. If this is not enough to pass our (five) 
% sample size checks, then increase the number of repetition to 2000, etc.
    
    pref.burnin = 150; % Number of iterations before we have burnin 
    pref.seed = rng(85029804, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = pref.burnin + 1; % Number of post-burnin iterations
    pref.repetitions = 100; % Number of repetitions of run
    
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
        
        % Decision rules: All-hunter
        pref2.rules = repmat( {'HUNTER'}, 1, pref2.N);
        
        % Repetitions
        data_mean_eccentricity_run   = NaN(pref.repetitions, pref.iterations);
        data_ENP_run                 = NaN(pref.repetitions, pref.iterations);
        data_mean_representation_run = NaN(pref.repetitions, pref.iterations);
        export_param_run             = NaN(pref.repetitions, 4);
        for rep=1:pref2.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, o_mean_representation] = ABM(pref3);
            
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
    % Does not apply. There is no prior expectation for ENP or mean eccentricty for the all-hunter model. 
    
    
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
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});    % Save file
    export_mean_representation = table(export_param_fmt(:,1), export_param_fmt(:,2), export_param_fmt(:,3), squeeze(est_mean_representation), squeeze(est_std_mean_representation), squeeze(est_se_mean_representation), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), SESD_ratio_mean_representation, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});    % Save file
    writetable(export_mean_eccentricity, strcat('data/MCP_hunter2_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/MCP_hunter2_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_mean_representation, strcat('data/MCP_hunter2_mean_representation_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');


