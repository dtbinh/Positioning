%% Monte Carlo Parameterization -- All-aggregator
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
% The All-aggregator model can be representent as a stationary Markov chain / 
% time-homogenous Markov chain since the transition probabilities are
% independent of time / iteration.
%
% The All-aggregator model does not contain random components (besides the 
% initial position of firms). 
%
% Additionally the All-aggregator model is in essence an implimatation of the 
% Lloyd's Algorithm. The continoues move to the centroid of the firm's 
% market will eventually lead to a Centroidal Voronoi Tessellation (CVT).
% The CVT is an optimal tessellation / partition. For a given set of 
% generating points / consumers / market there may exsist multiple CVT. 
%
% So the Markov chain is stationary and deterministic. Thus process 
% converges to a single state, but the process is nonergodic (so the 
% limiting distribution depends on the initial conditions).
%
% Thus we use the ensemble-average to calculate representative estimates.
%

clearvars;

%% 2. Determining Burn-In
% The process has burnin when there no longer is any change in variables. 
% Summary variables: mean eccentricity, effective number of firms (ENP)
    
    pref.seed = rng(74514339, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 100; % Number of iterations
    pref.repetitions = 1; % Number of repetitions of run

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    export_temp = NaN(pref.repetitions, 3, pref.runs); % The three extra coloumns are for: N, mu, n_ratio

    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Randomly sample parameters from uniform distribution of the range of the parameters
        % Number of firms
        pref.N = randi([2,12]); % N in [2,12]
        % Mean of subpopulation
        pref.mu = rand * 1.5; % mu in [0,1.5]
        % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
        pref.n_ratio = 1 + rand; % n_ratio in [1,2]
        
        % Decision rules: All-aggregator
        pref.rules = repmat( {'AGGREGATOR'}, 1, pref.N);
        
        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP] = ABM(pref);

            % Store summary variables from each run
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
            
            export_temp(rep,:,run) = [pref.N pref.mu pref.n_ratio];
        end

    end
    close(h);
    
    % Burnin at first occurrence of the final value / No change in variable
    burnin_mean_eccentricity = burnin(data_mean_eccentricity);
    burnin_ENP = burnin(data_ENP);
    
    % The maximum burn-in number of all test repetitions.
    max( burnin_ENP(:) )
    max( burnin_mean_eccentricity(:) )
    
    % Preliminary conclusion
    % Maximum burn-in for ENP is 69 [Run #472 (N=4, mu=1.0184, n_ratio=1.7948)].
    % Maximum burn-in for Mean Eccentricity is 70 [Run #472 (N=4, mu=1.0184, n_ratio=1.7948)].
    % So to be on the save side I will use a burn-in value of 100 for the all-aggregator model.
    
    
    
%% 3. Final model
% Using ensemble-average to calculate estimates. Using 100 burnin iterations.
% Starting with 1000 repetitions. If this is not enough to pass our (five) 
% sample size checks, then increase the number of repetition to 2000, etc.
    
    pref.burnin = 100; % Number of iterations before we have burnin 
    pref.seed = rng(39641227, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = pref.burnin + 1; % Number of post-burnin iterations
    pref.repetitions = 100; % Number of repetitions of run
    
    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = NaN(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = NaN(pref.repetitions, pref.iterations, pref.runs);
    export_temp = NaN(pref.repetitions, 4, pref.runs); % The four extra coloumns are for: repetition number, N, mu, n_ratio

    h = waitbar(0, 'Running...');
    for run=1:pref.runs
        waitbar(run/pref.runs);
        pref.run = run;

        % Randomly sample parameters from uniform distribution of the range of the parameters
        % Number of firms
        pref.N = randi([2,12]); % N in [2,12]
        % Mean of subpopulation
        pref.mu = rand * 1.5; % mu in [0,1.5]
        % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
        pref.n_ratio = 1 + rand; % n_ratio in [1,2]
        
        % Decision rules: All-aggregator
        pref.rules = repmat( {'AGGREGATOR'}, 1, pref.N);
        
        % Repetitions
        for rep=1:pref.repetitions
            pref.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP] = ABM(pref);
            
            % Store summary variables from each run and each repetition
            data_mean_eccentricity(rep,:,run) = o_mean_eccentricity';
            data_ENP(rep,:,run) = o_ENP';
            
            export_temp(rep,:,run) = [pref.N pref.mu pref.n_ratio rep ];
        end

    end
    close(h);
    
    % Post-burnin iterations
    data_burnin_mean_eccentricity = data_mean_eccentricity( :, pref.burnin+1:end, :);
    data_burnin_ENP = data_ENP( :, pref.burnin+1:end, :);
    
    % Ensemble average estimate
    est_mean_eccentricity = mean(data_burnin_mean_eccentricity, 1);
    est_ENP = mean(data_burnin_ENP, 1);
    
    % Ensemble average estimate standard deviation
    est_std_mean_eccentricity = std(data_burnin_mean_eccentricity, 0, 1);
    est_std_ENP = std(data_burnin_ENP, 0, 1);
    
    % Ensemble average estimate standard error
    est_se_mean_eccentricity = est_std_mean_eccentricity ./ sqrt(pref.repetitions);
    est_se_ENP = est_std_ENP ./ sqrt(pref.repetitions);

    
    
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
    % Have all summary variables been estimated with the same level of
    % precisions? This is (trivially) satisfied when the summary variables
    % have been estimated using the same number of post-burnin iterations,
    % since the SE/SD ratio simply returns one over the squareroot of the 
    % number of post-burnin iterations. 

    % Reshape the data to required format before exporting
    export_temp_fmt = reshape(permute(export_temp(1,:,:),[1 3 2]), [pref.runs, 4]);
    
    %%% 3.6 Export results
    % Format table before saving file
    export_mean_eccentricity = table(export_temp_fmt(:,1), export_temp_fmt(:,2), export_temp_fmt(:,3), squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), squeeze(est_se_mean_eccentricity), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), SESD_ratio_mean_eccentricity, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_ENP = table(export_temp_fmt(:,1), export_temp_fmt(:,2), export_temp_fmt(:,3), squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), SESD_ratio_ENP, ...
                      'VariableNames', {'N' 'mu' 'n_ratio' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});    % Save file
    writetable(export_mean_eccentricity, strcat('data/MCP_aggregator_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/MCP_aggregator_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');


