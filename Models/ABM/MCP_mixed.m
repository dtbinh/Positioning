%% Monte Carlo Parameterization -- Mixed decision rules + entry/exit
% version 0.01
% Jonas K. Sekamane
%
% Free paramter: a_f, tau 
%
% Endogenous number of firms.
% No polarization of subpopulation.
% Equal subpopulation size.
%
% Inspired in part by: 
%   Lever and Sergenti (2011), Laver, Sergenti and Schilperood (2012)



%% 1. Introduction
% Firms born with decision rule. Equal probabability of decision rule.
%
% The All-hunter model can be representent as a stationary Markov chain / 
% time-homogenous Markov chain since the transition probabilities are
% independent of time / iteration.
%
% Contains random components (decision rule, hunter direction).
%
% So the Markov chain is stationary (???) and stocastic, which implies that 
% the Markov chain is ergodic (the process tends to a unquie limiting 
% distribution that is independent of the initial conditions).
%
% Time-average do not give representative estimates. Using ensemble-average 
% to calculate estimates
%

clearvars;

%% 2. Determining Burn-In
% Graphically.    
    
    %%% 2.1 Idenfying runs that require most iteration
    % Idenfying from 500 runs with 200 iterations and 1 repetition.
    pref.seed = rng('default');
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 200; % Number of runs of the experiment
    pref.iterations = 200; % Number of iteration
    pref.repetitions = 1; % Number of repetitions of run
    pref.ruleset = cellstr([{'STICKER'}, {'HUNTER'}, {'AGGREGATOR'}]); % The set of possible decision rules
    pref.N = 1; % Number of initial firms
    pref.mu = 0; % Mean of subpopulation
    pref.n_ratio = 1; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.export_data = 0;
    pref.export_fig = 0;
    
    % Creating empty matrixes for summary variable 
    s_rule   = [1+length(pref.ruleset) pref.iterations pref.repetitions pref.runs];
    s_solo   = [1 pref.iterations pref.repetitions pref.runs];
    s_param  = [pref.repetitions 7 pref.runs];
    data_mean_eccentricity  = NaN(s_rule);
    data_ENP                = NaN(s_solo);
    data_mean_share         = NaN(s_rule);
    data_N                  = NaN(s_rule);
    data_mean_age_death     = NaN(s_rule);
    export_param            = NaN(s_param); % The seven extra coloumns are for: repetition number, N, mu, n_ratio, a_f, tau, psi

    poolobj = parpool('local',4)
    parfor_progress(pref.runs);
    parfor run=1:pref.runs
        pref2 = pref;
        pref2.run = run;

        % Randomly sample parameters from uniform distribution of the range of the parameters
        % Firm memory parameter (weight on past fitness).
        pref2.a_f = rand * 0.9; %a_f in [0.0,0.9]
        % The de facto survival threshold (market share/percentage).
        pref2.tau = 0.05 + rand * 0.25; % tau in [0.05, 0.30]
        % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
        pref2.psi = 10 + randi(15); % psi in [10, 25]
        
        
        % Repetitions
        data_mean_eccentricity_run  = NaN(s_rule(1:end-1));
        data_ENP_run                = NaN(s_solo(1:end-1));
        data_mean_share_run         = NaN(s_rule(1:end-1));
        data_N_run                  = NaN(s_rule(1:end-1));
        data_mean_age_death_run     = NaN(s_rule(1:end-1));
        export_param_run            = NaN(s_param(1:end-1));
        for rep=1:pref.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, o_mean_share, o_N, o_mean_age_death] = ABM_endo(pref3);

            % Store summary variables from each run
            data_mean_eccentricity_run(:,:,rep) = o_mean_eccentricity';
            data_ENP_run(:,:,rep) = o_ENP';
            data_mean_share_run(:,:,rep) = o_mean_share';
            data_N_run(:,:,rep) = o_N';
            data_mean_age_death_run(rep,:) = o_mean_age_death';
            
            export_param_run = [pref3.N pref3.mu pref3.n_ratio pref3.a_f pref3.tau pref3.psi rep];
        end
        data_mean_eccentricity(:,:,:,run)   = data_mean_eccentricity_run;
        data_ENP(:,:,:,run)                 = data_ENP_run;
        data_mean_share(:,:,:,run)          = data_mean_share_run;
        data_N(:,:,:,run)                   = data_N_run;
        data_mean_age_death(:,:,:,run)      = data_mean_age_death_run;
        export_param(:,:,run)               = export_param_run;

        parfor_progress;
    end
    parfor_progress(0);
    delete(poolobj)

    % Save summary variables
    % Reshape the data to required format before exporting
    export_param_fmt = reshape(permute(export_param,[1 3 2]), [pref.repetitions*pref.runs, 7]);
    data_mean_eccentricity_0 = NaN(s_rule([3 2 4]));
    data_mean_eccentricity_0(:,:,:) = permute(data_mean_eccentricity(1,:,:,:), [3 2 4 1]);
    export_mean_eccentricity    = [export_param_fmt reshape( data_mean_eccentricity_0, [pref.repetitions*pref.runs, pref.iterations] ) ];
    data_ENP_0 = NaN(s_solo([3 2 4]));
    data_ENP_0(:,:,:) = permute(data_ENP(1,:,:,:), [3 2 4 1]);
    export_ENP                  = [export_param_fmt reshape( data_ENP_0, [pref.repetitions*pref.runs, pref.iterations] ) ];
    % Export
    csvwrite(strcat('data/data_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), export_mean_eccentricity);
    csvwrite(strcat('data/data_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), export_ENP);

    % Identify runs that require most iterations visually using filters:
    % <https://github.com/jsekamane/filter-time>
    
    % Preliminary conclusion
    % Parameter settings that require most iterations
    % ENP: Low de facto survival threshold (tau 0.5-0.7).
    % mean eccentricity: Seem to converge very fast, however largest
    % variance for large de facto survival thresholds (tau 0.28-0.30).
    
    
    
    
    %%% 2.2 Burn-in
    
    % Preliminary conclusion
    % It seems as though the model has burnt in after 75 system
    % ticks/iterations. We will set the burnin period to 100 system ticks.

    
    
%% 3. Final model
% Using ensemble-average to calculate estimates. Using 100 burnin iterations.
% Starting with 100 repetitions. If this is not enough to pass our (five) 
% sample size checks, then increase the number of repetition to 1000, etc.
    
    pref.burnin = 100; % Number of iterations before we have burnin 
    pref.seed = rng(73475983, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 250; %Number of runs of the experiment
    pref.iterations = pref.burnin + 1; % Number of iteration
    pref.repetitions = 50; %Number of repetitions of run
    pref.ruleset = cellstr([{'STICKER'}, {'HUNTER'}, {'AGGREGATOR'}]); % The set of possible decision rules
    pref.N = 1; % Number of initial firms
    pref.mu = 0; % Mean of subpopulation
    pref.n_ratio = 1; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.export_data = 0;
    pref.export_fig = 0;
    
    % Creating empty matrixes for summary variable 
    s_rule   = [1+length(pref.ruleset) pref.iterations pref.repetitions pref.runs];
    s_solo   = [1 pref.iterations pref.repetitions pref.runs];
    s_param  = [pref.repetitions 7 pref.runs];
    data_mean_eccentricity  = NaN(s_rule);
    data_ENP                = NaN(s_solo);
    data_mean_share         = NaN(s_rule);
    data_N                  = NaN(s_rule);
    data_mean_age_death     = NaN(s_rule);
    export_param            = NaN(s_param); % The seven extra coloumns are for: repetition number, N, mu, n_ratio, a_f, tau, psi

    poolobj = parpool('local',4)
    parfor_progress(pref.runs);
    parfor run=1:pref.runs
        pref2 = pref;
        pref2.run = run;

        % Randomly sample parameters from uniform distribution of the range of the parameters
        % Firm memory parameter (weight on past fitness).
        pref2.a_f = rand * 0.9; %a_f in [0.0,0.9]
        % The de facto survival threshold (market share/percentage).
        pref2.tau = 0.05 + rand * 0.25; % tau in [0.05, 0.30]
        % Number of ticks per system ticks. The number of system tick must be integer, thus iterations/psi needs to be integer.
        pref2.psi = 10 + randi(15); % psi in [10, 25]
        
        
        % Repetitions
        data_mean_eccentricity_run  = NaN(s_rule(1:end-1));
        data_ENP_run                = NaN(s_solo(1:end-1));
        data_mean_share_run         = NaN(s_rule(1:end-1));
        data_N_run                  = NaN(s_rule(1:end-1));
        data_mean_age_death_run     = NaN(s_rule(1:end-1));
        export_param_run            = NaN(s_param(1:end-1));
        for rep=1:pref.repetitions
            pref3 = pref2;
            pref3.rep = rep;
            
            % Run ABM with parmeters
            [o_mean_eccentricity, o_ENP, o_mean_share, o_N, o_mean_age_death] = ABM_endo(pref3);
            
            % Store summary variables from each run
            data_mean_eccentricity_run(:,:,rep) = o_mean_eccentricity';
            data_ENP_run(:,:,rep) = o_ENP';
            data_mean_share_run(:,:,rep) = o_mean_share';
            data_N_run(:,:,rep) = o_N';
            data_mean_age_death_run(:,:,rep) = o_mean_age_death';
            
            export_param_run(rep,:) = [pref3.N pref3.mu pref3.n_ratio pref3.a_f pref3.tau pref3.psi rep];
        end
        data_mean_eccentricity(:,:,:,run)   = data_mean_eccentricity_run;
        data_ENP(:,:,:,run)                 = data_ENP_run;
        data_mean_share(:,:,:,run)          = data_mean_share_run;
        data_N(:,:,:,run)                   = data_N_run;
        data_mean_age_death(:,:,:,run)      = data_mean_age_death_run;
        export_param(:,:,run)               = export_param_run;

        parfor_progress;
    end
    parfor_progress(0);
    delete(poolobj)
    
    % Post-burnin iterations
    data_burnin_mean_eccentricity = data_mean_eccentricity( :, pref.burnin+1:end, :, :);
    data_burnin_ENP = data_ENP( :, pref.burnin+1:end, :, :);
    data_burnin_mean_share = data_mean_share( :, pref.burnin+1:end, :, :);
    data_burnin_N = data_N( :, pref.burnin+1:end, :, :);
    data_burnin_mean_age_death = data_mean_age_death( :, pref.burnin+1:end, :, :);
    
    % Ensemble average estimate
    est_mean_eccentricity = mean(data_burnin_mean_eccentricity, 3, 'omitnan');
    est_ENP = mean(data_burnin_ENP, 3);
    est_mean_share = mean(data_burnin_mean_share, 3, 'omitnan');
    est_N = mean(data_burnin_N, 3, 'omitnan');
    est_mean_age_death = mean(data_burnin_mean_age_death, 3, 'omitnan');
    
    % Ensemble average estimate standard deviation
    est_std_mean_eccentricity = std(data_burnin_mean_eccentricity, 0, 3, 'omitnan');
    est_std_ENP = std(data_burnin_ENP, 0, 3);
    est_std_mean_share = std(data_burnin_mean_share, 0, 3, 'omitnan');
    est_std_N = std(data_burnin_N, 0, 3, 'omitnan');
    est_std_mean_age_death = std(data_burnin_mean_age_death, 0, 3, 'omitnan');
    
    % Ensemble average estimate standard error
    est_se_mean_eccentricity = est_std_mean_eccentricity ./ sqrt(pref.repetitions);
    est_se_ENP = est_std_ENP ./ sqrt(pref.repetitions);
    est_se_mean_share = est_std_mean_share ./ sqrt(pref.repetitions);
    est_se_N = est_std_N ./ sqrt(pref.repetitions);
    est_se_mean_age_death = est_std_mean_age_death ./ sqrt(pref.repetitions);

    
    
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
    %SESD_ratio_mean_eccentricity = squeeze(est_se_mean_eccentricity) ./ squeeze(est_std_mean_eccentricity);
    %SESD_ratio_ENP = squeeze(est_se_ENP) ./ squeeze(est_std_ENP);
    % Have all summary variables been estimated with the same level of
    % precisions? This is (trivially) satisfied when the summary variables
    % have been estimated using the same number of post-burnin iterations,
    % since the SE/SD ratio simply returns one over the squareroot of the 
    % number of post-burnin iterations. 

    % Reshape the data to required format before exporting
    export_param_fmt_solo = reshape(permute(export_param(1,:,:),[1 3 2]), [pref.runs, 7]);
    export_param_fmt_rule = [ repelem(export_param_fmt_solo,1+length(pref.ruleset),1) repmat(0:length(pref.ruleset),1,size(export_param_fmt_solo,1))']; % Expand 'solo' so one line for each rule (rule=1,2,3...) plus line for total (rule=0).
    % Delete rep number (since we look at ensemble average).
    export_param_fmt_solo(:,7) = []; 
    export_param_fmt_rule(:,7) = [];
    % Convert to table
    export_param_fmt_solo_t = array2table(export_param_fmt_solo, ...
                                          'VariableNames', {'N' 'mu' 'n_ratio' 'a_f' 'tau' 'psi'});
    export_param_fmt_rule_t = array2table(export_param_fmt_rule, ...
                                          'VariableNames', {'N' 'mu' 'n_ratio' 'a_f' 'tau' 'psi' 'rule'});                                  
    
    %horzcat
    
    %%% 3.6 Export results
    % Format table before saving file
    l_rule = prod(s_rule([1 4]));
    
    export_mean_eccentricity_t = table(est_mean_eccentricity(:), est_std_mean_eccentricity(:), est_se_mean_eccentricity(:), NaN(l_rule,1), NaN(l_rule,1), NaN(l_rule,1), NaN(l_rule,1), NaN(l_rule,1), ...
                      'VariableNames', {'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_mean_eccentricity = horzcat(export_param_fmt_rule_t, export_mean_eccentricity_t);
    
    export_ENP_t = table(squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), NaN(pref.runs,1), ...
                      'VariableNames', {'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});    % Save file
    export_ENP = horzcat(export_param_fmt_solo_t, export_ENP_t);
    
    export_mean_share_t = table(est_mean_share(:), est_std_mean_share(:), est_se_mean_share(:), ...
                      'VariableNames', {'MeanEst' 'StdDev' 'StdError'});
    export_mean_share = horzcat(export_param_fmt_rule_t, export_mean_share_t);
    
    export_N_t = table(est_N(:), est_std_N(:), est_se_N(:), ...
                      'VariableNames', {'MeanEst' 'StdDev' 'StdError'});
    export_N = horzcat(export_param_fmt_rule_t, export_N_t);
    
    export_mean_age_death_t = table(est_mean_age_death(:), est_std_mean_age_death(:), est_se_mean_age_death(:), ...
                      'VariableNames', {'MeanEst' 'StdDev' 'StdError'});
    export_mean_age_death = horzcat(export_param_fmt_rule_t, export_mean_age_death_t);
    
    clear export_mean_eccentricity_t export_ENP_t export_mean_share_t export_N_t export_mean_age_death_t;
    
    time = char(pref.timestamp, 'yyyyMMdd_HHmmss');
    writetable(export_mean_eccentricity, strcat('data/MCP_mixed_', time, '_mean_eccentricity', '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/MCP_mixed_', time, '_ENP', '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_mean_share, strcat('data/MCP_mixed_', time, '_mean_share', '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_N, strcat('data/MCP_mixed_', time, '_N', '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_mean_age_death, strcat('data/MCP_mixed_', time, '_mean_age_death', '_i', num2str(pref.iterations), '_b', num2str(pref.burnin), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
