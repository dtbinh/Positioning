%% Grid Sweep -- All-sticker
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
% The All-stikcer model is a stationary Markov chain / time-homogenous 
% Markov chain. Furthermore the Markov chain is determenistic, the 
% sticker-firm does not change position. The transition probabilities are
% independent of time / iteration, and the probability is one that the 
% firms are in the same position in the next iteration.
%
% Thus we use the ensemble-average to calculate representative estimates.
%
% No need to determine burn-in since we know that the Markov process will
% not change state. So we just use the first iteration.
%

clearvars;


%% 2. Final model    
% Using ensemble-average to calculate estimates. Using a signle iteration.
% Starting with 1000 repetitions. If this is not enough to pass our (five) 
% sample size checks, then increase the number of repetition to 2000, etc.
    
    % Run through all firm sizes.
    test.N = 2:12;

    pref.seed = rng(76165364, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = length(test.N); % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 1; % Number of post-burnin iterations
    pref.repetitions = 1000; % Number of repetitions of run

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

        % Decision rules: All-sticker
        pref.rules = repmat( {'STICKER'}, 1, pref.N);

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
    
    
    % Ensemble average estimate
    est_mean_eccentricity = mean(data_mean_eccentricity, 1);
    est_ENP = mean(data_ENP, 1);
    
    % Ensemble average estimate standard deviation
    est_std_mean_eccentricity = std(data_mean_eccentricity, 0, 1);
    est_std_ENP = std(data_ENP, 0, 1);
    
    % Ensemble average estimate standard error
    est_se_mean_eccentricity = est_std_mean_eccentricity ./ sqrt(pref.repetitions);
    est_se_ENP = est_std_ENP ./ sqrt(pref.repetitions);
    
    
    
    %%% 3.1 Check 1 -- Mapping / R-Hat statistics
    % Not applicable. Because we don't use time averages (only 1 iteration).
    
    
    %%% 3.2 Check 2 -- Convergence / F-test p-value
    % The firm initial position is drawn uniformly from a cicle with mean
    % (0,0) and a radius of 3 std. dev. With no polarisation of the
    % subpopulations the total population mean is (0,0), thus we would
    % expect that the mean eccentricity from the all-sticker model is 1.5.
    % We insure a suffecient/representative sample size by requiring that 
    % the p-value of a F-test is greater than 0.1. H_0: mean eccentricity 
    % equal 1.5, H_A: different.
    %[h,p,~] = fishertest([ squeeze(est_mean_eccentricity) repmat(1.5, pref.runs, 1) ]);

    % TO-DO: Why not t-test ???
    [~,p] = ttest(data_mean_eccentricity(:,1,:), 1.5);
    pvalue_mean_eccentricity = squeeze(p);
    % No expecations for ENP
    pvalue_ENP = NaN(pref.runs, 1);
    % The p-value should be greater than 0.1.
    
    %%% 3.3 Check 3 -- Power zero
    % Calculating the power of the t-test where H_0: estimate is equal zero, H_A: different from zero
    power_zero_mean_eccentricity = powerzero(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.repetitions);
    power_zero_ENP = powerzero(squeeze(est_ENP), squeeze(est_std_ENP), pref.repetitions);
    % The power of the t-test should be at least 0.8
    
    
    %%% 3.4 Check 4 -- Power difference (grid sweeps)
    % Calculating the power of the two-sample t-test where H_0: estimate is equal estimate from adjacent grid, H_A: different
    power_diff_mean_eccentricity = powerdiff(squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), pref.repetitions);
    power_diff_ENP = powerdiff(squeeze(est_ENP), squeeze(est_std_ENP), pref.repetitions);
    % The power of the t-test should be at least 0.8 for ENP.
    % For mean eccentricity the power should be as low as possible, since
    % we expect that the estimates are 1.5. 
   
    
    %%% 3.5 Check 5 -- SE/SD ratio
    % 
    SESD_ratio_mean_eccentricity = squeeze(est_se_mean_eccentricity) ./ squeeze(est_std_mean_eccentricity);
    SESD_ratio_ENP = squeeze(est_se_ENP) ./ squeeze(est_std_ENP);
    % Have all summary variables been estimated with the same level of
    % precisions? This is (trivially) satisfied when the summary variables
    % have been estimated using the same number of repetitions,
    % since the SE/SD ratio simply returns one over the squareroot of the 
    % number of repetitions. 
   
   
   
    %%% 3.6 Export results
    % Format table before saving file
    export_mean_eccentricity = table(test.N', squeeze(est_mean_eccentricity), squeeze(est_std_mean_eccentricity), squeeze(est_se_mean_eccentricity), NaN(pref.runs,1), pvalue_mean_eccentricity, power_zero_mean_eccentricity, power_diff_mean_eccentricity, SESD_ratio_mean_eccentricity, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    export_ENP = table(test.N', squeeze(est_ENP), squeeze(est_std_ENP), squeeze(est_se_ENP), NaN(pref.runs,1), pvalue_ENP, power_zero_ENP, power_diff_ENP, SESD_ratio_ENP, ...
                      'VariableNames', {'N' 'MeanEst' 'StdDev' 'StdError' 'Check1_Rhat' 'Check2_Ftest' 'Check3_PowerZero' 'Check4_PowerDiff', 'Check5_SESD'});
    %export_ENP.Properties.Description = ['burnin ' num2str(pref.burnin) ' iterations ' num2str(pref.iterations)];
    % Save file
    writetable(export_mean_eccentricity, strcat('data/GS_sticker_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');
    writetable(export_ENP, strcat('data/GS_sticker_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_r', num2str(pref.repetitions), '.csv'),'Delimiter',',');

    % Preliminary conclusion
    %
    % With 100 repetitions the model almost satisfies check 2 for mean 
    % eccentricity (only for N=11 the p-value is below 0.1). Check 3 is
    % statisfied for both summary variables. Check 4 is not statisfied for
    % ENP, and the power is quite large for mean eccentricity. 
    %
    % With 1000 repetitions the model statisfies chech 2 for mean 
    % eccentricity (although the p-value for N=11 only just above 0.1).
    % Check 3 is statisfied. Check 4 is statisfied for ENP. For mean 
    % eccentricity check 4 is almost statisfied (the power is still quite 
    % high for N=11. It is 0.26).
    %
    % The ENP estimates are quite low compared with the results from 
    % Lever and Sergenti (2011).
 
    
    
