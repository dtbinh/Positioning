%% Monte Carlo Parameterization
% version 0.01
% Jonas K. Sekamane
%
% Exogenous number of firms (and thus decision rules).
%
% Inspired in part by: 
%   Lever and Sergenti (2011)


%% All-Hunter
%
% Time-homogenous Markov Chain:
% # The run enviroment: computational rules of the model, the run input parameters.
% # The run Markoc process: state space, transition matrix, initial state space vector
% # The run output variables
%
% Vector of variables: ???
% * coordinates of each firm
%
% ----------
%
% * Time-homogenous Markov Chain / Markov Chain with stationary transition probabilities.
% * Random component: the oppersite direction (range of 180 degrees) + starting position of firms ???
% * Small state space, or High off-diagonal transition probabilities.
%
% --> Ergodic process (Eisenhardt 1989; Ljungqvist and Sargent 2004)
% --> Stochastic processes for which a time average provides a representative estimate of ?
%
% * Repetition: 1.
% * Iteration: more than 1 (determined by diagnostic checks).
% * Burn-in: (determined by four steps below).
%
% # Specify Markov Chain representation; Specify vector of state space and summary variables.
% # Identify those runs that require most iterations to burn in.
% # Run several test repetitions of these; Use "second halves" to calculate R-hat (should be below 1.05) and calculate the estimate.
% # From estimate determine maximum burn-in of test repetition: within 1. std. dev. of estimate. 

clearvars;

%% 1. Determining burn in

    %%% 1.1 Markov Chain representation
    % Summary variables: mean eccentricity and ENP

    
    
    %%% 1.2 Idenfying runs that require most iteration
    % Idenfying from 500 runs with 200 iterations and 1 repetition.
    pref.seed = rng('default'); % Seed such that the randomly generated results are repeatable
    pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 500; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 200; % Number of iterations
    pref.repetitions = 1; % Number of repetitions of run

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = zeros(pref.runs, pref.iterations+3); % The three extra coloumns are for: N, mu, n_ratio
    data_ENP = zeros(pref.runs, pref.iterations+3);

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
        
        % Decision rules: All-hunter
        pref.rules = repmat( {'HUNTER'}, 1, pref.N);
        
        % Run ABM with parmeters
        [o_mean_eccentricity, o_ENP] = ABM(pref);

        % Store summary variables from each run
        data_mean_eccentricity(run,:) = [pref.N pref.mu pref.n_ratio o_mean_eccentricity'];
        data_ENP(run,:) = [pref.N pref.mu pref.n_ratio o_ENP'];

    end
    close(h);

    % Save summary variables
    csvwrite(strcat('data/data_mean_eccentricity_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), data_mean_eccentricity);
    csvwrite(strcat('data/data_ENP_', char(pref.timestamp, 'yyyyMMdd_HHmmss'),'.csv'), data_ENP);

    % Identify runs that require most iterations visually using filters:
    % <https://github.com/jsekamane/filter-time>
    
    % Preliminary conclusion
    % Parameter settings that require most iterations
    % mean eccentricity: Few firms (N: 2-3), and large difference in ideal  point of subpopulation (mu close to 1.5). No visible effect from relative size of population (n_ratio).
    % ENP: Many firms (N 10-12), and large difference in ideal point (mu close to 1.5). No visible effect from relative size of population (n_ratio).
    
 
    
    %%% 1.3 Test repetitions
    % Making 10 runs with follow parameters setings
    test.N       = [2 2 3 3 10 10 11 11 12 12];
    test.mu      = [1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5];
    test.n_ratio = [1 2 1 2 1 2 1 2 1 2];
    
    pref.seed = rng(85387475, 'twister'); % Seed such that the randomly generated results are repeatable (from random.org)
    pref.timestamp = datetime('now'); % Create time stamp so exported files don't overwrite exsisting files.
    pref.runs = 10; % Number of runs of the experiment
    pref.export_data = 0; % Exports the data
    pref.export_fig = 0; % Exports figures
    pref.iterations = 2000; % Number of iterations
    pref.repetitions = 50; % Number of repetitions of run

    % Creating empty matrixes for summary variable 
    data_mean_eccentricity = zeros(pref.repetitions, pref.iterations, pref.runs);
    data_ENP = zeros(pref.repetitions, pref.iterations, pref.runs);

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
            [o_mean_eccentricity, o_ENP] = ABM(pref);
            
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
    
    % R-hat statistics / potential scale reduction factor (PSRF) -- Brooks and Gelman (1998)    
    rhat_mean_eccentricity = rhat( data_2nd_mean_eccentricity );
    rhat_ENP =rhat( data_2nd_ENP );
    % Should be below 1.05

% repetition use same input parameter but different random seed.
% iterations is length of a repetition.
