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
% # Run several test repetitions of these; Use "second halves" to calculate R-hat (should be below 1.05) and estimate.
% # From estimate determine maximum burn-in of test repetition: within 1. std. dev. of estimate. 


%%% Monte Carlo parameterizations
% Rando sample parameters from uniform distribution of the range of the parameters

clearvars;

rng('default');
pref.timestamp = datetime('now'); % Number of runs of the experiment
pref.runs = 2; % Number of runs of the experiment
for run=1:pref.runs
    pref.run = run;
    
    % Number of firms
    pref.N = randi([2,12]); % N ? [2,12]
    % Mean of subpopulation
    pref.mu = rand * 1.5; % mu ? [0,1.5]
    % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    pref.n_ratio = 1 + rand; % n_ratio ? [1,2]
    
    %pref.seed = rng(21037635, 'twister'); % Seed such that the randomly generated results are repeatable
    pref.rules = repmat( {'HUNTER'}, 1, pref.N); % Decision rules: All-hunter
    pref.iterations = 24; % Number of iterations
    pref.repetitions = 1; % Number of repetitions of run
    
    ABM(pref);
    
end

% repetition use same input parameter but different random seed.
% iterations is length of a repetition.
