
%% ABM STATIC

% ALL-AGGREGATOR ? Determenistic THMC
% N = 12;
% mu = 1.5;
% n_ratio = 2;
% iterations = 100;


export1 = table((1:pref.iterations)', ENP, mean_eccentricity, mean_representation, ...
                      'VariableNames', {'iteration' 'ENP' 'mean_eccentricity' 'mean_representation'});
writetable(export1, strcat('data/All-aggregator_N12_mu15_nratio2_i100.csv'),'Delimiter',',');


% ALL-HUNTER - Stocastic THMC, time average representative
% N = 12;
% mu = 0;
% n_ratio = 1;
% iterations = 500;

export2 = table((1:pref.iterations)', ENP, mean_eccentricity, mean_representation, ...
                      'VariableNames', {'iteration' 'ENP' 'mean_eccentricity' 'mean_representation'});
writetable(export2, strcat('data/All-hunter_N12_mu0_nratio1_i500.csv'),'Delimiter',',');


% ALL-HUNTER - Stocastic THMC, time average not representative

pref.seed = rng('default'); % Seed such that the randomly generated results are repeatable
pref.timestamp = datetime('now');  % Create time stamp so exported files don't overwrite exsisting files.
pref.runs = 50; % Number of runs of the experiment
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


% ALL-HUNTER - Movement diagram
% N = 5;
% mu = 1.5;
% n_ratio = 2;
% iterations = 100;
xy_extend = reshape(permute(xy, [1 3 2]), [], 2);

export_xy = table(repelem(1:pref.iterations,pref.N)', repmat(1:pref.N,1,pref.iterations)', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_all-hunter_N5_mu15_nratio2_i100.csv'),'Delimiter',',');



% ALL-HUNTER - Eccentricity vs. Shares (long-run location behaviour vs oppertunistic short-run location behaviour)
% N = 5;
% mu = 0;
% n_ratio = 1;
% iterations = 1150; %150 burn-in iterations + 1000 post burn-in interations

burnin = 150;

% Test figure
figure(1);
clf reset; % Reset figure.
hold on;
for n=1:pref.N
    fig1_eccentricity = eccentricity(burnin+1:end,n);
    [fig1_eccentricity_sort, fig1sort] = sort(fig1_eccentricity);
    fig1_share = shares(burnin+1:end,n);
    %scatter(fig1_eccentricity, fig1_share, 'filled', 'MarkerFaceColor', color(n,:));
    %fitx = linspace(0,1.4,100);
    fity = smooth(fig1_eccentricity, fig1_share, 0.4,'rlowess');
    line(fig1_eccentricity_sort, fity(fig1sort), 'Color', color(n,:))
end
hold off;
xlim([0 1.4]); ylim([0 0.35]);

% Combine data
export = [repmat(1:pref.iterations,1,pref.N)' repelem(1:pref.N, 1, pref.iterations)' reshape(eccentricity, 1, [])' reshape(shares, 1, [])'];
% Remove burn-in iterations
postburnin = (export(:,1) > burnin);
export_burnin = array2table(export(postburnin,:), 'VariableNames', {'iteration','firm','eccentricity','share'});
% Export data
pref.timestamp = datetime('now');
writetable(export_burnin, strcat('data/All-hunter_', char(pref.timestamp, 'yyyyMMdd_HHmmss'), '_N5_mu0_n1_b150_i1150', '.csv'),'Delimiter',',');


%% ABM IND STATIC


% Four firms -- two configurations

% MAXCOV
% N = 4;
% mu = 1.5;
% n_ratio = 1.5;
% iterations = 1;
% psi = 200;
% rules = repmat({'MAXCOV'},1,pref.N);

burnin = 150;
xy_burnin = xy(:,:,(burnin+1):end);
xy_extend = reshape(permute(xy_burnin, [1 3 2]), [], 2);
export_xy = table(repelem(pref.iterations*(burnin+1):pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi-pref.iterations*(burnin))', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_maxcov_N4_mu15_nratio15_pbi49.csv'),'Delimiter',',');


% MAXCOV-INDUCTOR
% N = 4;
% mu = 1.5;
% n_ratio = 1.5;
% iterations = 1;
% psi = 1050;
% M = 100;
% a_a = 1-1/75;
% rules = repmat({'MAXCOV-INDUCTOR'},1,pref.N);

burnin = 1000;
xy_burnin = xy(:,:,(burnin+1):end);
xy_extend = reshape(permute(xy_burnin, [1 3 2]), [], 2);
export_xy = table(repelem(pref.iterations*(burnin+1):pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi-pref.iterations*(burnin))', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_maxcov-inductor_N4_mu15_nratio15_pbi49.csv'),'Delimiter',',');


% Twelve firms -- clustering

% MAXCOV - UNIMODAL
% N = 12;
% mu = 0;
% n_ratio = 1;
% iterations = 1;
% psi = 200;
% rules = repmat({'MAXCOV'},1,pref.N);

burnin = 150;
xy_burnin = xy(:,:,(burnin+1):end);
xy_extend = reshape(permute(xy_burnin, [1 3 2]), [], 2);
export_xy = table(repelem(pref.iterations*(burnin+1):pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi-pref.iterations*(burnin))', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_clustering_maxcov_N12_mu0_nratio1_pbi50.csv'),'Delimiter',',');

% MAXCOV - BIMODAL
% N = 12;
% mu = 1.5;
% n_ratio = 2;
% iterations = 1;
% psi = 200;
% rules = repmat({'MAXCOV'},1,pref.N);

burnin = 150;
xy_burnin = xy(:,:,(burnin+1):end);
xy_extend = reshape(permute(xy_burnin, [1 3 2]), [], 2);
export_xy = table(repelem(pref.iterations*(burnin+1):pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi-pref.iterations*(burnin))', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_clustering_maxcov_N12_mu15_nratio2_pbi50.csv'),'Delimiter',',');

% MAXCOV-INDUCTOR - UNIMODAL
% N = 12;
% mu = 0;
% n_ratio = 1;
% iterations = 1;
% psi = 1050;
% M = 100;
% a_a = 1-1/75;
% rules = repmat({'MAXCOV-INDUCTOR'},1,pref.N);

burnin = 1000;
xy_burnin = xy(:,:,(burnin+1):end);
xy_extend = reshape(permute(xy_burnin, [1 3 2]), [], 2);
export_xy = table(repelem(pref.iterations*(burnin+1):pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi-pref.iterations*(burnin))', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_clustering_maxcov-inductor_N12_mu0_nratio1_pbi50.csv'),'Delimiter',',');


% MAXCOV-INDUCTOR - BIMODAL
% N = 12;
% mu = 1.5;
% n_ratio = 2;
% iterations = 1;
% psi = 1050;
% M = 100;
% a_a = 1-1/75;
% rules = repmat({'MAXCOV-INDUCTOR'},1,pref.N);

burnin = 1000;
xy_burnin = xy(:,:,(burnin+1):end);
xy_extend = reshape(permute(xy_burnin, [1 3 2]), [], 2);
export_xy = table(repelem(pref.iterations*(burnin+1):pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi-pref.iterations*(burnin))', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_clustering_maxcov-inductor_N12_mu15_nratio2_pbi50.csv'),'Delimiter',',');



% Three firms -- long-run market share

% MAXCOV
% N = 3;
% mu = 1.5;
% n_ratio = 2;
% iterations = 1;
% psi = 250;
% rules = repmat({'MAXCOV'},1,pref.N);
ex_shares = reshape(shares, 1, [])';
matrix_shares = [repmat((1:pref.iterations*pref.psi)', pref.N, 1) repelem(1:pref.N, pref.iterations*pref.psi)' ex_shares];
csvwrite(strcat('data/three_maxcov_shares_N3_mu15_nratio2_i250.csv'), matrix_shares);


% MAXCOV-INDUCTOR
% N = 3;
% mu = 1.5;
% n_ratio = 2;
% iterations = 1;
% psi = 2500;
% rules = repmat({'MAXCOV-INDUCTOR'},1,pref.N);

ex_shares = reshape(shares, 1, [])';
matrix_shares = [repmat((1:pref.iterations*pref.psi)', pref.N, 1) repelem(1:pref.N, pref.iterations*pref.psi)' ex_shares];
csvwrite(strcat('data/three_maxcov-inductor_shares_N3_mu15_nratio2_i2500.csv'), matrix_shares);

export_cf_change = diff(cf(:,25,:,:), 1, 4);
export_cf_used = (export_cf_change > 0);

cf_change_1 = squeeze(export_cf_used(:,:,1,:));
cf_change_l_1 = reshape(cf_change_1, (pref.iterations*pref.psi-1)*pref.M, 1);
ex_cf_1 = [repelem((1:(pref.iterations*pref.psi-1))', pref.M) repmat((1:pref.M)', pref.iterations*pref.psi-1, 1) cf_change_l_1];
zero_idx = (ex_cf_1(:,3) == 0);
ex_cf_1(zero_idx,:) = [];
csvwrite(strcat('data/three_maxcov-inductor_cf_N3_mu15_nratio2_i2500_1.csv'), ex_cf_1); 

cf_change_2 = squeeze(export_cf_used(:,:,2,:));
cf_change_l_2 = reshape(cf_change_2, (pref.iterations*pref.psi-1)*pref.M, 1);
ex_cf_2 = [repelem((1:(pref.iterations*pref.psi-1))', pref.M) repmat((1:pref.M)', pref.iterations*pref.psi-1, 1) cf_change_l_2];
zero_idx = (ex_cf_2(:,3) == 0);
ex_cf_2(zero_idx,:) = [];
csvwrite(strcat('data/three_maxcov-inductor_cf_N3_mu15_nratio2_i2500_2.csv'), ex_cf_2); 

cf_change_3 = squeeze(export_cf_used(:,:,3,:));
cf_change_l_3 = reshape(cf_change_3, (pref.iterations*pref.psi-1)*pref.M, 1);
ex_cf_3 = [repelem((1:(pref.iterations*pref.psi-1))', pref.M) repmat((1:pref.M)', pref.iterations*pref.psi-1, 1) cf_change_l_3];
zero_idx = (ex_cf_3(:,3) == 0);
ex_cf_3(zero_idx,:) = [];
csvwrite(strcat('data/three_maxcov-inductor_cf_N3_mu15_nratio2_i2500_3.csv'), ex_cf_3); 


% MAXCOV-INDUCTOR-GA
% N = 3;
% mu = 1.5;
% n_ratio = 2;
% iterations = 50;
% psi = 50;
% M = 100;
% a_a = 1-1/75;
% C = 0.0005;
% crossover = 0.3;

ex_shares = reshape(shares, 1, [])';
matrix_shares = [repmat((1:pref.iterations*pref.psi)', pref.N, 1) repelem(1:pref.N, pref.iterations*pref.psi)' ex_shares];
csvwrite(strcat('data/three_maxcov-inductor-ga_shares_N3_mu15_nratio2_i2500.csv'), matrix_shares);

% MAT
save('data/three_maxcov-inductor-ga_xy_N3_mu15_nratio2_i2500', 'xy');
% CSV
ex_xy = reshape(permute(xy, [1 3 2]), pref.iterations*pref.psi*pref.N, 2);
export_xy = [repelem((1:pref.iterations*pref.psi)', pref.N) repmat((1:pref.N)', pref.iterations*pref.psi, 1) ex_xy];
csvwrite(strcat('data/three_maxcov-inductor-ga_xy_N3_mu15_nratio2_i2500.csv'), export_xy);



% Five firms -- accuracy

% MAXCOV-INDUCTOR-GA
% N = 5;
% mu = 1.3;
% n_ratio = 1.5;
% iterations = 50;
% psi = 50;
% M = 100;
% a_a = 1-1/75;
% C = 0.0005;
% crossover = 0.3;
% rules = repmat({'MAXCOV-INDUCTOR-GA'},1,pref.N);

xy_extend = reshape(permute(xy, [1 3 2]), [], 2);
export_xy = table(repelem(1:pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi)', xy_extend(:,1), xy_extend(:,2), ...
                      'VariableNames', {'iteration' 'firm' 'x' 'y'});
writetable(export_xy, strcat('data/xy_accuracy_maxcov-inductor-ga_N5_mu13_nratio15_psi50_i50.csv'),'Delimiter',',');

accuracy_i = squeeze(cf(:, 24, :, :));
accuracy_i_positive = accuracy_i;
accuracy_i_positive(accuracy_i_positive==0) = NaN; % Non-zero accuracy.
accuracy_i_positive_mean_firm = squeeze(mean(accuracy_i_positive, 1, 'omitnan'));  % Mean over all firm's condition/forcast rules.
accuracy_list = reshape(accuracy_i_positive_mean_firm, pref.iterations*pref.psi*pref.N, 1);
export_accuracy = table(repelem(1:pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi)', accuracy_list, ...
                      'VariableNames', {'iteration' 'firm' 'accuracy'});
writetable(export_accuracy, strcat('data/accuracy_maxcov-inductor-ga_N5_mu13_nratio15_psi50_i50.csv'),'Delimiter',',');

% MAXCOV-INDUCTOR-GA
% Using the same initial settings (xy_0 and cf_0). But re-runing model with
% rules = repmat({'MAXCOV-INDUCTOR'},1,pref.N);

accuracy_i = squeeze(cf(:, 24, :, :));
accuracy_i_positive = accuracy_i;
accuracy_i_positive(accuracy_i_positive==0) = NaN; % Non-zero accuracy.
accuracy_i_positive_mean_firm = squeeze(mean(accuracy_i_positive, 1, 'omitnan'));  % Mean over all firm's condition/forcast rules.
accuracy_list = reshape(accuracy_i_positive_mean_firm, pref.iterations*pref.psi*pref.N, 1);
export_accuracy = table(repelem(1:pref.iterations*pref.psi,pref.N)', repmat(1:pref.N,1,pref.iterations*pref.psi)', accuracy_list, ...
                      'VariableNames', {'iteration' 'firm' 'accuracy'});
writetable(export_accuracy, strcat('data/accuracy_maxcov-inductor_N5_mu13_nratio15_psi50_i50_same.csv'),'Delimiter',',');


