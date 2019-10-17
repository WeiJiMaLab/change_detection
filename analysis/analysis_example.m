% Generate model predictions for each model and set size. This could all be 
% done completely in parallel per set size and model. Note that model 15 
% (VP) takes quite a bit of memory and time to run for each set size.

fprintf('In this example analysis, the IL, SA, SR, EP, and VP\n') 
fprintf('models are fitted to the data of all subjects.\n')
fprintf('The entire analysis can take between 1-7 hours time, \n');
fprintf('depending on the speed of your computer, the number of \n')
fprintf('MC samples, and whether you parallelize the first "for" loop.\n\n');
fprintf('Press any key to start the analysis.\n');
pause()

num_samples = 500;
temp_subj_num = 1;

% set random seed
rng('shuffle')

for model_num = [5 4 17 16 15]
    fprintf('\nOrientation-change\n\n');
    for set_size_idx = 1:4
        % orientation-change
        is_color = 3;
        run_model_a_on_subj_b_with_c_samples_new(model_num, temp_subj_num, num_samples,is_color,set_size_idx);
    end
    fprintf('\nColor-change\n\n');
    for set_size_idx = 1:4
        % color-change
        is_color = 4;
        run_model_a_on_subj_b_with_c_samples_new(model_num, temp_subj_num, num_samples, is_color,set_size_idx)
    end
end

% Compute the log-likelihoods and maximum-likelihood fits to each subject. 
% Note that for model 15 (VP), this takes up about 8GB of RAM per subject. 
% There are 10 subjects for orientation-change data, 7 for color-change.
run_num = [];
for model_num = [5 4 17 16 15]
    fprintf('\nOrientation-change\n\n');
    for subj_num = 1:10
        is_color = 0;
        get_subj_LL_cluster_par(model_num,is_color,run_num,subj_num);
    end
    fprintf('\nColor-change\n\n');
    for subj_num = 1:7
        is_color = 1;
        get_subj_LL_cluster_par(model_num,is_color,run_num,subj_num);
    end
end

% Generate the psychometric curves for orientation and color
load ../data/subjCell
load ../data/subjColorCell
subjCell = revert_data(subjCell);
subjColorCell = revert_data(subjColorCell);

num_delta_bins = 11;
model_nums = [5 4 17 16 15];
is_color = 0;
[p_C_hat_mat HR FA] = compute_psych_curves_multiple(subjCell,num_delta_bins,model_nums,is_color)
is_color = 1;
figure
[p_C_hat_mat HR FA] = compute_psych_curves_multiple(subjColorCell,num_delta_bins,model_nums,is_color)

% Compute and plot BMC results
figure
draw_bars_new()
