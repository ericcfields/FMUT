%Run simulated normal data to check Type I error rate and power
%
%AUTHOR: Eric Fields
%VERSION DATE: 13 July 2017

clearvars;

%% Set-up

global VERBLEVEL
VERBLEVEL = 0;

%Design
n_electrodes = 1;
n_time_pts = 1;
n_subs = 16;
wg_design = [3 3 3];

%Effect
dims = [3 4 5];

%Parameters
n_exp = 1e3;
n_perm = 1e3;
alpha = 0.05;

%Add effects
main_effect = 0;
int_effect = 0;

%Pre-allocate results struct
test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), ...
                             'p', NaN(n_electrodes, n_time_pts), ...
                             'F_obs', NaN(n_electrodes, n_time_pts), ...
                             'Fmax_crit', NaN, ...
                             'df', [NaN, NaN], ...
                             'estimated_alpha', NaN, ...
                             'exact_test', NaN), ...
                       n_exp, 1);

                   
%% Simulate experiments

tic
parfor i = 1:n_exp
    
    %Simulate null data
    data = normrnd(0, 1, [n_electrodes, n_time_pts, wg_design, n_subs]);
    
    %Add effects
    if main_effect
        data(:, :, 1, :) = data(:, :, 1, :) + main_effect;
    end
    if int_effect
        if ndims(data) == 5
            data(:, :, 1, 1, :) = data(:, :, 1, 1, :) + int_effect;
            data(:, :, 1, 2, :) = data(:, :, 1, 2, :) - int_effect;
            data(:, :, 2, 1, :) = data(:, :, 2, 1, :) - int_effect;
            data(:, :, 2, 2, :) = data(:, :, 2, 2, :) + int_effect;
        elseif ndims(data) == 6
            data(:, :, 1, 1, 1, :) = data(:, :, 1, 1, 1, :) + int_effect;
            data(:, :, 1, 2, 1, :) = data(:, :, 1, 2, 1, :) - int_effect;
            data(:, :, 2, 1, 1, :) = data(:, :, 2, 1, 1, :) - int_effect;
            data(:, :, 2, 2, 1, :) = data(:, :, 2, 2, 1, :) + int_effect;
        end
    end
    
    %Calculate ANOVA
    test_results(i) = calc_Fmax(data, [], dims, n_perm, alpha);
    
end
toc


%% Report results

fprintf('\nRejection rate = %.3f\n', mean(mean([test_results(:).h])))
fprintf('Mean p = %.3f\n', mean(mean([test_results(:).p])))
fprintf('Mean F_obs = %.2f\n', mean(mean([test_results(:).F_obs])))
fprintf('Mean F_crit = %.2f\n', mean(mean([test_results(:).Fmax_crit])))
fprintf('Parametric F_crit = %.2f\n\n', finv(.95, test_results(1).df(1), test_results(1).df(2)))
