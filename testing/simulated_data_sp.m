%Run simulated normal data to check Type I error rate and power
%
%AUTHOR: Eric Fields
%VERSION DATE: 24 July 2017

clearvars;

%% Set-up

global VERBLEVEL
VERBLEVEL = 0;

%Design
n_electrodes = 1;
n_time_pts = 1;
wg_design = [3 3];
cond_subs = [8 8];
n_subs = sum(cond_subs);

%Effect
dims = [3, 4, 5];

%Parameters
n_exp = 1e3;
n_perm = 1e3;
alpha = 0.05;

%Add effects
wg_effect = 0;
bg_effect = 0;
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
for i = 1:n_exp
    
    %Simulate null data
    data = normrnd(0, 1, [n_electrodes, n_time_pts, wg_design, n_subs]);
    
    %Add effects
    if wg_effect
        data(:, :, 1, :) = data(:, :, 1, :) + wg_effect;
    end
    if bg_effect
        data(1:numel(data)/size(data,ndims(data))*cond_subs(1)) = data(1:numel(data)/size(data,ndims(data))*cond_subs(1)) + bg_effect;
    end
    if int_effect
        Asubs = 1:cond_subs(1);
        Bsubs = cond_subs(1)+1:cond_subs(1)+cond_subs(2);
        if ndims(data) == 4 && isequal(dims, [3, 4])
            data(:, :, 1, Asubs)  =  data(:, :, 1, Asubs) + int_effect;
            data(:, :, 1, Bsubs) =  data(:, :, 1, Bsubs) - int_effect;
            data(:, :, 2, Asubs)  =  data(:, :, 2, Asubs) - int_effect;
            data(:, :, 2, Bsubs) =  data(:, :, 2, Bsubs) + int_effect;
        elseif ndims(data) == 5
            data(:, :, 1, 1, Asubs) = data(:, :, 1, 1, Asubs) + int_effect;
            data(:, :, 1, 2, Asubs) = data(:, :, 1, 2, Asubs) - int_effect;
            data(:, :, 2, 1, Asubs) = data(:, :, 2, 1, Asubs) - int_effect;
            data(:, :, 2, 2, Asubs) = data(:, :, 2, 2, Asubs) + int_effect;
%             data(:, :, 1, 1, Bsubs) = data(:, :, 1, 1, Bsubs) - int_effect;
%             data(:, :, 1, 2, Bsubs) = data(:, :, 1, 2, Bsubs) + int_effect;
%             data(:, :, 2, 1, Bsubs) = data(:, :, 1, 1, Bsubs) + int_effect;
%             data(:, :, 2, 2, Bsubs) = data(:, :, 2, 2, Bsubs) - int_effect;
        end
    end
    
    %Calculate ANOVA
    %test_results(i) = calc_Fmax(data, cond_subs, dims, n_perm, alpha);
    test_results(i) = calc_param_ANOVA(data, cond_subs, dims, alpha, 'none', 'none');
    
end
toc


%% Report results

fprintf('\nRejection rate = %.3f\n', mean(mean([test_results(:).h])))
fprintf('Mean p = %.3f\n', mean(mean([test_results(:).p])))
fprintf('Mean F_obs = %.3f\n', mean(mean([test_results(:).F_obs])))
fprintf('Mean F_crit = %.3f\n', mean(mean([test_results(:).Fmax_crit])))
fprintf('Parametric F_crit = %.3f\n\n', finv(.95, test_results(1).df(1), test_results(1).df(2)))
