%Test step down procedures with fully simulated normal data
%
%Author: Eric Fields
%Version Date: 6 November 2017

clearvars;

%Basic parameters
n_exp = 1e3;
cond_subs = [];
n_electrodes = 1;
n_time_pts = 512;
dims = 3;
n_perm = 1e3;
alpha = 0.05;

%Step down procedure to use
step_down = false;

%Define effects
effect = 1;
effect_length = 100;

%Pre-allocate variables
test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                             'F_obs', NaN(n_electrodes, n_time_pts),  'Fmax_crit', NaN, ... 
                             'df', NaN(1, 2), 'estimated_alpha', NaN, 'exact_test', NaN), ...
                             n_exp, 1);
familywise_rej     = NaN(1, n_exp);
contrastwise_rej   = NaN(1, n_exp);
familywise_power   = NaN(1, n_exp);
contrastwise_power = NaN(1, n_exp);
FDR                = NaN(1, n_exp);

%Simulated experiments
tic
parfor exp = 1:n_exp
    
    %Simulate data
    data = normrnd(0, 1, [n_electrodes, n_time_pts, 2, 16]);
    if any(effect)
        data(:, 1:effect_length, 1, :) = data(:, 1:effect_length, 1, :) + effect;
    end
    
    %Conduct test
    test_results(exp) = calc_Fmax(data, cond_subs, dims, n_perm, alpha, step_down);
    
    %Summarize experiment results
    familywise_rej(exp) = any(test_results(exp).h(effect_length+1:end));
    contrastwise_rej(exp) = mean(test_results(exp).h(effect_length+1:end));
    if any(effect) && any(test_results(exp).h)
        familywise_power(exp) = any(test_results(exp).h);
        if any(test_results(exp).h)
            contrastwise_power(exp) = mean(test_results(exp).h(1:effect_length));
            FDR(exp) = sum(test_results(exp).h(effect_length+1:end)) / (sum(test_results(exp).h));
        else
            contrastwise_power(exp) = NaN;
            FDR(exp) = NaN;
        end
    end
    
end
toc

%Report results
fprintf('\nFamily-wise error rate = %f\n', mean(familywise_rej));
fprintf('Contrast-wise error rate = %f\n', mean(contrastwise_rej));
if any(effect)
    fprintf('Family-wise power = %f\n', mean(familywise_power));
    fprintf('Contrast-wise power = %f\n', nanmean(contrastwise_power));
    fprintf('FDR = %f\n', nanmean(FDR));
end
fprintf('\n');

beep;