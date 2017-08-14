clearvars;

n_exp = 500;
cond_subs = [];
n_electrodes = 1;
n_time_pts = 1000;
dims = 3;
n_perm = 1e3;
alpha = 0.05;
step_down = 3;

effect = 2;
effect_length = 100;

test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                             'F_obs', NaN(n_electrodes, n_time_pts),  'Fmax_crit', NaN, ... 
                             'df', NaN(1, 2), 'estimated_alpha', NaN, 'exact_test', NaN), ...
                             n_exp, 1);
familywise_rej = NaN(1, n_exp);
contrastwise_rej = NaN(1, n_exp);
tic
parfor exp = 1:n_exp

    data = normrnd(0, 1, [n_electrodes, n_time_pts, 2, 16]);
    
    if effect
        data(:, 1:effect_length, 1, :) = data(:, 1:effect_length, 1, :) + effect;
    end
    
    test_results(exp) = calc_Fmax(data, cond_subs, dims, n_perm, alpha, step_down);
    
    familywise_rej(exp) = any(test_results(exp).h(effect_length+1:end));
    contrastwise_rej(exp) = mean(test_results(exp).h(effect_length+1:end));
    if effect
        power(exp) = mean(test_results(exp).h(1:effect_length));
        FDR(exp) = sum(test_results(exp).h(effect_length+1:end)) / (sum(test_results(exp).h));
    end
    
end
toc

fprintf('\nFamily-wise error rate = %f\n', mean(familywise_rej));
fprintf('Contrast-wise error rate = %f\n', mean(contrastwise_rej));
if effect
    fprintf('Power = %f\n', mean(power));
    fprintf('FDR = %f\n', mean(FDR));
end
fprintf('\n');

beep;