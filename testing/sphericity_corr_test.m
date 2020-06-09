%Test Greenhouse-Geisser correction
%
%Author: Eric Fields
%Version Date: 9 June 2020

FMUT_dir = fileparts(fileparts(mfilename('fullpath')));

n_electrodes = 32;
n_time_pts = 40;
n_subs = 24;

%Random data for 3 x 2 ANOVA design
data = randn([n_electrodes, n_time_pts, 3, 2, n_subs]);

%Get epsilon estimates
epsilon_gg = estimate_epsilon(data, [], [3,4], 'gg');
epsilon_hf = estimate_epsilon(data, [], [3,4], 'hf');
epsilon_lb = estimate_epsilon(data, [], [3,4], 'lb');

%Calculate ANOVA
test_results_gg = calc_param_ANOVA(data, [], [3,4], 0.05, 'none', 'gg');
test_results_hf = calc_param_ANOVA(data, [], [3,4], 0.05, 'none', 'hf');
test_results_lb = calc_param_ANOVA(data, [], [3,4], 0.05, 'none', 'lb');

%Check results for a random electrode and time poin
e = randi([1,n_electrodes]);
t = randi([1,n_time_pts]);

fprintf('\nFMUT RESULTS\n');
fprintf('GG epsilon = %.4f\n', epsilon_gg(e,t));
fprintf('HF epsilon = %.4f\n', epsilon_hf(e,t));
fprintf('LB epsilon = %.4f\n', epsilon_lb(e,t));
fprintf('p(GG) = %.3f\n', test_results_gg.p(e,t));
fprintf('p(HF) = %.3f\n', test_results_hf.p(e,t));
fprintf('p(LB) = %.3f\n', test_results_lb.p(e,t));

%Create ANOVA table
oneway_data = reshape(data, n_electrodes, n_time_pts, [], n_subs);
T = [cell2table(cellfun(@(x) sprintf('S%d', x), num2cell(1:n_subs)', 'UniformOutput', false), 'VariableNames', {'sub'}) ...
     array2table(squeeze(oneway_data(e,t,:,:))')];
T.Properties.VariableNames = {'sub', 'A1B1', 'A2B1' ,'A3B1', 'A1B2', 'A2B2' ,'A3B2'};

%Output to csv for checking in jamovi/SPSS
writetable(T, fullfile(FMUT_dir ,'testing', 'outputs', 'gg_test.csv'));

%MATLAB repetaed measure ANOVA
withindesign = cell2table({'A1', 'B1'; 'A2', 'B1'; 'A3', 'B1'; ...
                           'A1', 'B2'; 'A2', 'B2'; 'A3', 'B2'}, ...
                           'VariableNames', {'A', 'B'});
rm = fitrm(T, 'A1B1-A3B2~1', 'WithinDesign', withindesign);
rm_anova_results = ranova(rm, 'WithinModel', 'A+B+A*B');
fprintf('\nMATLAB RM RESULTS\n');
disp(rm_anova_results);

%Check that MATLAB model gives the same result
assert(rm_anova_results{'(Intercept):A:B', 'pValueGG'} - test_results_gg.p(e,t) < 1e-7);
assert(rm_anova_results{'(Intercept):A:B', 'pValueHF'} - test_results_hf.p(e,t) < 1e-7);
assert(rm_anova_results{'(Intercept):A:B', 'pValueLB'} - test_results_lb.p(e,t) < 1e-7);

