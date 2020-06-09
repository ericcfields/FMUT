%Test Greenhouse-Geisser correction
%
%Author: Eric Fields
%Version Date: 9 June 2020

n_electrodes = 32;
n_time_pts = 40;
n_subs = 24;

%Random data for 3 x 2 ANOVA design
data = randn([n_electrodes, n_time_pts, 3, 2, n_subs]);

%Get epsilon estimates
epsilon = GG(data, [], [3,4]);

%Calculate ANOVA
test_results = calc_param_ANOVA(data, [], [3,4], 0.05, 'none', true);

%Check results for a random electrode and time poin
e = randi([1,n_electrodes]);
t = randi([1,n_time_pts]);

fprintf('\nFMUT RESULTS\n');
fprintf('epsilon = %.4f\n', epsilon(e,t));
fprintf('p = %.3f\n', test_results.p(e,t));

%Create ANOVA table
oneway_data = reshape(data, n_electrodes, n_time_pts, [], n_subs);
T = [cell2table(cellfun(@(x) sprintf('S%d', x), num2cell(1:n_subs)', 'UniformOutput', false), 'VariableNames', {'sub'}) ...
    array2table(squeeze(oneway_data(e,t,:,:))')];
T.Properties.VariableNames = {'sub', 'A1B1', 'A2B1' ,'A3B1', 'A1B2', 'A2B2' ,'A3B2'};

%Output to csv for checking in jamovi/SPSS
%writetable(T, 'gg_test.csv');

%MATLAB repetaed measure ANOVA
withindesign = cell2table({'A1', 'B1'; 'A2', 'B1'; 'A3', 'B1'; ...
                           'A1', 'B2'; 'A2', 'B2'; 'A3', 'B2'}, ...
                           'VariableNames', {'A', 'B'});
rm = fitrm(T, 'A1B1-A3B2~1', 'WithinDesign', withindesign);
rm_anova_results = ranova(rm, 'WithinModel', 'A+B+A*B');
fprintf('\nMATLAB RM RESULTS\n');
disp(rm_anova_results);

%Check that MATLAB model gives the same result
assert(rm_anova_results{'(Intercept):A:B', 'pValueGG'} - test_results.p(e,t) < 1e-7);
