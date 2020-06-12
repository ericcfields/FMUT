%For testing ANOVA output from various designs (with various reduction
%schemes) against stats software (e.g., jamovi, JASP, SPSS). Each section
%runs a test and outputs the data in a relevant csv for independent testing
%
%AUTHOR: Eric Fields
%VERSION DATE: 12 June 2020

%% Set-up

clearvars;
global VERBLEVEL;
VERBLEVEL = 0;

data = normrnd(0, 1, [1, 1, 3, 3, 2, 16]);
n_perm = 10;
alpha = .05;

%% AxB interaction

%FMUT
dims = [3, 4];
results = calc_Fmax(data, [], dims, n_perm, alpha);
fprintf('\nF_obs = %f\n', results.F_obs)
fprintf('df = [%d, %d]\n\n', results.df(1), results.df(2))

%Output
output_data = reshape(mean(data, 5), [9, 16])';
header = {};
for i = 1:9
    header = [header char(64+i)]; %#ok<AGROW>
end
T = cell2table(num2cell(output_data), 'VariableNames', header);
writetable(T, 'outputs/test_3x3_int.csv')

%% BxC interaction

%FMUT
dims = [4, 5];
results = calc_Fmax(data, [], dims, n_perm, alpha);
fprintf('\nF_obs = %f\n', results.F_obs)
fprintf('df = [%d, %d]\n\n', results.df(1), results.df(2))

%Output
output_data = reshape(mean(data, 3), [6, 16])';
header = {};
for i = 1:6
    header = [header char(64+i)]; %#ok<AGROW>
end
T = cell2table(num2cell(output_data), 'VariableNames', header);
writetable(T, 'outputs/test_3x2_int.csv')

%% AxBxC interaction

%FMUT
dims = [3, 4, 5];
results = calc_Fmax(data, [], dims, n_perm, alpha);
fprintf('\nF_obs = %f\n', results.F_obs)
fprintf('df = [%d, %d]\n\n', results.df(1), results.df(2))

%Output
output_data = reshape(data, [18, 16])';
header = {};
for i = 1:18
    header = [header char(64+i)]; %#ok<AGROW>
end
T = cell2table(num2cell(output_data), 'VariableNames', header);
writetable(T, 'outputs/test_3x3x2_int.csv')

