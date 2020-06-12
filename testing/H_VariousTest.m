%Test Various things
%
%AUTHOR: Eric Fields
%VERSION DATE: 12 June 2020

%Load a GND for testing
FMUT_dir = fileparts(fileparts(which('A_setupTest')));
load(fullfile(FMUT_dir, 'testing', 'data', 'EmProb_13subs_Test.GND'), '-mat');
load(fullfile(FMUT_dir, 'testing', 'data', 'Disflu_GroupLevel.GRP'), '-mat');


%% Check F-values calculated by approximate int code RB

%%% Two-way %%%

data = normrnd(0, 1, [32, 167, 2, 2, 16]);
dims = [3, 4];
n_perm = 10;

%F-test
F_obs = perm_rbANOVA(data, dims, n_perm, false);

%t-test
data = squeeze(data(:, :, 1, :, :) - data(:, :, 2, :, :));
data = squeeze(data(:, :, 1, :) - data(:, :, 2, :));
mn = mean(data, 3);
stderr = std(data, [], 3) / sqrt(size(data, 3));
t_vals = mn ./ stderr;

assert(all(F_obs(:) - t_vals(:).^2 < 1e-9))

%%% Three-way %%%

data = normrnd(0, 1, [32, 167, 2, 2, 2, 16]);
dims = [3, 4, 5];
n_perm = 10;

%F-test
F_obs = perm_rbANOVA(data, dims, n_perm, false);

%t-test
data = squeeze(data(:, :, 1, :, :, :) - data(:, :, 2, :, :, :));
data = squeeze(data(:, :, 1, :, :) - data(:, :, 2, :, :));
data = squeeze(data(:, :, 1, :) - data(:, :, 2, :));
mn = mean(data, 3);
stderr = std(data, [], 3) / sqrt(size(data, 3));
t_vals = mn ./ stderr;

assert(all(F_obs(:) - t_vals(:).^2 < 1e-9))


%% Check F-values calculated by approximate int code SP

%%% Two-way %%%

data = normrnd(0, 1, [32, 167, 2, 16]);
cond_subs = [8, 8];
dims = [3, 4];
n_perm = 10;

%F-test
F_obs = perm_spANOVA(data, cond_subs, dims, n_perm, false);

%t-test
data = reshape(data(:, :, 1, :) - data(:, :, 2, :), [32, 167, 16]);
data = permute(data, [3, 1, 2]);
[~, ~, ~, stats] = ttest2(data(1:8, :, :), data(9:16, :, :));

assert(all(F_obs(:) - stats.tstat(:).^2 < 1e-9))

%%% Three-way %%%

data = normrnd(0, 1, [32, 167, 2, 2, 16]);
cond_subs = [8, 8];
dims = [3, 4, 5];
n_perm = 10;

%F-test
F_obs = perm_spANOVA(data, cond_subs, dims, n_perm, false);

%t-test
data = reshape(data(:, :, 1, :, :) - data(:, :, 2, :, :), [32, 167, 2, 16]);
data = reshape(data(:, :, 1, :) - data(:, :, 2, :), [32, 167, 16]);
data = permute(data, [3, 1, 2]);
[~, ~, ~, stats] = ttest2(data(1:8, :, :), data(9:16, :, :));

assert(all(F_obs(:) - stats.tstat(:).^2 < 1e-9))


%% Output of functions is compatible and test get_mean_amplitude

time_wind = [300, 500];
n_perm = 100;
chan_hood = .9;

GND = FclustGND(GND, ...
               'bins',          [24, 26, 27, 29], ...
               'factor_names',  {'Probability', 'Emotion'}, ... 
               'factor_levels', [2, 2], ... 
               'time_wind', time_wind, ...
               'n_perm', n_perm, ...
               'chan_hood', chan_hood, ...
               'save_GND', 'no', ...
               'output_file', false, ...
               'plot_raster', 'no');
           
GND = FfdrGND(GND, ...
              'bins',          [24, 26, 27, 29], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [2, 2], ...
              'time_wind', time_wind, ...
              'output_file', false, ...
              'save_GND',      'no', ...
              'plot_raster',   'no');
          
GND = FmaxGND(GND, ...
              'bins',          [24, 26, 27, 29], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [2, 2], ... 
              'time_wind', time_wind, ...
              'n_perm', n_perm, ...
              'save_GND',      'no', ...
              'plot_raster',   'no');
          
GRP = FmaxGRP(GRP, ... 
              'bins',             4:7, ... 
              'bg_factor_name',   'reliability', ...
              'wg_factor_names',  {'expectedness', 'disfluency'}, ...
              'wg_factor_levels', [2, 2], ... 
              'time_wind', time_wind, ...
              'n_perm', n_perm, ...
              'save_GRP',         'no', ...
              'plot_raster',      'no');
          
GRP = FclustGRP(GRP, ... 
              'bins',             4:7, ... 
              'bg_factor_name',   'reliability', ...
              'wg_factor_names',  {'expectedness', 'disfluency'}, ...
              'wg_factor_levels', [2, 2], ... 
              'time_wind',        time_wind, ...
              'n_perm',           n_perm, ...
              'alpha',            0.05, ... 
              'output_file',      false, ...
              'chan_hood',        chan_hood, ...
              'save_GRP',         'no', ...
              'plot_raster',      'no');
          
GRP = FfdrGRP(GRP, ... 
              'bins',             4:7, ... 
              'bg_factor_name',   'reliability', ...
              'wg_factor_names',  {'expectedness', 'disfluency'}, ...
              'wg_factor_levels', [2, 2], ...
              'method',           'bh', ...
              'time_wind',        time_wind, ...
              'output_file',      false, ...
              'save_GRP',         'no', ...
              'plot_raster',      'no');

          
% Test get_mean_amplitude function

%With GND
for i = 1:length(GND.F_tests)
    mean_data_GND = get_mean_amplitude(GND, i, 'effect', 'Emotion', 'output_file', sprintf('outputs/test_get_mean_amplitude_GND_%s.csv', GND.F_tests(i).mult_comp_method));
end

%With GRP
for i = 1:length(GRP.F_tests)
    mean_data_GRP = get_mean_amplitude(GRP, i, 'effect', 'expectedness', 'output_file', sprintf('outputs/test_get_mean_amplitude_GRP_%s.csv', GRP.F_tests(i).mult_comp_method));
end
