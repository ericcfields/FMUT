%Test Various things
%AUTHOR: Eric Fields
%VERSION DATE: 13 July 2017

%Load a GND for testing
if ispc()
    load('R:\Public\GK_lab\Eric\FMUT_development\FMUT\testing\data\EmProb_13subs_Test.GND', '-mat');
    load('R:\Public\GK_lab\Eric\FMUT_development\FMUT\testing\data\Disflu_GroupLevel.GRP', '-mat');
elseif ismac()
    load('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT/testing/data/EmProb_13subs_Test.GND', '-mat');
    load('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT/testing/data/Disflu_GroupLevel.GRP', '-mat')
end

%% Check F-values calculated by approximate int code RB

%%% Two-way %%%

data = normrnd(0, 1, [32, 167, 2, 2, 16]);

%F-test
F_dist = perm_rbANOVA(data, 10);
F_obs = squeeze(F_dist(1, :, :));

%t-test
data = squeeze(data(:, :, 1, :, :) - data(:, :, 2, :, :));
data = squeeze(data(:, :, 1, :) - data(:, :, 2, :));
mn = mean(data, 3);
stderr = std(data, [], 3) / sqrt(size(data, 3));
t_vals = mn ./ stderr;

assert(all(F_obs(:) - t_vals(:).^2 < 1e-9))

%%% Three-way %%%

data = normrnd(0, 1, [32, 167, 2, 2, 2, 16]);

%F-test
F_dist = perm_rbANOVA(data, 10);
F_obs = squeeze(F_dist(1, :, :));

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
F_dist = perm_spANOVA(data, cond_subs, dims, n_perm);
F_obs = squeeze(F_dist(1, :, :));

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
F_dist = perm_spANOVA(data, cond_subs, dims, n_perm);
F_obs = squeeze(F_dist(1, :, :));

%t-test
data = reshape(data(:, :, 1, :, :) - data(:, :, 2, :, :), [32, 167, 2, 16]);
data = reshape(data(:, :, 1, :) - data(:, :, 2, :), [32, 167, 16]);
data = permute(data, [3, 1, 2]);
[~, ~, ~, stats] = ttest2(data(1:8, :, :), data(9:16, :, :));

assert(all(F_obs(:) - stats.tstat(:).^2 < 1e-9))


%% Output of functions is compatible

GND = FmaxGND(GND, ...
              'bins',          [24, 26, 27, 29], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [2, 2], ... 
              'n_perm',        10, ...
              'save_GND',      'no', ...
              'plot_raster',   'no');
          
GND = FclustGND(GND, ...
               'bins',          [24, 26, 27, 29], ...
               'factor_names',  {'Probability', 'Emotion'}, ... 
               'factor_levels', [2, 2], ... 
               'n_perm',        10, ...
               'save_GND',      'no', ...
               'plot_raster',   'no');
           
GND = FfdrGND(GND, ...
              'bins',          [24, 26, 27, 29], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [2, 2], ... 
              'save_GND',      'no', ...
              'plot_raster',   'no');
          
GRP = FmaxGRP(GRP, ... 
              'bins',             4:7, ... 
              'bg_factor_name',   'reliability', ...
              'wg_factor_names',  {'expectedness', 'disfluency'}, ...
              'wg_factor_levels', [2, 2], ... 
              'n_perm',           10, ...
              'save_GRP',         'no', ...
              'plot_raster',      'no');

