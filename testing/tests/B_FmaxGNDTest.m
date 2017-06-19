%Test FmaxGND function
%AUTHOR: Eric Fields
%VERSION DATE: 15 June 2017

%Load a GND for testing
if ispc()
    load('R:\Public\GK_lab\Eric\FMUT_development\FMUT_functions\testing\EmProb_13subs_Test.GND', '-mat');
elseif ismac()
    load('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT_functions/testing/EmProb_13subs_Test.GND', '-mat');
end

%Define some general variables
time_wind = [500, 800];
include_chans = {'Cz', 'CPz', 'Pz'};
[~, start_sample] = min(abs( GND.time_pts - time_wind(1) ));
[~, end_sample  ] = min(abs( GND.time_pts - time_wind(2) ));
electrodes = NaN(1, length(include_chans));
for c = 1:length(include_chans)
    electrodes(c) = find(strcmp(include_chans(c), {GND.chanlocs.labels}));
end

%% Exact interaction

GND = FmaxGND(GND, ...
              'bins',          [24, 26, 27, 29], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [2, 2], ... 
              'time_wind',     time_wind, ... 
              'include_chans', include_chans, ... 
              'n_perm',        1e4, ...
              'alpha',         0.05, ...
              'save_GND',      'no', ...
              'output_file',   fullfile('outputs', 'Fmax_test.xlsx'), ...
              'plot_raster',   'yes');

%F==t^2
assert(all(all(GND.F_tests(end).F_obs.ProbabilityXEmotion - GND.grands_t(electrodes, start_sample:end_sample, 51).^2 < 1e-4)));

%Check p-values with t-test
GND = tmaxGND(GND, 51, ...
              'time_wind', time_wind, ...
              'include_chans', include_chans, ...
              'n_perm', 1e4, ...
              'alpha', 0.05, ...
              'save_GND', 'no', ...
              'plot_gui', 'no', ...
              'plot_raster', 'no', ...
              'verblevel', 2);
ttest2xls(GND, length(GND.t_tests), fullfile('outputs', 'ttest2xls_test.xlsx'))


%Critical F == critical t^2
assert(GND.F_tests(end).F_crit.ProbabilityXEmotion - GND.t_tests(end).crit_t(1)^2 < 1);
%p-values are approximately the same
assert(all(all(GND.F_tests(end).adj_pval.ProbabilityXEmotion - GND.t_tests(end).adj_pval < .05)))


%% Oneway ANOVA

GND = FmaxGND(GND, ...
              'bins',          [24, 26, 27], ...
              'factor_names',  {'NEU_Probability'}, ... 
              'factor_levels', 3, ... 
              'time_wind',     time_wind, ... 
              'include_chans', include_chans, ... 
              'n_perm',        500, ...
              'alpha',         0.05, ...
              'save_GND',      'no', ...
              'output_file',   fullfile('outputs', 'Fmax_test_oneway.xlsx'), ...
              'plot_raster',   'yes');

%% Approximate interaction  

GND = FmaxGND(GND, ...
              'bins',          [24:29, 31, 32, 33], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [3, 3], ... 
              'time_wind',     time_wind, ... 
              'include_chans', include_chans, ... 
              'n_perm',        500, ...
              'alpha',         0.05, ...
              'save_GND',      'no', ...
              'plot_raster',   'no');

          
%% Mean time window

%F-test
GND = FmaxGND(GND, ...
              'bins',          [24, 26, 27, 29], ...
              'factor_names',  {'Probability', 'Emotion'}, ... 
              'factor_levels', [2, 2], ... 
              'time_wind',     time_wind, ...
              'mean_wind',     'yes', ...
              'include_chans', include_chans, ... 
              'n_perm',        1e4, ...
              'alpha',         0.05, ...
              'save_GND',      'no', ...
              'plot_raster',   'no');

%Calculate t-test on same data
GND = tmaxGND(GND, 51, ...
              'n_perm',       500, ...
              'time_wind',    time_wind, ...
              'mean_wind',    'yes', ...
              'include_chans', include_chans, ...
              'plot_raster',  'no', ...
              'plot_gui',     'no', ...
              'plot_mn_topo', 'no', ...
              'save_GND',     'no', ...
              'verblevel',    2);

%F==t^2
assert(all(GND.F_tests(end).F_obs.ProbabilityXEmotion - GND.t_tests(end).data_t.^2 < 1e-9))
