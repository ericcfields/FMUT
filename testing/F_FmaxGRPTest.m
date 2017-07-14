%Test FmaxGND function
%AUTHOR: Eric Fields
%VERSION DATE: 14 July 2017

global test_xls_output
if isempty(test_xls_output)
    test_xls_output = true;
end

%Load GRP
if ispc()
    load('R:\Public\GK_lab\Eric\FMUT_development\FMUT\testing\data\Disflu_GroupLevel.GRP', '-mat');
elseif ismac()
    load('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT/testing/data/Disflu_GroupLevel.GRP', '-mat')
end

%Define some general variables
time_wind = [300, 500];
include_chans = {'Fz', 'Cz', 'Pz'};
[~, start_sample] = min(abs( GRP.time_pts - time_wind(1) ));
[~, end_sample  ] = min(abs( GRP.time_pts - time_wind(2) ));
electrodes = NaN(1, length(include_chans));
for c = 1:length(include_chans)
    electrodes(c) = find(strcmp(include_chans(c), {GRP.chanlocs.labels}));
end
n_perm = 100;

%% Exact interaction

if test_xls_output
    output_file = fullfile('outputs', 'FmaxGRP_test.xlsx');
else
    output_file = false;
end 
GRP = FmaxGRP(GRP, ... 
              'bins',             4:7, ... 
              'bg_factor_name',   'reliability', ...
              'wg_factor_names',  {'expectedness', 'disfluency'}, ...
              'wg_factor_levels', [2, 2], ... 
              'time_wind',        time_wind, ...
              'include_chans',    include_chans, ... 
              'n_perm',           n_perm, ...
              'alpha',            0.05, ... 
              'output_file',      output_file, ...
              'save_GRP',         'no');

%F==t^2
assert(all(all(GRP.F_tests(end).F_obs.expectednessXdisfluencyXreliability - GRP.grands_t(electrodes, start_sample:end_sample, 78).^2 < 1e-4)));


% GND = tmaxGRP(GRP, 78, ...
%               'time_wind', time_wind, ...
%               'include_chans', include_chans, ...
%               'n_perm', 1e4, ...
%               'alpha', 0.05, ...
%               'save_GRP', 'no', ...
%               'plot_gui', 'no', ...
%               'plot_raster', 'no', ...
%               'verblevel', 2);
%           
% %Critical F == critical t^2
% assert(GRP.F_tests(end).F_crit.expectednessXdisfluencyXreliability - GND.t_tests(end).crit_t(1)^2 < 1);

          
%% One-way
%Completely randomized ANOVA

if test_xls_output
    output_file = fullfile('outputs', 'FmaxGRP_oneway_test.xlsx');
else
    output_file = false;
end 
GRP = FmaxGRP(GRP, ... 
              'bins',             4, ... 
              'bg_factor_name',   'reliability', ...
              'time_wind',        time_wind, ...
              'include_chans',    include_chans, ... 
              'n_perm',           n_perm, ...
              'alpha',            0.05, ... 
              'output_file',      output_file, ...
              'save_GRP',         'no');

%F==t^2
assert(all(all(GRP.F_tests(end).F_obs - GRP.grands_t(electrodes, start_sample:end_sample, 4).^2 < 1e-4)));

%% Aproximate interaction

GRP = FmaxGRP(GRP, ... 
              'bins',             [1, 4:9, 22:23],  ... 
              'bg_factor_name',   'reliability', ...
              'wg_factor_names',  {'expectedness', 'disfluency'}, ...
              'wg_factor_levels', [3, 3], ... 
              'time_wind',        time_wind, ...
              'include_chans',    include_chans, ... 
              'n_perm',           n_perm, ...
              'alpha',            0.05, ... 
              'save_GRP',         'no');
