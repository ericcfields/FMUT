%Function to conduct an Fmax permutation test for ANOVA designs with a between
%subjects factor
%
%EXAMPLE USAGE
% GRP = FmaxGRP(GRP, 'bins', 1:6, 'bg_factor_name', 'mood', ...
%               'wg_factor_names', {'probability', 'emotion'}, ...
%               'wg_factor_levels', [3, 2], 'time_wind', [500, 800], ...
%               'include_chans', {'Fz', 'Cz', 'Pz'}, 'n_perm', 1e4, ...
%               'alpha', 0.05);
%
%REQUIRED INPUTS
% GRP_or_fname      - A Mass Univariate Toolbox GRP struct or a string
%                     containing a filename of a GRP structure that 
%                     has been saved to disk (with full path if not in current
%                     working directory)
% bins              - Array with bins to use in ANOVA
%
%OPTIONAL INPUTS
% wg_factor_names   - cell array with names of within-subject factors in fastest 
%                     to slowest moving order within the bins provided
%                     {default: no within-subjects factors}
% wg_factor_levels  - number of levels in each within subject factorin fastest 
%                     to slowest moving order within the bins provided
%                     {default: no within-subjects factors}
% bg_factor_name    - A string specifying the name of the between-subjects
%                     factor {default: 'Group'}.
% time_wind      - 2D matrix of pairs of time values specifying the beginning 
%                  and end of the time windows in ms (e.g., 
%                  [200 400; 500, 800]). Every single time point in 
%                  the time window will be individually tested (i.e.,
%                  maximal temporal resolution). Note, boundaries of time 
%                  window(s) may not exactly correspond to desired time 
%                  window boundaries because of temporal digitization (i.e., 
%                  you only have samples every so many ms). 
%                  {default: 0 ms to the end of the epoch}
% mean_wind      - ['yes' or 'no'] If 'yes', the permutation test will be
%                  performed on the mean amplitude within each time window 
%                  specified by time_wind.  This sacrifices temporal 
%                  resolution to increase test power by reducing the number
%                  of comparisons.  If 'no', every single time point within
%                  time_wind's time windows will be tested individually.
%                  {default: 'no'}
% exclude_chans  - A cell array of channel labels to exclude from the
%                  permutation test (e.g., {'A2','VEOG','HEOG'}). This can 
%                  be used to exclude non-data channels (e.g. EOG channels) 
%                  or to increase test power by sacrificing spatial resolution
%                  (i.e., reducing the number of comparisons). Use headinfo.m 
%                  to see the channel labels stored in the GRP variable. You 
%                  cannot use both this option and 'include_chans' (below).
%                  {default: not used, all channels included in test}
% include_chans  - A cell array of channel labels to use in the permutation
%                  test (e.g., {'Fz','Cz','Pz'}). All other channels will
%                  be ignored. This option sacrifices spatial resolution to 
%                  increase test power by reducing the number of comparisons.
%                  Use headinfo.m to see the channel labels stored in the GRP
%                  variable. You cannot use both this option and 
%                  'exclude_chans' (above). 
%                  {default: not used, all channels included in test}
% n_perm         - number of permutations {default: 10,000}
% alpha          - A number between 0 and 1 specifying the family-wise 
%                  alpha level of the test. {default: 0.05}
% plot_raster    - ['yes' or 'no'] If 'yes', a two-dimensional (time x channel)
%                  binary "raster" diagram is created to illustrate the
%                  results of the permutation tests. This figure can be reproduced 
%                  with the function F_sig_raster.m. {default: 'yes'}
% save_GRP       - save GRP to disk, 'yes' or 'no' {default: user will be
%                  prompted}
% output_file    - Name of .xlsx file to output results. {default: no output}
% reproduce_test - [integer] The number of the permutation test stored in
%                  the GRP variable to reproduce.  For example, if 
%                  'reproduce_test' equals 2, the second F-test 
%                  stored in the GRP variable (i.e., GRP.F_tests{2}) will 
%                  be reproduced. Reproduction is accomplished by setting
%                  the random number generator used in the permutation test 
%                  to the same initial state it was in when the permutation 
%                  test was conducted. Obviously other options/inputs must
%                  also be the same to truly reproduce the test
% verblevel      - An integer specifiying the amount of information you want
%                  the Mass Univariate Toolbox to provide about what it is 
%                  doing during runtime.
%                      Options are:
%                        0 - quiet, only show errors, warnings, and EEGLAB reports
%                        1 - stuff anyone should probably know
%                        2 - stuff you should know the first time you start working
%                            with a data set {default value}
%                        3 - stuff that might help you debug (show all
%                            reports)
%
%OUTPUT
% GRP           - GRP struct, with results added in the F_tests field.
%
% optional additional output:
%  (Note: all optional outputs give information already contained in the
%   F_tests field of the GRP struct; they are simply available to make these
%   values more directly accesible and easier to work with)
%
% results       - The same struct added to the F_tests field, but assigned
%                 to its own variable; this might make it easier to do
%                 further operations.
% prm_pval      - The adj_pval field from the results struct. For a one-way
%                 ANOVA, this is an electrodes x time points array of
%                 p-values; for a multi-factor ANOVA, it is a struct with
%                 multiple such arrays
% F_obs         - The F_obs field from the results struct. For a one-way
%                 ANOVA, this is an electrodes x time points array of
%                 Fs; for a multi-factor ANOVA, it is a struct with
%                 multiple such arrays
% F_crit        - The F_crit field from the results struct. For a one-way
%                 ANOVA, this is a single number; for a multi-factor ANOVA, 
%                 it is a struct with F_crit for each effect.
%
%
%DESCRIPTION
%Main effects are calculated by permuting within each condition of the other
%factor(s). For a factor with two levels, this is equivalent to the tmaxGND 
%function (but will of course give results as an F-test).
%
%For interaction effects where more than one factor has more than two levels, 
%it is not possible to conduct a test that controls the Type I error exactly 
%at a specified level. For such designs, this function uses the permutation 
%of residuals method first described by Still & White (1981) and Freedman & Lane (1983). 
%The Type I error rate of this test is asymptotic to the nominal alpha as 
%sample size and/or signal to noise ratio increase.
%For designs where an exact test is possible, this function uses a
%restricted permutation method to conduct an exact test.
%
%See the FMUT documentation for more information:
%https://github.com/ericcfields/FMUT/wiki
%
%VERSION DATE: 14 July 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and SHOULD NOT be considered error free. 

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This function incorporates some code from the Mass Univariate Toolbox, 
%Copyright (c) 2015, David Groppe

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 7/11/17 - First version modified from FmaxGND
% 7/13/17 - updated to eliminated int_method
% 7/14/17 - Command window output moved to separate function


function [GRP, results, prm_pval, F_obs, F_crit] = FmaxGRP(GRP_or_fname, varargin)
    
    warning('You are using a beta version of this function. It needs further testing and should NOT be considered error free.');

    %% ~~~~~PARSE INPUT~~~~~

    global VERBLEVEL
    
    p=inputParser;
    p.addRequired('GRP_or_fname', @(x) ischar(x) || isstruct(x));
    p.addParameter('bins',          [],       @(x) isnumeric(x));
    p.addParameter('wg_factor_names',  '',       @(x) (ischar(x) || iscell(x)));
    p.addParameter('wg_factor_levels', [],       @(x) isnumeric(x));
    p.addParameter('bg_factor_name',   'Group',  @(x) ischar(x));
    p.addParameter('time_wind',     [],       @(x) (isnumeric(x) && size(x, 2)==2));
    p.addParameter('include_chans', [],       @(x) iscell(x));
    p.addParameter('exclude_chans', [],       @(x) iscell(x));
    p.addParameter('n_perm',        1e4,      @(x) isnumeric(x));
    p.addParameter('save_GRP',      'prompt', @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))) || islogical(x));
    p.addParameter('output_file',   false,    @(x) (ischar(x) || islogical(x)));
    p.addParameter('alpha',         0.05,     @(x) (isnumeric(x) && x<=1 && x>=0));
    p.addParameter('reproduce_test',false,    @(x) isnumeric(x));
    p.addParameter('mean_wind',     'no',     @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))));
    p.addParameter('verblevel',     [],       @(x) (isnumeric(x) && length(x)==1 && x>=0 && x<=3))
    p.addParameter('plot_raster',   'yes',    @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))));
    p.addParameter('time_block_dur', []);
    p.addParameter('plot_gui',       []);
    p.addParameter('plot_mn_topo',   []);
    
    p.parse(GRP_or_fname, varargin{:});
    
    if isempty(p.Results.verblevel),
        if isempty(VERBLEVEL),
            VERBLEVEL=2;
        end
    else
       VERBLEVEL=p.Results.verblevel; 
    end
    
    %Assign GRP
    if ischar(GRP_or_fname)
        load(GRP_or_fname, '-mat');
    elseif isstruct(GRP_or_fname)
        GRP = GRP_or_fname;
    else
        error('The GRP variable provided does not seem to be a valid GRP struct or filepath to a GRP struct.');
    end
    
    %Assign some variables for easier reference
    bins             = p.Results.bins;
    wg_factor_names  = p.Results.wg_factor_names;
    wg_factor_levels = p.Results.wg_factor_levels;
    time_wind        = p.Results.time_wind;
    n_perm           = p.Results.n_perm;
    alpha            = p.Results.alpha;
    
    %Check for required name-value inputs
    if isempty(bins)
        error('''bins'' is a required input. See >>help FmaxGRP.');
    end
    
    %Find id numbers for electrodes to use in analysis
    chan_labels = {GRP.chanlocs.labels};
    if ~isempty(p.Results.include_chans) && ~isempty(p.Results.exclude_chans)
        error('You cannot use BOTH ''include_chans'' and ''exclude_chans'' options.');
    elseif ~isempty(p.Results.include_chans)
        electrodes = NaN(1, length(p.Results.include_chans));
        for c = 1:length(p.Results.include_chans)
            if find(strcmp(p.Results.include_chans(c), chan_labels))
                electrodes(c) = find(strcmp(p.Results.include_chans(c), chan_labels));
            else
                error('Electrode ''%s'' does not exist.', p.Results.include_chans{c});
            end
        end
    elseif ~isempty(p.Results.exclude_chans)
        if ~all(ismember(p.Results.exclude_chans, chan_labels))
            missing_channels = p.Results.exclude_chans(~ismember(p.Results.exclude_chans, chan_labels));
            error([sprintf('The following channels appear in ''include_chans'' but do not appear in GRP.chanlocs.labels:\n') ... 
                   sprintf('%s ', missing_channels{:})])
        else
            electrodes = find(~ismember(chan_labels, p.Results.exclude_chans));
        end
    else
        electrodes = 1:length(GRP.chanlocs);
    end
    
    %MUT reatures not implemented here              
    if ~isempty(p.Results.time_block_dur)
        error('The ''time_block_dur'' option is not implemented for FmaxGRP. You''ll need to divide the time windows manually.');
    end              
    if ~isempty(p.Results.plot_gui)
        watchit('''plot_gui'' is not implemented for FmaxGRP.');
    end
    if ~isempty(p.Results.plot_mn_topo)
        watchit('''plot_mn_topo'' is not implemented for FmaxGRP.');
    end
    
    %Standardize formatting
    if ischar(wg_factor_names)
        wg_factor_names = {wg_factor_names};
    end
    time_wind = sort(time_wind, 2);
    time_wind = sort(time_wind, 1);
    
    %Set defaults for missing arguments
    if isempty(time_wind)
        time_wind = [0, GRP.time_pts(end)];
    end
    
    %Check for errors in input
    if ~isempty(wg_factor_levels)
        if length(wg_factor_names) ~= length(wg_factor_levels)
            error('The number of factors does not match in the ''wg_factor_names'' and ''wg_factor_levels'' inputs');
        end
        if isempty(wg_factor_names{1})
            error('''wg_factor_levels'' indicates a within-subjects factor, but no ''wg_factor_names'' input was given.')
        end
    elseif ~isempty(wg_factor_names{1})
        error('''wg_factor_names'' indicates a within-subjects factor, but no ''wg_factor_levels'' input was given.')
    end
    if sum(wg_factor_levels>2) > 2
        error('FmaxGRP cannot handle split plot designs with more than two within-subjects factors with more than two levels')
    end
    if prod(wg_factor_levels) ~= length(bins)
        error('Number of bins does not match the design specified by thte ''wg_factor_levels'' input.')
    end
    if ~isequal(reshape(time_wind', 1, []), unique(reshape(time_wind', 1, [])))
        error('When multiple time windows are provided, they cannot overlap.')
    end
    if alpha <= .01 && n_perm < 5000,
        watchit(sprintf('You are probably using too few permutations for an alpha level of %f.',alpha));
    elseif alpha <=.05 && n_perm < 1000,
        watchit(sprintf('You are probably using too few permutations for an alpha level of %f.',alpha));
    end
    if p.Results.reproduce_test
        if ~isfield(GRP, 'F_tests')
            error('You tried to reproduce test %d, but there are no results in GRP.F_tests.', p.Results.reproduce_test);
        elseif p.Results.reproduce_test > length(GRP.F_tests)
            error('You tried to reproduce test %d, but there are only %d tests in GRP.F_tests.', p.Results.reproduce_test, length(GRP.F_tests));
        end
    end

    
    %% ~~~~~ SET-UP ~~~~~

    %Get or set random # generator state
    if verLessThan('matlab','8.1')
        defaultStream=RandStream.getDefaultStream; 
    else
        defaultStream=RandStream.getGlobalStream;
    end
    if p.Results.reproduce_test
        seed_state = GRP.F_tests(p.Results.reproduce_test).seed_state;
        defaultStream.State = seed_state;
    else
        seed_state = defaultStream.State;
    end
    
    %Find time points or mean windows to use and extract the data for
    %analysis
    the_data = [];
    n_electrodes = length(electrodes);
    for g = 1:length(GRP.GND_fnames)
        
        %Load GND and check for errors
        load(GRP.GND_fnames{g}, '-mat')
        if ~exist(GRP.GND_fnames{g}, 'file')
            error('%s does not exist.', GRP.GND_fnames{g})
        elseif ~exist('GND', 'var')
            error('%s does not appear to contain a .GND variable', GRP.GND_fnames{g})
        end
        if ~all(all(GND.indiv_bin_ct(:, bins)))
            watchit(sprintf('Some subjects in\n%s\nappear to be missing data from bins used in this test!\nSee: GRP.indiv_bins_ct.', GRP.GND_fnames{g}));
        end
        
        cond_subs(1, g) = size(GND.indiv_erps, 4); %#ok<AGROW>
        
        %Get data (individual time points)
        if ~strcmpi(p.Results.mean_wind, 'yes') && ~strcmpi(p.Results.mean_wind, 'y')
            use_time_pts = [];
            for i = 1:size(time_wind, 1)
                [~, start_sample] = min(abs( GND.time_pts - time_wind(i, 1) ));
                [~, end_sample  ] = min(abs( GND.time_pts - time_wind(i, 2) ));
                time_wind(i, 1) = GND.time_pts(start_sample);
                time_wind(i, 2) = GND.time_pts(end_sample);
                use_time_pts = [use_time_pts start_sample:end_sample]; %#ok<AGROW>
                if VERBLEVEL && g == 1
                    if i == 1
                        fprintf('\nConducting Fmax permutation test from %d ms to %d ms', GND.time_pts(start_sample), GND.time_pts(end_sample));
                    else
                        fprintf(', %d ms to %d ms.', GND.time_pts(start_sample), GND.time_pts(end_sample));
                    end
                    if i == size(time_wind, 1)
                        fprintf('\n');
                    end
                end
            end
            n_time_pts = length(use_time_pts);
            the_data = cat(4, the_data, GND.indiv_erps(electrodes, use_time_pts, bins, :));
        
        %Get data (mean time window)
        else
            n_time_pts = size(time_wind, 1);
            use_time_pts = cell(1, n_time_pts);
            new_data = NaN(n_electrodes, n_time_pts, n_conds, cond_subs(g));
            for i = 1:size(time_wind, 1)
                [~, start_sample] = min(abs( GND.time_pts - time_wind(i, 1) ));
                [~, end_sample  ] = min(abs( GND.time_pts - time_wind(i, 2) ));
                time_wind(i, 1) = GND.time_pts(start_sample);
                time_wind(i, 2) = GND.time_pts(end_sample);
                use_time_pts{i} = start_sample:end_sample;
                new_data(:, i, :, :) = mean(GND.indiv_erps(electrodes, start_sample:end_sample, bins, :), 2);
                if VERBLEVEL && g == 1
                    if i == 1
                        fprintf('\nConducting Fmax permutation test in mean time windows %d-%d ms', GND.time_pts(start_sample), GND.time_pts(end_sample));
                    else
                        fprintf(', %d-%d ms', GND.time_pts(start_sample), GND.time_pts(end_sample));
                    end
                    if i == size(time_wind, 1)
                        fprintf('\n');
                    end
                end
            end
            the_data = cat(4, the_data, new_data);
        end
        
        clear GND
    
    end
    clear new_data
    
    %Report test information
    if VERBLEVEL
        fprintf('Number of channels: %d\n', size(the_data, 1));
        fprintf('Number of time points: %d\n', size(the_data, 2));
        fprintf('Total comparisons: %d\n', numel(the_data(:, :, 1, 1)));
        fprintf('Number of subjects: %d\n', size(the_data, 4));
    end
    
    %Divide the factors into separate dimensions for factorial ANOVA
    if length(wg_factor_levels) > 1
        the_data = reshape(the_data, [n_electrodes, n_time_pts, wg_factor_levels, sum(cond_subs)]);
    elseif isempty(wg_factor_levels)
        the_data = reshape(the_data, [n_electrodes, n_time_pts, sum(cond_subs)]);
    end
    
    %Figure out the effects we need to calculate
    if wg_factor_levels
        factor_names  = [wg_factor_names p.Results.bg_factor_name];
        factor_levels = [wg_factor_levels length(cond_subs)];
    else
        factor_names = {p.Results.bg_factor_name};
        factor_levels = length(cond_subs);
    end
    [effects, effects_labels] = get_effects(factor_names);
    
    
    %% ~~~~~ RUN PERMUTATION ANOVAS ~~~~~
    
    test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                                 'F_obs', NaN(n_electrodes, n_time_pts),  'Fmax_crit', NaN, ... 
                                 'df', NaN(1, 2), 'estimated_alpha', NaN, 'exact_test', NaN), ... 
                                 length(effects), 1);
    for i = 1:length(effects)
        if VERBLEVEL
            fprintf('\nCalculating %s effect\n', effects_labels{i});
        end
        %Calculate test
        test_results(i) = calc_Fmax(the_data, cond_subs, effects{i}+2, n_perm, alpha);       
    end
    

    %% ~~~~~ ADD RESULTS STRUCT TO GRP AND ASSIGN OTHER OUTPUT ~~~~~
    
    if (strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y'))
        use_time_pts = {use_time_pts};
    end
    %Create results struct with basic parameters
    results = struct('bins', bins, ...
                     'factors', {factor_names}, ...
                     'factor_levels', factor_levels, ...
                     'time_wind', time_wind, ...
                     'used_tpt_ids', use_time_pts, ...
                     'mean_wind', p.Results.mean_wind, ...
                     'include_chans', {{GRP.chanlocs(electrodes).labels}}, ...
                     'used_chan_ids', electrodes, ...
                     'mult_comp_method', 'Fmax perm test', ...
                     'n_perm', n_perm, ...
                     'desired_alphaORq', alpha, ...
                     'estimated_alpha', [], ...
                     'seed_state', seed_state, ...
                     'exact_test', [], ...
                     'null_test', [], ...
                     'adj_pval', [], ...
                     'F_obs', [], ...
                     'F_crit', [], ...
                     'df', [], ...              
                     'chan_hood', NaN, ...
                     'clust_info', NaN, ...
                     'fdr_rej', NaN);
    
    %Add statistical results
    assert(length(effects) == length(test_results));
    if length(effects) == 1
        results.null_test = test_results.h;
        results.adj_pval  = test_results.p;
        results.F_obs     = test_results.F_obs;
        results.F_crit    = test_results.Fmax_crit;
        results.df        = test_results.df;
        results.estimated_alpha = test_results.estimated_alpha;
        results.exact_test = test_results.exact_test;
    else
        for i = 1:length(effects)
            results.null_test.(effects_labels{i}) = test_results(i).h;
            results.adj_pval.(effects_labels{i})  = test_results(i).p;
            results.F_obs.(effects_labels{i})     = test_results(i).F_obs;
            results.F_crit.(effects_labels{i})    = test_results(i).Fmax_crit;
            results.df.(effects_labels{i})        = test_results(i).df;
            results.estimated_alpha.(effects_labels{i}) = test_results(i).estimated_alpha;
            results.exact_test.(effects_labels{i}) = test_results(i).exact_test;
        end
    end
                 
    %Add results struct to GRP
    if ~isfield(GRP, 'F_tests') || isempty(GRP.F_tests)
        GRP.F_tests = results;
    else
        GRP.F_tests(end+1) = results;
    end
    
    %Optional outputs
    if nargout > 2
        prm_pval = results.adj_pval;
        F_obs    = results.F_obs;
        F_crit   = results.F_crit;
    end
    

    %% ~~~~~ OUTPUT RESULTS ~~~~~
    
    %Output results to command window
    if VERBLEVEL
        report_results(GRP, length(GRP.F_tests))
    end
    
    %Plot results
    if ~strcmpi(p.Results.plot_raster, 'no') && ~strcmpi(p.Results.plot_raster, 'n')
        if VERBLEVEL
            fprintf('Generating raster plot:\n');
        end
        if length(effects_labels) == 1
            F_sig_raster(GRP, length(GRP.F_tests), 'use_color', 'rgb');
        else
            for i = 1:length(effects_labels)
                F_sig_raster(GRP, length(GRP.F_tests), 'effect', effects_labels{i}, 'use_color', 'rgb');
            end
        end
    end
    
    %Prompt user about saving GRP
    if ~strcmpi(p.Results.save_GRP, 'no') && ~strcmpi(p.Results.save_GRP, 'n')
        GRP = save_matmk(GRP);
    end
    
    %Output to spreadsheet if requested
    if p.Results.output_file
        if VERBLEVEL
            fprintf('\nWriting results to %s . . . ', p.Results.output_file)
        end
        Ftest2xls(GRP, length(GRP.F_tests), p.Results.output_file);
        if VERBLEVEL
            fprintf('DONE.\n\n')
        end
    end
    

end
