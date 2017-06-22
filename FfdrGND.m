%Function to conduct an ANOVA with FDR correction across electrodes and
%time points for one-way and factorial ANOVA
%
%EXAMPLE USAGE
%
% GND = FfdrGND(GND, 'bins', 1:6, 'factor_names', {'probability', 'emotion'}, ...
%                 'factor_levels', [3, 2], 'time_wind', [300, 900], ...
%                 'include_chans', {'Fz', 'Cz', 'Pz'}, 'n_perm', 1e4, ...
%                 'method', 'bh', 'q', 0.05) 
%
%
%REQUIRED INPUTS
% GND_or_fname   - A Mass Univariate Toolbox GND struct or a string
%                  containing a filename of a GND structure that 
%                  has been saved to disk (with full path if not in current
%                  working directory)
% bins           - array with bins to use in ANOVA
% factor_names   - cell array with names of factors in fastest to slowest
%                  moving order within the bins provided
% factor_levels  - number of factors in each level in the same order
%
%OPTIONAL INPUTS
% q              - A number between 0 and 1 specifying the family-wise
%                  q level of the test. q is the upper bound on the 
%                  expected proportion of rejected null hypotheses that are
%                  false rejections (i.e., the FDR). {default: 0.05}
% method         - ['bh', 'by', or 'bky'] The procedure used to control
%                  the FDR. 'bh' is the classic Benjamini & Hochberg (1995)
%                  procedure, which is guaranteed to control FDR when the 
%                  tests are independent or positively dependent (e.g., 
%                  positively correlated Gaussians). 'by' is a much more
%                  conservative version of 'bh' that always controls FDR
%                  (regardless of the dependency structure of the tests--
%                  Benjamini & Yekutieli, 2001). 'bky' is a "two-stage"
%                  version of 'bh' that is more powerful than 'bh' when a 
%                  lot of the null hypotheses are false (Benjamini, Krieger, &
%                  Yekutieli, 2006).  'bky' is guaranteed to control FDR when the
%                  tests are independent and tends to be slightly less
%                  powerful than 'bh' when few or no null hypothese are
%                  false. {default: 'bh'}
% time_wind      - 2D matrix of time values specifying the beginning 
%                  and end of the time windows in ms (e.g., 
%                  [500, 800]). Every single time point in 
%                  the time window will be individually tested (i.e.,
%                  maximal temporal resolution). Note, boundaries of time 
%                  window(s) may not exactly correspond to desired time 
%                  window boundaries because of temporal digitization (i.e., 
%                  you only have samples every so many ms). 
%                  {default: 0 ms to the end of the epoch}
% mean_wind      - ['yes' or 'no'] If 'yes', the test will be
%                  performed on the mean amplitude within the time window 
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
%                  to see the channel labels stored in the GND variable. You 
%                  cannot use both this option and 'include_chans' (below).
%                  {default: not used, all channels included in test}
% include_chans  - A cell array of channel labels to use in the permutation
%                  test (e.g., {'Fz','Cz','Pz'}). All other channels will
%                  be ignored. This option sacrifices spatial resolution to 
%                  increase test power by reducing the number of comparisons.
%                  Use headinfo.m to see the channel labels stored in the GND
%                  variable. You cannot use both this option and 
%                  'exclude_chans' (above). 
%                  {default: not used, all channels included in test}
% plot_raster    - ['yes' or 'no'] If 'yes', a two-dimensional (time x channel)
%                  binary "raster" diagram is created to illustrate the
%                  results of the tests. This figure can be reproduced 
%                  with the function F_sig_raster.m. {default: 'yes'}
% save_GND       - save GND to disk, 'yes' or 'no' {default: user will be
%                  prompted}
% output_file    - Name of .xlsx file to output results. Currently 
%                  only works on Windows. {default: no output}
% verblevel      - An integer specifiying the amount of information you want
%                  the Mass Univariate Toolbox to provide about what it is 
%                  doing during runtime. Note that little to no output is
%                  currently provided by this function.
%                      Options are:
%                        0 - quiet, only show errors, warnings, and EEGLAB reports
%                        1 - stuff anyone should probably know
%                        2 - stuff you should know the first time you start working
%                            with a data set {default value}
%                        3 - stuff that might help you debug (show all
%                            reports)
%
%
%OUTPUT
% GND           - GND struct, with results added in the F_tests field.
%
% optional additional output:
%  (Note: all optional outputs give information already contained in the
%   F_tests field of the GND struct; they are simply available to make these
%   values more directly accesible and easier to work with)
%
% results       - The same struct added to the F_tests field, but assigned
%                 to its own variable; this might make it easier to do
%                 further operations.
% adj_pval      - The adj_pval field from the results struct. For a one-way
%                 ANOVA, this is an electrodes x time points array of
%                 p-values; for a multi-factor ANOVA, it is a struct with
%                 multiple such arrays.
% F_obs         - The F_obs field from the results struct. For a one-way
%                 ANOVA, this is an electrodes x time points array of
%                 Fs; for a multi-factor ANOVA, it is a struct with
%                 multiple such arrays.
% F_crit        - The F-value corresponding to the FDR adjusted
%                 significance threshold
%
%FDR REFERENCES:
% Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. Journal of  the Royal Statistical Society. Series B (Methodological), 57(1), 289-300.
% Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006). Adaptive linear step-up procedures that control the false discovery rate. Biometrika, 93(3), 491-507. 
% Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. The Annals of Statistics, 29(4), 1165-1188. 
%
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%VERSION DATE: 21 June 2017
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
%
% 5/15/17         - First working version
% 6/2/17          - Electrode order in output when using include_chans now
%                   matches MUT functions
% 6/3/17          - Fixed problems with factorial ANOVA
% 6/12/17         - Added verblevel reports
% 6/14/17         - Updated error for incorrectly supplied electrode name;
%                   fixed command window output for one-way ANOVA
% 6/20/17         - Command window output for mean window analyses
% 6/21/17         - time_wind field of results struct is now accurate;
%                   changed used_tpt_ids field to cell array for mean window
%                   analyses

function [GND, results, adj_pval, F_obs, F_crit] = FfdrGND(GND_or_fname, varargin)

    warning('You are using a beta version of this function. It needs further testing and should NOT be considered error free.');

    %% ~~~~~PARSE INPUT~~~~~

    global VERBLEVEL
    
    p=inputParser;
    p.addRequired('GND_or_fname', @(x) ischar(x) || isstruct(x));
    p.addParameter('bins',          [],       @(x) isnumeric(x));
    p.addParameter('factor_names',  '',       @(x) (ischar(x) || iscell(x)));
    p.addParameter('factor_levels', '',       @(x) isnumeric(x));
    p.addParameter('time_wind',     [],       @(x) (isnumeric(x) && size(x, 2)==2));
    p.addParameter('include_chans', [],       @(x) iscell(x));
    p.addParameter('exclude_chans', [],       @(x) iscell(x));
    p.addParameter('save_GND',      'prompt', @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))) || islogical(x));
    p.addParameter('output_file',   false,    @(x) (ischar(x) || islogical(x)));
    p.addParameter('mean_wind',     'no',     @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))));
    p.addParameter('verblevel',     [],       @(x) (isnumeric(x) && length(x)==1 && x>=0 && x<=3))
    p.addParameter('plot_raster',   'yes',    @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))));
    p.addParameter('q',             0.05,     @(x) (isnumeric(x) && x<=1 && x>=0));
    p.addParameter('method',        'bh',     @(x) any(strcmpi(x, {'bh', 'by', 'bky'})));
    p.addParameter('time_block_dur', []);
    p.addParameter('plot_gui',       []);
    p.addParameter('plot_mn_topo',   []);
    
    p.parse(GND_or_fname, varargin{:});
    
    if isempty(p.Results.verblevel),
        if isempty(VERBLEVEL),
            VERBLEVEL=2;
        end
    else
       VERBLEVEL=p.Results.verblevel; 
    end
    
    %Assign GND
    if ischar(GND_or_fname)
        load(GND_or_fname, '-mat');
    elseif isstruct(GND_or_fname)
        GND = GND_or_fname;
    else
        error('The GND variable provided does not seem to be a valid GND struct or filepath to a GND struct');
    end
    
    %Assign some variables for easier reference
    bins          = p.Results.bins;
    factor_names  = p.Results.factor_names;
    factor_levels = p.Results.factor_levels;
    time_wind     = p.Results.time_wind;
    q             = p.Results.q;
    method        = p.Results.method;
    
    %Check for required name-value inputs
    if isempty(bins)
        error('''bins'' is a required input. See help FmaxGND');
    end
    if isempty(factor_names)
        error('''factor_names'' is a required input. See help FmaxGND');
    end
    if isempty(factor_levels)
        error('''factor_levels'' is a required input. See help FmaxGND');
    end
    
    %Find id numbers for electrodes to use in analysis
    chan_labels = {GND.chanlocs.labels};
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
            error([sprintf('The following channels appear in ''include_chans'' but do not appear in GND.chanlocs.labels:\n') ... 
                   sprintf('%s ', missing_channels{:})])
        else
            electrodes = find(~ismember(chan_labels, p.Results.exclude_chans));
        end
    else
        electrodes = 1:length(GND.chanlocs);
    end
    
    %MUT features not implemented here              
    if ~isempty(p.Results.time_block_dur)
        error('The ''time_block_dur'' option is not implemented for FmaxGND. You''ll need to divide the time windows manually');
    end              
    if ~isempty(p.Results.plot_gui)
        watchit('''plot_gui'' is not implemented for FmaxGND.');
    end
    if ~isempty(p.Results.plot_mn_topo)
        watchit('''plot_mn_topo'' is not implemented for FmaxGND.');
    end
    
    %Standardize formatting
    if ischar(factor_names)
        factor_names = {factor_names};
    end
    time_wind = sort(time_wind, 2);
    time_wind = sort(time_wind, 1);
    
    %Set defaults for missing arguments
    if isempty(time_wind)
        time_wind = [0, GND.time_pts(end)];
    end
    
    %Check for errors in input
    if length(factor_names) ~= length(factor_levels)
        error('The number of factors does not match in the ''factor_names'' and ''factor_levels'' inputs');
    end
    if length(factor_levels) > 2
        warning('This function has not been tested extensively with designs with more than two factors. Proceed with caution!');
    end
    if prod(factor_levels) ~= length(bins)
        error('Number of bins does not match design.')
    end
    if length(factor_levels)>3 && sum(factor_levels>2)>1
        error('FMUT has limited support for designs with more than three factors.')
    end
    if ~isequal(reshape(time_wind', 1, []), unique(reshape(time_wind', 1, [])))
        error('When multiple time windows are provided, they cannot overlap.')
    end
    if ~all(all(GND.indiv_bin_ct(:, bins)))
        watchit(sprintf('Some subjects appear to be missing data from bins used in this test!\nSee: GND.indiv_bins_ct'));
    end

    
    %% ~~~~~ SET-UP ~~~~~
    
    %Some useful numbers
    n_conds      = prod(factor_levels);
    n_subs       = size(GND.indiv_erps, 4);
    n_electrodes = length(electrodes); 

    %Find time points or mean windows to use and extract the data for
    %analysis
    if ~strcmpi(p.Results.mean_wind, 'yes') && ~strcmpi(p.Results.mean_wind, 'y')
        use_time_pts = [];
        for i = 1:size(time_wind, 1)
            [~, start_sample] = min(abs( GND.time_pts - time_wind(i, 1) ));
            [~, end_sample  ] = min(abs( GND.time_pts - time_wind(i, 2) ));
            time_wind(i, 1) = GND.time_pts(start_sample);
            time_wind(i, 2) = GND.time_pts(end_sample);
            use_time_pts = [use_time_pts start_sample:end_sample]; %#ok<AGROW>
            if VERBLEVEL
                if i == 1
                    fprintf('\nConducting Fmax permutation test from %d ms to %d ms', GND.time_pts(start_sample), GND.time_pts(end_sample));
                else
                    fprintf(', %d ms to %d ms', GND.time_pts(start_sample), GND.time_pts(end_sample));
                end
                if i == size(time_wind, 1)
                    fprintf('\n');
                end
            end
        end
        n_time_pts = length(use_time_pts);
        the_data = GND.indiv_erps(electrodes, use_time_pts, bins, :);
    else
        n_time_pts = size(time_wind, 1);
        use_time_pts = cell(n_time_pts, 1);
        the_data = NaN(n_electrodes, n_time_pts, n_conds, n_subs);
        for i = 1:size(time_wind, 1)
            [~, start_sample] = min(abs( GND.time_pts - time_wind(i, 1) ));
            [~, end_sample  ] = min(abs( GND.time_pts - time_wind(i, 2) ));
            time_wind(i, 1) = GND.time_pts(start_sample);
            time_wind(i, 2) = GND.time_pts(end_sample);
            use_time_pts{i} = start_sample:end_sample;
            the_data(:, i, :, :) = mean(GND.indiv_erps(electrodes, start_sample:end_sample, bins, :), 2);
            if VERBLEVEL
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
    end
    
    %Report test information
    if VERBLEVEL
        if strcmpi(method, 'bh')
            fprintf('FDR control procedure: Benjamini & Hochberg (independent or positive dependency)\n');
        elseif strcmpi(method, 'by')
            fprintf('FDR control procedure: Benjamini & Yekutieli (arbitrary dependency)\n');
        elseif strcmpi(method, 'bky')
            fprintf('FDR control procedure: Benjamini, Krieger, & Yekutieli (two-stage)\n');
        end
        fprintf('Number of channels: %d\n', size(the_data, 1));
        fprintf('Number of time points: %d\n', size(the_data, 2));
        fprintf('Total comparisons: %d\n', numel(the_data(:, :, 1, 1)));
        fprintf('Number of subjects: %d\n', size(the_data, 4));
    end
    
    %Divide the factors into separate dimensions for factorial ANOVA
    if length(factor_levels) > 1
        the_data = reshape(the_data,[n_electrodes, n_time_pts, factor_levels, n_subs]);
    end
    
    %Figure out the effects we need to calculate
    [effects, effects_labels] = get_effects(factor_names);
        
    
    %% ~~~~~ CALCULATE ANOVA AND FDR CORRECTION ~~~~~
    
    test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                                 'F_obs', NaN(n_electrodes, n_time_pts), 'F_crit', NaN, 'df', NaN(1, 2)), ...
                                 length(effects), 1);
    for i = 1:length(effects)
        if VERBLEVEL
            fprintf('\nCalculating %s effect\n', effects_labels{i});
        end
        test_results(i) = calc_param_ANOVA(the_data, effects{i}+2, .05);
        switch method
            case 'bh'
                [h, crit_p, ~, adj_p] = fdr_bh(test_results(i).p, q, 'pdep', 'no');
            case 'by'
                [h, crit_p, ~, adj_p] = fdr_bh(test_results(i).p, q, 'dep', 'no');
            case 'bky'
                [h, crit_p] = fdr_bky(test_results(i).p, q, 'no');
                adj_p = NaN;
        end
        test_results(i).p = adj_p;
        test_results(i).h = h;
        if crit_p == 0
            test_results(i).F_crit = NaN;
        else
            test_results(i).F_crit = finv(1-crit_p, test_results(i).df(1), test_results(i).df(2));
        end
    end
   

    %% ~~~~~ ADD RESULTS STRUCT TO GND AND ASSIGN OTHER OUTPUT ~~~~~
    
    if (strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y'))
        use_time_pts = {use_time_pts};
    end
    %Create results struct
    results = struct('bins', bins, ...
                     'factors', {factor_names}, ...
                     'factor_levels', factor_levels, ...
                     'time_wind', time_wind, ...
                     'used_tpt_ids', use_time_pts, ...
                     'mean_wind', p.Results.mean_wind, ...
                     'include_chans', {{GND.chanlocs(electrodes).labels}}, ...
                     'used_chan_ids', electrodes, ...
                     'mult_comp_method', method, ...
                     'interaction_method', NaN, ...
                     'n_perm', NaN, ...
                     'desired_alphaORq', q, ...
                     'estimated_alpha', NaN, ...
                     'seed_state', NaN, ...
                     'null_test', [], ...
                     'adj_pval', [], ...
                     'F_obs', [], ...
                     'F_crit', [], ...
                     'df', [], ...              
                     'chan_hood', NaN, ...
                     'clust_info', NaN, ...
                     'fdr_rej', []);
    
    %Add statistical results
    assert(length(effects) == length(test_results));
    if length(effects) == 1
        results.null_test       = test_results.h;
        results.fdr_rej         = test_results.h;
        results.adj_pval        = test_results.p;
        results.F_obs           = test_results.F_obs;
        results.df              = test_results.df;
        results.F_crit          = test_results.F_crit;
    else
        for i = 1:length(effects)
            results.null_test.(effects_labels{i})       = test_results(i).h;
            results.fdr_rej.(effects_labels{i})      = test_results(i).h;
            results.adj_pval.(effects_labels{i})        = test_results(i).p;
            results.F_obs.(effects_labels{i})           = test_results(i).F_obs;
            results.df.(effects_labels{i})              = test_results(i).df;
            results.F_crit.(effects_labels{i})          = test_results(i).F_crit;
        end
    end
                 
    %Add results struct to GND
    if ~isfield(GND, 'F_tests') || isempty(GND.F_tests)
        GND.F_tests = results;
    else
        GND.F_tests(end+1) = results;
    end
    
    %Optional outputs
    if nargout > 2
        adj_pval   = results.adj_pval;
        F_obs      = results.F_obs;
        F_crit     = results.F_crit;
    end
 
    %% ~~~~~ REPORT RESULTS TO COMMAND WINDOW ~~~~~
    
    if VERBLEVEL
        fprintf('\n##### RESULTS #####\n');
        if length(effects) == 1
                fprintf('\n%s effect\n', effects_labels{1});
                if any(results.null_test(:))
                    fprintf('Critical F-value: %.3f\n', results.F_crit);
                    fprintf('That corresponds to a test-wise alpha level of %.3f\n', ...
                            fcdf(results.F_crit, results.df(1), results.df(2), 'upper'));
                    fprintf('Total number of significant differences: %d\n', sum(results.null_test(:)));
                    fprintf('Estimated upper bound on expected number of false discoveries: %.1f\n', sum(results.null_test(:))*q);
                    if strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y')
                        for t = 1:size(time_wind, 1)
                            fprintf('Significant electrodes for time window %d - %d: ', time_wind(t, 1), time_wind(t, 2));
                            fprintf('%s ', results.include_chans{results.null_test(:, t)});
                            fprintf('\n');
                        end
                    else
                        fprintf('Electrodes and time points with significant effects:\n');
                        for t = 1:length(results.used_tpt_ids)
                            if any(results.null_test(:, t))
                                fprintf('%d ms, electrode(s): ', GND.time_pts(results.used_tpt_ids(t)));
                                fprintf('%s ', results.include_chans{results.null_test(:, t)});
                                fprintf('\n');
                            end
                        end
                    end
                    if ~strcmpi(method, 'bky')
                        fprintf('All significant corrected p-values are between %f and %f\n', ...
                                max(results.adj_pval(results.adj_pval <= .05)), ... 
                                min(results.adj_pval(:)));
                    end
                else
                    fprintf('NO significant time points or electrodes.\n\n');
                end
        else
            for i = 1:length(effects)
                fprintf('\n%s effect\n', effects_labels{i});
                if any(results.null_test.(effects_labels{i})(:))
                    fprintf('Critical F-value: %.3f\n', results.F_crit.(effects_labels{i}));
                    fprintf('That corresponds to a test-wise alpha level of %.3f\n', ...
                            fcdf(results.F_crit.(effects_labels{i}), results.df.(effects_labels{i})(1), results.df.(effects_labels{i})(2), 'upper'));
                    fprintf('Total number of significant differences: %d\n', sum(results.null_test.(effects_labels{i})(:)));
                    fprintf('Estimated upper bound on expected number of false discoveries: %.1f\n', sum(results.null_test.(effects_labels{i})(:))*q);
                    if strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y')
                        for t = 1:size(time_wind, 1)
                            fprintf('Significant electrodes for time windonw %d - %d: ', time_wind(t, 1), time_wind(t, 2));
                            fprintf('%s ', results.include_chans{results.null_test.(effects_labels{i})(:, t)});
                            fprintf('\n');
                        end
                    else
                        fprintf('Electrodes and time points with significant effects:\n');
                        for t = 1:length(results.used_tpt_ids)
                            if any(results.null_test.(effects_labels{i})(:, t))
                                fprintf('%d ms, electrode(s): ', GND.time_pts(results.used_tpt_ids(t)));
                                fprintf('%s ', results.include_chans{results.null_test.(effects_labels{i})(:, t)});
                                fprintf('\n');
                            end
                        end
                    end
                    if ~strcmpi(method, 'bky')
                        fprintf('All significant corrected p-values are between %f and %f\n', ...
                                max(results.adj_pval.(effects_labels{i})(results.adj_pval.(effects_labels{i}) <= .05)), ... 
                                min(results.adj_pval.(effects_labels{i})(:)));
                    end
                else
                    fprintf('NO significant time points or electrodes.\n\n');
                end
            end
        end
    end
    
    %% ~~~~~ PLOT & SAVE RESULTS TO DISK ~~~~~
    
    %Plot results
    if ~strcmpi(p.Results.plot_raster, 'no') && ~strcmpi(p.Results.plot_raster, 'n')
        if length(effects_labels) == 1
            F_sig_raster(GND, length(GND.F_tests), 'use_color', 'rgb');
        else
            for i = 1:length(effects_labels)
                F_sig_raster(GND, length(GND.F_tests), 'effect', effects_labels{i}, 'use_color', 'rgb');
            end
        end
    end
    
    %Prompt user about saving GND
    if ~strcmpi(p.Results.save_GND, 'no') && ~strcmpi(p.Results.save_GND, 'n')
        GND = save_matmk(GND);
    end
    
    %Output to spreadsheet if requested
    if p.Results.output_file
        if VERBLEVEL
            fprintf('\nWriting results to %s . . . ', p.Results.output_file)
        end
        Ftest2xls(GND, length(GND.F_tests), p.Results.output_file);
        if VERBLEVEL
            fprintf('DONE\n\n')
        end
    end

    
end
