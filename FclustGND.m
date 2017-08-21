%Function to conduct an cluster based permutation test for one-way and
%factorial within-subjects ANOVA
%
%EXAMPLE USAGE
% GND = FclustGND(GND, 'bins', 1:6, 'factor_names', {'probability', 'emotion'}, ...
%                 'factor_levels', [3, 2], 'time_wind', [500, 800], ...
%                 'include_chans', {'Fz', 'Cz', 'Pz'}, 'n_perm', 1e4, ...
%                 'chan_hood', 0.61, 'thresh_p', 0.05, 'alpha', 0.05);
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
% factor_levels  - number of factors in each level in fastest to slowest
%                  moving order within the bins provided
%
%OPTIONAL INPUTS
% chan_hood      - A scalar or a 2D symmetric binary matrix that indicates
%                  which channels are considered neighbors of other 
%                  channels. E.g., if chan_hood(2,10)=1, then Channel 2 
%                  and Channel 10 are nieghbors. You can produce a 
%                  chan_hood matrix using the function spatial_neighbors.m. 
%                  If a scalar is provided, then all electrodes within that 
%                  distance of a particular electrode are considered 
%                  neighbors. Note, EEGLAB's electrode coordinates assume
%                  the head has a radius of 1. See the help documentation 
%                  of the function spatial_neighbors to see how you could
%                  convert this distance threshold to centimeters. 
%                  {default: 0.61}
% head_radius    - The radius of the head in whatever units the Cartesian
%                  coordinates in GND.chanlocs are in. This is used to
%                  convert scalar values of chan_hood into centimeters.
%                  {default: []}
% thresh_p       - The test-wise p-value threshold for cluster inclusion. If
%                  a channel/time-point has a F-value that corresponds to an
%                  uncorrected p-value greater than thresh_p, it is assigned
%                  a p-value of 1 and not considered for clustering.
%                  {default: 0.05}
% time_wind      - 2D matrix of time values specifying the beginning 
%                  and end of the time windows in ms (e.g., 
%                  [500, 800]). Every single time point in 
%                  the time window will be individually tested (i.e.,
%                  maximal temporal resolution). Note, boundaries of time 
%                  window(s) may not exactly correspond to desired time 
%                  window boundaries because of temporal digitization (i.e., 
%                  you only have samples every so many ms). 
%                  {default: 0 ms to the end of the epoch}
% mean_wind      - ['yes' or 'no'] If 'yes', the permutation test will be
%                  performed on the mean amplitude within the time window 
%                  specified by time_wind.  If 'no', every single time point within
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
% n_perm         - number of permutations {default: 10,000}
% alpha          - A number between 0 and 1 specifying the family-wise 
%                  alpha level of the test. {default: 0.05}
% plot_raster    - ['yes' or 'no'] If 'yes', a two-dimensional (time x channel)
%                  binary "raster" diagram is created to illustrate the
%                  results of the permutation tests. This figure can be reproduced 
%                  with the function F_sig_raster.m. {default: 'yes'}
% save_GND       - save GND to disk, 'yes' or 'no' {default: user will be
%                  prompted}
% output_file    - Name of .xlsx file to output results. {default: no output}
% reproduce_test - [integer] The number of the permutation test stored in
%                  the GND variable to reproduce.  For example, if 
%                  'reproduce_test' equals 2, the second F-test 
%                  stored in the GND variable (i.e., GND.F_tests{2}) will 
%                  be reproduced. Reproduction is accomplished by setting
%                  the random number generator used in the permutation test 
%                  to the same initial state it was in when the permutation 
%                  test was conducted. Obviously other options/inputs must
%                  also be the same to truly reproduce the test.
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
% prm_pval      - The adj_pval field from the results struct. For a one-way
%                 ANOVA, this is an electrodes x time points array of
%                 p-values; for a multi-factor ANOVA, it is a struct with
%                 multiple such arrays.
% F_obs         - The F_obs field from the results struct. For a one-way
%                 ANOVA, this is an electrodes x time points array of
%                 Fs; for a multi-factor ANOVA, it is a struct with
%                 multiple such arrays.
% clust_info    - The clust_info field from the results struct. This is a
%                 struct that has sub-fields with information about each
%                 observed cluster and it's statistical significance
%
%
%DESCRIPTION
%
%
%See the FMUT documentation for more information:
%https://github.com/ericcfields/FMUT/wiki
%
%
%AUTHOR: Eric Fields
%VERSION DATE: 21 August 2017
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
% 11/28/16     - Moved calculations to sub function calc_F_clust_mass
% 12/8/16      - Added abilty to set global variable VERBLEVEL. Added
%                'mean_wind' option.
% 2/14/17      - Added cluster ID to spreadsheet output
% 3/31/17      - Moved spreadsheet output function to standalone
% 4/7/17       - Fixed error in checking interaction method
% 4/17/17      - Fixed inconsistency in used_tpt_id field with t-tests
%                results; changed desired_alpha to desired_alphaORq;
%                added plot_raster functionality; added estimated_alpha
%                to results
% 5/9/17       - Added informative error messages for incorrect
%                electrode names
% 6/2/17       - Electrode order in output when using include_chans now
%                matches MUT functions
% 6/12/17      - Added verblevel related reports
% 6/14/17      - Updated error for incorrectly supplied electrode name;
%                fixed command window output for one-way ANOVA
% 6/20/17      - Command window output for mean window analyses
% 6/21/17      - time_wind field of results struct is now accurate;
%                changed used_tpt_ids field to cell array for mean window
%                analyses
% 6/27/17      - More information in command window output
% 7/13/17      - int_method input eliminated; fixed spacing in command
%                window output
% 7/14/17      - Move command window output to separate function
% 7/15/17      - use_groups and group_n to F_tests

function [GND, results, prm_pval, F_obs, clust_info] = FclustGND(GND_or_fname, varargin)

    warning('You are using a beta version of FclustGND. Some bugs may remain and results should be interpreted with caution.');

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
    p.addParameter('n_perm',        1e4,      @(x) (isnumeric(x) && isscalar(x)));
    p.addParameter('save_GND',      'prompt', @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))) || islogical(x));
    p.addParameter('output_file',   false,    @(x) (ischar(x) || islogical(x)));
    p.addParameter('alpha',         0.05,     @(x) (isnumeric(x) && isscalar(x) && x<=1 && x>=0));
    p.addParameter('reproduce_test',false,    @(x) (isnumeric(x) && isscalar(x)));
    p.addParameter('mean_wind',     'no',     @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))));
    p.addParameter('verblevel',     [],       @(x) (isnumeric(x) && isscalar(x) && x>=0 && x<=3))
    p.addParameter('plot_raster',   'yes',    @(x) (any(strcmpi(x, {'yes', 'no', 'n', 'y'}))));
    p.addParameter('chan_hood',     0.61,     @(x) (isnumeric(x) && (isscalar(x) || (ismatrix(x) && size(x,1)==size(x,2)))));
    p.addParameter('head_radius',   [],       @(x) (isnumeric(x) && isscalar(x)));
    p.addParameter('thresh_p',      0.05,     @(x) (isnumeric(x) && isscalar(x) && x>=0 && x<=1));
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
        error('The GND variable provided does not seem to be a valid GND struct or filepath to a GND struct.');
    end
    
    %Assign some variables for easier reference
    bins          = p.Results.bins;
    factor_names  = p.Results.factor_names;
    factor_levels = p.Results.factor_levels;
    time_wind     = p.Results.time_wind;
    n_perm        = p.Results.n_perm;
    alpha         = p.Results.alpha;
    chan_hood     = p.Results.chan_hood;
    
    %Check for required name-value inputs
    if isempty(bins)
        error('''bins'' is a required input. See >>help FclustGND.');
    end
    if isempty(factor_names)
        error('''factor_names'' is a required input. See >>help FclustGND.');
    end
    if isempty(factor_levels)
        error('''factor_levels'' is a required input. See >>help FclustGND.');
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
            error([sprintf('The following channels appear in ''exclude_chans'' but do not appear in GND.chanlocs.labels:\n') ... 
                   sprintf('%s ', missing_channels{:})])
        else
            electrodes = find(~ismember(chan_labels, p.Results.exclude_chans));
        end
    else
        electrodes = 1:length(GND.chanlocs);
    end
    
    %MUT reatures not implemented here              
    if ~isempty(p.Results.time_block_dur)
        error('The ''time_block_dur'' option is not implemented for FclustGND. You''ll need to divide the time windows manually.');
    end              
    if ~isempty(p.Results.plot_gui)
        watchit('''plot_gui'' is not implemented for FclustGND.');
    end
    if ~isempty(p.Results.plot_mn_topo)
        watchit('''plot_mn_topo'' is not implemented for FclustGND.');
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
        error('The number of factors does not match in the ''factor_names'' and ''factor_levels'' inputs.');
    end
    if length(factor_levels) > 2
        warning('This function has not been tested extensively with designs with more than two factors. Proceed with caution!');
    end
    if sum(factor_levels>2) > 3
        error('Designs with more than three factors with more than two levels are not supported by FclustGND.')
    end
    if ~isequal(size(time_wind), [1, 2])
        error('''time_wind'' input must indicate a single time window with one starting point and one stopping point (e.g., [500, 800])');
    end
    if prod(factor_levels) ~= length(bins)
        error('Number of bins does not match design.')
    end
    if alpha <= .01 && n_perm < 5000,
        watchit(sprintf('You are probably using too few permutations for an alpha level of %f.',alpha));
    elseif alpha <=.05 && n_perm < 1000,
        watchit(sprintf('You are probably using too few permutations for an alpha level of %f.',alpha));
    end
    if ~all(all(GND.indiv_bin_ct(:, bins)))
        watchit(sprintf('Some subjects appear to be missing data from bins used in this test!\nSee: GND.indiv_bins_ct.'));
    end
    if p.Results.reproduce_test
        if ~isfield(GND, 'F_tests')
            error('You tried to reproduce test %d, but there are not results in GND.F_tests.', p.Results.reproduce_test);
        elseif p.Results.reproduce_test > length(GND.F_tests)
            error('You tried to reproduce test %d, but there are only %d tests in GND.F_tests.', p.Results.reproduce_test, length(GND.F_tests));
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
        seed_state = GND.F_tests(p.Results.reproduce_test).seed_state;
        defaultStream.State = seed_state;
    else
        seed_state = defaultStream.State;
    end
    
    %Some useful numbers
    n_subs       = size(GND.indiv_erps, 4);
    n_electrodes = length(electrodes); 

    %Find time points or mean windows to use and extract the data for
    %analysis
    if ~strcmpi(p.Results.mean_wind, 'yes') && ~strcmpi(p.Results.mean_wind, 'y')
        [~, start_sample] = min(abs( GND.time_pts - time_wind(1) ));
        [~, end_sample  ] = min(abs( GND.time_pts - time_wind(2) ));
        time_wind(1) = GND.time_pts(start_sample);
        time_wind(2) = GND.time_pts(end_sample);
        use_time_pts = start_sample:end_sample;
        n_time_pts = length(use_time_pts);
        the_data = GND.indiv_erps(electrodes, use_time_pts, bins, :);
        if VERBLEVEL
            fprintf('\nConducting cluster mass permutation test from %d ms to %d ms.\n', GND.time_pts(start_sample), GND.time_pts(end_sample));
        end
    else
        n_time_pts = 1;
        [~, start_sample] = min(abs( GND.time_pts - time_wind(1) ));
        [~, end_sample  ] = min(abs( GND.time_pts - time_wind(2) ));
        time_wind(1) = GND.time_pts(start_sample);
        time_wind(2) = GND.time_pts(end_sample);
        the_data = mean(GND.indiv_erps(electrodes, start_sample:end_sample, bins, :), 2);
        if VERBLEVEL
            fprintf('\nConducting cluster mass permutation test in the mean time window from %d ms to %d ms.\n', GND.time_pts(start_sample), GND.time_pts(end_sample));
        end
    end
    use_time_pts = start_sample:end_sample;
    
    if VERBLEVEL
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
    
    %Get chan_hood matrix if input was distance scalar
    if isscalar(chan_hood)
        if VERBLEVEL; fprintf('\n'); end;
        chan_hood = spatial_neighbors(GND.chanlocs(electrodes), chan_hood, p.Results.head_radius);
    end
        
    
    %% ~~~~~ RUN PERMUTATION ANOVAS ~~~~~
    
    test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                                 'F_obs', NaN(n_electrodes, n_time_pts), 'df', NaN(1, 2), 'clust_info', struct, ...
                                 'estimated_alpha', NaN, 'exact_test', NaN), ...
                                 length(effects), 1);
    for i = 1:length(effects)
        if VERBLEVEL
            fprintf('\nCalculating %s effect\n', effects_labels{i});
        end
        test_results(i) = calc_Fclust(the_data, [], effects{i}+2, n_perm, alpha, chan_hood, p.Results.thresh_p);
    end
   

    %% ~~~~~ ADD RESULTS STRUCT TO GND AND ASSIGN OTHER OUTPUT ~~~~~
    
    if (strcmpi(p.Results.mean_wind, 'yes') || strcmpi(p.Results.mean_wind, 'y'))
        use_time_pts = {{use_time_pts}};
    end
    %Create results struct
    results = struct('bins', bins, ...
                     'use_groups', NaN, ...
                     'group_n', n_subs, ...
                     'factors', {factor_names}, ...
                     'factor_levels', factor_levels, ...
                     'time_wind', time_wind, ...
                     'used_tpt_ids', use_time_pts, ...
                     'mean_wind', p.Results.mean_wind, ...
                     'include_chans', {{GND.chanlocs(electrodes).labels}}, ...
                     'used_chan_ids', electrodes, ...
                     'mult_comp_method', 'cluster mass perm test', ...
                     'n_perm', n_perm, ...
                     'desired_alphaORq', alpha, ...
                     'estimated_alpha', [], ...
                     'seed_state', seed_state, ...
                     'exact_test', [], ...
                     'null_test', [], ...
                     'adj_pval', [], ...
                     'F_obs', [], ...
                     'F_crit', NaN, ...
                     'df', [], ...              
                     'chan_hood', chan_hood, ...
                     'clust_info', [], ...
                     'fdr_rej', NaN);
    
    %Add statistical results
    assert(length(effects) == length(test_results));
    if length(effects) == 1
        results.null_test       = test_results.h;
        results.adj_pval        = test_results.p;
        results.F_obs           = test_results.F_obs;
        results.df              = test_results.df;
        results.clust_info      = test_results.clust_info;
        results.estimated_alpha = test_results.estimated_alpha;
        results.exact_test = test_results.exact_test;
    else
        for i = 1:length(effects)
            results.null_test.(effects_labels{i})       = test_results(i).h;
            results.adj_pval.(effects_labels{i})        = test_results(i).p;
            results.F_obs.(effects_labels{i})           = test_results(i).F_obs;
            results.df.(effects_labels{i})              = test_results(i).df;
            results.clust_info.(effects_labels{i})      = test_results(i).clust_info;
            results.estimated_alpha.(effects_labels{i}) = test_results(i).estimated_alpha;
            results.exact_test.(effects_labels{i}) = test_results(i).exact_test;
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
        prm_pval   = results.adj_pval;
        F_obs      = results.F_obs;
        clust_info = results.clust_info;
    end
 
    
    %% ~~~~~ OUTPUT RESULTS ~~~~~
    
    %Output results to command window
    if VERBLEVEL
        report_results(GND, length(GND.F_tests))
    end
    
    %Plot results
    if ~strcmpi(p.Results.plot_raster, 'no') && ~strcmpi(p.Results.plot_raster, 'n')
        if VERBLEVEL
            fprintf('Generating raster plot:\n');
        end
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
            fprintf('DONE.\n')
        end
    end
    
end
