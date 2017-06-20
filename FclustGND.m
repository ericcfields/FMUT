%Function to conduct an cluster based permutation test for one-way and
%factorial ANOVA
%
%
%EXAMPLE USAGE
% GND = FclustGND(GND, 'bins', 1:6, 'factor_names', {'probability', 'emotion'}, ...
%                 'factor_levels', [3, 2], 'time_wind', [500, 800], ...
%                 'include_chans', {'Fz', 'Cz', 'Pz'}, 'n_perm', 1e4, ...
%                 'chan_hood', 0.61, 'thresh_p', 0.05, 'alpha', 0.05, 
%                 'int_method', 'exact');
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
% int_method     - A string that should be either 'exact' or 'approximate'.
%                  If 'exact', the method of restricted permutations will
%                  be used to conduct a test that controls the Type I error
%                  rate at alpha (assuming enough permutations). 
%                  If 'approximate', the method of permutation of residuals 
%                  will be used to conduct a test with Type I error rate 
%                  asymptotic to alpha as noise decreases and/or number of 
%                  subjects increases. 
%                  See explanation and references below. {default:
%                  'exact' where possible, otherwise 'approximate'}
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
% n_perm         - number of permutations {default: 10,000}
% alpha          - A number between 0 and 1 specifying the family-wise 
%                  alpha level of the test. {default: 0.05}
% plot_raster    - ['yes' or 'no'] If 'yes', a two-dimensional (time x channel)
%                  binary "raster" diagram is created to illustrate the
%                  results of the permutation tests. This figure can be reproduced 
%                  with the function F_sig_raster.m. {default: 'yes'}
% save_GND       - save GND to disk, 'yes' or 'no' {default: user will be
%                  prompted}
% output_file    - Name of .xlsx file to output results. Currently 
%                  only works on Windows. {default: no output}
% reproduce_test - [integer] The number of the permutation test stored in
%                  the GND variable to reproduce.  For example, if 
%                  'reproduce_test' equals 2, the second F-test 
%                  stored in the GND variable (i.e., GND.F_tests{2}) will 
%                  be reproduced. Reproduction is accomplished by setting
%                  the random number generator used in the permutation test 
%                  to the same initial state it was in when the permutation 
%                  test was conducted.
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
%Main effects are calculated by permuting within each condition of the other
%factor(s). For a factor with two levels, this is similar to the cluster t-test 
%of the clustGND function. However, because F-values are always positive and are
%in squared units, the clusters and their relative sizes will not necessarily be 
%the same.
%
%For interaction effects where more than one factor has more than two levels, 
%it is not possible to conduct a test that controls the Type I error exactly 
%at a specified level. For such designs, this function uses the permutation 
%of residuals method first described by Still & White (1981) and Freedman & Lane (1983). 
%The Type I error rate of this test is asymptotic to the nominal alpha as 
%sample size and/or signal to noise ratio increase.
%For designs where an exact test is possible, this function can use a
%restricted permutation method to conduct an exact test. Optionally you
%can also use the approximate method for such cases. See below for references.
%NOTE: It is not recommended to use the approximate method for 2x2 designs.
%
%The statistic used is the cluster mass: i.e., the sum of all the F-values
%included in a given cluster.
%
%REFERENCES FOR PERMUTATION FACTORIAL ANOVA AND GLM
%Anderson, M. J. (2001). Permutation tests for univariate or multivariate analysis of variance and regression. Canadian Journal of Fisheries and Aquatic Sciences, 58(3), 626-639.
%Anderson, M., & Braak, C. T. (2003). Permutation tests for multi-factorial analysis of variance. Journal of statistical computation and simulation, 73(2), 85-113.
%Wheldon, M. C., Anderson, M. J., & Johnson, B. W. (2007). Identifying treatment effects in multi-channel measurements in electroencephalographic studies: Multivariate permutation tests and multiple comparisons. Australian & New Zealand Journal of Statistics, 49(4), 397-413. 
%Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E. (2014). Permutation inference for the general linear model. NeuroImage, 92, 381-397.
%
%REFERENCES FOR CLUSTER MASS METHOD
%Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., Taylor, E., & Brammer, M. J. (1999). Global, voxel, and cluster tests, by theory and permutation, for a difference between two groups of structural MR images of the brain. IEEE Transactions on Medical Imaging, 18(1), 32-42. 
%Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.
%Groppe, D. M., Urbach, T. P., & Kutas, M. (2011). Mass univariate analysis of event-related brain potentials/fields I: A critical tutorial review. Psychophysiology, 48(12), 1711-1725.
%
%
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%VERSION DATE: 20 June 2017
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
% 1.0.0        - Initial version modified from Fmax GND 2.5.2
% 1.1.0        - Added output to spreadsheet
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

function [GND, results, prm_pval, F_obs, clust_info] = FclustGND(GND_or_fname, varargin)

    warning('You are using a beta version of this function. It needs further testing and should NOT be considered error free.');

    %% ~~~~~PARSE INPUT~~~~~
    
    global VERBLEVEL
    
    %Assign GND
    if ischar(GND_or_fname)
        load(GND_or_fname, '-mat');
    elseif isstruct(GND_or_fname)
        GND = GND_or_fname;
    else
        error('The GND variable provided does not seem to be a valid GND struct or filepath to a GND struct');
    end
    
    %Assign name-value pair inputs to variables
    for i = 1:length(varargin)
        input = varargin{i};
        if ischar(input)
            switch input
                case 'bins'
                    bins = varargin{i+1};
                case 'factor_names'
                    factor_names = varargin{i+1};
                case 'factor_levels'
                    factor_levels = varargin{i+1};
                case 'time_wind'
                    time_wind = varargin{i+1};
                case 'include_chans'
                    chan_labels = {GND.chanlocs.labels};
                    include_chans = varargin{i+1};
                    electrodes = NaN(1, length(include_chans));
                    for c = 1:length(include_chans)
                        if find(strcmp(include_chans(c), chan_labels))
                            electrodes(c) = find(strcmp(include_chans(c), chan_labels));
                        else
                            error('Electrode ''%s'' does not exist.', include_chans{c});
                        end
                    end
                case 'exclude_chans'
                    chan_labels = {GND.chanlocs.labels};
                    exclude_chans = varargin{i+1};
                    electrodes = find(~ismember(chan_labels, exclude_chans));
                case 'n_perm'
                    n_perm = varargin{i+1};
                case 'save_GND'
                    save_GND = varargin{i+1};
                case 'output_file'
                    output_file = varargin{i+1};
                case 'alpha'
                    alpha = varargin{i+1};
                case 'reproduce_test'
                    rep_test = varargin{i+1};
                case 'int_method'
                    int_method = varargin{i+1};
                case 'chan_hood'
                    chan_hood = varargin{i+1};
                case 'thresh_p'
                    thresh_p = varargin{i+1};
                case 'head_radius'
                    head_radius = varargin{i+1};
                case 'mean_wind'
                    mean_wind = varargin{i+1};
                case 'verblevel'
                    VERBLEVEL = varargin{i+1};
                case 'plot_raster'
                    plot_raster = varargin{i+1};
                case 'time_block_dur'
                    error('The ''time_block_dur'' option is not implemented for FmaxGND. You''ll need to divide time windows manually');
                case 'plot_gui'
                    watchit('''plot_gui'' is not implemented for FclustGND.');
                case 'plot_mn_topo'
                    watchit('''plot_mn_topo'' is not implemented for FclustGND.');
            end
        end
    end
    
    %Set defaults for missing arguments
    if isempty(VERBLEVEL)
        VERBLEVEL = 2;
    end
    if ~exist('time_wind', 'var')
        time_wind = [0, GND.time_pts(end)];
    end
    if ~exist('electrodes', 'var')
        electrodes = 1:length(GND.chanlocs);
    end
    if ~exist('n_perm', 'var')
        n_perm = 1e4;
    end
    if ~exist('save_GND', 'var')
        save_GND = 'prompt';
    end
    if ~exist('output_file', 'var')
        output_file = false;
    end
    if ~exist('alpha', 'var')
        alpha = 0.05;
    end
    if ~exist('rep_test', 'var')
        rep_test = false;
    end
    if ~exist('int_method', 'var')
        if length(factor_levels) == 1
            int_method = 'none';
        elseif sum(factor_levels > 2) <= 1
            int_method = 'exact';
            fprintf('\n\nUsing restricted permutations to conduct an exact test of the interaction effect.\nSee help FclustGND for more information.\n\n')
        else
            int_method = 'approx';
            fprintf('\n\nAn exact test of the interaction is not possible for this design.\nUsing permutation of residuals method to conduct an approximate test.\nSee help FclustGND for more information.\n\n')
        end
    end
    if ~exist('chan_hood', 'var')
        chan_hood = 0.61;
    end
    if ~exist('thresh_p', 'var')
        thresh_p = 0.05;
    end
    if ~exist('head_radius', 'var')
        head_radius = [];
    end
    if ~exist('mean_wind', 'var')
        mean_wind = 'no';
    end
    if ~exist('plot_raster', 'var')
        plot_raster = 'yes';
    end
    
    %Standardize formatting
    if ischar(factor_names)
        factor_names = {factor_names};
    end
    if strcmpi(int_method, 'approximate')
        int_method = 'approx';
    end
    time_wind = sort(time_wind);
    
    %Check for errors and problems in input
    if ~exist('bins', 'var')
        error('''bins'' is a required input. See help FclustGND');
    end
    if ~exist('factor_names', 'var')
        error('''factor_names'' is a required input. See help FclustGND');
    end
    if ~exist('factor_levels', 'var')
        error('''factor_levels'' is a required input. See help FclustGND');
    end
    if exist('include_chans', 'var') && exist('exclude_chans', 'var')
        error('You cannot use BOTH ''include_chans'' and ''exclude_chans'' options.');
    end
    if exist('exclude_chans', 'var') && ~all(ismember(exclude_chans, chan_labels))
        missing_channels = exclude_chans(~ismember(exclude_chans, chan_labels));
        error([sprintf('The following channels appear in ''exclude_chans'' but do not appear in GND.chanlocs.labels:\n') ... 
               sprintf('%s ', missing_channels{:})])
    end
    if length(factor_names) ~= length(factor_levels)
        error('The number of factors does not match in the ''factor_names'' and ''factor_levels'' inputs');
    end
    if length(factor_levels) > 2
        warning('This function has not been tested extensively with designs with more than two factors. Proceed with caution!');
    end
    if length(factor_levels) > 3 && strcmpi(int_method, 'approx')
        error('FclustGND does not currently support designs with more than three factors using the approximate method of calculating interaction effects.');
    end
    if strcmpi(int_method, 'exact') && sum(factor_levels > 2) > 1
        error('An exact test of the interaction is not possible if more than one factor has more than two levels. See help FclustGND for more information.');
    end
    if isequal(factor_levels, [2, 2]) && strcmpi(int_method, 'approx')
        button = questdlg('WARNING: The type I error rate is not well-controlled by the approximate method of calculating the interaction for a 2x2 design. Are you sure you want to proceed?', 'WARNING');
        if ~strcmp(button, 'Yes')
            return;
        end
    end
    if ~strcmpi(int_method, 'approx') && ~strcmpi(int_method, 'exact') && ~strcmpi(int_method, 'none')
        error('The int_method argument must be either ''approximate'' or ''exact''.')
    end
    if prod(factor_levels) ~= length(bins)
        error('Number of bins does not match design.')
    end
    if ~isequal(size(time_wind), [1, 2])
        error('''time_wind'' input must indicate a single time window with one starting point and one stopping point (e.g., [500, 800])');
    end
    if alpha <= .01 && n_perm < 5000,
        watchit(sprintf('You are probably using too few permutations for an alpha level of %f.',alpha));
    elseif alpha <=.05 && n_perm < 1000,
        watchit(sprintf('You are probably using too few permutations for an alpha level of %f.',alpha));
    end
    if ~all(all(GND.indiv_bin_ct(:, bins)))
        watchit(sprintf('Some subjects appear to be missing data from bins used in this test!\nSee: GND.indiv_bins_ct'));
    end

    
    %% ~~~~~ SET-UP ~~~~~

    %Get or set random # generator state
    if verLessThan('matlab','8.1')
        defaultStream=RandStream.getDefaultStream; 
    else
        defaultStream=RandStream.getGlobalStream;
    end
    if rep_test
        seed_state = GND.F_tests{rep_test}.seed_state;
        defaultStream.State = seed_state;
    else
        seed_state = defaultStream.State;
    end
    
    %Some useful numbers
    n_subs       = size(GND.indiv_erps, 4);
    n_electrodes = length(electrodes); 

    %Find time points or mean windows to use and extract the data for
    %analysis
    if ~strcmpi(mean_wind, 'yes') && ~strcmpi(mean_wind, 'y')
        [~, start_sample] = min(abs( GND.time_pts - time_wind(1) ));
        [~, end_sample  ] = min(abs( GND.time_pts - time_wind(2) ));
        use_time_pts = start_sample:end_sample;
        n_time_pts = length(use_time_pts);
        the_data = GND.indiv_erps(electrodes, use_time_pts, bins, :);
        if VERBLEVEL
            fprintf('\nConducting cluster mass permutation test from %d ms to %d ms\n', GND.time_pts(start_sample), GND.time_pts(end_sample));
        end
    else
        n_time_pts = 1;
        [~, start_sample] = min(abs( GND.time_pts - time_wind(1) ));
        [~, end_sample  ] = min(abs( GND.time_pts - time_wind(2) ));
        the_data = mean(GND.indiv_erps(electrodes, start_sample:end_sample, bins, :), 2);
        if VERBLEVEL
            fprintf('\nConducting cluster mass permutation test in the mean time window from %d ms to %d ms\n', GND.time_pts(start_sample), GND.time_pts(end_sample));
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
        chan_hood = spatial_neighbors(GND.chanlocs(electrodes), chan_hood, head_radius);
    end
        
    
    %% ~~~~~ RUN PERMUTATION ANOVAS ~~~~~
    
    test_results = repmat(struct('h', NaN(n_electrodes, n_time_pts), 'p', NaN(n_electrodes, n_time_pts), ... 
                                 'F_obs', NaN(n_electrodes, n_time_pts), 'df', NaN(1, 2), 'clust_info', struct, ...
                                 'estimated_alpha', NaN), ...
                                 length(effects), 1);
    for i = 1:length(effects)
        if VERBLEVEL
            fprintf('\nCalculating %s effect\n', effects_labels{i});
        end
        test_results(i) = calc_F_clust_mass(the_data, effects{i}+2, n_perm, alpha, int_method, chan_hood, thresh_p);
    end
   

    %% ~~~~~ ADD RESULTS STRUCT TO GND AND ASSIGN OTHER OUTPUT ~~~~~
    
    %Create results struct
    results = struct('bins', bins, ...
                     'factors', {factor_names}, ...
                     'factor_levels', factor_levels, ...
                     'time_wind', time_wind, ...
                     'used_tpt_ids', use_time_pts, ...
                     'mean_wind', mean_wind, ...
                     'include_chans', {{GND.chanlocs(electrodes).labels}}, ...
                     'used_chan_ids', electrodes, ...
                     'mult_comp_method', 'cluster mass perm test', ...
                     'interaction_method', int_method, ...
                     'n_perm', n_perm, ...
                     'desired_alphaORq', alpha, ...
                     'estimated_alpha', [], ...
                     'seed_state', seed_state, ...
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
    else
        for i = 1:length(effects)
            results.null_test.(effects_labels{i})       = test_results(i).h;
            results.adj_pval.(effects_labels{i})        = test_results(i).p;
            results.F_obs.(effects_labels{i})           = test_results(i).F_obs;
            results.df.(effects_labels{i})              = test_results(i).df;
            results.clust_info.(effects_labels{i})      = test_results(i).clust_info;
        end
        results.estimated_alpha = test_results(1).estimated_alpha;
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
    
    %% ~~~~~ REPORT RESULTS TO COMMAND WINDOW ~~~~~
    
    if VERBLEVEL
        fprintf('\n##### RESULTS #####\n');
        if length(effects) == 1
                fprintf('\n%s effect\n', effects_labels{1});
                fprintf('# of clusters found: %d\n', length(results.clust_info.null_test));
                fprintf('# of significant clusters found: %d\n', sum(results.clust_info.null_test));
                if sum(results.clust_info.null_test)
                    fprintf('Significant cluster p-values range from %f to %f\n\n', ...
                            max(results.clust_info.pval(results.clust_info.null_test)), ...
                            min(results.clust_info.pval(results.clust_info.null_test)));
                else
                    fprintf('All p-values >= %f\n\n', min(results.clust_info.pval));
                end
        else
            for i = 1:length(effects)
                fprintf('\n%s effect\n', effects_labels{i});
                fprintf('# of clusters found: %d\n', length(results.clust_info.(effects_labels{i}).null_test));
                fprintf('# of significant clusters found: %d\n', sum(results.clust_info.(effects_labels{i}).null_test));
                if sum(results.clust_info.(effects_labels{i}).null_test)
                    fprintf('Significant cluster p-values range from %f to %f\n\n', ...
                            max(results.clust_info.(effects_labels{i}).pval(results.clust_info.(effects_labels{i}).null_test)), ...
                            min(results.clust_info.(effects_labels{i}).pval(results.clust_info.(effects_labels{i}).null_test)));
                else
                    fprintf('All p-values >= %f\n\n', min(results.clust_info.(effects_labels{i}).pval));
                end
            end
        end
    end
 
    
    %% ~~~~~ PLOT & SAVE RESULTS TO DISK ~~~~~
    
    %Plot results
    if ~strcmpi(plot_raster, 'no') && ~strcmpi(plot_raster, 'n')
        if length(effects_labels) == 1
            F_sig_raster(GND, length(GND.F_tests), 'use_color', 'rgb');
        else
            for i = 1:length(effects_labels)
                F_sig_raster(GND, length(GND.F_tests), 'effect', effects_labels{i}, 'use_color', 'rgb');
            end
        end
    end
    
    %Prompt user about saving GND
    if ~strcmpi(save_GND, 'no')
        GND = save_matmk(GND);
    end
    
    %Output to spreadsheet if requested
    if output_file
        if VERBLEVEL
            fprintf('\nWriting results to spreadsheet . . . \n')
        end
        Ftest2xls(GND, length(GND.F_tests), output_file);
        if VERBLEVEL
            fprintf('DONE\n')
        end
    end
    
end
