function [GND, results] = FmaxGND(GND_or_fname, varargin)
%Function to conduct an Fmax permutation test for a two-way factorial
%design.
%
%The main effects permute the data within each condition of the other
%factor(s). For a factor with two levels, this is equivalent to the tmaxGND 
%function (but will of course give results as an F-test).
%
%For interaction effects where more than one factor has more than two levels, 
%it is not possible to conduct a test that controls the Type I error exactly 
%at a specified level. For such designs, this function uses the permutation 
%of residuals method first described by Still & White (1981) and Freedman & Lane (1983) 
%The Type I error rate of this test is asymptotic to the nominal alpha as 
%ample size and/or signal to noise ratio increase.
%For designs where an exact test is possible, this function can use a
%restricted permutation method to conduct an exact test. Optionally you
%can also use the approximate method for such tests (which may have more
%power). See below for references.
%NOTE: It is not recommended to use the approximate method for 2x2 designs.
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
% time_wind      - 2D matrix of pairs of time values specifying the beginning 
%                  and end of the time windows in ms (e.g., 
%                  [500, 800]). Every single time point in 
%                  the time window will be individually tested (i.e.,
%                  maximal temporal resolution). Note, boundaries of time 
%                  window(s) may not exactly correspond to desired time 
%                  window boundaries because of temporal digitization (i.e., 
%                  you only have samples every so many ms). {default: 0 ms 
%                  to the end of the epoch}
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
% save_GND       - save GND to disk, 'yes' or 'no' {default: user will be
%                  prompted}
% output_file    - Name of .xlsx file to output results. Currently only works
%                  with MATLAB 2013b and later, on Windows, for 2x2 desings.          
%                  {default: no output}
% reproduce_test - [integer] The number of the permutation test stored in
%                  the GND variable to reproduce.  For example, if 
%                  'reproduce_test' equals 2, the second F-test 
%                  stored in the GND variable (i.e., GND.F_tests{2}) will 
%                  be reproduced. Reproduction is accomplished by setting
%                  the random number generator used in the permutation test 
%                  to the same initial state it was in when the permutation 
%                  test was conducted.
% int_method     - A string that should be either 'exact' or 'approximate'.
%                  If 'exact', the method of restricted permutations will
%                  be used to conduct a test that controls the Type I error
%                  rate at alpha (assuming enough permutations). 
%                  If 'approximate', the method of permutation of residuals 
%                  will be used to conduct a test with Type I error rate 
%                  asymptotic to alpha as noise decreases and/or number of 
%                  subjects increases. 
%                  See explanation above and references below. {default:
%                  'exact' where possible, otherwise 'approximate'}
%
%OUTPUT
% GND           - GND struct, with results added in the F_tests field.
% results       - The same struct added to the F_tests field, but assigned
%                 to it's own variable; this might make it easier to do
%                 further operations.
%
%EXAMPLE USAGE
% GND = Fmax_factorial_GND(GND, 'bins', 1:6, 'factor_names', {'probability', 'emotion'}, ...
%                          'factor_levels', [3, 2], 'time_wind', [500, 800], ...
%                          'include_chans', {'Fz', 'Cz', 'Pz'}, 'n_perm', 1e4, ...
%                          'alpha', 0.05, 'int_method', 'exact');
%
%
%References for permutation factorial ANOVA and GLM:
%Still, A. W., & White, A. P. (1981). The approximate randomization test as an alternative to the F test in analysis of variance. British Journal of Mathematical and Statistical Psychology, 34(2), 243-252.
%Freedman, D., & Lane, D. (1983). A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics, 1(4), 292-298.
%Anderson, M. J. (2001). Permutation tests for univariate or multivariate analysis of variance and regression. Canadian Journal of Fisheries and Aquatic Sciences, 58(3), 626-639.
%Anderson, M., & Braak, C. T. (2003). Permutation tests for multi-factorial analysis of variance. Journal of statistical computation and simulation, 73(2), 85-113.
%Wheldon, M. C., Anderson, M. J., & Johnson, B. W. (2007). Identifying treatment effects in multi-channel measurements in electroencephalographic studies: Multivariate permutation tests and multiple comparisons. Australian & New Zealand Journal of Statistics, 49(4), 397-413. 
%Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E. (2014). Permutation inference for the general linear model. NeuroImage, 92, 381-397. 
%
%VERSION: 3.0.0 alpha
%version date: 10 November 2016
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. It needs additional testing and SHOULD NOT be
%considered error free. 

%Copyright (c) 2016, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.
%This function incorporates small amounts of code from the Mass Univariate
%Toolbox, Copyright (c) 2015, David Groppe: https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox/blob/master/LICENSE

%%%%%%%%%%% REVISION LOG %%%%%%%%%%%
%
% --- As Fmax_factorial_GND ---
%(two-way ANOVA only)
%
% 0.1.0 beta     - First working version. Takes mean over electrodes rather than conducting test at each. 
% 1.0.0 beta     - Function now finds Fmax across time AND electrodes
% 1.1.0 beta     - 1) Added alpha level as optional input. 2) Changed electrodes
%                  input to match Mass Univariate Toolbox functions (i.e.,: 
%                  'include_chans' and 'exclude_chans'). 3) Added ability to
%                  use filepath to GND variable rather than GND variable in
%                  workspace.
% 1.2.0 beta     - vectorized ANOVA calculations, significantly increasing
%                  speed
% 1.3.0 beta     - added option to output result struct directly; fixed bug in
%                  spreadsheet output--now uses correct electrode name for sheets
% 1.4.0 beta     - added reproduce test option
%   1.4.1 beta   - fixed some of the code for testing the function
% 1.5.0 beta     - added requested alpha level to output
%
% --- As FmaxGND ---
%(generalized for different numbers of factors)
%
% 2.0.0 beta     - This version still only supports two-way designs.
%                  int_method parameter added with exact tests of interaction 
%                  now available for designs that allow it
% 2.1.0 beta     - Now using averaging and one-way ANOVA to calculate main
%                  effects. This simplifies the code and runs faster
% 2.2.0 beta     - Now can handle one-way or two-way ANOVAs. Many aspects
%                  of the code generalized to handle more than two-factors, 
%                  but not all. Output to spreadsheet still only works for
%                  two-way
% 2.3.0 beta     - Generalized main_effect and exact_interaction_effect 
%                  functions to handle any number of factors
%   2.3.1 beta   - Re-factored some code to use precalculated
%                  num_electrodes, num_time_pts, etc. More readable and
%                  probably slightly faster
%   2.3.2 beta   - Moved spreadsheet output to sub-function
%   2.3.3 beta   - Fixed backward incompatibilty in get_effects function
%                  core functionality now works at least back to MATLAB
%                  2012. Added assertion to make sure results in h and
%                  p-values match
%   2.3.4 beta   - Now use the same function for data reduction for both
%                  main effects and exact interactions. Fixed error that
%                  caused save_GND default to be 'no' rather than 'prompt'.
%                  A few other small changes.
% 3.0.0 alpha    - Generalizing approximate interaction to any number of
%                  factors


    %% ~~~~~PARSE INPUT~~~~~
    
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
                    electrodes = find(ismember(chan_labels, include_chans));
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
            end
        end
    end
    
    %Set defaults for missing arguments
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
        elseif sum(factor_levels > 2) < 2
            int_method = 'exact';
            fprintf('\n\nUsing restricted permutations to conduct an exact test of the interaction effect.\nSee help FmaxGND for mor information.\n\n')
        else
            int_method = 'approx';
            fprintf('\n\nAn exact test of the interaction is not possible for this design.\nUsing permutation of residuals method to conduct and approximate test.\nSee help FmaxGND for mor information.\n\n')
        end
    end
    
    %Standardize formatting
    if ischar(factor_names)
        factor_names = {factor_names};
    end
    if strcmpi(int_method, 'approximate')
        int_method = 'approx';
    end
    
    %Check for errors and problems in input
    if exist('include_chans', 'var') && exist('exclude_chans', 'var')
        error('You cannot use BOTH ''include_chans'' and ''exclude_chans'' options.');
    end
    if length(factor_names) ~= length(factor_levels)
        error('The number of factors does not match in the factor_names and factor_level inputs');
    end
    if length(factor_levels) > 2
        fprintf('\n\nWARNING: This function has not been tested extensively with designs with more than two factors. Proceed with caution!\n\n');
    end
    if strcmpi(int_method, 'exact') && sum(factor_levels > 2) > 1
        error('An exact test of the interaction is not possible if more than one factor has more than two levels. See help FmaxGND for more information.');
    end
    if isequal(factor_levels, [2, 2]) && strcmpi(int_method, 'approx')
        button = questdlg('WARNING: The type I error rate is not well-controlled by the approximate method of calculating the interaction for a 2x2 design. Are you sure you want to proceed?', 'WARNING');
        if ~strcmp(button, 'Yes')
            return;
        end
    end
    if ~strcmpi(int_method, 'approx') && ~strcmpi(int_method, 'exact')
        error('The int_method argument must be either ''approximate'' or ''exact''')
    end
    if prod(factor_levels) ~= length(bins)
        error('Number of bins does not match design.')
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

    %Convert times in ms to samples
    [~, start_sample] = min(abs( GND.time_pts - time_wind(1) ));
    [~, end_sample  ] = min(abs( GND.time_pts - time_wind(2) ));
    
    %Some useful numbers
    num_time_pts   = end_sample - start_sample + 1;
    num_subs       = size(GND.indiv_erps, 4);
    num_electrodes = length(electrodes);

    %Extract relevant data based on parameters above
    the_data = GND.indiv_erps(electrodes, start_sample:end_sample, bins, :);
    %Divide the factors into separate dimensions for factorial ANOVA
    if length(factor_levels) > 1
        the_data = reshape(the_data,[num_electrodes, num_time_pts, factor_levels, num_subs]);
    end
    
    %Figure out the effects we need to calculate
    [effects, effects_labels] = get_effects(factor_names);
    
    
    %% ~~~~~ RUN PERMUTATION ANOVAS ~~~~~
    

    test_results = cell(length(effects), 1);
    for i = 1:length(effects)
        effect = effects{i};
        if length(effects) == 1
            test_results{i} = oneway(the_data, n_perm, alpha);
        elseif length(effect) > 1 && strcmpi(int_method, 'approx')
            test_results{i} = approx_interaction_effect(the_data, n_perm, alpha);
        else
            test_results{i} = exact_factorial_effect(the_data, effect+2, n_perm, alpha);
        end
    end

    
    %% ~~~~~ ADD RESULTS STRUCT TO GND ~~~~~
    
    %Create results struct with basic parameters
    results = struct('bins', bins, ...
                     'factors', {factor_names}, ...
                     'factor_levels', factor_levels, ...
                     'time_wind', time_wind, ...
                     'used_tpt_ids', GND.time_pts(start_sample:end_sample)', ...
                     'include_chans', {{GND.chanlocs(electrodes).labels}}, ...
                     'used_chan_ids', electrodes, ...
                     'mult_comp_method', 'Fmax perm test', ...
                     'interaction_method', int_method, ...
                     'n_perm', n_perm, ...
                     'desired_alpha', alpha, ...
                     'seed_state', seed_state);
    
    %Add statistical results
    assert(length(effects) == length(test_results));
    for i = 1:length(effects)
        results.(['h_' effects_labels{i}])     = test_results{i}{1};
        results.(['p_' effects_labels{i}])     = test_results{i}{2};
        results.(['Fobs_' effects_labels{i}])  = test_results{i}{3};
        results.(['Fcrit_' effects_labels{i}]) = test_results{i}{4};
    end
                 
    %Add results struct to GND
    if ~isfield(GND, 'F_tests')
        GND.F_tests = {results};
    else
        GND.F_tests{end+1} = results;
    end
    
    %% ~~~~~ SAVE RESULTS TO DISK ~~~~~
    
    %Prompt user about saving GND
    if ~strcmpi(save_GND, 'no')
        GND = save_matmk(GND);
    end
    
    %Output to spreadsheet if requested
    if output_file
        save_to_spreadsheet(GND, test_results, factor_names, output_file);
    end

end




%% ########################################################################
%#########################    SUB-FUNCTIONS    ############################
%##########################################################################

%% 
function [effects, effects_labels] = get_effects(factor_names)

    effects = {};
    for i = 1:length(factor_names)
        effects = [effects; num2cell(nchoosek(1:length(factor_names), i), 2)]; %#ok<AGROW>
    end
    
    effects_labels = cell(size(effects));
    for i = 1:length(effects)
        if verLessThan('matlab','8.1')
            %Workaround for the fact that strjoin doesn't exist in older versions of MATLAB
            effects_labels{i} = '';
            for j = 1:length(effects{i})
                if j < length(effects{i})
                    effects_labels{i} = [effects_labels{i} factor_names{effects{i}(j)} 'X'];
                else
                    effects_labels{i} = [effects_labels{i} factor_names{effects{i}(j)}];
                end
            end
        else
            effects_labels{i} = strjoin(factor_names(effects{i}), 'X');            
        end
    end
    
end


%% 
function test_results = oneway(data, n_perm, alpha)
    
    %Make sure there's only one factor
    assert(ndims(data) == 4);
    
    %Some useful numbers
    num_electrodes = size(data, 1);
    num_time_pts   = size(data, 2);
    num_conds      = size(data, 3);
    num_subs       = size(data, 4);
    
    %Calculate degrees of freedom
    %(Always the same, so no point calculating in the loop)
    dfA   = num_conds - 1;
    dfBL  = num_subs - 1;
    dfRES = dfA * dfBL;
    
    %Perform n_perm permutations
    F_dist = NaN(n_perm, num_electrodes, num_time_pts);
    for i = 1:n_perm
        
        %Permute the data
        if i ==1
            perm_data = data;
        else
            for n = 1:num_subs
                perm_data(:, :, :, n) = data(:, :, randperm(size(data, 3)), n);
            end
        end
        
        %Calculate sums of squares
        SSyint = (sum(sum(perm_data, 3), 4).^2) / (num_conds * num_subs);
        %SSTO   = sum(sum(perm_data.^2, 3), 4) - SSyint;
        SSA    = (sum(sum(perm_data, 4).^2, 3) / num_subs) - SSyint;
        SSBL   = (sum(sum(perm_data, 3).^2, 4) / num_conds) - SSyint;
        SSRES  = sum(sum(perm_data.^2, 3), 4) - SSA - SSBL - SSyint;
        %assert(all(abs(SSTO(:) - (SSA(:) + SSBL(:) + SSRES(:))) < 1e-9));
        
        %Calculate F
        SSA(SSA < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSA/dfA) ./ (SSRES/dfRES);
        
    end
    
    %Calculate Fobs, Fmax distribution, and Fmax critical value
    F_obs = reshape(F_dist(1, :, :), [num_electrodes, num_time_pts]);
    Fmax_dist = max(max(F_dist, [], 2), [], 3);
    Fmax_dist = sort(Fmax_dist);
    Fmax_crit = Fmax_dist(ceil((1-alpha) * length(Fmax_dist)));

    %Null hypothesis test
    h = F_obs > Fmax_crit;
    
    %Calculate p-value
    p = NaN(num_electrodes, num_time_pts);
    for e = 1:num_electrodes
        for t = 1:num_time_pts
            p(e,t) = mean(Fmax_dist >= F_obs(e,t));
        end
    end
    assert(isequal(h, p<=alpha));
    
    test_results = {h, p, F_obs, Fmax_crit, Fmax_dist};
end


%%
function test_results = exact_factorial_effect(data, dims, n_perm, alpha)
%Use restricted permutation method to calculate an exact test of main
%effects and interactions in factorial designs
% 1. Reduce data to one-way design by:
%    a. averaging across any factor not involved in effect
%    b. for interactions, subtracting out all factors but one (assumes
%       factors have two levels)
% 2. Calculate one-way permutation ANOVA from reduced data
    
    %Average across factors not involved in this effect
    if length(dims) < ndims(data) - 3
        %Put the factors to average across as the initial dimensions
        num_dims = ndims(data);
        all_dims = 1:ndims(data);
        reorder = [all_dims(~ismember(all_dims, [1,2,dims,num_dims])), all_dims(ismember(all_dims, [1,2,dims,num_dims]))];
        reduced_data = permute(data, reorder);
        %Reduce all the factors to average across to a single dimension
        dim_sizes = size(reduced_data);
        num_dims_to_avg = ndims(data) - length(dims) - 3;
        reduced_data = reshape(reduced_data, [prod(dim_sizes(1:num_dims_to_avg)), dim_sizes((num_dims_to_avg+1):end)]);
        %Take the mean across that dimension
        reduced_data = mean(reduced_data, 1);
        %Get rid of the extra singleton dimension
        reduced_data = reshape(reduced_data, dim_sizes((num_dims_to_avg+1):end));
    else
        reduced_data = data;
    end
    assert(ndims(reduced_data) == length(dims) + 3)
    
    %For interactions, reduce data to one-way design via subtraction
    if length(dims) > 1
        %Get/check design structure
        dim_sizes = size(reduced_data);
        factor_levels = dim_sizes(3:(ndims(reduced_data)-1));
        assert(sum(factor_levels > 2) < 2);
        %Put the factors to subtract across as the initial dimensions
        [~, factor_order] = sort(factor_levels);
        reorder = [factor_order(1:end-1)+2, 1, 2, factor_order(end)+2 ndims(reduced_data)];
        reduced_data = permute(reduced_data, reorder);
        %Flatten data and subtract until the data is reduced to the right size
        reduced_data = reshape(reduced_data, 1, []);
        while length(reduced_data) > size(data,1)*size(data,2)*max(factor_levels)*size(data,ndims(data))
            reduced_data = reduced_data(1:2:length(reduced_data)) - reduced_data(2:2:length(reduced_data));
        end
        %Put back in regular format
        reduced_data = reshape(reduced_data, [size(data,1), size(data,2), max(factor_levels), size(data,ndims(data))]);
    end
    
    assert(ndims(reduced_data) == 4); %data has been succesfully reuduced to one-way

    %Calculate a one-way ANOVA from this reduced data
    test_results = oneway(reduced_data, n_perm, alpha);

end


%%
function save_to_spreadsheet(GND, test_results, factor_names, output_file)
%Output results to spreadsheet
    
    %THIS NEEDS A LOT OF WORK!
    
    results = GND.F_tests{end};
    
    %Check against all the current limitations . . .
    if verLessThan('matlab','8.2')
        watchit('Output to spreadsheet currently only supported for MATLAB 2013b and later. No output file written.');
        return;
    elseif ~ispc()
        watchit('Output to spreadsheet currently only supported for Windows. No output file written.');
        return;
    elseif length(results.factor_levels) ~= 2
        watchit('Output to spreadsheet currently only supported for two-way ANOVA. No output file written.')
        return;
    end
    
    warning('off','MATLAB:xlswrite:AddSheet')

    %Summary
    test_summary = {'Fmax permutation test', ''; ...
                    'GND: ', [GND.filepath GND.filename]; ...
                    'Factors: ', sprintf('%s, %s', results.factors{1}, results.factors{2}); ...
                    'Factors Levels: ', sprintf('%d, %d', results.factor_levels(1), results.factor_levels(2)); ...
                    'Number of permutations = ', results.n_perm};
    xlswrite(output_file, test_summary,'test_summary');


    %Results
    num_electrodes = length(results.include_chans);
    num_time_pts   = length(results.used_tpt_ids);
    for e = 1:num_electrodes

        output_table = table;
        warning('off','MATLAB:table:RowsAddedNewVars')

        for t = 1:num_time_pts
            output_table( sprintf('Fobs_%s', factor_names{1}),                      sprintf('t_%d', results.used_tpt_ids(t)))  =  {test_results{1}{3}(e,t)};
            output_table( sprintf('Fcrit_%s', factor_names{1}),                     sprintf('t_%d', results.used_tpt_ids(t)))  =  test_results{1}(4);
            output_table( sprintf('p_%s', factor_names{1}),                         sprintf('t_%d', results.used_tpt_ids(t)))  =  {test_results{1}{2}(e,t)};

            output_table( sprintf('Fobs_%s', factor_names{2}),                      sprintf('t_%d', results.used_tpt_ids(t)))  =  {test_results{2}{3}(e,t)};
            output_table( sprintf('Fcrit_%s', factor_names{2}),                     sprintf('t_%d', results.used_tpt_ids(t)))  =  test_results{2}(4);
            output_table( sprintf('p_%s', factor_names{2}),                         sprintf('t_%d', results.used_tpt_ids(t)))  =  {test_results{2}{2}(e,t)};

            output_table( sprintf('Fobs_%sX%s', factor_names{1}, factor_names{2}),  sprintf('t_%d', results.used_tpt_ids(t)))  =  {test_results{3}{3}(e,t)};
            output_table( sprintf('Fcrit_%sX%s', factor_names{1}, factor_names{2}), sprintf('t_%d', results.used_tpt_ids(t)))  =  test_results{3}(4);
            output_table( sprintf('p_%sX%s', factor_names{1}, factor_names{2}),     sprintf('t_%d', results.used_tpt_ids(t)))  =  {test_results{3}{2}(e,t)};         
        end

        writetable(output_table, output_file, 'WriteRowNames', true, 'Sheet', results.include_chans{e});

    end


end

%%
function test_results = approx_interaction_effect(data, n_perm, alpha)
%Use permutation of residuals method to calculate an approximate test of
%an interaction effect in a factorial design
% 1. Subtract main effects within each subject from all data points to obtain permutation residuals
% 2. Randomly permute all data points within each subject across all time
%    points
% 3. For each permutation, perform factorial ANOVA and save the largest F for the 
%    interaction effect across all time points (Fmax)
% 4. Use distribution of Fmax and F ratios obtained above to reject or fail to reject null for each
%    time point

    %Some useful numbers
    n_electrodes = size(data, 1);
    n_time_pts   = size(data, 2);
    n_subs       = size(data, ndims(data));
    
    %Subtract main effects within each subject so that the data is 
    %exchangeable under the null hypothesis for the interaction
    int_res = NaN(size(data));
    for e = 1:n_electrodes
        for t = 1:n_time_pts
            for p = 1:size(data, 3)
                for q = 1:size(data, 4)
                    for n = 1:n_subs
                        int_res(e,t,p,q,n) = data(e,t,p,q,n) - mean(squeeze(data(e,t,p,:,n))) - mean(squeeze(data(e,t,:,q,n))) + mean(mean(squeeze(data(e,t,:,:,n))));
                    end
                end
            end
        end
    end
    assert(abs(mean(int_res(:))) < 1e-9)
        
    
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    %                     GENERALIZING INT RES
    
    %??????????????????????????????????????
    
    %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    %Calculate degrees of freedom
    dim_sizes = size(data);
    dfnum = prod(dim_sizes(3:end-1)-1);
    dfdenom = prod(dim_sizes(3:end)-1);

    %Re-arrange data for permutation
    flat_data = reshape(int_res, n_electrodes, n_time_pts, [], n_subs);

    %Perform n_perm permutations
    F_dist = NaN(n_perm, n_electrodes, n_time_pts);
    flat_perm_data = NaN(size(flat_data));
    for i = 1:n_perm

        %Permute the data
        if i == 1
            perm_data = int_res;
        else
            for s = 1:n_subs
                flat_perm_data(:, :, :, s) = flat_data(:, :, randperm(size(flat_data, 3)), s);
            end
            perm_data = reshape(flat_perm_data, size(data));
        end
        
        %Calculate F at each time point and electrode combination
        
        %Numerator SS
        SSlower = 0;
        for d = 3:ndims(data)
            reduced_data = reshape(mean(perm_data, d), n_electrodes, n_time_pts, []);
            SSlower = SSlower + (sum(reduced_data.^2, 3) - (sum(reduced_data, 3).^2)/size(reduced_data, 3))*size(data, d);
        end
        collapsed_data = reshape(perm_data, n_electrodes, n_time_pts, []);
        SSdenom = sum(collapsed_data.^2, 3) - (sum(collapsed_data, 3).^2)/size(collapsed_data, 3) - SSlower;
        
        %Denominator SS
        SSlower = 0;
        num_data = mean(perm_data, ndims(perm_data));
        for d = 3:ndims(num_data)
            reduced_data = reshape(mean(num_data, d), n_electrodes, n_time_pts, []);
            SSlower = SSlower + (sum(reduced_data.^2, 3) - (sum(reduced_data, 3).^2)/size(reduced_data, 3))*size(num_data, d)*size(data, ndims(data));
        end
        collapsed_data = reshape(num_data, n_electrodes, n_time_pts, []);
        SSnum = sum(collapsed_data.^2, 3)*size(data, ndims(data)) - (sum(collapsed_data, 3).^2)/size(collapsed_data, 3) - SSlower;
        
        %Calculate F
        SSnum(SSnum < 1e-12) = 0; %Eliminates large F values that result from floating point error 
        F_dist(i, :, :) = (SSnum/dfnum) ./ (SSdenom/dfdenom);

    end
    
    %Calculate Fobs, Fmax distribution, and Fmax critical value
    F_obs = reshape(F_dist(1,:,:), n_electrodes, n_time_pts);
    Fmax_dist = max(max(F_dist, [], 2), [], 3);
    Fmax_dist = sort(Fmax_dist);
    Fmax_crit = Fmax_dist(ceil((1-alpha) * length(Fmax_dist)));

    %Null hypothesis test
    h = F_obs > Fmax_crit;
    
    %Calculate p-value
    p = NaN(n_electrodes, n_time_pts);
    for e = 1:size(data, 1)
        for t = 1:size(data, 2)
            p(e,t) = mean(Fmax_dist >= F_obs(e,t));
        end
    end
    assert(isequal(h, p<=alpha));
    
    test_results = {h, p, F_obs, Fmax_crit, Fmax_dist};

end
