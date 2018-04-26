%Function to get effect data for each subject based on a mass univariate
%test. Returns an array with a voltage value for each subject and each bin
%averaged across the locations that were significnat in the mass univariate
%analysis.
%
%EXAMPLE USAGE
% >> sub_data = get_sub_effects1(GND, 2, 'effect', 'Congruency')
%
%REQUIRED INPUTS
% GND_or_fname       - A Mass Univariate Toolbox GND or GRP struct or a string
%                      containing a filename of a GND or GRP structure that 
%                      has been saved to disk (with full path if not in current
%                      working directory)
%  test_id           - [integer] The index # of the F-test results 
%                      stored in the GND/GRP/specGND variable that you wish to visualize.  
%                      To see what test results are available, look at 
%                      the "F_tests" field of your variable (e.g., GND.F_tests)
%
%OPTIONAL INPUTS:
%  effect           - The effect within the test to plot. Required for
%                     factorial ANOVA.
%  clust_id         - For cluster mass tests, which cluster to average
%                     across {default: all significant clusters}
%  bins             - The bins to get data from. {default: all bins used in
%                     the F-test}
%  output_file      - Name of .xlsx file to output results. {default: no output}
%
%OUTPUT
%  sub_data        - a subjects x bins matrix of the mean amplitude across
%                    all locations significant in the F-test
%
%AUTHOR: Eric Fields
%VERSION DATE: 24 April 2018
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2018, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function sub_data = get_sub_effects(GND_or_fname, test_id, varargin)

    %% PARSE INPUT
    
    p=inputParser;
    p.addRequired('GND_or_fname', @(x) ischar(x) || isstruct(x));
    p.addRequired('test_id', @(x) isnumeric(x) && length(x)==1);
    p.addParameter('effect', [], @(x) ischar(x));
    p.addParameter('clust_id', [], @(x) isnumeric(x));
    p.addParameter('bins', [], @(x) isnumeric(x));
    p.addParameter('output_file', false);
    p.parse(GND_or_fname, test_id, varargin{:});
    
    %Assign GND
    if ischar(GND_or_fname) && exists(GND_or_fname, 'file')
        load(GND_or_fname, '-mat'); %#ok<LOAD>
        if exists('GRP', 'var')
            GNDorGRP = GRP;
            clear GRP;
        elseif exists('GND', 'var')
            GNDorGRP = GND; %#ok<NODEF>
        else
            error('Did not find a GND or GRP variable in %s', GND_or_fname);
        end
    elseif isstruct(GND_or_fname)
        GNDorGRP = GND_or_fname;
    else
        error('The GND variable provided does not seem to be a valid GND struct or filepath to a GND struct.');
    end
    
    if isfield(GNDorGRP, 'GND_fnames')
        GNDs = {};
        for i = 1:length(GNDorGRP.GND_fnames)
            load(GNDorGRP.GND_fnames{i}, '-mat')
            GNDs{i} = GND; %#ok<AGROW>
        end
        clear GND
    else
        GNDs = {GNDorGRP};
    end
    
    %Assign variables
    effect = p.Results.effect;
    clust_id = p.Results.clust_id;
    bins = p.Results.bins;
    output_file = p.Results.output_file;
    results = GNDorGRP.F_tests(test_id);
    
    %Assign defaults
    if isempty(bins)
        bins = results.bins;
    end
    
    %Check for errors in input
    if length(results.factor_levels) > 1 && isempty(effect)
        error('You must specify which effect you want data from. See >> help get_sub_effects');
    end
    if ~isempty(clust_id) && ~strcmpi(results.mult_comp_method, 'cluster mass perm test')
        error('Test %d in the GND is not a cluster mass test', test_id);
    elseif clust_id > length(results.clust_info.(effect).null_test)
        error('There is no cluster %d for the %s effect', clust_id, effect);
    end

    %% GET SUBJECT DATA

    %useful numbers
    n_electrodes = length(GNDorGRP.chanlocs);
    n_time_pts = length(GNDorGRP.time_pts);
    n_bins = length(bins);
    
    %Find locations that are significant
    sig_locs = zeros(n_electrodes, n_time_pts);
    if isempty(clust_id)
        sig_locs(results.used_chan_ids, results.used_tpt_ids) = results.null_test.(effect);
    else
        sig_locs(results.used_chan_ids, results.used_tpt_ids) = ismember(results.clust_info.(effect).clust_ids, clust_id);
    end
    sig_locs = logical(sig_locs);

    %Extract mean amplitudes
    for i = 1:length(GNDs)
        GND = GNDs{i};
        if i == 1
            if length(GNDs) > 1
                sub_data = ['subject' 'group' {GND.bin_info(bins).bindesc}];
            else
                sub_data = ['subject' {GND.bin_info(bins).bindesc}];
            end
        end
        n_subs = size(GND.indiv_erps, 4);
        group_data = NaN(n_subs, n_bins);
        for s = 1:n_subs
            for b = 1:n_bins
                bin = bins(b);
                data = GND.indiv_erps(:, :, bin, s);
                group_data(s, b) = mean(data(sig_locs));
            end
        end
        assert(~any(isnan(group_data(:))));
        if length(GNDs) > 1
            sub_data = [sub_data; GND.indiv_subnames' repmat(GNDorGRP.group_desc(i), [n_subs, 1]) num2cell(group_data)]; %#ok<AGROW>
        else
            sub_data = [sub_data; GND.indiv_subnames' num2cell(group_data)]; %#ok<AGROW>
        end
    end

    %% OUTPUT
    
    if output_file
        xlswrite(output_file, sub_data);
    end
    
end
