%Function to get amplitude averaged across the time points and electrodes
%that were significant in a mass univariate test. Mean amplitude is
%returned for each subject and bin involved in the test.
%
%Note that care should be taken with regard to analyses conducted on such
%mean amplitudes: if analyses are not statistically indpendent of the
%original mass univariate test, results will be biased (see Kriegeskorte et
%al., Nature Neuroscience, 2009).
%
%EXAMPLE USAGE
% >> mean_amp_data = get_mean_amplitude(GND, 2, 'effect', 'Congruency');
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
%  effect           - Name of the effect of interest (e.g., 'Probability or
%                     'ProbabilityXEmotion'). Required for factorial ANOVA.
%  clust_id         - For cluster mass tests, which cluster to average
%                     across {default: all significant clusters}
%  bins             - The bins to get data from. {default: all bins used in
%                     the F-test}
%  output_file      - Name of .csv file to output results. {default: no output}
%
%OUTPUT
%  mean_amp_data    - struct with the requested data as well as bin and
%                     subject information
%
%AUTHOR: Eric Fields
%VERSION DATE: 14 March 2019
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 

%Copyright (c) 2018, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function mean_amp_data = get_mean_amplitude(GND_or_fname, test_id, varargin)

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
    
    %Load all GNDs in a GRP variable
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
    if isempty(effect)
        null_test = results.null_test;
        clust_info = results.clust_info;
    else
        null_test = results.null_test.(effect);
        clust_info = results.null_test.(effect);
    end
    
    %Assign defaults
    if isempty(bins)
        bins = results.bins;
    end
    
    %Check for errors in input
    if length(results.factor_levels) > 1 && isempty(effect)
        error('You must specify which effect you want data from. See >> help get_sub_effects');
    end
    if ~isempty(clust_id)
        if ~strcmpi(results.mult_comp_method, 'cluster mass perm test')
            error('Test %d in the GND is not a cluster mass test', test_id);
        elseif clust_id > length(clust_info.null_test)
            error('There is no cluster %d for the %s effect', clust_id, effect);
        end
    end

    %% GET SUBJECT DATA
    
    %Initialize output
    mean_amp_data = struct;
    mean_amp_data.subjects = {};
    mean_amp_data.bins = bins;
    mean_amp_data.bindesc = {GNDorGRP.bin_info(bins).bindesc};
    mean_amp_data.data = [];

    %useful numbers
    n_electrodes = length(GNDorGRP.chanlocs);
    n_time_pts   = length(GNDorGRP.time_pts);
    n_bins       = length(bins);
    
    %Find locations that are significant
    sig_locs = zeros(n_electrodes, n_time_pts);
    if isempty(clust_id)
        sig_locs(results.used_chan_ids, results.used_tpt_ids) = null_test;
    else
        sig_locs(results.used_chan_ids, results.used_tpt_ids) = ismember(clust_info.clust_ids, clust_id);
    end
    sig_locs = logical(sig_locs);

    %Extract mean amplitudes
    for i = 1:length(GNDs)
        GND = GNDs{i};
        n_subs = size(GND.indiv_erps, 4);
        if length(GNDs) == 1
            mean_amp_data.subjects = [mean_amp_data.subjects; GND.indiv_subnames'];
        else
            mean_amp_data.subjects = [mean_amp_data.subjects; repmat(GNDorGRP.group_desc(i), [n_subs,1]) GND.indiv_subnames'];
        end
        group_data = NaN(n_subs, n_bins);
        for s = 1:n_subs
            for b = 1:n_bins
                bin = bins(b);
                data = GND.indiv_erps(:, :, bin, s);
                group_data(s, b) = mean(data(sig_locs));
            end
        end
        assert(~any(isnan(group_data(:))));
        mean_amp_data.data = [mean_amp_data.data; group_data];
    end

    %% OUTPUT
    
    if output_file
        f_out = fopen(output_file, 'w');
        if length(GNDs) == 1
            fprintf(f_out, 'subject_ids,%s\n', strjoin(mean_amp_data.bindesc, ',')); %header
            for s = 1:size(mean_amp_data.data, 1)
                data_str = num2str(mean_amp_data.data(s,:),'%f,');
                data_str = data_str(1:end-1);
                fprintf(f_out, '%s,%s\n', mean_amp_data.subjects{s}, data_str);
            end
        else
            fprintf(f_out,'subject_ids,group,%s\n', strjoin(mean_amp_data.bindesc, ',')); %header
            for s = 1:size(mean_amp_data.data, 1)
                data_str = num2str(mean_amp_data.data(s,:),'%f,');
                data_str = data_str(1:end-1);
                fprintf(f_out, '%s,%s,%s\n', mean_amp_data.subjects{s,2}, mean_amp_data.subjects{s,1}, data_str);
            end
        end
        fclose(f_out);
    end
    
end
