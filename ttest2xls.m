%Output results of Mass Univariate Toolbox t-test to a spreadsheet.
%
%REQUIRED INPUTS
% GND            - A GND variable with t-test results
% test_id        - The test number within the t_tests field of the GND
%                  struct
% output_fname   - The filename for the spreadsheet that will be saved. If
%                  you don't want to save in the current working directory, 
%                  include a full filepath
%OPTIONAL INPUT
% format_output  - A boolean specifying whether to apply formatting to the 
%                  spreadsheet output. {default: true}
%
%VERSION DATE: 20 June 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 3/31/17   - First version
% 4/7/17    - Fixed small error
% 5/9/17    - Added cluster summary output; fixed error in output of time
%             point ids
% 5/15/17   - Updated to work with FDR results
% 5/17/17   - Added ability to use Python to format spreadsheet; added #
%             subjects to test summary sheet
% 5/24/17   - Updated for new xls formatting function; added optional argument
%             to specify whether to format or not
% 6/15/17   - Updated to work with xlwrite and improved output of critical
%             values
% 6/20/17   - Output for mean window analyses
function ttest2xls(GND, test_id, output_fname, format_output)
    
    %% Set-up

	%Set formatting option if no input
    if nargin < 4
        if ispc()
            format_output = true;
        else
            format_output = false;
        end
    end
    if format_output && ~ispc()
        watchit(sprintf('Spreadsheet formatting on non-Windows systems is buggy.\nSee the FMUT documentation for an explanation and possible workaround.'))
    end
    
    %Define function for writing to spreadsheet
    if ispc()
        writexls = @xlswrite;
    else
        % Add Java POI Libs to matlab javapath
        javaaddpath(fullfile(fileparts(which('Ftest2xls')), 'poi_library/poi-3.8-20120326.jar'));
        javaaddpath(fullfile(fileparts(which('Ftest2xls')), 'poi_library/poi-ooxml-3.8-20120326.jar'));
        javaaddpath(fullfile(fileparts(which('Ftest2xls')), 'poi_library/poi-ooxml-schemas-3.8-20120326.jar'));
        javaaddpath(fullfile(fileparts(which('Ftest2xls')), 'poi_library/xmlbeans-2.3.0.jar'));
        javaaddpath(fullfile(fileparts(which('Ftest2xls')), 'poi_library/dom4j-1.6.1.jar'));
        javaaddpath(fullfile(fileparts(which('Ftest2xls')), 'poi_library/stax-api-1.0.1.jar'));
        writexls = @xlwrite;
    end
    
    %Make sure we're not adding sheets to existing file
    if exist(output_fname, 'file')
        user_resp = questdlg(sprintf('WARNING: %s already exists. Do you want to overwrite it?', output_fname), 'WARNING');
        if strcmp(user_resp, 'No')
            return;
        else
            delete(output_fname)
        end
    end
    
    %Create assign t-tests results for easier reference
    results = GND.t_tests(test_id);
    
    warning('off','MATLAB:xlswrite:AddSheet')
    
    %% Test summary

    summary = {'Study', GND.exp_desc; ...
               'GND', [GND.filepath GND.filename]; ...
               'Bin', results.bin; ...
               'Bin Description', GND.bin_info(results.bin).bindesc; ...
               'Time Window', sprintf('%d - %d', round(results.time_wind(:))); ...
               'Mean window', results.mean_wind; ...
               'Electrodes', [sprintf('%s, ', results.include_chans{1:end-1}), results.include_chans{end}]; ...
               'Multiple comparisons correction method', results.mult_comp_method; ...
               'Number of permutations', results.n_perm; ...
               'Alpha or q(FDR)', results.desired_alphaORq; ...
               '# subjects', results.df+1};
    if ~strcmpi(results.mult_comp_method, 'cluster mass perm test')
        summary(end+1, :) = {'t critical value', results.crit_t(2)};
    end
    writexls(output_fname, summary, 'test summary');
    
    %% Cluster summary
    
    if strcmpi(results.mult_comp_method, 'cluster mass perm test')
        if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
            test_tobs = results.data_t;
        else
            test_tobs = GND.grands_t(results.used_chan_ids, results.used_tpt_ids, results.bin);
        end
        clust_sum = cell(11, 2);
        clust_sum{1, 1} = GND.bin_info(results.bin).bindesc;
        clust_sum{3, 1} = 'POSITIVE CLUSTERS';
        clust_sum{3, 4} = 'NEGATIVE CLUSTERS';
        for clust_type = 1:2
            row = 5;
            if clust_type == 1
                col = 1;
                num_clusters = length(results.clust_info.pos_clust_pval);
            elseif clust_type == 2
                col = 4;
                num_clusters = length(results.clust_info.neg_clust_pval);
            end
            for cluster = 1:num_clusters
                %labels
                clust_sum(row:row+8, col) = {sprintf('CLUSTER %d', cluster); ... 
                                             'cluster mass'; 'p-value'; 'spatial extent'; ...
                                             'temporal extent'; 'spatial peak'; 'temporal peak'; ...
                                             'spatial mass peak'; 'temporal mass peak'};
                %assign cluster mass and p-value;
                %get array of t-values in cluster
                if clust_type == 1
                    clust_sum{row+1, col+1} = results.clust_info.pos_clust_mass(cluster); 
                    clust_sum{row+2, col+1} = results.clust_info.pos_clust_pval(cluster);
                    clust_tobs = test_tobs;
                    clust_tobs(results.clust_info.pos_clust_ids ~= cluster) = 0;
                elseif clust_type == 2
                    clust_sum{row+1, col+1} = results.clust_info.neg_clust_mass(cluster); 
                    clust_sum{row+2, col+1} = results.clust_info.neg_clust_pval(cluster);
                    clust_tobs = test_tobs;
                    clust_tobs(results.clust_info.neg_clust_ids ~= cluster) = 0;
                end
                %spatial extent
                clust_sum{row+3, col+1} = sprintf('%s, ', results.include_chans{any(clust_tobs, 2)});
                %temporal extent
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+4, col+1} = sprintf('Mean window: %d-%d', results.time_wind(1), results.time_wind(2));
                else
                    clust_sum{row+4, col+1} = sprintf('%d - %d', ... 
                                                      GND.time_pts(min(results.used_tpt_ids(any(clust_tobs, 1)))), ... 
                                                      GND.time_pts(max(results.used_tpt_ids(any(clust_tobs, 1)))));
                end
                %Spatial and temporal peak
                [max_elec, max_timept] = find(abs(clust_tobs) == max(abs(clust_tobs(:)))); %find location of max F in cluster
                clust_sum{row+5, col+1} = results.include_chans{max_elec};
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+6, col+1} = sprintf('Mean window: %d-%d', results.time_wind(1), results.time_wind(2));
                else
                    clust_sum{row+6, col+1} = GND.time_pts(results.used_tpt_ids(max_timept));
                end
                %Spatial and temporal center (collapsed across the other
                %dimension)
                [~, max_elec_clust] = max(abs(sum(clust_tobs, 2)));
                clust_sum{row+7, col+1} = results.include_chans{max_elec_clust};
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+8, col+1} = sprintf('Mean window: %d-%d', results.time_wind(1), results.time_wind(2));
                else
                    [~, max_time_clust] = max(abs(sum(clust_tobs, 1)));
                    clust_sum{row+8, col+1} = GND.time_pts(results.used_tpt_ids(max_time_clust));
                end
                row = row+10;
            end
        end
        writexls(output_fname, clust_sum, 'cluster summary');
    end
    
    %% cluster IDs, t obs, p-values for each effect
    
    %Create headers
    chan_header = [' '; results.include_chans'];
    if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
        time_header = cell(1, size(results.time_wind, 1));
        for t = 1:size(results.time_wind, 1)
            time_header{t} = sprintf('%d-%d', results.time_wind(t,1), results.time_wind(t,2));
        end 
    else
        time_header = num2cell(GND.time_pts(results.used_tpt_ids));
    end
        
    %Cluster_ids
    if strcmpi(results.mult_comp_method, 'cluster mass perm test')
        pos_clust_ids = [chan_header, [time_header; num2cell(results.clust_info.pos_clust_ids)]];
        neg_clust_ids = [chan_header, [time_header; num2cell(results.clust_info.neg_clust_ids)]];
        writexls(output_fname, pos_clust_ids, 'pos_clust_IDs');
        writexls(output_fname, neg_clust_ids, 'neg_clust_IDs');
    end

    %t observed
    if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
        t_obs = num2cell(results.data_t);
    else
        t_obs = num2cell(GND.grands_t(results.used_chan_ids, results.used_tpt_ids, results.bin));
    end
    t_obs_table = [chan_header, [time_header; t_obs]];
    writexls(output_fname, t_obs_table, 't_obs');

    %p-values
    if ~strcmpi(results.mult_comp_method, 'bky')
        adj_pvals = [chan_header, [time_header; num2cell(results.adj_pval)]];
        writexls(output_fname, adj_pvals, 'adj_pvals');
    end
    
    %Try to format
    if format_output
		try
			format_xls(output_fname)
		catch ME
			watchit(sprintf('Could not format spreadsheet because of the following error:\n%s', ME.message))
		end
    end
	
end
