%Output results of F-test to a spreadsheet.
%
%REQUIRED INPUTS
% GND            - A GND variable with F-test results
% test_id        - The test number within the F_tests field of the GND
%                  struct
% output_fname   - The filename for the spreadsheet that will be saved. If
%                  you don't want to save in the current working directory, 
%                  include a full filepath
%OPTIONAL INPUT
% format_output     - A boolean specifying whether to apply formatting to the 
%                  spreadsheet output. {default: true}
%
%VERSION DATE: 20 June 2017
%AUTHOR: Eric Fields, Tufts University (Eric.Fields@tufts.edu)
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 3/31/17   - Moved to separate function from FmaxGND and FclustGND
% 4/7/17    - Now works with one factor F-test; fixed problem with sheetnames
%             being too long due to long effect names; fixed a few other small
%             errors
% 5/9/17    - Added cluster summary sheet
% 5/15/17   - Updated work with FDR results
% 5/17/17   - Added ability to use Python to format spreadsheet; added #
%             subjects to test summary sheet
% 5/24/17   - Updated for new xls formatting function; added optional argument
%             to specify whether to format or not
% 6/15/17   - Updated to use xlwrite and better output of critical values
% 6/20/17   - Output for mean window analyses

function Ftest2xls(GND, test_id, output_fname, format_output)
    
    %% Set-up

    %Set formatting option if no input
	if nargin < 4 %#ok<ALIGN>
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
    
    %Create some variables for easier reference
    results = GND.F_tests(test_id);
    [~, effects_labels] = get_effects(results.factors);
    if length(effects_labels) ==1
        n_subs = results.df(2)/results.df(1) + 1;
    else
        n_subs = results.df.(effects_labels{1})(2)/results.df.(effects_labels{1})(1) + 1;
    end
    
    warning('off','MATLAB:xlswrite:AddSheet')

    %% Summary sheet
    
    fact_levels = sprintf('%d X ', results.factor_levels);
    fact_levels = fact_levels(1:end-2);
    summary = {'Study', GND.exp_desc; ...
               'GND', [GND.filepath GND.filename]; ...
               'Bins', sprintf('%d ', results.bins); ...
               'Factors', [sprintf('%s X ', results.factors{1:end-1}), results.factors{end}]; ...
               'Factor_levels', fact_levels; ...
               'Time Window', sprintf('%d-%d ', results.time_wind'); ...
               'Mean window', results.mean_wind; ...
               'Electrodes', [sprintf('%s, ', results.include_chans{1:end-1}), results.include_chans{end}]; ...
               'Multiple comparisons correction method', results.mult_comp_method; ...
               'Interaction method', results.interaction_method; ...
               'Number of permutations', results.n_perm; ...
               'Alpha or q(FDR)', results.desired_alphaORq; ...
               '# subjects', n_subs};
    if ~strcmpi(results.mult_comp_method, 'cluster mass perm test') 
        if ~isstruct(results.F_crit)
            summary(end+1, :) = {'F critical value', results.F_crit};
        else
            for i = 1:length(effects_labels)
                summary(end+1, :) = {sprintf('%s F critical value', effects_labels{i}), ...
                                     results.F_crit.(effects_labels{i})}; %#ok<AGROW>
            end
        end
    end
    writexls(output_fname, summary, 'test summary');
    
    %% Cluster summary sheet
    
    if strcmpi(results.mult_comp_method, 'cluster mass perm test')
        clust_sum = cell(11, 3*length(effects_labels)-1);
        col = 1;
        for i = 1:length(effects_labels)
            effect = effects_labels{i};
            clust_sum{1, col} = effect;
            row = 3;
            if length(effects_labels) == 1
                num_clusters = length(results.clust_info.pval);
            else
                num_clusters = length(results.clust_info.(effect).pval);
            end
            for cluster = 1:num_clusters
                %labels
                clust_sum(row:row+8, col) = {sprintf('CLUSTER %d', cluster); ... 
                                             'cluster mass'; 'p-value'; 'spatial extent'; ...
                                             'temporal extent'; 'spatial peak'; 'temporal peak'; ...
                                             'spatial mass peak'; 'temporal mass peak'};
                %assign cluster mass and p-value;
                %get array of F-values in cluster
                if length(effects_labels) == 1
                    clust_sum{row+1, col+1} = results.clust_info.clust_mass(cluster); 
                    clust_sum{row+2, col+1} = results.clust_info.pval(cluster);
                    clust_Fobs = results.F_obs;
                    clust_Fobs(results.clust_info.clust_ids ~= cluster) = 0;
                else
                    clust_sum{row+1, col+1} = results.clust_info.(effect).clust_mass(cluster); 
                    clust_sum{row+2, col+1} = results.clust_info.(effect).pval(cluster);
                    clust_Fobs = results.F_obs.(effect);
                    clust_Fobs(results.clust_info.(effect).clust_ids ~= cluster) = 0;
                end
                %spatial extent
                clust_sum{row+3, col+1} = sprintf('%s, ', results.include_chans{any(clust_Fobs, 2)});
                %temporal extent
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+4, col+1} = sprintf('Mean window: %d-%d', results.time_wind(1), results.time_wind(2));
                else
                    clust_sum{row+4, col+1} = sprintf('%d - %d', ...
                                                      GND.time_pts(min(results.used_tpt_ids(any(clust_Fobs, 1)))), ... 
                                                      GND.time_pts(max(results.used_tpt_ids(any(clust_Fobs, 1)))));
                end
                %Spatial and temporal peak
                [max_elec, max_timept] = find(clust_Fobs == max(clust_Fobs(:))); %find location of max F in cluster
                clust_sum{row+5, col+1} = results.include_chans{max_elec};
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+6, col+1} = sprintf('Mean window: %d-%d', results.time_wind(1), results.time_wind(2));
                else
                    clust_sum{row+6, col+1} = GND.time_pts(results.used_tpt_ids(max_timept));
                end
                %Spatial and temporal center (collapsed across the other
                %dimension)
                [~, max_elec_clust] = max(sum(clust_Fobs, 2));
                clust_sum{row+7, col+1} = results.include_chans{max_elec_clust};
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+8, col+1} = sprintf('Mean window: %d-%d', results.time_wind(1), results.time_wind(2));
                else
                    [~, max_time_clust] = max(sum(clust_Fobs, 1));
                    clust_sum{row+8, col+1} = GND.time_pts(results.used_tpt_ids(max_time_clust));
                end
                row = row+10;
            end
            col = col+3;
        end
        writexls(output_fname, clust_sum, 'cluster summary');
    end
    
    %% Results: Cluster IDs, F obs, p-values
    
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
    
    %Write sheets
    if length(effects_labels) == 1
        
        %Sheet names are limited to 31 characters
        if length(effects_labels{1}) > 19
            effects_labels{1} = effects_labels{1}(1:19);
        end
        
        %Cluster_ids
        if strcmpi(results.mult_comp_method, 'cluster mass perm test')
            %cluster ids
            clust_ids = [chan_header, [time_header; num2cell(results.clust_info.clust_ids)]];
            writexls(output_fname, clust_ids, sprintf('%s_clust_IDs', effects_labels{1}));
        end

        %F_obs
        F_obs_table = [chan_header, [time_header; num2cell(results.F_obs)]];
        writexls(output_fname, F_obs_table, sprintf('%s_F_obs', effects_labels{1}));

        %p-values
        if ~strcmpi(results.mult_comp_method, 'bky')
            adj_pvals = [chan_header, [time_header; num2cell(results.adj_pval)]];
            writexls(output_fname, adj_pvals, sprintf('%s_adj_pvals', effects_labels{1}));
        end
        
    else
        for i = 1:length(effects_labels)
            
            %Sheet names are limited to 31 characters
            if length(effects_labels{i}) > 19
                effects_labels{i} = effects_labels{i}(1:19);
            end

            %Cluster_ids
            if strcmpi(results.mult_comp_method, 'cluster mass perm test')
                clust_ids = [chan_header, [time_header; num2cell(results.clust_info.(effects_labels{i}).clust_ids)]];
                writexls(output_fname, clust_ids, sprintf('%s_clust_IDs', effects_labels{i}));
            end

            %F_obs
            F_obs_table = [chan_header, [time_header; num2cell(results.F_obs.(effects_labels{i}))]];
            writexls(output_fname, F_obs_table, sprintf('%s_F_obs', effects_labels{i}));

            %p-values
            if ~strcmpi(results.mult_comp_method, 'bky')
                adj_pvals = [chan_header, [time_header; num2cell(results.adj_pval.(effects_labels{i}))]];
                writexls(output_fname, adj_pvals, sprintf('%s_adj_pvals', effects_labels{i}));
            end

        end
    end
    
    %Try to format
	if format_output
		try
			format_xls(output_fname)
		catch ME
			watchit(sprintf('Could not format spreadsheet because of the following error:\n%s\n', ME.message))
		end
	end
    
end
