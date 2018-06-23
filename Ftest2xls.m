%Output results of F-test to a spreadsheet.
%
%EXAMPLE USAGE
% >> Ftest2xls(GND, 1, 'example.xlsx')
%
%REQUIRED INPUTS
% GND            - A GND or GRP variable with F-test results
% test_id        - The test number within the F_tests field of the GND
%                  struct
% output_fname   - The filename for the spreadsheet that will be saved. If
%                  you don't want to save in the current working directory, 
%                  include a full filepath
%OPTIONAL INPUT
% format_output  - A boolean specifying whether to apply formatting to the 
%                  spreadsheet output. {default: true}
%
%VERSION DATE: 23 June 2018
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.


function Ftest2xls(GND, test_id, output_fname, format_output)
    
    %% Set-up
    
    global VERBLEVEL;
    if isempty(VERBLEVEL)
        VERBLEVEL = 2;
    end

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
    
    %Define function for writing to spreadsheet and update Java class path
    %if necessary
    if ispc()
        writexls = @xlswrite;
    else 
        %Get full path for POI .jar files not on static Java class path
        [~, missing_poi_files] = get_poi_paths();
        %Add POI file to dynamic Java class path if they aren't on the
        %static path
        if ~isempty(missing_poi_files)
            if VERBLEVEL
                watchit(sprintf(['POI libraries are not on the static Java class path.\n', ...
                                 'This can delete global variables and create other problems.\n', ...
                                 'You can avoid this by running add_poi_path.']));
            end
            %Add each .jar file to the dynamic path
            add_poi_path('-dynamic');
        end
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
    n_subs = sum(results.group_n);
    
    warning('off', 'MATLAB:xlswrite:AddSheet')

    %% Summary sheet
    
    fact_levels = sprintf('%d X ', results.factor_levels);
    fact_levels = fact_levels(1:end-2);
    if iscell(results.use_groups)
        use_groups = [sprintf('%s, ', results.use_groups{1:end-1}), results.use_groups{end}];
    elseif isnan(results.use_groups)
        use_groups = 'N/A (no between-subjects factor)';
    else
        error('Cannot interpret F_tests.use_groups');
    end
    summary = {'Study', GND.exp_desc; ...
               'GND', [GND.filepath GND.filename]; ...
               'Groups', use_groups; ...
               '# Subjects in group', [sprintf('%d  ', results.group_n(1:end-1)), num2str(results.group_n(end))];
               'Bins', sprintf('%d ', results.bins); ...
               'Factors', [sprintf('%s X ', results.factors{1:end-1}), results.factors{end}]; ...
               'Factor_levels', fact_levels; ...
               'Time Window', sprintf('%.0f - %.0f ', results.time_wind'); ...
               'Mean window', results.mean_wind; ...
               'Electrodes', [sprintf('%s, ', results.include_chans{1:end-1}), results.include_chans{end}]; ...
               'Multiple comparisons correction method', results.mult_comp_method; ...
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
                    clust_sum{row+4, col+1} = sprintf('Mean window: %.0f - %.0f', results.time_wind(1), results.time_wind(2));
                else
                    clust_sum{row+4, col+1} = sprintf('%.0f - %.0f', ...
                                                      GND.time_pts(min(results.used_tpt_ids(any(clust_Fobs, 1)))), ... 
                                                      GND.time_pts(max(results.used_tpt_ids(any(clust_Fobs, 1)))));
                end
                %Spatial and temporal peak
                [max_elec, max_timept] = find(clust_Fobs == max(clust_Fobs(:))); %find location of max F in cluster
                clust_sum{row+5, col+1} = results.include_chans{max_elec};
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+6, col+1} = sprintf('Mean window: %.0f - %.0f', results.time_wind(1), results.time_wind(2));
                else
                    clust_sum{row+6, col+1} = sprintf('%.0f', GND.time_pts(results.used_tpt_ids(max_timept)));
                end
                %Spatial and temporal center (collapsed across the other
                %dimension)
                [~, max_elec_clust] = max(sum(clust_Fobs, 2));
                clust_sum{row+7, col+1} = results.include_chans{max_elec_clust};
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    clust_sum{row+8, col+1} = sprintf('Mean window: %.0f - %.0f', results.time_wind(1), results.time_wind(2));
                else
                    [~, max_time_clust] = max(sum(clust_Fobs, 1));
                    clust_sum{row+8, col+1} = sprintf('%0.f', GND.time_pts(results.used_tpt_ids(max_time_clust)));
                end
                row = row+10;
            end
            col = col+3;
        end
        writexls(output_fname, clust_sum, 'cluster summary');
    end
    
    %% Cluster IDs, F obs, p-values for each effect
    
    %Create headers
    chan_header = [' '; results.include_chans'];
    if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
        time_header = cell(1, size(results.time_wind, 1));
        for t = 1:size(results.time_wind, 1)
            time_header{t} = sprintf('%.0f - %.0f', results.time_wind(t,1), results.time_wind(t,2));
        end 
    else
        time_header = num2cell(GND.time_pts(results.used_tpt_ids));
    end
    
    %Write sheets
    if length(effects_labels) == 1
        
        %Sheet names are limited to 31 characters
        if length(effects_labels{1}) > 19
            sheet_label = effects_labels{1}(1:19);
        else
            sheet_label = effects_labels{1};
        end
        
        %Cluster_ids
        if strcmpi(results.mult_comp_method, 'cluster mass perm test')
            %cluster ids
            clust_ids = [chan_header, [time_header; num2cell(results.clust_info.clust_ids)]];
            writexls(output_fname, clust_ids, sprintf('%s_clust_IDs', sheet_label));
        end

        %F_obs
        F_obs_table = [chan_header, [time_header; num2cell(results.F_obs)]];
        writexls(output_fname, F_obs_table, sprintf('%s_F_obs', sheet_label));

        %p-values
        if ~strcmpi(results.mult_comp_method, 'bky')
            adj_pvals = [chan_header, [time_header; num2cell(results.adj_pval)]];
            writexls(output_fname, adj_pvals, sprintf('%s_adj_pvals', sheet_label));
        end
        
    else
        for i = 1:length(effects_labels)
            
            %Sheet names are limited to 31 characters
            if length(effects_labels{i}) > 19
                sheet_label = effects_labels{i}(1:19);
            else
                sheet_label = effects_labels{i};
            end

            %Cluster_ids
            if strcmpi(results.mult_comp_method, 'cluster mass perm test')
                clust_ids = [chan_header, [time_header; num2cell(results.clust_info.(effects_labels{i}).clust_ids)]];
                writexls(output_fname, clust_ids, sprintf('%s_clust_IDs', sheet_label));
            end

            %F_obs
            F_obs_table = [chan_header, [time_header; num2cell(results.F_obs.(effects_labels{i}))]];
            writexls(output_fname, F_obs_table, sprintf('%s_F_obs', sheet_label));

            %p-values
            if ~strcmpi(results.mult_comp_method, 'bky')
                adj_pvals = [chan_header, [time_header; num2cell(results.adj_pval.(effects_labels{i}))]];
                writexls(output_fname, adj_pvals, sprintf('%s_adj_pvals', sheet_label));
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
