%Summarize results of F-test at the command window
%
%EXAMPLE USAGE
% >> report_results(GND, 1)
%
%REQUIRED INPUTS
% GND        - A GND or GRP variable with F-test results
% test_id    - The test number of the results to report within the GND or
%              GRP variable
% 
%VERSION DATE: 14 July 2014
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 7/14/17  - Moved to separate function

function report_results(GND, test_id)

    switch GND.F_tests(test_id).mult_comp_method
        case 'Fmax perm test'
            report_Fmax(GND, test_id)
        case 'cluster mass perm test'
            report_clust(GND, test_id)
        otherwise
            if any(strcmpi(GND.F_tests(test_id), {'bh', 'by', 'bky'}))
                report_fdr(GND, test_id)
            end
    end
            
end


function report_Fmax(GND, test_id)
    
    results = GND.F_tests(test_id);
    [effects, effects_labels] = get_effects(results.factors);
    if isstruct(results.F_obs)
        assert(all(strcmpi(fieldnames(results.F_obs), effects_labels)))
    end
    
    fprintf('\n##### RESULTS #####\n');
    for i = 1:length(effects)
        fprintf('\n%s effect\n', effects_labels{i});
        if length(effects) == 1
            if any(results.null_test(:))
                fprintf('Critical F-value: %.3f\n', results.F_crit);
                fprintf('That corresponds to a test-wise alpha level of %.3f.\n', ...
                        fcdf(results.F_crit, results.df(1), results.df(2), 'upper'));
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    for t = 1:size(results.time_wind, 1)
                        fprintf('Significant electrodes for time window %d - %d: ', results.time_wind(t, 1), results.time_wind(t, 2));
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
                fprintf('All significant corrected p-values are between %f and %f.\n', ... 
                         max(results.adj_pval(results.adj_pval <= .05)), ... 
                         min(results.adj_pval(:)));
            else
                fprintf('NO significant time points or electrodes.\n');
            end
        else
            if any(results.null_test.(effects_labels{i})(:))
                fprintf('Critical F-value: %.3f\n', results.F_crit.(effects_labels{i}));
                fprintf('That corresponds to a test-wise alpha level of %.3f.\n', ...
                        fcdf(results.F_crit.(effects_labels{i}), results.df.(effects_labels{i})(1), results.df.(effects_labels{i})(2), 'upper'));
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    for t = 1:size(results.time_wind, 1)
                        fprintf('Significant electrodes for time windonw %d - %d: ', results.time_wind(t, 1), results.time_wind(t, 2));
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
                fprintf('All significant corrected p-values are between %f and %f.\n', ... 
                         max(results.adj_pval.(effects_labels{i})(results.adj_pval.(effects_labels{i}) <= .05)), ... 
                         min(results.adj_pval.(effects_labels{i})(:)));
            else
                fprintf('NO significant time points or electrodes.\n');
            end
        end
    end

end


function report_clust(GND, test_id)

    results = GND.F_tests(test_id);
    [effects, effects_labels] = get_effects(results.factors);
    if isstruct(results.F_obs)
        assert(all(strcmpi(fieldnames(results.F_obs), effects_labels)))
    end

    fprintf('\n##### RESULTS #####\n\n');
    if length(effects) == 1
            fprintf('%s effect\n', effects_labels{1});
            fprintf('# of clusters found: %d\n', length(results.clust_info.null_test));
            fprintf('# of significant clusters found: %d\n', sum(results.clust_info.null_test));
            if sum(results.clust_info.null_test)
                fprintf('Significant cluster p-values range from %f to %f.\n', ...
                        max(results.clust_info.pval(results.clust_info.null_test)), ...
                        min(results.clust_info.pval(results.clust_info.null_test)));
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    for t = 1:size(results.time_wind, 1)
                        fprintf('Electrodes in a significant cluster for time window %d - %d: ', results.time_wind(t, 1), results.time_wind(t, 2));
                        fprintf('%s ', results.include_chans{results.null_test(:, t)});
                        fprintf('\n');
                    end
                    fprintf('\n');
                else
                    fprintf('Electrodes and time points included in a significant cluster:\n');
                    for t = 1:length(results.used_tpt_ids)
                        if any(results.null_test(:, t))
                            fprintf('%d ms, electrode(s): ', GND.time_pts(results.used_tpt_ids(t)));
                            fprintf('%s ', results.include_chans{results.null_test(:, t)});
                            fprintf('\n');
                        end
                    end
                    fprintf('\n');
                end
            else
                fprintf('All p-values >= %f\n\n', min(results.clust_info.pval));
            end
    else
        for i = 1:length(effects)
            fprintf('%s effect\n', effects_labels{i});
            fprintf('# of clusters found: %d\n', length(results.clust_info.(effects_labels{i}).null_test));
            fprintf('# of significant clusters found: %d\n', sum(results.clust_info.(effects_labels{i}).null_test));
            if sum(results.clust_info.(effects_labels{i}).null_test)
                fprintf('Significant cluster p-values range from %f to %f.\n', ...
                        max(results.clust_info.(effects_labels{i}).pval(results.clust_info.(effects_labels{i}).null_test)), ...
                        min(results.clust_info.(effects_labels{i}).pval(results.clust_info.(effects_labels{i}).null_test)));
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    for t = 1:size(results.time_wind, 1)
                        fprintf('Electrodes in a significant cluster for time window %d - %d: ', results.time_wind(t, 1), results.time_wind(t, 2));
                        fprintf('%s ', results.include_chans{results.null_test.(effects_labels{i})(:, t)});
                        fprintf('\n');
                    end
                    fprintf('\n');
                else
                    fprintf('Electrodes and time points included in a significant cluster:\n');
                    for t = 1:length(results.used_tpt_ids)
                        if any(results.null_test.(effects_labels{i})(:, t))
                            fprintf('%d ms, electrode(s): ', GND.time_pts(results.used_tpt_ids(t)));
                            fprintf('%s ', results.include_chans{results.null_test.(effects_labels{i})(:, t)});
                            fprintf('\n');
                        end
                    end
                    fprintf('\n');
                end
            else
                fprintf('All p-values >= %f.\n\n', min(results.clust_info.(effects_labels{i}).pval));
            end
        end
    end

end


function report_fdr(GND, test_id)

    results = GND.F_tests(test_id);
    [effects, effects_labels] = get_effects(results.factors);
    if isstruct(results.F_obs)
        assert(all(strcmpi(fieldnames(results.F_obs), effects_labels)))
    end

    fprintf('\n##### RESULTS #####\n');
    if length(effects) == 1
            fprintf('\n%s effect\n', effects_labels{1});
            if any(results.null_test(:))
                fprintf('Critical F-value: %.3f\n', results.F_crit);
                fprintf('That corresponds to a test-wise alpha level of %.3f.\n', ...
                        fcdf(results.F_crit, results.df(1), results.df(2), 'upper'));
                fprintf('Total number of significant differences: %d\n', sum(results.null_test(:)));
                fprintf('Estimated upper bound on expected number of false discoveries: %.1f.\n', sum(results.null_test(:))*q);
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    for t = 1:size(results.time_wind, 1)
                        fprintf('Significant electrodes for time window %d - %d: ', results.time_wind(t, 1), results.time_wind(t, 2));
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
                    fprintf('All significant corrected p-values are between %f and %f.\n', ...
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
                fprintf('That corresponds to a test-wise alpha level of %.3f.\n', ...
                        fcdf(results.F_crit.(effects_labels{i}), results.df.(effects_labels{i})(1), results.df.(effects_labels{i})(2), 'upper'));
                fprintf('Total number of significant differences: %d\n', sum(results.null_test.(effects_labels{i})(:)));
                fprintf('Estimated upper bound on expected number of false discoveries: %.1f\n', sum(results.null_test.(effects_labels{i})(:))*q);
                if strcmpi(results.mean_wind, 'yes') || strcmpi(results.mean_wind, 'y')
                    for t = 1:size(results.time_wind, 1)
                        fprintf('Significant electrodes for time windonw %d - %d: ', results.time_wind(t, 1), results.time_wind(t, 2));
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
                    fprintf('All significant corrected p-values are between %f and %f.\n', ...
                            max(results.adj_pval.(effects_labels{i})(results.adj_pval.(effects_labels{i}) <= .05)), ... 
                            min(results.adj_pval.(effects_labels{i})(:)));
                end
            else
                fprintf('NO significant time points or electrodes.\n\n');
            end
        end
    end
        
end
